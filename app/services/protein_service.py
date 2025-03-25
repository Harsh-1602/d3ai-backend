import logging
import uuid
import httpx
from datetime import datetime
from typing import List, Optional, Dict, Any, Tuple
from app.db.mongodb import db
from app.schemas.docking import Protein, ProteinCreate
from app.schemas.drug import Drug
from app.core.config import settings
import asyncio

logger = logging.getLogger(__name__)

class ProteinService:
    """Service for protein-related operations."""
    
    async def get_protein_by_id(self, protein_id: str) -> Optional[Protein]:
        """Get a specific protein by ID."""
        logger.info(f"Getting protein by ID: {protein_id}")
        
        # Query database
        proteins_collection = db["proteins"]
        doc = await proteins_collection.find_one({"_id": protein_id})
        
        if not doc:
            return None
        
        # Get disease name
        disease_name = None
        if doc.get("disease_id"):
            diseases_collection = db["diseases"]
            disease_doc = await diseases_collection.find_one({"_id": doc.get("disease_id")})
            if disease_doc:
                disease_name = disease_doc.get("name")
        
        # Create Protein object
        protein = Protein(
            id=str(doc.get("_id")),
            name=doc.get("name"),
            pdb_id=doc.get("pdb_id"),
            sequence=doc.get("sequence"),
            disease_id=doc.get("disease_id"),
            description=doc.get("description"),
            binding_site=doc.get("binding_site"),
            created_at=doc.get("created_at", datetime.now()),
            updated_at=doc.get("updated_at", datetime.now()),
            disease_name=disease_name
        )
        
        return protein
    
    async def create_protein(self, protein: ProteinCreate) -> Protein:
        """Create a new protein."""
        logger.info(f"Creating new protein: {protein.name}")
        
        # Create document
        protein_id = str(uuid.uuid4())
        now = datetime.now()
        
        protein_doc = {
            "_id": protein_id,
            "name": protein.name,
            "pdb_id": protein.pdb_id,
            "sequence": protein.sequence,
            "disease_id": protein.disease_id,
            "description": protein.description,
            "binding_site": protein.binding_site,
            "created_at": now,
            "updated_at": now
        }
        
        # Insert into database
        proteins_collection = db["proteins"]
        await proteins_collection.insert_one(protein_doc)
        
        # Get disease name
        disease_name = None
        if protein.disease_id:
            diseases_collection = db["diseases"]
            disease_doc = await diseases_collection.find_one({"_id": protein.disease_id})
            if disease_doc:
                disease_name = disease_doc.get("name")
        
        # Return created protein
        return Protein(
            id=protein_id,
            name=protein.name,
            pdb_id=protein.pdb_id,
            sequence=protein.sequence,
            disease_id=protein.disease_id,
            description=protein.description,
            binding_site=protein.binding_site,
            created_at=now,
            updated_at=now,
            disease_name=disease_name
        )
    
    async def get_chembl_target_id(self, uniprot_id: str) -> Optional[str]:
        """Fetch ChEMBL Target ID using UniProt ID."""
        if not uniprot_id:
            logger.warning("Cannot fetch ChEMBL target ID: No UniProt ID provided")
            return None
        
        logger.info(f"Getting ChEMBL target ID for UniProt ID: {uniprot_id}")
        
        # Use our enhanced get_protein_external_links method which includes multiple fallbacks
        external_links = await self.get_protein_external_links(uniprot_id)
        
        # Check if we found a ChEMBL ID in the external links
        if 'ChEMBL' in external_links:
            chembl_id = external_links['ChEMBL']
            logger.info(f"Found ChEMBL ID {chembl_id} for UniProt ID {uniprot_id}")
            return chembl_id
            
        logger.warning(f"No ChEMBL ID found for UniProt ID {uniprot_id} after all attempts")
        return None
    
    async def get_drugs_by_target_id(self, target_id: str, limit: int = 20) -> List[Dict[str, Any]]:
        """Get drugs that interact with a specific target (protein) from ChEMBL."""
        logger.info(f"Getting drugs for ChEMBL target ID: {target_id}")
        
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}"
        headers = {"Accept": "application/json"}
        
        try:
            async with httpx.AsyncClient() as client:
                response = await client.get(url, headers=headers)
                
                if response.status_code == 200:
                    data = response.json()
                    mechanisms = data.get("activities", [])
                    
                    # Get molecule details for each mechanism
                    drugs = []
                    for mechanism in mechanisms:
                        molecule_chembl_id = mechanism.get("molecule_chembl_id")
                        if not molecule_chembl_id:
                            continue
                            
                        molecule_url = f"{settings.CHEMBL_API_BASE_URL}/molecule/{molecule_chembl_id}.json"
                        molecule_response = await client.get(molecule_url, headers=headers)
                        
                        if molecule_response.status_code == 200:
                            molecule_data = molecule_response.json()
                            
                            # Extract SMILES string
                            smiles = molecule_data.get("molecule_structures", {}).get("canonical_smiles", "")
                            
                            # Create drug object
                            drug = {
                                "name": molecule_data.get("pref_name", "Unknown"),
                                "molecule_chembl_id": molecule_chembl_id,
                                "smile": smiles,
                                "mechanism_of_action": mechanism.get("mechanism_of_action", ""),
                                "target_chembl_id": target_id,
                                "inchi_key": molecule_data.get("molecule_structures", {}).get("standard_inchi_key", ""),
                                "molecular_formula": molecule_data.get("molecule_properties", {}).get("full_molformula", ""),
                                "molecular_weight": molecule_data.get("molecule_properties", {}).get("full_mwt", ""),
                                "activity_comments": mechanism.get("action_type", ""),
                                "alogp": molecule_data.get("molecule_properties", {}).get("alogp", ""),
                                "aromatic_rings": molecule_data.get("molecule_properties", {}).get("aromatic_rings", ""),
                                "cx_logd": molecule_data.get("molecule_properties", {}).get("cx_logd", ""),
                                "cx_logp": molecule_data.get("molecule_properties", {}).get("cx_logp", ""),
                                "cx_most_apka": molecule_data.get("molecule_properties", {}).get("cx_most_apka", ""),
                                "cx_most_bpka": molecule_data.get("molecule_properties", {}).get("cx_most_bpka", ""),
                                "hba": molecule_data.get("molecule_properties", {}).get("hba", ""),
                                "hba_lipinski": molecule_data.get("molecule_properties", {}).get("hba_lipinski", ""),
                                "hbd": molecule_data.get("molecule_properties", {}).get("hbd", ""),
                                "hbd_lipinski": molecule_data.get("molecule_properties", {}).get("hbd_lipinski", ""),
                                "heavy_atoms": molecule_data.get("molecule_properties", {}).get("heavy_atoms", ""),
                                "molecular_species": molecule_data.get("molecule_properties", {}).get("molecular_species", ""),
                                "mw_freebase": molecule_data.get("molecule_properties", {}).get("mw_freebase", ""),
                                "mw_monoisotopic": molecule_data.get("molecule_properties", {}).get("mw_monoisotopic", ""),
                                "np_likeness_score": molecule_data.get("molecule_properties", {}).get("np_likeness_score", ""),
                                "num_lipinski_ro5_violations": molecule_data.get("molecule_properties", {}).get("num_lipinski_ro5_violations", ""),
                                "num_ro5_violations": molecule_data.get("molecule_properties", {}).get("num_ro5_violations", ""),
                                "psa": molecule_data.get("molecule_properties", {}).get("psa", ""),
                                "qed_weighted": molecule_data.get("molecule_properties", {}).get("qed_weighted", ""),
                                "ro3_pass": molecule_data.get("molecule_properties", {}).get("ro3_pass", ""),
                                "rtb": molecule_data.get("molecule_properties", {}).get("rtb", "")
                            }
                            drugs.append(drug)
                    
                    return drugs
        except Exception as e:
            logger.error(f"Error fetching drugs for target: {e}")
            
        return []
    
    async def get_proteins_by_disease_id(self, disease_id: str) -> List[Protein]:
        """Get proteins associated with a specific disease."""
        logger.info(f"Getting proteins for disease ID: {disease_id}")
        
        # Query database
        proteins_collection = db["proteins"]
        cursor = proteins_collection.find({"disease_id": disease_id})
        
        # Get disease name
        disease_name = None
        diseases_collection = db["diseases"]
        disease_doc = await diseases_collection.find_one({"_id": disease_id})
        if disease_doc:
            disease_name = disease_doc.get("name")
        
        # Convert to Protein objects
        proteins = []
        async for doc in cursor:
            protein = Protein(
                id=str(doc.get("_id")),
                name=doc.get("name"),
                pdb_id=doc.get("pdb_id"),
                sequence=doc.get("sequence"),
                disease_id=doc.get("disease_id"),
                description=doc.get("description"),
                binding_site=doc.get("binding_site"),
                created_at=doc.get("created_at", datetime.now()),
                updated_at=doc.get("updated_at", datetime.now()),
                disease_name=disease_name
            )
            proteins.append(protein)
        
        return proteins
    
    async def get_drugs_by_disease_proteins(self, disease_id: str, limit_per_protein: int = 10) -> List[Dict[str, Any]]:
        """
        Get drugs associated with proteins related to a disease.
        
        This method:
        1. Gets proteins associated with the disease
        2. For each protein, gets the ChEMBL target ID
        3. For each target, finds associated drugs/molecules
        4. Returns the drugs with their SMILES representations
        """
        logger.info(f"Getting drugs for proteins related to disease ID: {disease_id}")
        
        # Get proteins for the disease
        proteins = await self.get_proteins_by_disease_id(disease_id)
        
        if not proteins:
            logger.warning(f"No proteins found for disease ID: {disease_id}")
            return []
        
        # Get drugs for each protein
        all_drugs = []
        protein_drug_mapping = {}
        
        for protein in proteins:
            # Get UniProt ID (this would need to be added to your protein schema or fetched from external source)
            uniprot_id = protein.pdb_id  # This is a simplification; in reality, you'd need to map PDB ID to UniProt ID
            
            # Get ChEMBL target ID
            target_id = await self.get_chembl_target_id(uniprot_id)
            
            if not target_id:
                logger.warning(f"No ChEMBL target ID found for protein: {protein.name}")
                continue
            
            # Get drugs for the target
            drugs = await self.get_drugs_by_target_id(target_id, limit=limit_per_protein)
            
            if drugs:
                # Store mapping of protein to drugs
                protein_drug_mapping[protein.name] = drugs
                
                # Add drugs to overall list
                all_drugs.extend(drugs)
        
        # Remove duplicates based on molecule_chembl_id
        unique_drugs = {}
        for drug in all_drugs:
            if drug["molecule_chembl_id"] not in unique_drugs:
                unique_drugs[drug["molecule_chembl_id"]] = drug
        
        # Add protein association information
        for drug_id, drug in unique_drugs.items():
            associated_proteins = []
            for protein_name, protein_drugs in protein_drug_mapping.items():
                if any(d["molecule_chembl_id"] == drug_id for d in protein_drugs):
                    associated_proteins.append(protein_name)
            drug["associated_proteins"] = associated_proteins
        
        return list(unique_drugs.values())
    
    async def get_protein_info_by_disease(self, disease_name: str, skip: int = 0, limit: int = 100, disease_ids: List[str] = None) -> List[Protein]:
        """
        Get protein information for proteins associated with a disease.
        
        This method:
        1. Searches for disease using UniProt disease search API (if disease_ids not provided)
        2. Gets disease ID from search results (if disease_ids not provided)
        3. Gets proteins associated with disease ID(s)
        4. Returns formatted protein information
        
        Args:
            disease_name: Name of the disease to search for (used if disease_ids not provided)
            skip: Number of results to skip (for pagination)
            limit: Maximum number of results to return
            disease_ids: Optional list of disease IDs to use instead of searching
            
        Returns:
            List of Protein objects
        """
        logger.info(f"Getting protein info for disease: {disease_name}, disease_ids: {disease_ids}")
        
        try:
            all_proteins = []
            
            # If disease_ids provided, use them directly
            if disease_ids and len(disease_ids) > 0:
                for disease_id in disease_ids:
                    proteins = await self._get_proteins_for_disease_id(disease_id, disease_name)
                    all_proteins.extend(proteins)
            else:
                # First search for disease to get disease ID
                disease_search_url = "https://rest.uniprot.org/diseases/search"
                disease_params = {
                    "query": disease_name,
                    "format": "json"
                }
                
                async with httpx.AsyncClient() as client:
                    # Get disease search results
                    disease_response = await client.get(disease_search_url, params=disease_params)
                    if disease_response.status_code != 200:
                        logger.error(f"Error searching disease: {disease_response.status_code}")
                        return []
                        
                    disease_data = disease_response.json()
                    disease_results = disease_data.get("results", [])
                    
                    if not disease_results:
                        logger.warning(f"No disease found for name: {disease_name}")
                        return []
                    
                    # Get proteins for ALL matching diseases instead of just the first one
                    logger.info(f"Found {len(disease_results)} disease variants for '{disease_name}'")
                    for disease_result in disease_results:
                        disease_id = disease_result.get("id")
                        if not disease_id:
                            continue
                            
                        disease_variant_name = disease_result.get("name", disease_name)
                        logger.info(f"Getting proteins for disease variant: {disease_variant_name} (ID: {disease_id})")
                        
                        # Get proteins for this disease ID
                        proteins = await self._get_proteins_for_disease_id(disease_id, disease_variant_name)
                        logger.info(f"Found {len(proteins)} proteins for disease variant {disease_variant_name}")
                        all_proteins.extend(proteins)
            
            # Apply pagination
            start = skip
            end = skip + limit
            logger.info(f"Found total of {len(all_proteins)} proteins across all disease variants, returning {min(limit, len(all_proteins)-skip)} after pagination")
            return all_proteins[start:end]
            
        except Exception as e:
            logger.error(f"Error fetching protein info: {e}")
            return []
            
        return []

    async def _get_proteins_for_disease_id(self, disease_id: str, disease_name: str) -> List[Protein]:
        """Helper method to get proteins for a specific disease ID"""
        logger.info(f"Getting proteins for disease ID: {disease_id}")
        
        try:
            async with httpx.AsyncClient() as client:
                # Now get proteins associated with this disease ID
                protein_search_url = "https://rest.uniprot.org/uniprotkb/search"
                protein_params = {
                    "query": f"cc_disease:{disease_id}",
                    "format": "json"
                }
                
                protein_response = await client.get(protein_search_url, params=protein_params)
                if protein_response.status_code != 200:
                    logger.error(f"Error getting proteins: {protein_response.status_code}")
                    return []
                
                protein_data = protein_response.json()
                protein_results = protein_data.get("results", [])
                
                proteins = []
                for result in protein_results:
                    # Get protein name from recommended name
                    protein_desc = result.get("proteinDescription", {})
                    recommended_name = protein_desc.get("recommendedName", {})
                    name = recommended_name.get("fullName", {}).get("value", "Unknown")
                    
                    # Get UniProt ID from primaryAccession
                    uniprot_id = result.get("primaryAccession", "")
                    
                    # Get PDB IDs from database references
                    pdb_refs = [ref.get("id") for ref in result.get("dbReferences", []) 
                              if ref.get("type") == "PDB"]
                    pdb_id = pdb_refs[0] if pdb_refs else None
                    
                    # Get disease description from comments
                    disease_comments = [c for c in result.get("comments", []) 
                                     if c.get("commentType") == "DISEASE"]
                    description = ""
                    if disease_comments:
                        texts = disease_comments[0].get("texts", [])
                        if texts:
                            description = texts[0].get("value", "").split(".")[0]
                    
                    # Create protein object
                    now = datetime.now()
                    protein = Protein(
                        id=str(uuid.uuid4()),
                        name=name,
                        pdb_id=pdb_id,
                        uniprot_id=uniprot_id,
                        sequence=result.get("sequence", {}).get("value", ""),
                        disease_id=disease_id,
                        description=description,
                        binding_site={},  # Not available from UniProt API
                        created_at=now,
                        updated_at=now,
                        disease_name=disease_name
                    )
                    proteins.append(protein)
                
                return proteins
        except Exception as e:
            logger.error(f"Error getting proteins for disease ID {disease_id}: {e}")
            return []

    async def get_all_drug_smiles_for_disease(self, disease_name: str, skip: int = 0, limit: int = 100) -> Tuple[List[Dict[str, Any]], List[str]]:
        """
        Get all drugs and their SMILES representations for a disease.
        
        This method:
        1. Gets proteins associated with the disease
        2. For each protein, gets the ChEMBL target ID
        3. For each target, finds associated drugs/molecules
        4. Returns the drugs with their SMILES representations
        
        Returns:
            Tuple containing:
            - List of drug dictionaries with full information
            - List of SMILES strings only (for easy processing)
        """
        logger.info(f"Getting all drug SMILES for disease: {disease_name}")
        
        # Get proteins for the disease
        proteins = await self.get_protein_info_by_disease(disease_name, skip, limit)
        
        if not proteins:
            logger.warning(f"No proteins found for disease: {disease_name}")
            return [], []
        
        # Get drugs for each protein
        all_drugs = []
        
        for protein in proteins:
            # Get ChEMBL target ID using the UniProt ID
            target_id = await self.get_chembl_target_id(protein.uniprot_id)
            
            if not target_id:
                continue
            
            # Get drugs for the target
            drugs = await self.get_drugs_by_target_id(target_id)
            logger.info(f"Found {len(drugs)} drugs for protein {protein.name}")
            all_drugs.extend(drugs)
            
        # Remove duplicates
        unique_drugs = {}
        for drug in all_drugs:
            if drug["molecule_chembl_id"] not in unique_drugs:
                unique_drugs[drug["molecule_chembl_id"]] = drug
        
        # Extract SMILES strings
        smiles_list = [drug["smile"] for drug in unique_drugs.values() if drug["smile"]]
        
        return list(unique_drugs.values()), smiles_list

    async def get_drugs_for_protein(self, protein: Protein) -> Dict[str, Any]:
        """Get drugs for a single protein."""
        logger.info(f"Getting drugs for protein: {protein.name}")
        
        try:
            # Get ChEMBL target ID using the UniProt ID
            target_id = await self.get_chembl_target_id(protein.uniprot_id)
            if not target_id:
                logger.warning(f"No ChEMBL target found for protein: {protein.name} (UniProt ID: {protein.uniprot_id})")
                return {}
            
            # Get drugs for the target
            drugs = await self.get_drugs_by_target_id(target_id)
            
            if not drugs:
                logger.warning(f"No drugs found for protein: {protein.name}")
                return {}
            
            # Create drug dictionary
            drug_dict = {
                "name": protein.name,
                "molecule_chembl_id": drugs[0]["molecule_chembl_id"],
                "smile": drugs[0]["smile"],
                "mechanism_of_action": drugs[0]["mechanism_of_action"],
                "target_chembl_id": drugs[0]["target_chembl_id"],
                "inchi_key": drugs[0]["inchi_key"],
                "molecular_formula": drugs[0]["molecular_formula"],
                "molecular_weight": drugs[0]["molecular_weight"],
                "activity_comments": drugs[0]["activity_comments"],
                "alogp": drugs[0]["alogp"],
                "aromatic_rings": drugs[0]["aromatic_rings"],
                "cx_logd": drugs[0]["cx_logd"],
                "cx_logp": drugs[0]["cx_logp"],
                "cx_most_apka": drugs[0]["cx_most_apka"],
                "cx_most_bpka": drugs[0]["cx_most_bpka"],
                "hba": drugs[0]["hba"],
                "hba_lipinski": drugs[0]["hba_lipinski"],
                "hbd": drugs[0]["hbd"],
                "hbd_lipinski": drugs[0]["hbd_lipinski"],
                "heavy_atoms": drugs[0]["heavy_atoms"],
                "molecular_species": drugs[0]["molecular_species"],
                "mw_freebase": drugs[0]["mw_freebase"],
                "mw_monoisotopic": drugs[0]["mw_monoisotopic"],
                "np_likeness_score": drugs[0]["np_likeness_score"],
                "num_lipinski_ro5_violations": drugs[0]["num_lipinski_ro5_violations"],
                "num_ro5_violations": drugs[0]["num_ro5_violations"],
                "psa": drugs[0]["psa"],
                "qed_weighted": drugs[0]["qed_weighted"],
                "ro3_pass": drugs[0]["ro3_pass"],
                "rtb": drugs[0]["rtb"],
                "associated_proteins": drugs[0]["associated_proteins"],
                "timestamp": datetime.now().isoformat()
            }
            
            return drug_dict
        except Exception as e:
            logger.error(f"Error fetching drugs for protein: {e}")
            return {}

    async def get_drugs_for_proteins(self, proteins: List[Protein]) -> List[Dict[str, Any]]:
        """Get drugs for multiple proteins concurrently."""
        logger.info(f"Getting drugs for {len(proteins)} proteins")
        
        tasks = [self.get_drugs_for_protein(protein) for protein in proteins]
        results = await asyncio.gather(*tasks)
        
        return results

    async def clear_protein_drug_cache(self):
        """Clear the protein drug cache."""
        logger.info("Clearing protein drug cache")
        # Implementation of clearing the cache
        pass 

    async def get_diseases_by_name(self, disease_name: str) -> List[Dict[str, Any]]:
        """
        Get all diseases matching the search query.
        
        This method:
        1. Searches for diseases using UniProt disease search API
        2. Returns all matching diseases with their details
        """
        logger.info(f"Getting diseases matching query: {disease_name}")
        
        try:
            # Search for diseases
            disease_search_url = "https://rest.uniprot.org/diseases/search"
            disease_params = {
                "query": disease_name,
                "format": "json"
            }
            
            async with httpx.AsyncClient() as client:
                # Get disease search results
                disease_response = await client.get(disease_search_url, params=disease_params)
                if disease_response.status_code != 200:
                    logger.error(f"Error searching diseases: {disease_response.status_code}")
                    return []
                    
                disease_data = disease_response.json()
                disease_results = disease_data.get("results", [])
                
                if not disease_results:
                    logger.warning(f"No diseases found for query: {disease_name}")
                    return []
                
                # Return all matching diseases with their details
                diseases = []
                for disease in disease_results:
                    disease_info = {
                        "id": disease.get("id"),
                        "name": disease.get("name"),
                        "acronym": disease.get("acronym", ""),
                        "definition": disease.get("definition", ""),
                        "alternativeNames": disease.get("alternativeNames", []),
                        "reviewedProteinCount": disease.get("statistics", {}).get("reviewedProteinCount", 0),
                        "keywords": [kw.get("name") for kw in disease.get("keywords", [])]
                    }
                    diseases.append(disease_info)
                
                return diseases
                
        except Exception as e:
            logger.error(f"Error fetching diseases: {e}")
            return [] 

    async def get_protein_external_links(self, uniprot_id: str) -> Dict[str, str]:
        """
        Fetch external links/cross-references for a protein using its UniProt ID
        
        Args:
            uniprot_id: The UniProt ID of the protein
            
        Returns:
            A dictionary mapping database names to their IDs
        """
        if not uniprot_id:
            logger.warning("Cannot fetch external links: No UniProt ID provided")
            return {}
            
        logger.info(f"Getting external links for protein with UniProt ID: {uniprot_id}")
        
        # Initialize the external_links dictionary
        external_links = {'UniProt': uniprot_id}
        
        # First try UniProt API
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        params = {"format": "json"}
        
        try:
            logger.info(f"Making request to UniProt API: {url}")
            async with httpx.AsyncClient() as client:
                response = await client.get(url, params=params)
                print(f"UniProt API response: {response.json()}")
                logger.info(f"UniProt API response status: {response.status_code}")
                
                if response.status_code != 200:
                    logger.error(f"Error fetching external links for protein {uniprot_id}: Status {response.status_code}")
                    try:
                        error_text = await response.text()
                        logger.error(f"Error response: {error_text[:200]}")
                    except Exception as e:
                        logger.error(f"Could not get error text: {e}")
                else:
                    data = response.json()
                    logger.info(f"Successfully retrieved UniProt data for {uniprot_id}")
                    
                    # Extract all cross-references
                    cross_refs = data.get('uniProtKBCrossReferences', [])
                    logger.info(f"Found {len(cross_refs)} cross-references")
                    
                    for ref in cross_refs:
                        db = ref.get('database')
                        db_id = ref.get('id')
                        if db and db_id:
                            external_links[db] = db_id
                            logger.debug(f"Found cross-reference: {db} = {db_id}")
                
                # Check if we found a ChEMBL ID
                if 'ChEMBL' in external_links:
                    logger.info(f"Found ChEMBL ID: {external_links['ChEMBL']} from UniProt for UniProt ID: {uniprot_id}")
                else:
                    logger.warning(f"No ChEMBL ID found in UniProt cross-references for {uniprot_id}")
        
        except Exception as e:
            logger.error(f"Error fetching external links from UniProt API for protein {uniprot_id}: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
        
        # If no ChEMBL ID found, try ChEMBL API directly
        if 'ChEMBL' not in external_links:
            logger.info(f"Trying ChEMBL API directly for UniProt ID: {uniprot_id}")
            try:
                chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/target.json?uniprot={uniprot_id}"
                headers = {"Accept": "application/json"}
                
                async with httpx.AsyncClient() as client:
                    chembl_response = await client.get(chembl_url, headers=headers)
                    
                    if chembl_response.status_code == 200:
                        chembl_data = chembl_response.json()
                        targets = chembl_data.get("targets", [])
                        
                        if targets:
                            chembl_id = targets[0]["target_chembl_id"]
                            external_links['ChEMBL'] = chembl_id
                            logger.info(f"Found ChEMBL ID: {chembl_id} directly from ChEMBL API for UniProt ID: {uniprot_id}")
                        else:
                            logger.warning(f"No targets found in ChEMBL API for UniProt ID: {uniprot_id}")
                    else:
                        logger.error(f"ChEMBL API error: {chembl_response.status_code}")
            except Exception as e:
                logger.error(f"Error fetching ChEMBL ID from ChEMBL API: {e}")
                
        # If still no ChEMBL ID found, check hardcoded mappings for common proteins
        if 'ChEMBL' not in external_links:
            # Common known mappings for important proteins
            # These are manually curated mappings for proteins that might not have
            # proper cross-references in UniProt or ChEMBL
            hardcoded_mappings = {
                'P05067': 'CHEMBL2487',  # Amyloid beta precursor protein (APP)
                'P37840': 'CHEMBL4296002',  # Alpha-synuclein
                'P10636': 'CHEMBL2364663',  # Tau protein
                'P49841': 'CHEMBL2095194',  # GSK3B
                'P06213': 'CHEMBL1824',  # Insulin receptor
                'P14672': 'CHEMBL2074',  # GLUT4
                'P43220': 'CHEMBL1784',  # GLP-1 receptor
                'P56817': 'CHEMBL4822',  # BACE1
                'O60260': 'CHEMBL4303387',  # Parkin
                'Q99497': 'CHEMBL2334209',  # DJ-1
                'P01375': 'CHEMBL4679',  # TNF-alpha
                'P01584': 'CHEMBL4679',  # IL-1B
            }
            
            if uniprot_id in hardcoded_mappings:
                chembl_id = hardcoded_mappings[uniprot_id]
                external_links['ChEMBL'] = chembl_id
                logger.info(f"Using hardcoded ChEMBL ID: {chembl_id} for UniProt ID: {uniprot_id}")
        
        # Log a list of all found databases for debugging
        if external_links:
            logger.info(f"External links found for {uniprot_id}: {', '.join(external_links.keys())}")
        else:
            logger.warning(f"No external links found for {uniprot_id}")
            
        return external_links 