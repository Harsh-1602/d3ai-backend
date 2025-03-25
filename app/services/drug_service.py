import logging
import uuid
import httpx
from datetime import datetime
from typing import List, Optional, Dict, Any
from app.db.mongodb import db
from app.schemas.drug import (
    DrugCreate, Drug, DrugSearchResult, DrugGenerationResult
)
from app.core.config import settings

logger = logging.getLogger(__name__)

class DrugService:
    """Service for drug-related operations."""
    
    async def search_drugs(
        self, 
        disease: Optional[str] = None, 
        keyword: Optional[str] = None, 
        source: str = "all", 
        limit: int = 20
    ) -> DrugSearchResult:
        """
        Search for drug candidates.
        
        This method searches for drug candidates based on disease name or keywords.
        It can search across multiple databases including ChEMBL, PubChem, ZINC, and research papers.
        """
        logger.info(f"Searching drugs with disease={disease}, keyword={keyword}, source={source}")
        
        # Determine the search query
        query = disease or keyword or ""
        
        # Initialize results
        all_drugs = []
        
        # Search in database first
        if source in ["all", "database"]:
            db_drugs = await self._search_drugs_in_database(disease, keyword, limit)
            all_drugs.extend(db_drugs)
        
        # Search in external sources if needed
        remaining_limit = limit - len(all_drugs)
        if remaining_limit > 0:
            if source in ["all", "chembl"] and disease:
                chembl_drugs = await self._search_drugs_in_chembl(disease, remaining_limit)
                all_drugs.extend(chembl_drugs)
            
            remaining_limit = limit - len(all_drugs)
            if remaining_limit > 0 and source in ["all", "pubchem"] and keyword:
                pubchem_drugs = await self._search_drugs_in_pubchem(keyword, remaining_limit)
                all_drugs.extend(pubchem_drugs)
            
            remaining_limit = limit - len(all_drugs)
            if remaining_limit > 0 and source in ["all", "papers"] and disease:
                paper_drugs = await self._search_drugs_in_papers(disease, remaining_limit)
                all_drugs.extend(paper_drugs)
        
        # Ensure we don't exceed the limit
        all_drugs = all_drugs[:limit]
        
        # Return the search result
        return DrugSearchResult(
            query=query,
            total_results=len(all_drugs),
            drugs=all_drugs,
            search_source=source
        )
    
    async def _search_drugs_in_database(
        self, 
        disease: Optional[str] = None, 
        keyword: Optional[str] = None, 
        limit: int = 20
    ) -> List[Drug]:
        """Search for drugs in the database."""
        # Build query
        query = {}
        
        if disease:
            # First, find the disease ID
            diseases_collection = db["diseases"]
            disease_doc = await diseases_collection.find_one({"name": {"$regex": disease, "$options": "i"}})
            
            if disease_doc:
                query["disease_id"] = str(disease_doc.get("_id"))
        
        if keyword:
            query["$or"] = [
                {"name": {"$regex": keyword, "$options": "i"}},
                {"description": {"$regex": keyword, "$options": "i"}},
                {"mechanism_of_action": {"$regex": keyword, "$options": "i"}}
            ]
        
        # Query database
        drugs_collection = db["drugs"]
        cursor = drugs_collection.find(query).limit(limit)
        
        # Convert to Drug objects
        drugs = []
        async for doc in cursor:
            # Get disease name
            disease_name = None
            if doc.get("disease_id"):
                diseases_collection = db["diseases"]
                disease_doc = await diseases_collection.find_one({"_id": doc.get("disease_id")})
                if disease_doc:
                    disease_name = disease_doc.get("name")
            
            # Create Drug object
            drug = Drug(
                id=str(doc.get("_id")),
                name=doc.get("name"),
                smile=doc.get("smile"),
                disease_id=doc.get("disease_id"),
                description=doc.get("description"),
                brand_names=doc.get("brand_names", []),
                drug_class=doc.get("drug_class"),
                mechanism_of_action=doc.get("mechanism_of_action"),
                research_paper_url=doc.get("research_paper_url"),
                created_at=doc.get("created_at", datetime.now()),
                updated_at=doc.get("updated_at", datetime.now()),
                disease_name=disease_name,
                properties={},  # TODO: Fetch actual properties
                similar_drugs=[]  # TODO: Fetch similar drugs
            )
            drugs.append(drug)
        
        return drugs
    
    async def _search_drugs_in_chembl(self, disease: str, limit: int = 10) -> List[Drug]:
        """Search for drugs in ChEMBL database."""
        logger.info(f"Searching drugs in ChEMBL for disease: {disease}")
        
        try:
            # Search for disease in ChEMBL
            async with httpx.AsyncClient() as client:
                # First, search for the disease
                disease_url = f"{settings.CHEMBL_API_BASE_URL}/disease.json?search={disease}&limit=1"
                disease_response = await client.get(disease_url)
                disease_data = disease_response.json()
                
                if not disease_data.get("diseases") or len(disease_data["diseases"]) == 0:
                    return []
                
                disease_chembl_id = disease_data["diseases"][0]["disease_chembl_id"]
                
                # Then, get drugs for this disease
                drugs_url = f"{settings.CHEMBL_API_BASE_URL}/mechanism.json?target_chembl_id={disease_chembl_id}&limit={limit}"
                drugs_response = await client.get(drugs_url)
                drugs_data = drugs_response.json()
                
                if not drugs_data.get("mechanisms") or len(drugs_data["mechanisms"]) == 0:
                    return []
                
                # Convert to Drug objects
                drugs = []
                for mechanism in drugs_data["mechanisms"]:
                    # Get molecule details
                    molecule_chembl_id = mechanism.get("molecule_chembl_id")
                    if not molecule_chembl_id:
                        continue
                    
                    molecule_url = f"{settings.CHEMBL_API_BASE_URL}/molecule/{molecule_chembl_id}.json"
                    molecule_response = await client.get(molecule_url)
                    molecule_data = molecule_response.json()
                    
                    if not molecule_data:
                        continue
                    
                    # Create a new drug ID
                    drug_id = str(uuid.uuid4())
                    
                    # Create Drug object
                    drug = Drug(
                        id=drug_id,
                        name=molecule_data.get("pref_name", "Unknown"),
                        smile=molecule_data.get("molecule_structures", {}).get("canonical_smiles", ""),
                        disease_id="",  # We don't have a disease ID in our database yet
                        description=mechanism.get("mechanism_of_action", ""),
                        brand_names=molecule_data.get("molecule_synonyms", []),
                        drug_class=molecule_data.get("molecule_type", ""),
                        mechanism_of_action=mechanism.get("mechanism_of_action", ""),
                        research_paper_url="",
                        created_at=datetime.now(),
                        updated_at=datetime.now(),
                        disease_name=disease,
                        properties={},
                        similar_drugs=[]
                    )
                    drugs.append(drug)
                
                return drugs
                
        except Exception as e:
            logger.error(f"Error searching drugs in ChEMBL: {e}")
            return []
    
    async def _search_drugs_in_pubchem(self, keyword: str, limit: int = 10) -> List[Drug]:
        """Search for drugs in PubChem database."""
        logger.info(f"Searching drugs in PubChem for keyword: {keyword}")
        
        try:
            # Search for compounds in PubChem
            async with httpx.AsyncClient() as client:
                search_url = f"{settings.PUBCHEM_API_BASE_URL}/compound/name/{keyword}/JSON?MaxRecords={limit}"
                response = await client.get(search_url)
                data = response.json()
                
                if not data.get("PC_Compounds"):
                    return []
                
                # Convert to Drug objects
                drugs = []
                for compound in data["PC_Compounds"]:
                    # Extract compound information
                    cid = compound.get("id", {}).get("id", {}).get("cid")
                    if not cid:
                        continue
                    
                    # Get detailed information
                    detail_url = f"{settings.PUBCHEM_API_BASE_URL}/compound/cid/{cid}/JSON"
                    detail_response = await client.get(detail_url)
                    detail_data = detail_response.json()
                    
                    if not detail_data.get("PC_Compounds"):
                        continue
                    
                    detail = detail_data["PC_Compounds"][0]
                    
                    # Extract SMILES
                    smile = ""
                    for prop in detail.get("props", []):
                        if prop.get("urn", {}).get("label") == "SMILES":
                            smile = prop.get("value", {}).get("sval", "")
                            break
                    
                    # Extract name
                    name = keyword
                    for prop in detail.get("props", []):
                        if prop.get("urn", {}).get("label") == "IUPAC Name":
                            name = prop.get("value", {}).get("sval", keyword)
                            break
                    
                    # Create a new drug ID
                    drug_id = str(uuid.uuid4())
                    
                    # Create Drug object
                    drug = Drug(
                        id=drug_id,
                        name=name,
                        smile=smile,
                        disease_id="",  # We don't have a disease ID in our database yet
                        description="",
                        brand_names=[],
                        drug_class="",
                        mechanism_of_action="",
                        research_paper_url="",
                        created_at=datetime.now(),
                        updated_at=datetime.now(),
                        disease_name="",
                        properties={},
                        similar_drugs=[]
                    )
                    drugs.append(drug)
                
                return drugs
                
        except Exception as e:
            logger.error(f"Error searching drugs in PubChem: {e}")
            return []
    
    async def _search_drugs_in_papers(self, disease: str, limit: int = 10) -> List[Drug]:
        """Search for drugs in research papers."""
        logger.info(f"Searching drugs in research papers for disease: {disease}")
        
        # TODO: Implement actual search in vector database
        # For now, return an empty list
        return []
    
    async def get_drugs_by_disease(
        self, 
        disease_id: str, 
        skip: int = 0, 
        limit: int = 100
    ) -> List[Drug]:
        """Get drug candidates for a specific disease."""
        logger.info(f"Getting drugs for disease ID: {disease_id}")
        
        # Query database
        drugs_collection = db["drugs"]
        cursor = drugs_collection.find({"disease_id": disease_id}).skip(skip).limit(limit)
        
        # Get disease name
        disease_name = None
        diseases_collection = db["diseases"]
        disease_doc = await diseases_collection.find_one({"_id": disease_id})
        if disease_doc:
            disease_name = disease_doc.get("name")
        
        # Convert to Drug objects
        drugs = []
        async for doc in cursor:
            # Create Drug object
            drug = Drug(
                id=str(doc.get("_id")),
                name=doc.get("name"),
                smile=doc.get("smile"),
                disease_id=doc.get("disease_id"),
                description=doc.get("description"),
                brand_names=doc.get("brand_names", []),
                drug_class=doc.get("drug_class"),
                mechanism_of_action=doc.get("mechanism_of_action"),
                research_paper_url=doc.get("research_paper_url"),
                created_at=doc.get("created_at", datetime.now()),
                updated_at=doc.get("updated_at", datetime.now()),
                disease_name=disease_name,
                properties={},  # TODO: Fetch actual properties
                similar_drugs=[]  # TODO: Fetch similar drugs
            )
            drugs.append(drug)
        
        return drugs
    
    async def generate_drug_candidates(
        self,
        disease_id: str,
        seed_drug_ids: List[str],
        count: int = 10,
        diversity_factor: float = 0.7
    ) -> DrugGenerationResult:
        """Generate new drug candidates."""
        logger.info(f"Generating drug candidates for disease ID: {disease_id}")
        
        # Get disease name
        disease_name = "Unknown Disease"
        diseases_collection = db["diseases"]
        disease_doc = await diseases_collection.find_one({"_id": disease_id})
        if disease_doc:
            disease_name = disease_doc.get("name")
        
        # Get seed drugs
        seed_drugs = []
        seed_smiles = []
        
        if seed_drug_ids:
            drugs_collection = db["drugs"]
            for drug_id in seed_drug_ids:
                drug_doc = await drugs_collection.find_one({"_id": drug_id})
                if drug_doc:
                    seed_drugs.append(drug_doc.get("name"))
                    seed_smiles.append(drug_doc.get("smile"))
        else:
            # If no seed drugs provided, get some from the database
            drugs_collection = db["drugs"]
            cursor = drugs_collection.find({"disease_id": disease_id}).limit(3)
            async for doc in cursor:
                seed_drugs.append(doc.get("name"))
                seed_smiles.append(doc.get("smile"))
        
        # TODO: Implement actual molecule generation using MegaMolBERT or other models
        # For now, we'll return mock data
        
        # Create mock generated molecules
        generated_molecules = []
        for i in range(count):
            molecule = {
                "id": str(uuid.uuid4()),
                "name": f"Generated Molecule {i+1}",
                "smile": seed_smiles[0] if seed_smiles else "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",
                "similarity": 0.8 - (i * 0.05),  # Decreasing similarity
                "properties": {
                    "bioavailability": 0.7 + (i * 0.02),
                    "toxicity": 0.3 - (i * 0.01),
                    "solubility": 0.6 + (i * 0.03)
                }
            }
            generated_molecules.append(molecule)
        
        # Return the result
        return DrugGenerationResult(
            disease_id=disease_id,
            disease_name=disease_name,
            seed_drugs=seed_drugs,
            generated_molecules=generated_molecules
        )
    
    async def get_drug_by_id(self, drug_id: str) -> Optional[Drug]:
        """Get a specific drug by ID."""
        logger.info(f"Getting drug by ID: {drug_id}")
        
        # Query database
        drugs_collection = db["drugs"]
        doc = await drugs_collection.find_one({"_id": drug_id})
        
        if not doc:
            return None
        
        # Get disease name
        disease_name = None
        if doc.get("disease_id"):
            diseases_collection = db["diseases"]
            disease_doc = await diseases_collection.find_one({"_id": doc.get("disease_id")})
            if disease_doc:
                disease_name = disease_doc.get("name")
        
        # Create Drug object
        drug = Drug(
            id=str(doc.get("_id")),
            name=doc.get("name"),
            smile=doc.get("smile"),
            disease_id=doc.get("disease_id"),
            description=doc.get("description"),
            brand_names=doc.get("brand_names", []),
            drug_class=doc.get("drug_class"),
            mechanism_of_action=doc.get("mechanism_of_action"),
            research_paper_url=doc.get("research_paper_url"),
            created_at=doc.get("created_at", datetime.now()),
            updated_at=doc.get("updated_at", datetime.now()),
            disease_name=disease_name,
            properties={},  # TODO: Fetch actual properties
            similar_drugs=[]  # TODO: Fetch similar drugs
        )
        
        return drug
    
    async def create_drug(self, drug: DrugCreate) -> Drug:
        """Create a new drug."""
        logger.info(f"Creating new drug: {drug.name}")
        
        # Create document
        drug_id = str(uuid.uuid4())
        now = datetime.now()
        
        drug_doc = {
            "_id": drug_id,
            "name": drug.name,
            "smile": drug.smile,
            "disease_id": drug.disease_id,
            "description": drug.description,
            "brand_names": drug.brand_names,
            "drug_class": drug.drug_class,
            "mechanism_of_action": drug.mechanism_of_action,
            "research_paper_url": drug.research_paper_url,
            "created_at": now,
            "updated_at": now
        }
        
        # Insert into database
        drugs_collection = db["drugs"]
        await drugs_collection.insert_one(drug_doc)
        
        # Get disease name
        disease_name = None
        if drug.disease_id:
            diseases_collection = db["diseases"]
            disease_doc = await diseases_collection.find_one({"_id": drug.disease_id})
            if disease_doc:
                disease_name = disease_doc.get("name")
        
        # Return created drug
        return Drug(
            id=drug_id,
            name=drug.name,
            smile=drug.smile,
            disease_id=drug.disease_id,
            description=drug.description,
            brand_names=drug.brand_names,
            drug_class=drug.drug_class,
            mechanism_of_action=drug.mechanism_of_action,
            research_paper_url=drug.research_paper_url,
            created_at=now,
            updated_at=now,
            disease_name=disease_name,
            properties={},
            similar_drugs=[]
        ) 