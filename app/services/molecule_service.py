import logging
from datetime import datetime
import uuid
from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import base64
import io
from bson.objectid import ObjectId
from pymongo import MongoClient
from pymongo.collection import Collection
import requests
import json
import random
import urllib.parse
import aiohttp

from app.core.config import settings
from app.schemas.molecule import Molecule, MoleculeCreate, MoleculeValidationResult, BulkMoleculeValidationResult, NvidiaGenMolRequest
from app.utils.db import get_database

logger = logging.getLogger(__name__)

class MoleculeService:
    def __init__(self):
        """Initialize MoleculeService with database connection."""
        self.db = get_database()
        self.collection: Collection = self.db["molecules"]

    async def generate_molecules(
        self, 
        seed_smiles: str, 
        count: int = 10, 
        diversity_factor: float = 0.5,
        similarity_threshold: float = 0.4
    ) -> List[Molecule]:
        """
        Generate new molecules from seed SMILES.
        
        Args:
            seed_smiles: SMILES string to use as a seed
            count: Number of molecules to generate
            diversity_factor: Factor to control molecular diversity (0-1)
            similarity_threshold: Minimum similarity to seed molecule (0-1)
            
        Returns:
            List of generated molecules
        """
        logger.info(f"Generating {count} molecules from seed SMILES: {seed_smiles}")
        
        # For now, return dummy data since we don't have an actual generative model
        # In a real implementation, this would use an AI model to generate molecules
        
        molecules = []
        # Create some simple variations of the seed molecule for demonstration
        mol = Chem.MolFromSmiles(seed_smiles)
        
        if mol is None:
            logger.error(f"Invalid seed SMILES: {seed_smiles}")
            return []
        
        # Add the seed molecule first
        molecules.append({
            "_id": str(uuid.uuid4()),
            "smile": seed_smiles,
            "name": f"Generated from {seed_smiles[:10]}...",
            "description": "Seed molecule",
            "parent_drug_id": None,
            "generation_method": "manual",
            "is_valid": True,
            "created_at": datetime.now(),
            "updated_at": datetime.now()
        })
        
        # Create simple variations for demonstration
        # In a real implementation, this would use an AI model
        for i in range(1, min(count, 5)):
            # Create a simple variation by adding methyl groups
            # This is just for demonstration - not a real molecular generator
            rdkit_mol = Chem.MolFromSmiles(seed_smiles)
            if rdkit_mol:
                # Try different modifications based on the index
                if i % 3 == 0:
                    # Add OH group to first atom if possible
                    rdkit_mol = Chem.AddHs(rdkit_mol)
                    AllChem.Compute2DCoords(rdkit_mol)
                    modified_smile = Chem.MolToSmiles(rdkit_mol)
                elif i % 3 == 1:
                    # Try to replace a hydrogen with methyl if possible
                    modified_smile = seed_smiles + "C"  # Naive modification, not chemically valid
                else:
                    # Just add a carbon
                    modified_smile = seed_smiles + ".C"
            else:
                modified_smile = seed_smiles
                
            molecules.append({
                "_id": str(uuid.uuid4()),
                "smile": modified_smile,
                "name": f"Variant {i} of {seed_smiles[:10]}...",
                "description": f"Generated molecule variant {i}",
                "parent_drug_id": None,
                "generation_method": "ai_assisted",
                "is_valid": True,
                "created_at": datetime.now(),
                "updated_at": datetime.now()
            })
        
        return molecules

    async def generate_molecules_nvidia_genmol(self, request: NvidiaGenMolRequest) -> List[Molecule]:
        """
        Generate new molecules using NVIDIA GenMol API.
        
        Args:
            request: NvidiaGenMolRequest with parameters for generation
            
        Returns:
            List of generated molecules
        """
        logger.info(f"Generating molecules using NVIDIA GenMol API from SMILES: {request.smiles}")
        
        # Create the seed molecule for similarity comparison
        seed_mol = Chem.MolFromSmiles(request.smiles)
        if seed_mol:
            # Calculate fingerprint for the seed molecule for similarity comparison
            seed_fp = AllChem.GetMorganFingerprintAsBitVect(seed_mol, 2, nBits=2048)
        
        # NVIDIA GenMol API endpoint
        # Note: API endpoints may change based on API version
        api_urls = [
            "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate",
             # Another alternative
        ]
        
        # Prepare the payload
        # NVIDIA may expect different formats depending on the API version
        payload = {
            "smiles": request.smiles,
            "num_molecules": request.num_molecules,
            "temperature": request.temperature,
            "noise": request.noise,
            "step_size": request.step_size,
            "scoring": request.scoring,
            "unique": request.unique
        }
        
        # Also try alternative format
       
        
        # API headers
        headers = {
            "accept": "application/json",
            "content-type": "application/json",
            "authorization": f"Bearer {settings.NVIDIA_API_KEY}"
        }
        
        last_error = None
        for url in api_urls:
            for current_payload in [payload]:
                try:
                    # Make the API call
                    add_string = ".[*{20-20}]"
                    current_payload["smiles"] = f"{current_payload['smiles']}{add_string}"
                    logger.info(f"Trying NVIDIA GenMol API with URL: {url}")
                    response = requests.post(url, json=current_payload, headers=headers, timeout=30)
            
                    response.raise_for_status()  # Raise exception for non-200 status codes
                    
                    # Parse the response
                    result = response.json()
                    
                    # Create molecule objects from the API response
                    molecules = []
                    
                    # Handle response based on format - different API versions may have different response structures
                    if "status" in result and "molecules" in result:
                        # Standard format
                        for idx, mol_data in enumerate(result["molecules"]):
                            if "smiles" in mol_data:
                                mol_smiles = mol_data["smiles"]
                                # Generate a name based on the index
                                name = f"GenMol Generated {idx+1}"
                                
                                # Extract properties if available
                                properties = {}
                                if "score" in mol_data:
                                    properties[request.scoring] = mol_data["score"]
                                
                                # Generate a unique ID
                                molecule_id = str(uuid.uuid4())
                                
                                # Generate structure image URL from PubChem
                                structure_img = await self.get_pubchem_image(mol_smiles)
                                
                                # Add image to properties
                                if structure_img:
                                    if not properties:
                                        properties = {}
                                    properties["image"] = structure_img
                                
                                # Calculate properties using RDKit
                                mol = Chem.MolFromSmiles(mol_smiles)
                                if mol:
                                    try:
                                        # Add computed properties using RDKit
                                        from rdkit.Chem import Descriptors, Lipinski, DataStructs
                                        properties.update({
                                            "molecular_weight": round(Descriptors.MolWt(mol), 2),
                                            "logP": round(Descriptors.MolLogP(mol), 2),
                                            "num_h_donors": Lipinski.NumHDonors(mol),
                                            "num_h_acceptors": Lipinski.NumHAcceptors(mol),
                                            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                                            "tpsa": round(Descriptors.TPSA(mol), 2)
                                        })
                                        
                                        # Calculate similarity to seed molecule
                                        if seed_mol:
                                            mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                                            similarity = DataStructs.TanimotoSimilarity(seed_fp, mol_fp)
                                            properties["similarity_to_seed"] = round(similarity, 3)
                                            
                                    except Exception as e:
                                        logger.error(f"Error calculating molecular properties: {e}")
                                
                                # Create a molecule object
                                molecule = {
                                    "_id": molecule_id,
                                    "id": molecule_id,
                                    "smile": mol_smiles,
                                    "name": name,
                                    "description": f"Generated using NVIDIA GenMol API with {request.scoring} scoring",
                                    "parent_drug_id": None,
                                    "generation_method": "nvidia_genmol",
                                    "is_valid": True,
                                    "created_at": datetime.now(),
                                    "updated_at": datetime.now(),
                                    "properties": properties,
                                    "structure_img": structure_img
                                }
                                
                                molecules.append(molecule)
                    print(f"Molecules: {molecules}")
                    if molecules:
                        logger.info(f"Successfully generated {len(molecules)} molecules using NVIDIA GenMol API at {url}")
                        return molecules
                
                except Exception as e:
                    logger.warning(f"Error generating molecules with NVIDIA GenMol API at {url}: {e}")
                    last_error = e
                    continue  # Try next URL or payload format
        
        # If we reach here, all attempts failed
        logger.error(f"All attempts to generate molecules with NVIDIA GenMol API failed. Last error: {last_error}")
        
        # Generate dummy data as fallback
        logger.info("Generating dummy data as fallback")
        molecules = []
        # for i in range(request.num_molecules):
        #     molecules.append({
        #         "_id": str(uuid.uuid4()),
        #         "smile": request.smiles,  # Use the input SMILES as a placeholder
        #         "name": f"GenMol Fallback {i+1}",
        #         "description": "Dummy data (NVIDIA GenMol API unavailable)",
        #         "parent_drug_id": None,
        #         "generation_method": "dummy_data",
        #         "is_valid": True,
        #         "created_at": datetime.now(),
        #         "updated_at": datetime.now(),
        #         "properties": {
        #             request.scoring: random.random()  # Random score between 0 and 1
        #         }
        #     })
        
        return molecules

    async def validate_molecule(self, smile: str) -> MoleculeValidationResult:
        """
        Validate a molecule structure.
        
        Args:
            smile: SMILES string to validate
            
        Returns:
            Validation result with status and errors if any
        """
        logger.info(f"Validating molecule: {smile}")
        
        result = {
            "is_valid": False,
            "errors": [],
            "warnings": [],
            "smile": smile
        }
        
        try:
            mol = Chem.MolFromSmiles(smile)
            if mol is None:
                result["errors"].append("Invalid SMILES string")
                return result
                
            # Check for other potential issues
            mol = Chem.AddHs(mol)
            
            # Basic checks
            if mol.GetNumAtoms() == 0:
                result["errors"].append("Molecule has no atoms")
            elif mol.GetNumAtoms() > 100:
                result["warnings"].append("Molecule is very large (>100 atoms)")
                
            # Check for undesirable substructures (example)
            if Chem.MolFromSmarts("C(=O)N") and mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)N")):
                result["warnings"].append("Contains amide group")
            
            # If no errors, it's valid
            if not result["errors"]:
                result["is_valid"] = True
                
            return result
            
        except Exception as e:
            logger.error(f"Error validating SMILES {smile}: {e}")
            result["errors"].append(f"Validation error: {str(e)}")
            return result

    async def validate_molecules_bulk(self, smiles_list: List[str]) -> BulkMoleculeValidationResult:
        """
        Validate multiple molecule structures.
        
        Args:
            smiles_list: List of SMILES strings to validate
            
        Returns:
            Bulk validation result with status and errors for each molecule
        """
        logger.info(f"Validating {len(smiles_list)} molecules")
        
        results = []
        
        for smile in smiles_list:
            result = await self.validate_molecule(smile)
            results.append(result)
            
        return {
            "total": len(results),
            "valid_count": sum(1 for r in results if r["is_valid"]),
            "invalid_count": sum(1 for r in results if not r["is_valid"]),
            "results": results
        }

    async def get_molecules_by_drug(self, drug_id: str, skip: int = 0, limit: int = 100) -> List[Molecule]:
        """
        Get molecules associated with a specific drug.
        
        Args:
            drug_id: The ID of the drug
            skip: Number of records to skip for pagination
            limit: Maximum number of records to return
            
        Returns:
            List of molecules associated with the drug
        """
        logger.info(f"Getting molecules for drug: {drug_id}")
        
        try:
            # Query the database for molecules with this parent_drug_id
            cursor = self.collection.find({"parent_drug_id": drug_id}).skip(skip).limit(limit)
            
            molecules = []
            async for doc in cursor:
                # Convert ObjectId to string for serialization
                doc["_id"] = str(doc["_id"])
                molecules.append(doc)
                
            return molecules
            
        except Exception as e:
            logger.error(f"Error fetching molecules for drug {drug_id}: {e}")
            return []

    async def get_molecule_by_id(self, molecule_id: str) -> Optional[Molecule]:
        """
        Get detailed information about a specific molecule.
        
        Args:
            molecule_id: The ID of the molecule
            
        Returns:
            Molecule object if found, None otherwise
        """
        logger.info(f"Getting molecule with ID: {molecule_id}")
        
        try:
            # Try to convert string ID to ObjectId if it's in that format
            query = {"_id": molecule_id}
            try:
                if ObjectId.is_valid(molecule_id):
                    query = {"_id": ObjectId(molecule_id)}
            except:
                pass
                
            # Query the database
            doc = await self.collection.find_one(query)
            
            if doc:
                # Convert ObjectId to string for serialization
                doc["_id"] = str(doc["_id"])
                return doc
                
            return None
            
        except Exception as e:
            logger.error(f"Error fetching molecule {molecule_id}: {e}")
            return None

    async def create_molecule(self, molecule: MoleculeCreate) -> Molecule:
        """
        Add a new molecule to the database.
        
        Args:
            molecule: Molecule data to add
            
        Returns:
            Added molecule with generated ID
        """
        logger.info(f"Creating molecule: {molecule.smile}")
        
        try:
            # Validate the molecule first
            validation = await self.validate_molecule(molecule.smile)
            
            # Create a new molecule document
            new_molecule = {
                "_id": str(uuid.uuid4()),
                "smile": molecule.smile,
                "name": molecule.name or f"Molecule {molecule.smile[:10]}...",
                "description": molecule.description or "",
                "parent_drug_id": molecule.parent_drug_id,
                "generation_method": molecule.generation_method or "manual",
                "is_valid": validation["is_valid"],
                "created_at": datetime.now(),
                "updated_at": datetime.now()
            }
            
            # Insert into the database
            result = await self.collection.insert_one(new_molecule)
            
            # Return the created molecule
            new_molecule["_id"] = str(new_molecule["_id"])
            return new_molecule
            
        except Exception as e:
            logger.error(f"Error creating molecule: {e}")
            # Re-raise the exception for the API to handle
            raise 

    async def get_pubchem_image(self, smiles: str) -> str:
        """Fetch molecule image from PubChem and return as base64 encoded string"""
        logger.info(f"Fetching PubChem image for SMILES: {smiles[:30]}...")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles)}/PNG"
        print(f"URL: {url}")
        try:
            async with aiohttp.ClientSession() as session:
                logger.debug(f"Sending request to PubChem: {url}")
                async with session.get(url, timeout=10) as response:
                    if response.status == 200:
                        image_bytes = await response.read()
                        image_size = len(image_bytes)
                        logger.info(f"Successfully fetched PubChem image ({image_size} bytes)")
                        # Return just the base64 encoded string without the prefix
                        # The frontend will add the prefix as needed
                        return base64.b64encode(image_bytes).decode('utf-8')
                    else:
                        error_text = await response.text()
                        logger.warning(f"Failed to get PubChem image: status={response.status}, response={error_text[:100]}")
                        return ""
        except aiohttp.ClientError as e:
            logger.error(f"Network error fetching PubChem image: {str(e)}")
            return ""
        except Exception as e:
            logger.error(f"Unexpected error fetching PubChem image: {str(e)}")
            return "" 