import logging
from datetime import datetime
import uuid
from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from bson.objectid import ObjectId
from pymongo import MongoClient
from pymongo.collection import Collection

from app.core.config import settings
from app.schemas.molecule import Molecule, MoleculeCreate, MoleculeValidationResult, BulkMoleculeValidationResult
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