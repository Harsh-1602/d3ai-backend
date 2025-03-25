import logging
from typing import Dict, Any, List, Optional
from datetime import datetime
from bson.objectid import ObjectId
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Crippen
from rdkit.Chem import Draw
import numpy as np

from app.utils.db import get_database
from app.schemas.property import PropertyCreate, Property, PropertyFilterRequest

logger = logging.getLogger(__name__)

class PropertyService:
    def __init__(self):
        """Initialize PropertyService."""
        self.db = get_database()
        self.collection = self.db["properties"]
        self.available_properties = {
            "molecular_weight": Descriptors.ExactMolWt,
            "logp": Descriptors.MolLogP,
            "hbd": Descriptors.NumHDonors,
            "hba": Descriptors.NumHAcceptors,
            "tpsa": Descriptors.TPSA,
            "rotatable_bonds": Descriptors.NumRotatableBonds,
            "aromatic_rings": lambda m: len(list(m.GetAromaticRings())),
            "heavy_atoms": Descriptors.HeavyAtomCount,
            "formal_charge": lambda m: Chem.GetFormalCharge(m)
        }

    async def calculate_properties(self, smiles: str) -> Dict[str, Any]:
        """
        Calculate molecular properties for a given SMILES string.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            Dictionary of calculated properties
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": "Invalid SMILES string"}

            # Basic RDKit properties
            properties = {}
            for name, descriptor in self.available_properties.items():
                try:
                    value = descriptor(mol)
                    properties[name] = float(value) if isinstance(value, (int, float, np.float64)) else value
                except Exception as e:
                    logger.error(f"Error calculating {name}: {e}")
                    properties[name] = None

            # Calculate additional properties from schema
            properties.update({
                "bioavailability": self._estimate_bioavailability(mol),
                "toxicity": self._estimate_toxicity(mol),
                "solubility": self._estimate_solubility(mol),
                "druglikeness": self._estimate_druglikeness(mol),
                "lipophilicity": Crippen.MolLogP(mol) if mol else None,
                "half_life": self._estimate_half_life(mol),
                "clearance_rate": self._estimate_clearance_rate(mol)
            })

            # Add Lipinski's Rule of 5 analysis
            properties["lipinski"] = self._check_lipinski(properties)
            
            return properties
            
        except Exception as e:
            logger.error(f"Error calculating properties for {smiles}: {e}")
            return {"error": str(e)}

    def _estimate_bioavailability(self, mol) -> float:
        """Estimate bioavailability score (0-1)."""
        if mol is None:
            return 0.0
        # Simple estimation based on Lipinski's rules
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        score = 1.0
        if mw > 500: score *= 0.8
        if logp > 5: score *= 0.8
        if hbd > 5: score *= 0.8
        if hba > 10: score *= 0.8
        return score

    def _estimate_toxicity(self, mol) -> float:
        """Estimate toxicity score (0-1), where 0 is least toxic."""
        if mol is None:
            return 1.0
        # Simple estimation based on structural features
        score = 0.0
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[N+]")): score += 0.2
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[S+]")): score += 0.2
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[#15]")): score += 0.3
        if Descriptors.ExactMolWt(mol) > 800: score += 0.3
        return min(score, 1.0)

    def _estimate_solubility(self, mol) -> float:
        """Estimate solubility score (0-1)."""
        if mol is None:
            return 0.0
        # Simple estimation based on properties
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        score = 1.0
        if logp > 5: score *= 0.7
        if tpsa < 20: score *= 0.7
        return score

    def _estimate_druglikeness(self, mol) -> float:
        """Estimate druglikeness score (0-1)."""
        if mol is None:
            return 0.0
        # Based on Lipinski's Rule of 5 and additional drug-like properties
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        
        score = 1.0
        if mw < 160 or mw > 500: score *= 0.8
        if logp < -0.4 or logp > 5.6: score *= 0.8
        if hbd > 5: score *= 0.8
        if hba > 10: score *= 0.8
        if tpsa > 140: score *= 0.8
        if rotatable > 10: score *= 0.8
        return score

    def _estimate_half_life(self, mol) -> float:
        """Estimate half-life in hours."""
        if mol is None:
            return 0.0
        # Simplified estimation based on molecular properties
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        # Very simple approximation
        return max(1.0, min(24.0, (mw/500) * (abs(logp)/5) * 12))

    def _estimate_clearance_rate(self, mol) -> float:
        """Estimate clearance rate."""
        if mol is None:
            return 0.0
        # Simplified estimation, inverse relationship with half-life
        half_life = self._estimate_half_life(mol)
        return 0.693 / half_life if half_life > 0 else 0.0

    async def calculate_properties_bulk(self, smiles_list: List[str]) -> List[Dict[str, Any]]:
        """
        Calculate properties for multiple SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            List of property dictionaries
        """
        results = []
        for smiles in smiles_list:
            properties = await self.calculate_properties(smiles)
            results.append({
                "smiles": smiles,
                "properties": properties
            })
        return results

    async def create_property(self, property_data: PropertyCreate) -> Property:
        """
        Create a new property record in the database.
        
        Args:
            property_data: Property data to create
            
        Returns:
            Created property record
        """
        property_dict = property_data.model_dump()
        property_dict["_id"] = str(ObjectId())
        property_dict["created_at"] = datetime.utcnow()
        property_dict["updated_at"] = datetime.utcnow()
        
        await self.collection.insert_one(property_dict)
        return Property(**property_dict)

    async def get_property(self, property_id: str) -> Optional[Property]:
        """
        Get a property record by ID.
        
        Args:
            property_id: ID of the property record
            
        Returns:
            Property record if found
        """
        result = await self.collection.find_one({"_id": property_id})
        if result:
            return Property(**result)
        return None

    async def filter_molecules(self, filter_request: PropertyFilterRequest) -> List[str]:
        """
        Filter molecules based on property criteria.
        
        Args:
            filter_request: Filter criteria
            
        Returns:
            List of molecule IDs that pass the filters
        """
        query = {"molecule_id": {"$in": filter_request.molecule_ids}}
        
        # Add property filters
        if filter_request.bioavailability_min is not None:
            query["bioavailability"] = {"$gte": filter_request.bioavailability_min}
        if filter_request.toxicity_max is not None:
            query["toxicity"] = {"$lte": filter_request.toxicity_max}
        if filter_request.solubility_min is not None:
            query["solubility"] = {"$gte": filter_request.solubility_min}
        if filter_request.druglikeness_min is not None:
            query["druglikeness"] = {"$gte": filter_request.druglikeness_min}
            
        # Add custom filters
        for prop, conditions in filter_request.custom_filters.items():
            for op, value in conditions.items():
                if op == "min":
                    query[prop] = query.get(prop, {})
                    query[prop]["$gte"] = value
                elif op == "max":
                    query[prop] = query.get(prop, {})
                    query[prop]["$lte"] = value
                    
        cursor = self.collection.find(query)
        results = []
        async for doc in cursor:
            results.append(doc["molecule_id"])
            
        return results

    def _check_lipinski(self, properties: Dict[str, Any]) -> Dict[str, bool]:
        """
        Check Lipinski's Rule of 5 compliance.
        
        Args:
            properties: Dictionary of molecular properties
            
        Returns:
            Dictionary indicating compliance with each rule
        """
        return {
            "molecular_weight_ok": properties.get("molecular_weight", 0) <= 500,
            "logp_ok": -0.4 <= properties.get("logp", 0) <= 5.6,
            "hbd_ok": properties.get("hbd", 0) <= 5,
            "hba_ok": properties.get("hba", 0) <= 10,
            "passes_all": all([
                properties.get("molecular_weight", 0) <= 500,
                -0.4 <= properties.get("logp", 0) <= 5.6,
                properties.get("hbd", 0) <= 5,
                properties.get("hba", 0) <= 10
            ])
        }

    async def get_property_ranges(self, smiles_list: List[str]) -> Dict[str, Dict[str, float]]:
        """
        Calculate property ranges for a set of molecules.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            Dictionary of property ranges with min, max, mean values
        """
        all_properties = []
        for smiles in smiles_list:
            props = await self.calculate_properties(smiles)
            if "error" not in props:
                all_properties.append(props)

        if not all_properties:
            return {}

        ranges = {}
        for prop in self.available_properties.keys():
            values = [p[prop] for p in all_properties if prop in p and p[prop] is not None]
            if values:
                ranges[prop] = {
                    "min": float(min(values)),
                    "max": float(max(values)),
                    "mean": float(np.mean(values)),
                    "std": float(np.std(values))
                }

        return ranges 