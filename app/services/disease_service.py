import logging
import uuid
from datetime import datetime
from typing import List, Optional, Dict, Any
from app.db.mongodb import db
from app.schemas.disease import (
    DiseaseCreate, Disease, DiseaseInDB, DiseaseIdentificationResult, DiseasePrediction
)
from app.core.config import settings

logger = logging.getLogger(__name__)

class DiseaseService:
    """Service for disease-related operations."""
    
    async def identify_disease_from_symptoms(self, symptoms: List[str]) -> DiseaseIdentificationResult:
        """
        Identify a disease based on provided symptoms.
        
        This method uses a machine learning model to predict the most likely disease
        based on the provided symptoms.
        """
        logger.info(f"Identifying disease from symptoms: {symptoms}")
        
        # TODO: Implement actual disease identification using ML model
        # For now, we'll use a simple mock implementation
        
        # Get diseases from database that have the most matching symptoms
        diseases_collection = db["diseases"]
        cursor = diseases_collection.find({"symptoms": {"$in": symptoms}})
        diseases = []
        
        for doc in await cursor.to_list(length=10):
            # Calculate a simple confidence score based on symptom overlap
            disease_symptoms = set(doc.get("symptoms", []))
            input_symptoms = set(symptoms)
            matching_symptoms = disease_symptoms.intersection(input_symptoms)
            
            if matching_symptoms:
                confidence = len(matching_symptoms) / max(len(disease_symptoms), len(input_symptoms))
                
                # Create a Disease object
                disease = Disease(
                    id=str(doc.get("_id")),
                    name=doc.get("name"),
                    description=doc.get("description"),
                    symptoms=doc.get("symptoms", []),
                    causes=doc.get("causes", []),
                    icd_code=doc.get("icd_code"),
                    category=doc.get("category"),
                    created_at=doc.get("created_at", datetime.now()),
                    updated_at=doc.get("updated_at", datetime.now()),
                    associated_drugs=[],
                    research_papers_count=0
                )
                
                diseases.append((disease, confidence))
        
        # Sort by confidence score
        diseases.sort(key=lambda x: x[1], reverse=True)
        
        # If no diseases found, return a default response
        if not diseases:
            # Create a mock disease for demonstration
            mock_disease = Disease(
                id="unknown",
                name="Unknown Disease",
                description="Could not identify a specific disease based on the provided symptoms.",
                symptoms=symptoms,
                created_at=datetime.now(),
                updated_at=datetime.now(),
                associated_drugs=[],
                research_papers_count=0
            )
            return DiseaseIdentificationResult(
                disease=mock_disease,
                confidence=0.0,
                alternative_diseases=[]
            )
        
        # Return the result
        main_disease, main_confidence = diseases[0]
        alternative_diseases = [d for d, _ in diseases[1:5]]  # Take up to 4 alternatives
        
        return DiseaseIdentificationResult(
            disease=main_disease,
            confidence=main_confidence,
            alternative_diseases=alternative_diseases
        )
    
    async def get_diseases(
        self, 
        name: Optional[str] = None, 
        symptom: Optional[str] = None, 
        skip: int = 0, 
        limit: int = 100
    ) -> List[Disease]:
        """Get a list of diseases with optional filtering."""
        logger.info(f"Getting diseases with filters: name={name}, symptom={symptom}")
        
        # Build query
        query = {}
        if name:
            query["name"] = {"$regex": name, "$options": "i"}  # Case-insensitive search
        if symptom:
            query["symptoms"] = {"$regex": symptom, "$options": "i"}
        
        # Query database
        diseases_collection = db["diseases"]
        cursor = diseases_collection.find(query).skip(skip).limit(limit)
        
        # Convert to Disease objects
        diseases = []
        async for doc in cursor:
            # Get associated drugs
            drugs_collection = db["drugs"]
            drugs_cursor = drugs_collection.find({"disease_id": str(doc.get("_id"))})
            associated_drugs = [drug.get("name") async for drug in drugs_cursor]
            
            # Create Disease object
            disease = Disease(
                id=str(doc.get("_id")),
                name=doc.get("name"),
                description=doc.get("description"),
                symptoms=doc.get("symptoms", []),
                causes=doc.get("causes", []),
                icd_code=doc.get("icd_code"),
                category=doc.get("category"),
                created_at=doc.get("created_at", datetime.now()),
                updated_at=doc.get("updated_at", datetime.now()),
                associated_drugs=associated_drugs,
                research_papers_count=0  # TODO: Implement actual count
            )
            diseases.append(disease)
        
        return diseases
    
    async def get_disease_by_id(self, disease_id: str) -> Optional[Disease]:
        """Get a specific disease by ID."""
        logger.info(f"Getting disease by ID: {disease_id}")
        
        # Query database
        diseases_collection = db["diseases"]
        doc = await diseases_collection.find_one({"_id": disease_id})
        
        if not doc:
            return None
        
        # Get associated drugs
        drugs_collection = db["drugs"]
        drugs_cursor = drugs_collection.find({"disease_id": disease_id})
        associated_drugs = [drug.get("name") async for drug in drugs_cursor]
        
        # Create Disease object
        disease = Disease(
            id=str(doc.get("_id")),
            name=doc.get("name"),
            description=doc.get("description"),
            symptoms=doc.get("symptoms", []),
            causes=doc.get("causes", []),
            icd_code=doc.get("icd_code"),
            category=doc.get("category"),
            created_at=doc.get("created_at", datetime.now()),
            updated_at=doc.get("updated_at", datetime.now()),
            associated_drugs=associated_drugs,
            research_papers_count=0  # TODO: Implement actual count
        )
        
        return disease
    
    async def create_disease(self, disease: DiseaseCreate) -> Disease:
        """Create a new disease."""
        logger.info(f"Creating new disease: {disease.name}")
        
        # Create document
        disease_id = str(uuid.uuid4())
        now = datetime.now()
        
        disease_doc = {
            "_id": disease_id,
            "name": disease.name,
            "description": disease.description,
            "symptoms": disease.symptoms,
            "causes": disease.causes,
            "icd_code": disease.icd_code,
            "category": disease.category,
            "created_at": now,
            "updated_at": now
        }
        
        # Insert into database
        diseases_collection = db["diseases"]
        await diseases_collection.insert_one(disease_doc)
        
        # Return created disease
        return Disease(
            id=disease_id,
            name=disease.name,
            description=disease.description,
            symptoms=disease.symptoms,
            causes=disease.causes,
            icd_code=disease.icd_code,
            category=disease.category,
            created_at=now,
            updated_at=now,
            associated_drugs=[],
            research_papers_count=0
        )
    
    async def search_diseases(self, query: str, limit: int = 10) -> List[DiseasePrediction]:
        """Search for diseases by name or symptoms."""
        logger.info(f"Searching diseases with query: {query}")
        
        # Build search query
        search_query = {
            "$or": [
                {"name": {"$regex": query, "$options": "i"}},
                {"symptoms": {"$regex": query, "$options": "i"}},
                {"description": {"$regex": query, "$options": "i"}}
            ]
        }
        
        # Query database
        diseases_collection = db["diseases"]
        cursor = diseases_collection.find(search_query).limit(limit)
        
        # Convert to DiseasePrediction objects
        predictions = []
        async for doc in cursor:
            # Calculate confidence and matching symptoms
            disease_symptoms = doc.get("symptoms", [])
            matching_symptoms = [s for s in disease_symptoms if query.lower() in s.lower()]
            not_matching_symptoms = [s for s in disease_symptoms if s not in matching_symptoms]
            
            # Simple confidence calculation
            confidence = 0.5  # Base confidence
            if matching_symptoms:
                confidence += 0.5 * (len(matching_symptoms) / len(disease_symptoms))
            
            # Create DiseasePrediction object
            prediction = DiseasePrediction(
                disease_name=doc.get("name"),
                confidence=min(confidence, 1.0),  # Cap at 1.0
                symptoms_matched=matching_symptoms,
                symptoms_not_matched=not_matching_symptoms
            )
            predictions.append(prediction)
        
        # Sort by confidence
        predictions.sort(key=lambda x: x.confidence, reverse=True)
        
        return predictions
    
    async def suggest_diseases(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Suggest diseases based on a partial query string.
        
        This method provides autocomplete functionality for disease names
        as the user types in the search box.
        
        Args:
            query: Partial query string entered by the user
            limit: Maximum number of suggestions to return
            
        Returns:
            List of disease suggestions with id and name
        """
        logger.info(f"Suggesting diseases for query: {query}")
        
        if not query or len(query) < 2:
            return []
            
        # Create regex for case-insensitive partial matching
        regex_pattern = {"$regex": f".*{query}.*", "$options": "i"}
        
        # Query database for diseases that match the pattern
        diseases_collection = db["diseases"]
        cursor = diseases_collection.find(
            {"$or": [
                {"name": regex_pattern},
                {"aliases": regex_pattern}
            ]},
            projection={"_id": 1, "name": 1, "aliases": 1}
        ).limit(limit)
        
        # Convert to list and format results
        suggestions = []
        async for doc in cursor:
            suggestions.append({
                "id": doc["_id"],
                "name": doc["name"],
                "aliases": doc.get("aliases", [])
            })
            
        return suggestions 