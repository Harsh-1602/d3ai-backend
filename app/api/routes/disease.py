from fastapi import APIRouter, Depends, HTTPException, Body, Query
from fastapi.responses import JSONResponse
from typing import List, Optional, Dict, Any
from app.schemas.disease import (
    SymptomInput, Disease, DiseaseCreate, DiseaseIdentificationResult, DiseasePrediction
)
from app.services.disease_service import DiseaseService
from app.core.config import settings
import logging

logger = logging.getLogger(__name__)

router = APIRouter(prefix=f"{settings.API_V1_STR}/diseases", tags=["Diseases"])

@router.post("/identify", response_model=DiseaseIdentificationResult)
async def identify_disease(symptoms: SymptomInput):
    """
    Identify a disease based on provided symptoms.
    
    This endpoint takes a list of symptoms and uses AI to predict the most likely disease.
    It returns the predicted disease along with confidence scores and alternative options.
    """
    try:
        disease_service = DiseaseService()
        return await disease_service.identify_disease_from_symptoms(symptoms.symptoms)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error identifying disease: {str(e)}")

@router.get("/", response_model=List[Disease])
async def get_diseases(
    name: Optional[str] = Query(None, description="Filter by disease name"),
    symptom: Optional[str] = Query(None, description="Filter by symptom"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get a list of diseases.
    
    This endpoint returns a list of diseases from the database.
    You can filter by name, symptom, or get all diseases with pagination.
    """
    try:
        disease_service = DiseaseService()
        return await disease_service.get_diseases(name=name, symptom=symptom, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching diseases: {str(e)}")

@router.get("/{disease_id}", response_model=Disease)
async def get_disease(disease_id: str):
    """
    Get a specific disease by ID.
    
    This endpoint returns detailed information about a specific disease.
    """
    try:
        disease_service = DiseaseService()
        disease = await disease_service.get_disease_by_id(disease_id)
        if not disease:
            raise HTTPException(status_code=404, detail=f"Disease with ID {disease_id} not found")
        return disease
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching disease: {str(e)}")

@router.post("/", response_model=Disease)
async def create_disease(disease: DiseaseCreate = Body(...)):
    """
    Create a new disease.
    
    This endpoint allows adding a new disease to the database.
    """
    try:
        disease_service = DiseaseService()
        return await disease_service.create_disease(disease)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating disease: {str(e)}")

@router.get("/search/{query}", response_model=List[DiseasePrediction])
async def search_diseases(
    query: str,
    limit: int = Query(10, description="Maximum number of results to return")
):
    """
    Search for diseases by name or symptoms.
    
    This endpoint searches for diseases that match the provided query string.
    """
    try:
        disease_service = DiseaseService()
        return await disease_service.search_diseases(query, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error searching diseases: {str(e)}")

@router.get("/suggest/{query}")
async def suggest_diseases(
    query: str,
    limit: int = Query(10, description="Maximum number of suggestions to return")
):
    """
    Suggest diseases based on a partial query string.
    
    This endpoint provides autocomplete functionality for disease names
    as the user types in the search box.
    """
    try:
        logger.info(f"Disease suggestion request for: {query}")
        # Instead of using the database, use a static list of diseases
        # This is a temporary solution until we integrate with a more comprehensive source
        results = await suggest_diseases_static(query, limit)
        logger.info(f"Found {len(results)} suggestions for '{query}'")
        return JSONResponse(content=results)
    except Exception as e:
        logger.error(f"Error suggesting diseases: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Error suggesting diseases: {str(e)}")
        
async def suggest_diseases_static(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """
    Suggest diseases from a static list without using the database.
    """
    # List of common diseases with their IDs and aliases
    common_diseases = [
        {"id": "d001", "name": "Diabetes Mellitus", "aliases": ["Type 2 Diabetes", "Diabetes"]},
        {"id": "d002", "name": "Hypertension", "aliases": ["High Blood Pressure", "HTN"]},
        {"id": "d003", "name": "Asthma", "aliases": ["Bronchial Asthma"]},
        {"id": "d004", "name": "Alzheimer's Disease", "aliases": ["Alzheimer Disease", "Senile Dementia"]},
        {"id": "d005", "name": "Parkinson's Disease", "aliases": ["Parkinson Disease"]},
        {"id": "d006", "name": "Coronary Heart Disease", "aliases": ["Coronary Artery Disease", "CAD", "Ischemic Heart Disease"]},
        {"id": "d007", "name": "Cancer", "aliases": ["Malignant Neoplasm", "Malignancy"]},
        {"id": "d008", "name": "Lung Cancer", "aliases": ["Pulmonary Cancer", "Bronchogenic Carcinoma"]},
        {"id": "d009", "name": "Breast Cancer", "aliases": ["Mammary Cancer"]},
        {"id": "d010", "name": "Prostate Cancer", "aliases": ["Prostatic Carcinoma"]},
        {"id": "d011", "name": "Colorectal Cancer", "aliases": ["Colon Cancer", "Rectal Cancer"]},
        {"id": "d012", "name": "Stroke", "aliases": ["Cerebrovascular Accident", "CVA", "Brain Attack"]},
        {"id": "d013", "name": "Chronic Obstructive Pulmonary Disease", "aliases": ["COPD", "Emphysema", "Chronic Bronchitis"]},
        {"id": "d014", "name": "Arthritis", "aliases": ["Joint Inflammation"]},
        {"id": "d015", "name": "Rheumatoid Arthritis", "aliases": ["RA"]},
        {"id": "d016", "name": "Osteoarthritis", "aliases": ["Degenerative Joint Disease", "DJD"]},
        {"id": "d017", "name": "Chronic Kidney Disease", "aliases": ["CKD", "Renal Failure"]},
        {"id": "d018", "name": "Multiple Sclerosis", "aliases": ["MS"]},
        {"id": "d019", "name": "Schizophrenia", "aliases": []},
        {"id": "d020", "name": "Bipolar Disorder", "aliases": ["Manic Depression"]},
        {"id": "d021", "name": "Depression", "aliases": ["Major Depressive Disorder", "MDD", "Clinical Depression"]},
        {"id": "d022", "name": "Anxiety Disorders", "aliases": ["Generalized Anxiety Disorder", "GAD"]},
        {"id": "d023", "name": "HIV/AIDS", "aliases": ["Human Immunodeficiency Virus Infection", "Acquired Immunodeficiency Syndrome"]},
        {"id": "d024", "name": "Tuberculosis", "aliases": ["TB"]},
        {"id": "d025", "name": "Malaria", "aliases": []},
        {"id": "d026", "name": "Influenza", "aliases": ["Flu"]},
        {"id": "d027", "name": "Pneumonia", "aliases": ["Lung Infection"]},
        {"id": "d028", "name": "Hepatitis", "aliases": ["Liver Inflammation"]},
        {"id": "d029", "name": "Cirrhosis", "aliases": ["Liver Cirrhosis"]},
        {"id": "d030", "name": "Gastroesophageal Reflux Disease", "aliases": ["GERD", "Acid Reflux"]},
    ]
    
    if not query or len(query.strip()) < 2:
        return []
    
    # Case-insensitive search in disease names and aliases
    query = query.lower().strip()
    results = []
    
    for disease in common_diseases:
        # Check if query matches the disease name
        if query in disease["name"].lower():
            results.append(disease)
            continue
        
        # Check if query matches any of the aliases
        for alias in disease["aliases"]:
            if query in alias.lower():
                results.append(disease)
                break
    
    # Sort results by relevance (exact match first, then starts with, then contains)
    def sort_key(disease):
        name_lower = disease["name"].lower()
        if name_lower == query:
            return 0
        elif name_lower.startswith(query):
            return 1
        else:
            return 2
    
    results.sort(key=sort_key)
    
    # Limit the number of results
    return results[:limit] 