from fastapi import APIRouter, Depends, HTTPException, Body, Query, Path
from typing import List, Optional
from app.schemas.drug import (
    Drug, DrugCreate, DrugSearchResult, DrugGenerationRequest, DrugGenerationResult
)
from app.services.drug_service import DrugService
from app.core.config import settings

router = APIRouter(prefix=f"{settings.API_V1_STR}/drug-discovery", tags=["Drug Discovery"])

@router.get("/search", response_model=DrugSearchResult)
async def search_drugs(
    disease: Optional[str] = Query(None, description="Disease name to search for drugs"),
    keyword: Optional[str] = Query(None, description="Keyword to search for drugs"),
    source: Optional[str] = Query(None, description="Source to search (chembl, pubchem, zinc, papers, all)"),
    limit: int = Query(20, description="Maximum number of results to return")
):
    """
    Search for drug candidates.
    
    This endpoint searches for drug candidates based on disease name or keywords.
    It can search across multiple databases including ChEMBL, PubChem, ZINC, and research papers.
    """
    if not disease and not keyword:
        raise HTTPException(status_code=400, detail="Either 'disease' or 'keyword' parameter is required")
    
    try:
        drug_service = DrugService()
        return await drug_service.search_drugs(
            disease=disease, 
            keyword=keyword, 
            source=source or "all", 
            limit=limit
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error searching drugs: {str(e)}")

@router.get("/by-disease/{disease_id}", response_model=List[Drug])
async def get_drugs_by_disease(
    disease_id: str = Path(..., description="Disease ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get drug candidates for a specific disease.
    
    This endpoint returns a list of drug candidates that are associated with a specific disease.
    """
    try:
        drug_service = DrugService()
        return await drug_service.get_drugs_by_disease(disease_id, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching drugs: {str(e)}")

@router.post("/generate", response_model=DrugGenerationResult)
async def generate_drug_candidates(request: DrugGenerationRequest = Body(...)):
    """
    Generate new drug candidates.
    
    This endpoint generates new drug candidates based on seed drugs for a specific disease.
    It uses generative models to create structurally valid and diverse molecules.
    """
    try:
        drug_service = DrugService()
        return await drug_service.generate_drug_candidates(
            disease_id=request.disease_id,
            seed_drug_ids=request.seed_drug_ids,
            count=request.count,
            diversity_factor=request.diversity_factor
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating drug candidates: {str(e)}")

@router.get("/{drug_id}", response_model=Drug)
async def get_drug(drug_id: str = Path(..., description="Drug ID")):
    """
    Get detailed information about a specific drug.
    
    This endpoint returns detailed information about a specific drug candidate,
    including its properties, similar drugs, and associated molecules.
    """
    try:
        drug_service = DrugService()
        drug = await drug_service.get_drug_by_id(drug_id)
        if not drug:
            raise HTTPException(status_code=404, detail=f"Drug with ID {drug_id} not found")
        return drug
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching drug: {str(e)}")

@router.post("/", response_model=Drug)
async def create_drug(drug: DrugCreate = Body(...)):
    """
    Add a new drug to the database.
    
    This endpoint allows adding a new drug candidate to the database.
    """
    try:
        drug_service = DrugService()
        return await drug_service.create_drug(drug)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating drug: {str(e)}") 