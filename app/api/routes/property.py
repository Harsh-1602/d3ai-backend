from fastapi import APIRouter, Depends, HTTPException, Body, Query, Path
from typing import List, Optional, Union
from app.schemas.property import (
    Property, PropertyCreate, PropertyPredictionRequest, PropertyPredictionResult,
    BulkPropertyPredictionResult, PropertyFilterRequest, PropertyFilterResult
)
from app.services.property_service import PropertyService
from app.core.config import settings

router = APIRouter(prefix=f"{settings.API_V1_STR}/properties", tags=["Molecular Properties"])

@router.post("/predict", response_model=Union[PropertyPredictionResult, BulkPropertyPredictionResult])
async def predict_properties(request: PropertyPredictionRequest = Body(...)):
    """
    Predict properties for molecules.
    
    This endpoint predicts molecular properties for one or more SMILES strings.
    Properties include bioavailability, toxicity, solubility, and more.
    """
    try:
        property_service = PropertyService()
        
        # Check if we're dealing with a single SMILE or multiple
        if isinstance(request.smiles, str):
            return await property_service.predict_properties(request.smiles)
        else:
            return await property_service.predict_properties_bulk(request.smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error predicting properties: {str(e)}")

@router.post("/filter", response_model=PropertyFilterResult)
async def filter_molecules_by_properties(request: PropertyFilterRequest = Body(...)):
    """
    Filter molecules by properties.
    
    This endpoint filters a list of molecules based on their properties.
    You can specify minimum/maximum thresholds for various properties.
    """
    try:
        property_service = PropertyService()
        return await property_service.filter_molecules_by_properties(
            molecule_ids=request.molecule_ids,
            bioavailability_min=request.bioavailability_min,
            toxicity_max=request.toxicity_max,
            solubility_min=request.solubility_min,
            druglikeness_min=request.druglikeness_min,
            custom_filters=request.custom_filters
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error filtering molecules: {str(e)}")

@router.get("/molecule/{molecule_id}", response_model=Property)
async def get_molecule_properties(molecule_id: str = Path(..., description="Molecule ID")):
    """
    Get properties for a specific molecule.
    
    This endpoint returns the properties of a specific molecule.
    """
    try:
        property_service = PropertyService()
        properties = await property_service.get_properties_by_molecule_id(molecule_id)
        if not properties:
            raise HTTPException(status_code=404, detail=f"Properties for molecule with ID {molecule_id} not found")
        return properties
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching properties: {str(e)}")

@router.post("/", response_model=Property)
async def create_property_record(property_data: PropertyCreate = Body(...)):
    """
    Add a new property record to the database.
    
    This endpoint allows adding a new property record for a molecule to the database.
    """
    try:
        property_service = PropertyService()
        return await property_service.create_property_record(property_data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating property record: {str(e)}")

@router.get("/compare", response_model=List[Property])
async def compare_molecule_properties(
    molecule_ids: str = Query(..., description="Comma-separated list of molecule IDs to compare")
):
    """
    Compare properties of multiple molecules.
    
    This endpoint returns the properties of multiple molecules for comparison.
    """
    try:
        property_service = PropertyService()
        molecule_id_list = [id.strip() for id in molecule_ids.split(",")]
        return await property_service.get_properties_for_molecules(molecule_id_list)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error comparing properties: {str(e)}") 