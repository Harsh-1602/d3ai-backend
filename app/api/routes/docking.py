from fastapi import APIRouter, Depends, HTTPException, Body, Query, Path
from typing import List, Optional
from app.schemas.docking import (
    Protein, ProteinCreate, Docking, DockingCreate, DockingRequest, 
    DockingResult, BulkDockingResult
)
from app.services.docking_service import DockingService
from app.core.config import settings

router = APIRouter(prefix=f"{settings.API_V1_STR}/docking", tags=["Molecular Docking"])

@router.post("/dock", response_model=BulkDockingResult)
async def dock_molecules(request: DockingRequest = Body(...)):
    """
    Perform molecular docking.
    
    This endpoint performs docking of molecules with a protein target.
    It predicts binding affinity and binding poses for each molecule.
    """
    try:
        docking_service = DockingService()
        return await docking_service.dock_molecules(
            molecule_ids=request.molecule_ids,
            protein_id=request.protein_id,
            exhaustiveness=request.exhaustiveness,
            num_modes=request.num_modes
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error performing docking: {str(e)}")

@router.get("/proteins", response_model=List[Protein])
async def get_proteins(
    disease_id: Optional[str] = Query(None, description="Filter by disease ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get a list of proteins.
    
    This endpoint returns a list of proteins from the database.
    You can filter by disease ID or get all proteins with pagination.
    """
    try:
        docking_service = DockingService()
        return await docking_service.get_proteins(disease_id=disease_id, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching proteins: {str(e)}")

@router.get("/proteins/{protein_id}", response_model=Protein)
async def get_protein(protein_id: str = Path(..., description="Protein ID")):
    """
    Get a specific protein by ID.
    
    This endpoint returns detailed information about a specific protein.
    """
    try:
        docking_service = DockingService()
        protein = await docking_service.get_protein_by_id(protein_id)
        if not protein:
            raise HTTPException(status_code=404, detail=f"Protein with ID {protein_id} not found")
        return protein
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching protein: {str(e)}")

@router.post("/proteins", response_model=Protein)
async def create_protein(protein: ProteinCreate = Body(...)):
    """
    Add a new protein to the database.
    
    This endpoint allows adding a new protein to the database.
    """
    try:
        docking_service = DockingService()
        return await docking_service.create_protein(protein)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating protein: {str(e)}")

@router.get("/results", response_model=List[Docking])
async def get_docking_results(
    molecule_id: Optional[str] = Query(None, description="Filter by molecule ID"),
    protein_id: Optional[str] = Query(None, description="Filter by protein ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get docking results.
    
    This endpoint returns docking results from the database.
    You can filter by molecule ID, protein ID, or get all results with pagination.
    """
    try:
        docking_service = DockingService()
        return await docking_service.get_docking_results(
            molecule_id=molecule_id, 
            protein_id=protein_id, 
            skip=skip, 
            limit=limit
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching docking results: {str(e)}")

@router.get("/results/{docking_id}", response_model=Docking)
async def get_docking_result(docking_id: str = Path(..., description="Docking ID")):
    """
    Get a specific docking result by ID.
    
    This endpoint returns detailed information about a specific docking result.
    """
    try:
        docking_service = DockingService()
        docking = await docking_service.get_docking_result_by_id(docking_id)
        if not docking:
            raise HTTPException(status_code=404, detail=f"Docking result with ID {docking_id} not found")
        return docking
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching docking result: {str(e)}")

@router.post("/results", response_model=Docking)
async def create_docking_result(docking: DockingCreate = Body(...)):
    """
    Add a new docking result to the database.
    
    This endpoint allows adding a new docking result to the database.
    """
    try:
        docking_service = DockingService()
        return await docking_service.create_docking_result(docking)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating docking result: {str(e)}") 