from fastapi import APIRouter, Depends, HTTPException, Body, Query, Path
from typing import List, Optional, Dict, Any
from app.schemas.molecule import (
    Molecule, MoleculeCreate, MoleculeGenerationRequest, 
    MoleculeValidationResult, BulkMoleculeValidationRequest, BulkMoleculeValidationResult
)
from app.services.molecule_service import MoleculeService
from app.utils.molecular_visualization import (
    smiles_to_2d_image, smiles_to_3d_coordinates, smiles_to_molecule_properties,
    compare_molecules_similarity, create_molecule_grid_image
)
from app.core.config import settings

router = APIRouter(prefix=f"{settings.API_V1_STR}/molecules", tags=["Molecules"])

@router.post("/generate", response_model=List[Molecule])
async def generate_molecules(request: MoleculeGenerationRequest = Body(...)):
    """
    Generate new molecules from seed SMILES.
    
    This endpoint generates new molecules based on provided seed SMILES strings.
    It uses generative models to create structurally valid and diverse molecules.
    """
    try:
        molecule_service = MoleculeService()
        return await molecule_service.generate_molecules(
            seed_smiles=request.seed_smiles,
            count=request.count,
            diversity_factor=request.diversity_factor,
            similarity_threshold=request.similarity_threshold
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating molecules: {str(e)}")

@router.post("/validate", response_model=MoleculeValidationResult)
async def validate_molecule(smile: str = Body(..., embed=True)):
    """
    Validate a molecule structure.
    
    This endpoint checks if a SMILES string represents a valid molecular structure.
    It returns validation status and errors if any.
    """
    try:
        molecule_service = MoleculeService()
        return await molecule_service.validate_molecule(smile)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error validating molecule: {str(e)}")

@router.post("/validate-bulk", response_model=BulkMoleculeValidationResult)
async def validate_molecules_bulk(request: BulkMoleculeValidationRequest = Body(...)):
    """
    Validate multiple molecule structures.
    
    This endpoint checks if multiple SMILES strings represent valid molecular structures.
    It returns validation status and errors for each molecule.
    """
    try:
        molecule_service = MoleculeService()
        return await molecule_service.validate_molecules_bulk(request.smiles)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error validating molecules: {str(e)}")

@router.get("/by-drug/{drug_id}", response_model=List[Molecule])
async def get_molecules_by_drug(
    drug_id: str = Path(..., description="Drug ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get molecules associated with a specific drug.
    
    This endpoint returns a list of molecules that are associated with a specific drug.
    """
    try:
        molecule_service = MoleculeService()
        return await molecule_service.get_molecules_by_drug(drug_id, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching molecules: {str(e)}")

@router.get("/{molecule_id}", response_model=Molecule)
async def get_molecule(molecule_id: str = Path(..., description="Molecule ID")):
    """
    Get detailed information about a specific molecule.
    
    This endpoint returns detailed information about a specific molecule,
    including its properties and visualization data.
    """
    try:
        molecule_service = MoleculeService()
        molecule = await molecule_service.get_molecule_by_id(molecule_id)
        if not molecule:
            raise HTTPException(status_code=404, detail=f"Molecule with ID {molecule_id} not found")
        return molecule
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching molecule: {str(e)}")

@router.post("/", response_model=Molecule)
async def create_molecule(molecule: MoleculeCreate = Body(...)):
    """
    Add a new molecule to the database.
    
    This endpoint allows adding a new molecule to the database.
    """
    try:
        molecule_service = MoleculeService()
        return await molecule_service.create_molecule(molecule)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating molecule: {str(e)}")

@router.post("/visualize/2d", response_model=Dict[str, Any])
async def visualize_molecule_2d(
    smile: str = Body(..., embed=True),
    width: int = Body(300, embed=True),
    height: int = Body(200, embed=True)
):
    """
    Create a 2D visualization of a molecule from its SMILES string.
    
    This endpoint converts a SMILES string to a 2D molecular image in base64 format.
    """
    try:
        # Create 2D visualization
        image = smiles_to_2d_image(smile, width, height)
        
        if not image:
            raise HTTPException(status_code=400, detail="Failed to create 2D visualization")
        
        # Get basic properties
        properties = smiles_to_molecule_properties(smile)
        
        return {
            "image": image,
            "properties": properties,
            "smile": smile
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating visualization: {str(e)}")

@router.post("/visualize/3d", response_model=Dict[str, Any])
async def visualize_molecule_3d(smile: str = Body(..., embed=True)):
    """
    Get 3D coordinates of a molecule from its SMILES string for visualization.
    
    This endpoint converts a SMILES string to 3D coordinates that can be used for visualization.
    """
    try:
        # Create 3D coordinates
        coordinates = smiles_to_3d_coordinates(smile)
        
        if not coordinates:
            raise HTTPException(status_code=400, detail="Failed to create 3D coordinates")
        
        # Get basic properties
        properties = smiles_to_molecule_properties(smile)
        
        return {
            "coordinates": coordinates,
            "properties": properties,
            "smile": smile
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating 3D coordinates: {str(e)}")

@router.post("/visualize/grid", response_model=Dict[str, Any])
async def visualize_molecules_grid(
    smiles: List[str] = Body(..., embed=True),
    legends: Optional[List[str]] = Body(None, embed=True),
    mols_per_row: int = Body(4, embed=True)
):
    """
    Create a grid visualization of multiple molecules from their SMILES strings.
    
    This endpoint creates a grid image of multiple molecules in base64 format.
    """
    try:
        # Validate input
        if not smiles:
            raise HTTPException(status_code=400, detail="No SMILES provided")
        
        if len(smiles) > 20:
            raise HTTPException(status_code=400, detail="Maximum 20 molecules allowed")
        
        # Create grid image
        image = create_molecule_grid_image(smiles, legends, mols_per_row)
        
        if not image:
            raise HTTPException(status_code=400, detail="Failed to create grid visualization")
        
        return {
            "image": image,
            "count": len(smiles)
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating grid visualization: {str(e)}")

@router.post("/compare", response_model=Dict[str, Any])
async def compare_molecules(
    smile1: str = Body(..., embed=True),
    smile2: str = Body(..., embed=True)
):
    """
    Compare two molecules based on their SMILES strings.
    
    This endpoint calculates the similarity between two molecules and returns visualizations for both.
    """
    try:
        # Calculate similarity
        similarity = compare_molecules_similarity(smile1, smile2)
        
        # Create visualizations
        image1 = smiles_to_2d_image(smile1)
        image2 = smiles_to_2d_image(smile2)
        
        if not image1 or not image2:
            raise HTTPException(status_code=400, detail="Failed to create visualizations")
        
        # Get properties
        properties1 = smiles_to_molecule_properties(smile1)
        properties2 = smiles_to_molecule_properties(smile2)
        
        return {
            "similarity": similarity,
            "molecule1": {
                "smile": smile1,
                "image": image1,
                "properties": properties1
            },
            "molecule2": {
                "smile": smile2,
                "image": image2,
                "properties": properties2
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error comparing molecules: {str(e)}")

@router.post("/properties", response_model=Dict[str, Any])
async def get_molecule_properties_from_smile(smile: str = Body(..., embed=True)):
    """
    Get properties of a molecule from its SMILES string.
    
    This endpoint calculates various molecular properties based on the SMILES string.
    """
    try:
        # Get properties
        properties = smiles_to_molecule_properties(smile)
        
        if not properties["is_valid"]:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        return {
            "smile": smile,
            "properties": properties
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error calculating properties: {str(e)}") 