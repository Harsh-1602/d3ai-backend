from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime

class ProteinBase(BaseModel):
    """Base schema for protein data."""
    name: str = Field(..., description="Name of the protein")
    pdb_id: Optional[str] = Field(None, description="PDB ID of the protein")
    uniprot_id: Optional[str] = Field(None, description="UniProt ID of the protein")
    sequence: Optional[str] = Field(None, description="Amino acid sequence of the protein")
    disease_id: Optional[str] = Field(None, description="ID of the associated disease")
    description: Optional[str] = Field(None, description="Description of the protein")
    binding_site: Optional[Dict[str, Any]] = Field(None, description="Binding site information")

class ProteinCreate(ProteinBase):
    """Schema for creating a new protein."""
    pass

class ProteinInDB(ProteinBase):
    """Schema for a protein that exists in the database."""
    id: str = Field(..., description="Unique identifier for the protein")
    created_at: datetime = Field(..., description="When the protein was first added to the database")
    updated_at: datetime = Field(..., description="When the protein was last updated")
    
    class Config:
        from_attributes = True

class Protein(ProteinInDB):
    """Schema for protein response."""
    disease_name: Optional[str] = Field(None, description="Name of the associated disease")
    chembl_id: Optional[str] = Field(None, description="Target ChEMBL ID of the protein")
    
class DockingBase(BaseModel):
    """Base schema for docking data."""
    molecule_id: str = Field(..., description="ID of the molecule")
    protein_id: str = Field(..., description="ID of the protein")
    binding_affinity: float = Field(..., description="Binding affinity score (lower is better)")
    binding_pose: Optional[Dict[str, Any]] = Field(default_factory=dict, description="3D coordinates of the binding pose")
    interactions: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Molecular interactions")

class DockingCreate(DockingBase):
    """Schema for creating a new docking record."""
    pass

class DockingInDB(DockingBase):
    """Schema for a docking record that exists in the database."""
    id: str = Field(..., description="Unique identifier for the docking record")
    created_at: datetime = Field(..., description="When the docking record was first added to the database")
    updated_at: datetime = Field(..., description="When the docking record was last updated")
    
    class Config:
        from_attributes = True

class Docking(DockingInDB):
    """Schema for docking response."""
    molecule_smile: Optional[str] = Field(None, description="SMILE string of the molecule")
    molecule_name: Optional[str] = Field(None, description="Name of the molecule")
    protein_name: Optional[str] = Field(None, description="Name of the protein")
    
class DockingRequest(BaseModel):
    """Schema for requesting docking."""
    molecule_ids: List[str] = Field(..., description="IDs of the molecules to dock")
    protein_id: str = Field(..., description="ID of the protein to dock with")
    exhaustiveness: Optional[int] = Field(8, description="Exhaustiveness of the docking search")
    num_modes: Optional[int] = Field(9, description="Number of binding modes to generate")

# New schemas for NVIDIA DiffDock API
class DiffDockRequest(BaseModel):
    """Schema for requesting docking via NVIDIA DiffDock API."""
    protein: str = Field(..., description="Protein structure data in PDB format")
    ligand: str = Field(..., description="Ligand data as SMILES string or in SDF format")
    ligand_file_type: str = Field("smiles", description="Type of ligand data: 'smiles' or 'sdf'")
    num_poses: Optional[int] = Field(2, description="Number of binding poses to generate")
    time_divisions: Optional[int] = Field(20, description="Time divisions for the optimization process")
    steps: Optional[int] = Field(18, description="Number of denoising steps")
    save_trajectory: Optional[bool] = Field(False, description="Whether to save the trajectory")
    is_staged: Optional[bool] = Field(False, description="Whether to use staged sampling")

class DiffDockResponse(BaseModel):
    """Schema for the result of docking via NVIDIA DiffDock API."""
    ligand_positions: List[str] = Field(..., description="SDF content for each predicted ligand pose")
    position_confidence: List[float] = Field(..., description="Confidence score for each predicted pose (lower is better)")
    status: str = Field(..., description="Status of the docking operation")
    message: Optional[str] = Field(None, description="Error or additional message")

class VisualizationRequest(BaseModel):
    """Schema for requesting 3D visualization of docking results."""
    protein: str = Field(..., description="Protein structure data in PDB format")
    ligand_poses: List[str] = Field(..., description="SDF content for each predicted ligand pose")
    confidence_scores: List[float] = Field(..., description="Confidence score for each predicted pose")

class VisualizationResponse(BaseModel):
    """Schema for the result of visualization generation."""
    viewer_html: str = Field(..., description="HTML content for the 3D viewer")
    status: str = Field(..., description="Status of the visualization generation")
    message: Optional[str] = Field(None, description="Error or additional message")
    
class DockingResult(BaseModel):
    """Schema for the result of docking."""
    molecule_id: str = Field(..., description="ID of the molecule")
    protein_id: str = Field(..., description="ID of the protein")
    binding_affinity: float = Field(..., description="Binding affinity score (lower is better)")
    binding_pose: Dict[str, Any] = Field(..., description="3D coordinates of the binding pose")
    interactions: Dict[str, Any] = Field(..., description="Molecular interactions")
    
class BulkDockingResult(BaseModel):
    """Schema for the result of bulk docking."""
    total: int = Field(..., description="Total number of docking jobs")
    completed: int = Field(..., description="Number of completed docking jobs")
    results: List[DockingResult] = Field(..., description="Docking results for each molecule")
    protein_name: str = Field(..., description="Name of the protein")
    best_molecule_id: str = Field(..., description="ID of the molecule with the best binding affinity") 