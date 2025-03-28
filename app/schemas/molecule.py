from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime

class MoleculeBase(BaseModel):
    """Base schema for molecule data."""
    smile: str = Field(..., description="SMILE string representation of the molecule")
    parent_drug_id: Optional[str] = Field(None, description="ID of the parent drug (if generated)")
    name: Optional[str] = Field(None, description="Name of the molecule (if available)")
    description: Optional[str] = Field(None, description="Description of the molecule")
    generation_method: Optional[str] = Field(None, description="Method used to generate the molecule")
    is_valid: bool = Field(True, description="Whether the molecule is structurally valid")

class MoleculeCreate(MoleculeBase):
    """Schema for creating a new molecule."""
    pass

class MoleculeInDB(MoleculeBase):
    """Schema for a molecule that exists in the database."""
    id: str = Field(..., description="Unique identifier for the molecule")
    created_at: datetime = Field(..., description="When the molecule was first added to the database")
    updated_at: datetime = Field(..., description="When the molecule was last updated")
    
    class Config:
        from_attributes = True

class Molecule(MoleculeInDB):
    """Schema for molecule response that includes properties."""
    properties: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Properties of the molecule")
    visualization_data: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Data for 3D visualization")
    similarity_to_parent: Optional[float] = Field(None, description="Similarity to parent drug (0-1)")
    
class MoleculeGenerationRequest(BaseModel):
    """Schema for requesting molecule generation."""
    seed_smiles: List[str] = Field(..., description="List of SMILE strings to use as seeds")
    count: int = Field(10, description="Number of molecules to generate")
    diversity_factor: float = Field(0.7, description="Diversity factor for generation (0-1)")
    similarity_threshold: float = Field(0.6, description="Minimum similarity to parent molecules (0-1)")
    
class MoleculeValidationResult(BaseModel):
    """Schema for the result of molecule validation."""
    smile: str = Field(..., description="SMILE string of the molecule")
    is_valid: bool = Field(..., description="Whether the molecule is structurally valid")
    validation_errors: Optional[List[str]] = Field(default_factory=list, description="List of validation errors")
    visualization_data: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Data for 2D visualization")
    
class BulkMoleculeValidationRequest(BaseModel):
    """Schema for requesting bulk molecule validation."""
    smiles: List[str] = Field(..., description="List of SMILE strings to validate")
    
class BulkMoleculeValidationResult(BaseModel):
    """Schema for the result of bulk molecule validation."""
    total: int = Field(..., description="Total number of molecules processed")
    valid_count: int = Field(..., description="Number of valid molecules")
    results: List[MoleculeValidationResult] = Field(..., description="Validation results for each molecule")

class NvidiaGenMolRequest(BaseModel):
    """Schema for requesting molecule generation from NVIDIA GenMol API."""
    smiles: str = Field(..., description="SMILES string to use as seed")
    num_molecules: int = Field(30, description="Number of molecules to generate")
    temperature: str = Field("1", description="Temperature parameter for generation")
    noise: str = Field("1", description="Noise parameter for generation")
    step_size: int = Field(1, description="Step size parameter")
    scoring: str = Field("QED", description="Scoring function to use (QED, logP, etc.)")
    unique: bool = Field(False, description="Whether to generate unique molecules") 