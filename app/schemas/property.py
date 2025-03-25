from typing import List, Optional, Dict, Any, Union
from pydantic import BaseModel, Field
from datetime import datetime

class PropertyBase(BaseModel):
    """Base schema for molecular property data."""
    molecule_id: str = Field(..., description="ID of the molecule")
    bioavailability: Optional[float] = Field(None, description="Bioavailability score (0-1)")
    toxicity: Optional[float] = Field(None, description="Toxicity score (0-1)")
    solubility: Optional[float] = Field(None, description="Solubility score (0-1)")
    druglikeness: Optional[float] = Field(None, description="Druglikeness score (0-1)")
    lipophilicity: Optional[float] = Field(None, description="Lipophilicity score (0-1)")
    half_life: Optional[float] = Field(None, description="Half-life in hours")
    clearance_rate: Optional[float] = Field(None, description="Clearance rate")
    additional_properties: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Additional properties")

class PropertyCreate(PropertyBase):
    """Schema for creating a new property record."""
    pass

class PropertyInDB(PropertyBase):
    """Schema for a property record that exists in the database."""
    id: str = Field(..., description="Unique identifier for the property record")
    created_at: datetime = Field(..., description="When the property record was first added to the database")
    updated_at: datetime = Field(..., description="When the property record was last updated")
    
    class Config:
        from_attributes = True

class Property(PropertyInDB):
    """Schema for property response."""
    molecule_smile: Optional[str] = Field(None, description="SMILE string of the molecule")
    molecule_name: Optional[str] = Field(None, description="Name of the molecule")
    
class PropertyPredictionRequest(BaseModel):
    """Schema for requesting property prediction."""
    smiles: Union[str, List[str]] = Field(..., description="SMILE string(s) to predict properties for")
    
class PropertyPredictionResult(BaseModel):
    """Schema for the result of property prediction."""
    smile: str = Field(..., description="SMILE string of the molecule")
    properties: Dict[str, Any] = Field(..., description="Predicted properties")
    
class BulkPropertyPredictionResult(BaseModel):
    """Schema for the result of bulk property prediction."""
    total: int = Field(..., description="Total number of molecules processed")
    results: List[PropertyPredictionResult] = Field(..., description="Prediction results for each molecule")
    
class PropertyFilterRequest(BaseModel):
    """Schema for filtering molecules by properties."""
    molecule_ids: List[str] = Field(..., description="List of molecule IDs to filter")
    bioavailability_min: Optional[float] = Field(None, description="Minimum bioavailability score (0-1)")
    toxicity_max: Optional[float] = Field(None, description="Maximum toxicity score (0-1)")
    solubility_min: Optional[float] = Field(None, description="Minimum solubility score (0-1)")
    druglikeness_min: Optional[float] = Field(None, description="Minimum druglikeness score (0-1)")
    custom_filters: Optional[Dict[str, Dict[str, float]]] = Field(default_factory=dict, description="Custom property filters")
    
class PropertyFilterResult(BaseModel):
    """Schema for the result of property filtering."""
    total: int = Field(..., description="Total number of molecules")
    filtered_count: int = Field(..., description="Number of molecules that passed the filter")
    filtered_molecules: List[str] = Field(..., description="IDs of molecules that passed the filter") 