from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime

class DrugBase(BaseModel):
    """Base schema for drug data."""
    name: str = Field(..., description="Generic name of the drug")
    smile: str = Field(..., description="SMILE string representation of the drug molecule")
    disease_id: str = Field(..., description="ID of the disease this drug is associated with")
    description: Optional[str] = Field(None, description="Description of the drug")
    brand_names: Optional[List[str]] = Field(default_factory=list, description="Commercial/brand names of the drug")
    drug_class: Optional[str] = Field(None, description="Class of the drug (e.g., antibiotic, antiviral)")
    mechanism_of_action: Optional[str] = Field(None, description="How the drug works")
    research_paper_url: Optional[str] = Field(None, description="URL to the research paper describing the drug")

class DrugCreate(DrugBase):
    """Schema for creating a new drug."""
    pass

class DrugInDB(DrugBase):
    """Schema for a drug that exists in the database."""
    id: str = Field(..., description="Unique identifier for the drug")
    created_at: datetime = Field(..., description="When the drug was first added to the database")
    updated_at: datetime = Field(..., description="When the drug was last updated")
    
    class Config:
        from_attributes = True

class Drug(DrugInDB):
    """Schema for drug response that includes associated molecules."""
    disease_name: Optional[str] = Field(None, description="Name of the disease this drug treats")
    properties: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Key properties of the drug")
    similar_drugs: Optional[List[str]] = Field(default_factory=list, description="Similar drugs")
    
class DrugSearchResult(BaseModel):
    """Schema for the result of drug search."""
    query: str = Field(..., description="The disease or keyword that was searched")
    total_results: int = Field(..., description="Total number of results found")
    drugs: List[Drug] = Field(..., description="List of drugs found")
    search_source: str = Field(..., description="Source of the search results (e.g., ChEMBL, PubChem, research papers)")
    
class DrugGenerationRequest(BaseModel):
    """Schema for requesting drug generation."""
    disease_id: str = Field(..., description="ID of the disease to generate drugs for")
    seed_drug_ids: Optional[List[str]] = Field(default_factory=list, description="List of drug IDs to use as seeds")
    count: int = Field(10, description="Number of drug candidates to generate")
    diversity_factor: float = Field(0.7, description="Diversity factor for generation (0-1)")
    
class DrugGenerationResult(BaseModel):
    """Schema for the result of drug generation."""
    disease_id: str = Field(..., description="ID of the disease")
    disease_name: str = Field(..., description="Name of the disease")
    seed_drugs: List[str] = Field(..., description="Names of seed drugs used")
    generated_molecules: List[Dict[str, Any]] = Field(..., description="Generated drug molecules") 