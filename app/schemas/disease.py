from typing import List, Optional
from pydantic import BaseModel, Field
from datetime import datetime

class SymptomInput(BaseModel):
    """Schema for symptom input."""
    symptoms: List[str] = Field(..., description="List of symptoms experienced by the patient")

class DiseaseBase(BaseModel):
    """Base schema for disease data."""
    name: str = Field(..., description="Name of the disease")
    description: Optional[str] = Field(None, description="Description of the disease")
    symptoms: List[str] = Field(default_factory=list, description="Common symptoms of the disease")
    causes: Optional[List[str]] = Field(default_factory=list, description="Causes of the disease")
    icd_code: Optional[str] = Field(None, description="ICD-10 code for the disease")
    category: Optional[str] = Field(None, description="Category of the disease (e.g., infectious, genetic, etc.)")

class DiseaseCreate(DiseaseBase):
    """Schema for creating a new disease."""
    pass

class DiseaseInDB(DiseaseBase):
    """Schema for a disease that exists in the database."""
    id: str = Field(..., description="Unique identifier for the disease")
    created_at: datetime = Field(..., description="When the disease was first added to the database")
    updated_at: datetime = Field(..., description="When the disease was last updated")
    
    class Config:
        from_attributes = True

class Disease(DiseaseInDB):
    """Schema for disease response that includes associated drugs."""
    associated_drugs: Optional[List[str]] = Field(default_factory=list, description="Drug names known to treat this disease")
    research_papers_count: Optional[int] = Field(0, description="Number of research papers related to this disease")

class DiseaseIdentificationResult(BaseModel):
    """Schema for the result of disease identification from symptoms."""
    disease: Disease
    confidence: float = Field(..., description="Confidence score of the disease identification (0-1)")
    alternative_diseases: List[Disease] = Field(default_factory=list, description="Alternative possible diseases based on symptoms")
    
class DiseasePrediction(BaseModel):
    """Schema for disease prediction from symptoms."""
    disease_name: str = Field(..., description="Predicted disease name")
    confidence: float = Field(..., description="Confidence score (0-1)")
    symptoms_matched: List[str] = Field(..., description="Symptoms that matched the disease")
    symptoms_not_matched: List[str] = Field(..., description="Symptoms that didn't match the disease") 