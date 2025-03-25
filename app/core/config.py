import os
from pydantic import Field
from pydantic_settings import BaseSettings
from dotenv import load_dotenv

# Load environment variables from .env file if it exists
load_dotenv()

class Settings(BaseSettings):
    # API settings
    API_V1_STR: str = "/api/v1"
    PROJECT_NAME: str = "D3AI - Disease-based Drug Discovery"
    
    # Database settings
    MONGODB_URL: str = Field(default=os.getenv("MONGODB_URL", "mongodb://localhost:27017/"))
    DATABASE_NAME: str = Field(default=os.getenv("DATABASE_NAME", "d3ai"))
    
    # Service keys
    CHEMBL_API_BASE_URL: str = "https://www.ebi.ac.uk/chembl/api/data"
    PUBCHEM_API_BASE_URL: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    ZINC_API_BASE_URL: str = "http://zinc.docking.org/substances"
    
    # AI Models
    MEGAMOLBERT_MODEL_PATH: str = Field(default=os.getenv("MEGAMOLBERT_MODEL_PATH", "./models/megamolbert"))
    NVIDIA_DOCKING_API_KEY: str = Field(default=os.getenv("NVIDIA_DOCKING_API_KEY", ""))
    NVIDIA_DOCKING_API_URL: str = Field(default=os.getenv("NVIDIA_DOCKING_API_URL", ""))
    
    # Vector Database settings
    VECTOR_DB_URL: str = Field(default=os.getenv("VECTOR_DB_URL", "localhost:6333"))
    VECTOR_DB_COLLECTION: str = Field(default=os.getenv("VECTOR_DB_COLLECTION", "drug_papers"))
    
    # Frontend URL for CORS
    FRONTEND_URL: str = Field(default=os.getenv("FRONTEND_URL", "http://localhost:3000"))
    
    # Environment
    ENVIRONMENT: str = Field(default=os.getenv("ENVIRONMENT", "development"))
    
    # Logging
    LOG_LEVEL: str = Field(default=os.getenv("LOG_LEVEL", "INFO"))
    
    # Search settings
    MAX_SEARCH_RESULTS: int = 50
    PAPER_SEARCH_LIMIT: int = 100
    
    # Cache settings
    CACHE_TTL: int = 3600  # seconds
    
    # Molecular properties
    DEFAULT_BIOAVAILABILITY_THRESHOLD: float = 0.7
    DEFAULT_TOXICITY_THRESHOLD: float = 0.3
    
    class Config:
        case_sensitive = True
        env_file = ".env"

settings = Settings() 