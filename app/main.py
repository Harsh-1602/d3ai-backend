from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.routes import disease, drug_discovery, molecule, property, docking, protein
import os

app = FastAPI(
    title="D3AI - Disease-based Drug Discovery API",
    description="API for drug discovery based on disease identification and molecular analysis",
    version="1.0.0",
)

# Configure CORS
origins = [
    "http://localhost:3000",  # React app (local)
    "http://localhost:8000",  # FastAPI itself (local)
    "http://127.0.0.1:3000",
    "http://127.0.0.1:8000",
    os.getenv("FRONTEND_URL", ""),  # Deployed frontend URL from environment variable
    "https://d3ai-frontend.vercel.app",  # Add your actual frontend URL here
    # Add other production URLs as needed
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
)

# Include routers
app.include_router(disease.router)
app.include_router(drug_discovery.router)
app.include_router(molecule.router)
app.include_router(property.router)
app.include_router(docking.router)
app.include_router(protein.router)

@app.get("/")
async def root():
    return {
        "message": "Welcome to D3AI API",
        "docs": "/docs",
        "version": "1.0.0"
    } 