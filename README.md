# D3AI Backend - Disease-based Drug Discovery API

This is the backend API for the D3AI (Disease-based Drug Discovery AI) platform, which helps scientists screen drug candidates faster by leveraging AI and machine learning.

## Overview

The D3AI platform streamlines the drug discovery process by:

1. Identifying diseases from symptoms or direct input
2. Searching for potential drug candidates from research papers and databases
3. Generating new drug molecules from seed molecules
4. Validating molecular structures
5. Predicting molecular properties (bioavailability, toxicity, etc.)
6. Visualizing protein-drug interactions (docking)

## Technology Stack

- **FastAPI**: Modern, high-performance web framework for building APIs
- **MongoDB**: NoSQL database for storing disease, drug, and molecule data
- **RDKit**: Cheminformatics library for molecular operations
- **PyTorch & Transformers**: For AI models like MegaMolBERT
- **Autodock Vina**: For molecular docking simulations

## Project Structure

```
d3ai-backend/
├── app/
│   ├── api/
│   │   └── routes/         # API endpoints
│   ├── core/               # Core configuration
│   ├── db/                 # Database connections
│   ├── models/             # AI models
│   ├── schemas/            # Pydantic schemas
│   ├── services/           # Business logic
│   └── main.py             # Application entry point
├── tests/                  # Unit and integration tests
├── .env                    # Environment variables (not in repo)
└── requirements.txt        # Python dependencies
```

## API Endpoints

The API is organized into the following sections:

- `/api/v1/diseases`: Disease identification and management
- `/api/v1/drug-discovery`: Drug search and generation
- `/api/v1/molecules`: Molecule generation and validation
- `/api/v1/properties`: Molecular property prediction
- `/api/v1/docking`: Protein-ligand docking

## Getting Started

### Prerequisites

- Python 3.8+
- MongoDB
- RDKit

### Installation

1. Clone the repository
2. Create a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```
4. Create a `.env` file with the following variables:
   ```
   MONGODB_URL=mongodb://localhost:27017/
   DATABASE_NAME=d3ai
   MEGAMOLBERT_MODEL_PATH=./models/megamolbert
   ```

### Running the API

```
uvicorn app.main:app --reload
```

The API will be available at http://localhost:8000, and the interactive documentation at http://localhost:8000/docs.

## Development

### Adding New Features

1. Define schemas in `app/schemas/`
2. Implement business logic in `app/services/`
3. Create API endpoints in `app/api/routes/`
4. Update the main application in `app/main.py` if needed

### Testing

```
pytest
```

## License

This project is licensed under the MIT License - see the LICENSE file for details. 