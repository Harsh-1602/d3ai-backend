"""
Database initialization script for D3AI.

This script creates the necessary collections in MongoDB and seeds them with initial data.
"""

import os
import asyncio
import sys
from datetime import datetime
import uuid
from pymongo import MongoClient
from dotenv import load_dotenv

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load environment variables
load_dotenv()

# MongoDB connection details
MONGODB_URL = os.getenv("MONGODB_URL", "mongodb://localhost:27017/")
DATABASE_NAME = os.getenv("DATABASE_NAME", "d3ai")

# Initial data
DISEASES = [
    {
        "_id": str(uuid.uuid4()),
        "name": "Diabetes Mellitus",
        "description": "A metabolic disorder characterized by high blood sugar levels over a prolonged period.",
        "symptoms": ["excessive thirst", "frequent urination", "increased hunger", "fatigue", "blurred vision", "slow healing sores"],
        "causes": ["insulin resistance", "insufficient insulin production", "genetic factors", "obesity"],
        "icd_code": "E11",
        "category": "Metabolic",
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    },
    {
        "_id": str(uuid.uuid4()),
        "name": "Hypertension",
        "description": "A condition in which the force of the blood against the artery walls is too high.",
        "symptoms": ["headaches", "shortness of breath", "nosebleeds", "flushing", "dizziness", "chest pain"],
        "causes": ["obesity", "excessive salt intake", "stress", "alcohol consumption", "genetic factors"],
        "icd_code": "I10",
        "category": "Cardiovascular",
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    },
    {
        "_id": str(uuid.uuid4()),
        "name": "Alzheimer's Disease",
        "description": "A progressive neurologic disorder that causes brain cells to die and the brain to shrink (atrophy).",
        "symptoms": ["memory loss", "confusion", "difficulty speaking", "anxiety", "mood swings", "apathy"],
        "causes": ["genetic mutations", "age", "family history", "head injuries", "cardiovascular disease"],
        "icd_code": "G30",
        "category": "Neurological",
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    }
]

DRUGS = [
    {
        "_id": str(uuid.uuid4()),
        "name": "Metformin",
        "smile": "CN(C)C(=N)NC(=N)N",
        "disease_id": None,  # Will be set to Diabetes ID
        "description": "Metformin is a biguanide used to treat type 2 diabetes mellitus.",
        "brand_names": ["Glucophage", "Fortamet", "Glumetza", "Riomet"],
        "drug_class": "Biguanide",
        "mechanism_of_action": "Decreases hepatic glucose production and intestinal glucose absorption",
        "research_paper_url": "https://pubmed.ncbi.nlm.nih.gov/12186605/",
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    },
    {
        "_id": str(uuid.uuid4()),
        "name": "Lisinopril",
        "smile": "CCOC(=O)C(CCc1ccccc1)NC(C(=O)N1CCCC1C(=O)O)C",
        "disease_id": None,  # Will be set to Hypertension ID
        "description": "Lisinopril is an ACE inhibitor used to treat hypertension and heart failure.",
        "brand_names": ["Prinivil", "Zestril"],
        "drug_class": "ACE Inhibitor",
        "mechanism_of_action": "Inhibits angiotensin-converting enzyme, which decreases the production of angiotensin II",
        "research_paper_url": "https://pubmed.ncbi.nlm.nih.gov/2568760/",
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    }
]

PROTEINS = [
    {
        "_id": str(uuid.uuid4()),
        "name": "Insulin receptor",
        "pdb_id": "P06213",  # This is actually the UniProt ID
        "sequence": "MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNYIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICKSHGCTAEGLCCHSECLGNCSQPDDPTKCVACRNFYLDGRCVETCPPPYYHFQDWRCVNFSFCQDLHHKCKNSRRQGCHQYVIHNNKCIPECPSGYTMNSSNLLCTPCLGPCPKVCHLLEGEKTIDSVTSAQELRGCTVINGSLIINIRGGNNLAAELEANLGLIEEISGYLKIRRSYALVSLSFFRKLRLIRGETLEIGNYSFYALDNQNLRQLWDWSKHNLTITQGKLFFHYNPKLCLSEIHKMEEVSGTKGRQERNDIALKTNGDQASCENELLKFSYIRTSFDKILLRWEPYWPPDFRDLLGFMLFYKEAPYQNVTEFDGQDACGSNSWTVVDIDPPLRSNDPKSQNHPGWLMRGLKPWTQYAIFVKTLVTFSDERRTYGAKSDIIYVQTDATNPSVPLDPISVSNSSSQIILKWKPPSDPNGNITHYLVFWERQAEDSELFELDYCLKGLKLPSRTWSPPFESEDSQKHNQSEYEDSAGECCSCPKTDSQILKELEESSFRKTFEDYLHNVVFVPRKTSSGTGAEDPRPSRKRRSLGDVGNVTVAVPTVAAFPNTSSTSVPTSPEEHRPFEKVVNKESLVISGLRHFTGYRIELQACNQDTPEERCSVAAYVSARTMPEAKADDIVGPVTHEIFENNVVHLMWQEPKEPNGLIVLYEVSYRRYGDEELHLCVSRKHFALERGCRLRGLSPGNYSVRIRATSLAGNGSWTEPTYFYVTDYLDVPSNIAKIIIGPLIFVFLFSVVIGSIYLFLRKRQPDGPLGPLYASSNPEYLSASDVFPCSVYVPDEWEVSREKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNESASLRERIEFLNEASVMKGFTCHHVVRLLGVVSKGQPTLVVMELMAHGDLKSYLRSLRPEAENNPGRPPPTLQEMIQMAAEIADGMAYLNAKKFVHRDLAARNCMVAHDFTVKIGDFGMTRDIYETDYYRKGGKGLLPVRWMAPESLKDGVFTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQVLKFVMDGGYLDQPDNCPERVTDLMRMCWQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFFHSEENKAPESEELEMEFEDMENVPLDRSSHCQREEAGGRDGGSSLGFKRSYEEHIPYTHMNGGKKNGRILTLPRSNPS",
        "disease_id": None,  # Will be set to Diabetes ID
        "description": "The insulin receptor is a transmembrane receptor that is activated by insulin and IGF-I, IGF-II.",
        "binding_site": {
            "residues": ["Phe714", "Val715", "Asp716", "Arg717", "Leu718"],
            "coordinates": {"x": 10.5, "y": 12.3, "z": 15.7}
        },
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    },
    {
        "_id": str(uuid.uuid4()),
        "name": "Angiotensin-converting enzyme",
        "pdb_id": "P12821",  # This is actually the UniProt ID
        "sequence": "MGAASGRRGPGLLLPLPLLLLLPPQPALALDPGLQPGNFSADEAGAQLFAQSYNSSAEQVLFQSVAASWAHDTNITAENARRQEEAALLSQEFAEAWGQKAKELYEPIWQNFTDPQLRRIIGAVRTLGSANLPLAKRQQYNALLSNMSRIYSTAKVCLPNKTATCWSLDPDLTNILASSRSYAMLLFAWEGWHNAAGIPLKPLYEDFTALSNEAYKQDGFTDTGAYWRSWYNSPTEFDNSYFVQYMDRIGVLALPEFSNSYLRPRPL",
        "disease_id": None,  # Will be set to Hypertension ID
        "description": "The angiotensin-converting enzyme is an enzyme involved in the synthesis of angiotensin II, a potent vasoconstrictor.",
        "binding_site": {
            "residues": ["His383", "His387", "Glu411", "His523", "Tyr520"],
            "coordinates": {"x": 15.3, "y": 22.1, "z": 8.4}
        },
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    }
]

MOLECULES = [
    {
        "_id": str(uuid.uuid4()),
        "smile": "CN(C)C(=N)NC(=N)N",
        "parent_drug_id": None,  # Will be set to Metformin ID
        "name": "Metformin",
        "description": "Biguanide antihyperglycemic agent",
        "generation_method": "Manual",
        "is_valid": True,
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    },
    {
        "_id": str(uuid.uuid4()),
        "smile": "CCOC(=O)C(CCc1ccccc1)NC(C(=O)N1CCCC1C(=O)O)C",
        "parent_drug_id": None,  # Will be set to Lisinopril ID
        "name": "Lisinopril",
        "description": "ACE inhibitor",
        "generation_method": "Manual",
        "is_valid": True,
        "created_at": datetime.now(),
        "updated_at": datetime.now()
    }
]

async def init_database():
    """Initialize the MongoDB database with collections and sample data."""
    print(f"Connecting to MongoDB at {MONGODB_URL}")
    client = MongoClient(MONGODB_URL)
    db = client[DATABASE_NAME]
    
    # Create collections if they don't exist
    collections = ["diseases", "drugs", "molecules", "proteins", "properties", "docking_results"]
    existing_collections = db.list_collection_names()
    
    for collection in collections:
        if collection not in existing_collections:
            db.create_collection(collection)
            print(f"Created collection: {collection}")
    
    # Insert initial data: Diseases
    print("Inserting sample diseases...")
    diseases_collection = db["diseases"]
    if diseases_collection.count_documents({}) == 0:
        result = diseases_collection.insert_many(DISEASES)
        print(f"Inserted {len(result.inserted_ids)} diseases")
    else:
        print("Diseases collection already has data. Skipping...")

    # Update disease_id references in other collections
    diabetes_id = diseases_collection.find_one({"name": "Diabetes Mellitus"})["_id"]
    hypertension_id = diseases_collection.find_one({"name": "Hypertension"})["_id"]
    
    # Update Drugs
    DRUGS[0]["disease_id"] = diabetes_id  # Metformin -> Diabetes
    DRUGS[1]["disease_id"] = hypertension_id  # Lisinopril -> Hypertension
    
    # Update Proteins
    PROTEINS[0]["disease_id"] = diabetes_id  # Insulin receptor -> Diabetes
    PROTEINS[1]["disease_id"] = hypertension_id  # ACE -> Hypertension
    
    # Insert Drugs
    print("Inserting sample drugs...")
    drugs_collection = db["drugs"]
    if drugs_collection.count_documents({}) == 0:
        result = drugs_collection.insert_many(DRUGS)
        print(f"Inserted {len(result.inserted_ids)} drugs")
    else:
        print("Drugs collection already has data. Skipping...")
    
    # Update drug_id references
    metformin_id = drugs_collection.find_one({"name": "Metformin"})["_id"]
    lisinopril_id = drugs_collection.find_one({"name": "Lisinopril"})["_id"]
    
    # Update Molecules
    MOLECULES[0]["parent_drug_id"] = metformin_id
    MOLECULES[1]["parent_drug_id"] = lisinopril_id
    
    # Insert Proteins
    print("Inserting sample proteins...")
    proteins_collection = db["proteins"]
    if proteins_collection.count_documents({}) == 0:
        result = proteins_collection.insert_many(PROTEINS)
        print(f"Inserted {len(result.inserted_ids)} proteins")
    else:
        print("Proteins collection already has data. Skipping...")
    
    # Insert Molecules
    print("Inserting sample molecules...")
    molecules_collection = db["molecules"]
    if molecules_collection.count_documents({}) == 0:
        result = molecules_collection.insert_many(MOLECULES)
        print(f"Inserted {len(result.inserted_ids)} molecules")
    else:
        print("Molecules collection already has data. Skipping...")
    
    print("Database initialization complete!")

if __name__ == "__main__":
    asyncio.run(init_database())