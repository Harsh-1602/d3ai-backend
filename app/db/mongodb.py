import logging
from pymongo import MongoClient
from pymongo.database import Database
from app.core.config import settings
import time

logger = logging.getLogger(__name__)

class MongoDB:
    """MongoDB database connection manager."""
    
    client: MongoClient = None
    db: Database = None
    
    @classmethod
    def connect_to_mongodb(cls):
        """Connect to MongoDB database with retry logic."""
        if cls.client is None:
            max_retries = 3
            retries = 0
            while retries < max_retries:
                try:
                    logger.info(f"Connecting to MongoDB at {settings.MONGODB_URL}")
                    cls.client = MongoClient(
                        settings.MONGODB_URL,
                        serverSelectionTimeoutMS=5000,
                        connectTimeoutMS=10000
                    )
                    # Force a command to test the connection
                    cls.client.admin.command('ping')
                    cls.db = cls.client[settings.DATABASE_NAME]
                    logger.info(f"Connected to MongoDB at {settings.MONGODB_URL}")
                    
                    # Create collections if they don't exist
                    cls._create_collections()
                    
                    return cls.db
                except Exception as e:
                    retries += 1
                    logger.error(f"Failed to connect to MongoDB (attempt {retries}/{max_retries}): {e}")
                    if retries >= max_retries:
                        raise
                    time.sleep(1)  # Wait before retrying
        return cls.db
    
    @classmethod
    def close_mongodb_connection(cls):
        """Close MongoDB connection."""
        if cls.client is not None:
            cls.client.close()
            cls.client = None
            cls.db = None
            logger.info("Closed connection to MongoDB")
    
    @classmethod
    def _create_collections(cls):
        """Create necessary collections if they don't exist."""
        collections = ["diseases", "drugs", "molecules", "properties", "docking_results"]
        existing_collections = cls.db.list_collection_names()
        
        for collection in collections:
            if collection not in existing_collections:
                cls.db.create_collection(collection)
                logger.info(f"Created collection: {collection}")
                
                # Create indexes for each collection
                if collection == "diseases":
                    cls.db[collection].create_index("name", unique=True)
                    cls.db[collection].create_index("symptoms")
                elif collection == "drugs":
                    cls.db[collection].create_index("name", unique=True)
                    cls.db[collection].create_index("disease_id")
                elif collection == "molecules":
                    cls.db[collection].create_index("smile", unique=True)
                    cls.db[collection].create_index("parent_drug_id")
                elif collection == "properties":
                    cls.db[collection].create_index("molecule_id", unique=True)
                elif collection == "docking_results":
                    cls.db[collection].create_index([("molecule_id", 1), ("protein_id", 1)], unique=True)

# Global instance
try:
    db = MongoDB.connect_to_mongodb()
except Exception as e:
    logger.error(f"Failed to connect to MongoDB at startup: {e}")
    db = None  # Will be initialized on first request 