import os
from pymongo import MongoClient
from motor.motor_asyncio import AsyncIOMotorClient
from app.core.config import settings

# Global variables for database connections
_sync_client = None
_async_client = None

def get_sync_client() -> MongoClient:
    """
    Get a synchronous MongoDB client.
    
    Returns:
        MongoDB client
    """
    global _sync_client
    
    if _sync_client is None:
        mongodb_url = settings.MONGODB_URL
        _sync_client = MongoClient(mongodb_url)
        
    return _sync_client

def get_async_client() -> AsyncIOMotorClient:
    """
    Get an asynchronous MongoDB client.
    
    Returns:
        Async MongoDB client
    """
    global _async_client
    
    if _async_client is None:
        mongodb_url = settings.MONGODB_URL
        _async_client = AsyncIOMotorClient(mongodb_url)
        
    return _async_client

def get_database(db_name: str = None):
    """
    Get the MongoDB database.
    
    Args:
        db_name: Name of the database to use, defaults to DATABASE_NAME from settings
        
    Returns:
        MongoDB database
    """
    if db_name is None:
        db_name = settings.DATABASE_NAME
        
    client = get_async_client()
    return client[db_name]

def close_connections():
    """Close all database connections."""
    global _sync_client, _async_client
    
    if _sync_client is not None:
        _sync_client.close()
        _sync_client = None
        
    if _async_client is not None:
        _async_client.close()
        _async_client = None 