import uvicorn
import os
from dotenv import load_dotenv
from app.utils.logging_config import setup_logging

# Load environment variables
load_dotenv()

if __name__ == "__main__":
    # Set up logging
    log_level = os.getenv("LOG_LEVEL", "INFO")
    setup_logging(log_level)
    
    # Get port from environment variable or use default
    port = int(os.getenv("PORT", 8000))
    
    # Run the application
    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=port,
        reload=True,
        log_level=log_level.lower()
    ) 