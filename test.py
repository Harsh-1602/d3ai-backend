import requests
# URL encode the SMILES string
url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/C/PNG"

# Use GET request since we're retrieving an image
response = requests.get(url)

if response.status_code == 200:
    # Get the image data as bytes
    image_data = response.content
    
    # Convert to base64 for displaying/storing
    import base64
    image_base64 = base64.b64encode(image_data).decode('utf-8')
    
    print("Image data retrieved successfully")
    print(f"Image size: {len(image_data)} bytes")
    print(f"Base64 encoded image: {image_base64[:50]}...")
else:
    print(f"Failed to get image: {response.status_code}")
    print(response.text)