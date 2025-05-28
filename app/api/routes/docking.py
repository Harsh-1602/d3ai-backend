from fastapi import APIRouter, Depends, HTTPException, Body, Query, Path
from typing import List, Optional
import requests
import json
import os
from pathlib import Path as FilePath
from app.schemas.docking import (
    Protein, ProteinCreate, Docking, DockingCreate, DockingRequest, 
    DockingResult, BulkDockingResult, DiffDockRequest, DiffDockResponse,
    VisualizationRequest, VisualizationResponse
)
from app.services.docking_service import DockingService
from app.core.config import settings

router = APIRouter(prefix=f"{settings.API_V1_STR}/docking", tags=["Molecular Docking"])

@router.post("/diffdock", response_model=DiffDockResponse)
async def perform_diffdock(request: DiffDockRequest = Body(...)):
    """
    Perform molecular docking using NVIDIA's DiffDock API.
    
    This endpoint performs docking of a ligand with a protein target using deep learning.
    It predicts binding poses and confidence scores for each pose.
    """
    try:
        # Convert SMILES to SDF if ligand_file_type is 'smiles'
        ligand_data = request.ligand
        ligand_file_type = "sdf"  # Always use SDF for NVIDIA API
        
        if request.ligand_file_type.lower() == "smiles":
            try:
                # Convert SMILES to SDF using RDKit
                from rdkit import Chem
                from rdkit.Chem import AllChem
                
                mol = Chem.MolFromSmiles(request.ligand)
                if not mol:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Invalid SMILES string: {request.ligand}"
                    )
                
                # Add hydrogen atoms and generate 3D coordinates
                mol = Chem.AddHs(mol)
                success = AllChem.EmbedMolecule(mol, randomSeed=42)
                if success == -1:
                    raise HTTPException(
                        status_code=400,
                        detail="Failed to generate 3D coordinates for the molecule"
                    )
                
                # Optimize the structure
                AllChem.MMFFOptimizeMolecule(mol)
                
                # Convert to SDF
                sdf_stream = Chem.SDWriter(FilePath("temp.sdf"))
                sdf_stream.write(mol)
                sdf_stream.close()
                
                # Read the SDF file
                with open("temp.sdf", "r") as f:
                    ligand_data = f.read()
                
                # Remove temporary file
                os.remove("temp.sdf")
                
                print(f"Successfully converted SMILES to SDF: {request.ligand}")
            except ImportError:
                raise HTTPException(
                    status_code=500,
                    detail="RDKit is not installed. Cannot convert SMILES to SDF."
                )
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Error converting SMILES to SDF: {str(e)}"
                )
        
        # NVIDIA API call
        payload = {
            "ligand": ligand_data,
            "ligand_file_type": ligand_file_type,
            "protein": request.protein,
            "num_poses": request.num_poses,
            "time_divisions": request.time_divisions,
            "steps": request.steps,
            "save_trajectory": request.save_trajectory,
            "is_staged": request.is_staged
        }

        headers = {
            "accept": "application/json",
            "content-type": "application/json",
            "Authorization": f"Bearer {settings.NVIDIA_API_KEY}"
        }

        url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
        
        print(f"Sending DiffDock request with {len(payload['protein'])} bytes of protein data and {len(payload['ligand'])} bytes of ligand data")
        response = requests.post(url, json=payload, headers=headers)
        
        # Check if the request was successful
        if response.status_code != 200:
            error_detail = f"NVIDIA API Error: {response.text}"
            print(f"DiffDock API error: {error_detail}")
            raise HTTPException(
                status_code=response.status_code,
                detail=error_detail
            )
        
        result = response.json()
        
        return DiffDockResponse(
            ligand_positions=result.get("ligand_positions", []),
            position_confidence=result.get("position_confidence", []),
            status="success"
        )
    except HTTPException:
        raise
    except Exception as e:
        error_msg = f"Error performing DiffDock: {str(e)}"
        print(error_msg)
        raise HTTPException(status_code=500, detail=error_msg)


@router.post("/visualization", response_model=VisualizationResponse)
async def generate_visualization(request: VisualizationRequest = Body(...)):
    """
    Generate HTML visualization for protein-ligand docking.
    
    This endpoint takes protein data and ligand poses to generate an interactive 3D visualization.
    """
    try:
        # Path to the template
        template_path = FilePath("app/templates/viewer_template.html")
        
        # Create templates directory if it doesn't exist
        os.makedirs(template_path.parent, exist_ok=True)
        
        # Check if the template exists, if not, create it
        if not template_path.exists():
            with open(template_path, "w") as f:
                f.write(VIEWER_TEMPLATE)
        
        # Read the template file
        with open(template_path, "r") as f:
            template = f.read()
        
        # Create ligand data for the visualization
        ligand_data = []
        for i, (pose, confidence) in enumerate(zip(request.ligand_poses, request.confidence_scores)):
            ligand_data.append({
                "path": f"rank{i+1}_confidence_{confidence}.sdf",
                "name": f"rank{i+1}_confidence_{confidence}",
                "rank": i+1,
                "confidence": confidence,
                "content": pose
            })
        
        # Replace placeholders in the template
        html = template.replace(
            "const receptorStr = \"\";", 
            f"const receptorStr = `{request.protein}`;"
        ).replace(
            "const ligandData = [];",
            f"const ligandData = {json.dumps(ligand_data)};"
        )
        
        return VisualizationResponse(
            viewer_html=html,
            status="success"
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating visualization: {str(e)}")


@router.post("/dock", response_model=BulkDockingResult)
async def dock_molecules(request: DockingRequest = Body(...)):
    """
    Perform molecular docking.
    
    This endpoint performs docking of molecules with a protein target.
    It predicts binding affinity and binding poses for each molecule.
    """
    try:
        docking_service = DockingService()
        return await docking_service.dock_molecules(
            molecule_ids=request.molecule_ids,
            protein_id=request.protein_id,
            exhaustiveness=request.exhaustiveness,
            num_modes=request.num_modes
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error performing docking: {str(e)}")

@router.get("/proteins", response_model=List[Protein])
async def get_proteins(
    disease_id: Optional[str] = Query(None, description="Filter by disease ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get a list of proteins.
    
    This endpoint returns a list of proteins from the database.
    You can filter by disease ID or get all proteins with pagination.
    """
    try:
        docking_service = DockingService()
        return await docking_service.get_proteins(disease_id=disease_id, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching proteins: {str(e)}")

@router.get("/proteins/{protein_id}", response_model=Protein)
async def get_protein(protein_id: str = Path(..., description="Protein ID")):
    """
    Get a specific protein by ID.
    
    This endpoint returns detailed information about a specific protein.
    """
    try:
        docking_service = DockingService()
        protein = await docking_service.get_protein_by_id(protein_id)
        if not protein:
            raise HTTPException(status_code=404, detail=f"Protein with ID {protein_id} not found")
        return protein
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching protein: {str(e)}")

@router.post("/proteins", response_model=Protein)
async def create_protein(protein: ProteinCreate = Body(...)):
    """
    Add a new protein to the database.
    
    This endpoint allows adding a new protein to the database.
    """
    try:
        docking_service = DockingService()
        return await docking_service.create_protein(protein)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating protein: {str(e)}")

@router.get("/results", response_model=List[Docking])
async def get_docking_results(
    molecule_id: Optional[str] = Query(None, description="Filter by molecule ID"),
    protein_id: Optional[str] = Query(None, description="Filter by protein ID"),
    skip: int = 0,
    limit: int = 100
):
    """
    Get docking results.
    
    This endpoint returns docking results from the database.
    You can filter by molecule ID, protein ID, or get all results with pagination.
    """
    try:
        docking_service = DockingService()
        return await docking_service.get_docking_results(
            molecule_id=molecule_id, 
            protein_id=protein_id, 
            skip=skip, 
            limit=limit
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching docking results: {str(e)}")

@router.get("/results/{docking_id}", response_model=Docking)
async def get_docking_result(docking_id: str = Path(..., description="Docking ID")):
    """
    Get a specific docking result by ID.
    
    This endpoint returns detailed information about a specific docking result.
    """
    try:
        docking_service = DockingService()
        docking = await docking_service.get_docking_result_by_id(docking_id)
        if not docking:
            raise HTTPException(status_code=404, detail=f"Docking result with ID {docking_id} not found")
        return docking
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error fetching docking result: {str(e)}")

@router.post("/results", response_model=Docking)
async def create_docking_result(docking: DockingCreate = Body(...)):
    """
    Add a new docking result to the database.
    
    This endpoint allows adding a new docking result to the database.
    """
    try:
        docking_service = DockingService()
        return await docking_service.create_docking_result(docking)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating docking result: {str(e)}")

# HTML template for the 3D visualization
VIEWER_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
    <title>DiffDock Results Viewer</title>
    <script src="https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js"></script>
    <style>
        body { margin: 0; padding: 0; overflow: hidden; }
        #viewport { width: 100vw; height: 100vh; }
        #controls { 
            position: absolute; 
            top: 10px; 
            left: 10px; 
            background: rgba(255,255,255,0.8); 
            padding: 10px; 
            border-radius: 5px;
            max-height: 85vh;
            overflow-y: auto;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
            z-index: 100;
            min-width: 180px;
            max-width: 200px;
            font-size: 0.9rem;
            transition: transform 0.3s ease, opacity 0.3s ease;
        }
        #controls.minimized {
            transform: translateX(-95%);
            opacity: 0.6;
        }
        #controls.minimized:hover {
            opacity: 0.8;
        }
        #controls.minimized .control-panel,
        #controls.minimized .debug-info {
            display: none;
        }
        #controls.hidden {
            display: none;
        }
        .control-panel {
            margin-bottom: 10px;
        }
        h3 {
            margin-top: 0;
            margin-bottom: 8px;
            font-size: 1.1rem;
        }
        .header-container {
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        #toggle-controls {
            background: #2196F3;
            color: white;
            border: none;
            border-radius: 3px;
            width: 26px;
            height: 26px;
            padding: 2px;
            cursor: pointer;
            font-weight: bold;
            font-size: 14px;
            display: flex;
            align-items: center;
            justify-content: center;
        }
        #toggle-controls:hover {
            background: #0b7dda;
        }
        #close-controls {
            background: #f44336;
            color: white;
            border: none;
            border-radius: 3px;
            width: 26px;
            height: 26px;
            padding: 2px;
            cursor: pointer;
            font-weight: bold;
            font-size: 14px;
            display: flex;
            align-items: center;
            justify-content: center;
            margin-left: 5px;
        }
        #close-controls:hover {
            background: #d32f2f;
        }
        #reopen-button {
            position: absolute;
            top: 10px;
            left: 10px;
            background: #4CAF50;
            color: white;
            border: none;
            border-radius: 5px;
            padding: 8px 12px;
            cursor: pointer;
            font-weight: bold;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
            z-index: 99;
            display: none;
        }
        #reopen-button:hover {
            background: #45a049;
        }
        #toggle-icon {
            transition: transform 0.3s ease;
        }
        .minimized #toggle-icon {
            transform: rotate(180deg);
        }
        button {
            margin: 3px;
            padding: 6px 10px;
            cursor: pointer;
            border: none;
            background: #4CAF50;
            color: white;
            border-radius: 3px;
            font-weight: bold;
            font-size: 0.8rem;
        }
        button:hover {
            background: #45a049;
        }
        button.disabled {
            background: #cccccc;
            cursor: not-allowed;
        }
        .debug-info {
            font-family: monospace;
            font-size: 11px;
            color: #666;
            margin-top: 10px;
            border-top: 1px solid #ddd;
            padding-top: 8px;
        }
        #debug-log {
            max-height: 150px;
            overflow-y: auto;
            background: #f5f5f5;
            padding: 5px;
            border-radius: 3px;
            margin-top: 8px;
            font-size: 10px;
            white-space: pre-wrap;
        }
        .ligand-button {
            display: block;
            width: 100%;
            margin: 3px 0;
            padding: 6px;
            cursor: pointer;
            border: 1px solid #ddd;
            background: #f0f0f0;
            border-radius: 3px;
            text-align: left;
            font-size: 0.8rem;
        }
        .ligand-button:hover {
            background: #e0e0e0;
        }
        .ligand-button.active {
            background: #4CAF50;
            color: white;
        }
        .confidence {
            font-size: 0.7em;
            display: inline-block;
            margin-left: 3px;
        }
        .control-buttons {
            display: flex;
            gap: 3px;
            flex-wrap: wrap;
            justify-content: space-between;
        }
        .view-options button {
            flex: 1;
            min-width: 70px;
            background: #2196F3;
            font-size: 0.75rem;
            padding: 5px 8px;
        }
        .view-options button:hover {
            background: #0b7dda;
        }
        .minimized-hint {
            display: none;
            writing-mode: vertical-rl;
            text-orientation: mixed;
            transform: rotate(180deg);
            position: absolute;
            right: 5px;
            top: 50%;
            transform: translateY(-50%) rotate(180deg);
            color: #333;
            font-weight: bold;
        }
        #controls.minimized .minimized-hint {
            display: block;
        }
        .controls-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
        }
        .header-buttons {
            display: flex;
            align-items: center;
        }
    </style>
</head>
<body>
    <div id="viewport"></div>
    <button id="reopen-button">Show Controls</button>
    # <div id="controls">
    #     <div class="controls-header">
    #         <h3>DiffDock</h3>
    #         <div class="header-buttons">
    #             <button id="toggle-controls" title="Minimize/Expand">◀</button>
    #             <button id="close-controls" title="Close">✕</button>
    #         </div>
    #     </div>
    #     <span class="minimized-hint">DiffDock</span>
    #     <div class="control-panel">
    #         <div class="control-buttons">
    #             <button onclick="showAllLigands()" id="showAllBtn">Show All</button>
    #             <button onclick="hideAllLigands()" id="hideAllBtn">Hide All</button>
    #         </div>
            
    #         <div class="view-options control-buttons">
    #             <button onclick="resetView()">Reset</button>
    #             <button onclick="toggleSpin()" id="spinBtn">Spin</button>
    #             <button onclick="toggleProtein()" id="proteinBtn">Protein</button>
    #         </div>
    #     </div>
        
    #     <div class="control-panel">
    #         <h3>Poses</h3>
    #         <div id="ligand-buttons"></div>
    #     </div>
        
    #     <div class="debug-info">
    #         <div>Poses: <span id="pose-count">0</span></div>
    #         <div id="debug-log" style="display: none;"></div>
    #     </div>
    # </div>

    <script>
        // Initialize debugLog
        const debugLog = document.getElementById('debug-log');
        function log(message) {
            console.log(message);
            const logEntry = document.createElement('div');
            logEntry.textContent = message;
            debugLog.appendChild(logEntry);
            debugLog.scrollTop = debugLog.scrollHeight;
        }

        // Add toggle controls functionality
        document.getElementById('toggle-controls').addEventListener('click', function() {
            const controls = document.getElementById('controls');
            controls.classList.toggle('minimized');
            
            // Save state to localStorage
            const isMinimized = controls.classList.contains('minimized');
            localStorage.setItem('diffDockControlsMinimized', isMinimized);
            
            log("Controls panel " + (isMinimized ? "minimized" : "expanded"));
        });
        
        // Add close button functionality
        document.getElementById('close-controls').addEventListener('click', function() {
            const controls = document.getElementById('controls');
            controls.classList.add('hidden');
            document.getElementById('reopen-button').style.display = 'block';
            
            // Save state to localStorage
            localStorage.setItem('diffDockControlsHidden', 'true');
            
            log("Controls panel closed");
        });
        
        // Add reopen button functionality
        document.getElementById('reopen-button').addEventListener('click', function() {
            const controls = document.getElementById('controls');
            controls.classList.remove('hidden');
            this.style.display = 'none';
            
            // Save state to localStorage
            localStorage.setItem('diffDockControlsHidden', 'false');
            
            log("Controls panel reopened");
        });

        // Check saved state on load
        document.addEventListener('DOMContentLoaded', function() {
            const isMinimized = localStorage.getItem('diffDockControlsMinimized') === 'true';
            const isHidden = localStorage.getItem('diffDockControlsHidden') === 'true';
            
            if (isMinimized) {
                document.getElementById('controls').classList.add('minimized');
            }
            
            if (isHidden) {
                document.getElementById('controls').classList.add('hidden');
                document.getElementById('reopen-button').style.display = 'block';
            }
        });
        
        function displayDebugInfo() {
            const debugLog = document.getElementById('debug-log');
            if (debugLog.style.display === 'none') {
                debugLog.style.display = 'block';
            } else {
                debugLog.style.display = 'none';
            }
        }

        // Initialize NGL viewer
        const stage = new NGL.Stage("viewport");
        stage.setParameters({ backgroundColor: "white" });

        // Store loaded components
        let proteinComponent = null;
        const ligands = {};
        const activeLigands = new Set();
        let isSpinning = false;
        
        // Set receptor string - will be replaced by the server
        const receptorStr = "";
        
        // Set ligand data - will be replaced by the server
        const ligandData = [];
        
        // Update pose count
        document.getElementById('pose-count').textContent = ligandData.length;
        
        // Load receptor if available
        if (receptorStr.trim() !== "") {
            const receptorBlob = new Blob([receptorStr], {type: 'text/plain'});
            log("Creating receptor blob: " + receptorBlob.size + " bytes");
            
            stage.loadFile(receptorBlob, { ext: "pdb" })
                .then(function(component) {
                    proteinComponent = component;
                    
                    component.addRepresentation("cartoon", { 
                        color: "chainid",
                        opacity: 0.7
                    });
                    
                    component.addRepresentation("licorice", { 
                        sele: "hetero and not (water or ion)", 
                        multipleBond: true 
                    });
                    
                    log("Protein loaded successfully");
                    stage.autoView();
                })
                .catch(function(error) {
                    log("Error loading protein: " + error);
                });
        } else {
            log("No protein data available");
        }

        // Function to create color based on confidence
        function getColorForConfidence(confidence) {
            // Simple gradient from red (0.0) to green (1.0)
            const minConfidence = Math.min(...ligandData.map(l => l.confidence));
            const maxConfidence = Math.max(...ligandData.map(l => l.confidence));
            
            // Normalize confidence to 0-1 range
            const normalizedValue = (confidence - minConfidence) / (maxConfidence - minConfidence);
            
            // Create color string in hex format
            let r, g;
            if (normalizedValue < 0.5) {
                // Red to Yellow
                r = 255;
                g = Math.round(255 * 2.0 * normalizedValue);
                return "#" + r.toString(16).padStart(2, '0') + g.toString(16).padStart(2, '0') + "00";
            } else {
                // Yellow to Green
                r = Math.round(255 * (1.0 - 2.0 * (normalizedValue - 0.5)));
                g = 255;
                return "#" + r.toString(16).padStart(2, '0') + g.toString(16).padStart(2, '0') + "00";
            }
        }

        // Function to load a ligand
        function loadLigand(ligandInfo) {
            log("Loading ligand: " + ligandInfo.name);
            
            if (ligands[ligandInfo.name]) {
                // Already loaded, just make visible
                log("Ligand already loaded, making visible: " + ligandInfo.name);
                ligands[ligandInfo.name].setVisibility(true);
                activeLigands.add(ligandInfo.name);
                updateButtons();
                return Promise.resolve();
            }
            
            try {
                const content = ligandInfo.content;
                log("Ligand content length: " + content.length + " bytes");
                
                if (content.length < 10) {
                    log("ERROR: Ligand content too short for " + ligandInfo.name);
                    return Promise.reject("Content too short");
                }
                
                const ligandBlob = new Blob([content], {type: 'chemical/x-sdf'});
                log("Created ligand blob: " + ligandBlob.size + " bytes");
                
                return stage.loadFile(ligandBlob, { ext: "sdf", name: ligandInfo.name })
                    .then(function(component) {
                        // Get color based on confidence
                        const colorScheme = getColorForConfidence(ligandInfo.confidence);
                        log("Using color " + colorScheme + " for ligand " + ligandInfo.name);
                        
                        component.addRepresentation("licorice", { 
                            multipleBond: true,
                            radius: 0.3,
                            scale: 1.5,
                            colorValue: colorScheme
                        });
                        
                        // Store the component
                        ligands[ligandInfo.name] = component;
                        activeLigands.add(ligandInfo.name);
                        log("Successfully loaded ligand " + ligandInfo.name);
                        
                        updateButtons();
                        
                        // Center on ligand
                        if (Object.keys(ligands).length === 1) {
                            setTimeout(() => {
                                stage.autoView();
                                log("Centered on first ligand");
                            }, 100);
                        }
                        
                        return component;
                    })
                    .catch(function(error) {
                        log("ERROR loading ligand " + ligandInfo.name + ": " + error);
                        return Promise.reject(error);
                    });
            } catch (e) {
                log("ERROR creating ligand blob for " + ligandInfo.name + ": " + e);
                return Promise.reject(e);
            }
        }

        // Function to toggle a ligand
        function toggleLigand(ligandName) {
            log("Toggling ligand: " + ligandName);
            
            if (activeLigands.has(ligandName)) {
                // Hide this ligand
                if (ligands[ligandName]) {
                    ligands[ligandName].setVisibility(false);
                    log("Hide ligand: " + ligandName);
                }
                activeLigands.delete(ligandName);
            } else {
                // Show this ligand
                const ligandInfo = ligandData.find(l => l.name === ligandName);
                if (ligandInfo) {
                    log("Show ligand: " + ligandName);
                    loadLigand(ligandInfo);
                }
            }
            updateButtons();
        }

        // Show all ligands
        function showAllLigands() {
            log("Show all ligands");
            
            const promises = ligandData.map(ligandInfo => {
                return loadLigand(ligandInfo).catch(e => {
                    log("Error in showAllLigands: " + e);
                    return null;
                });
            });
            
            Promise.all(promises).then(() => {
                log("All ligands loaded");
                stage.autoView();
            });
        }

        // Hide all ligands
        function hideAllLigands() {
            log("Hide all ligands");
            
            Object.keys(ligands).forEach(ligandName => {
                if (ligands[ligandName]) {
                    ligands[ligandName].setVisibility(false);
                }
            });
            activeLigands.clear();
            updateButtons();
        }
        
        // Reset view
        function resetView() {
            log("Reset view");
            stage.autoView();
        }
        
        // Toggle spin
        function toggleSpin() {
            isSpinning = !isSpinning;
            stage.setSpin(isSpinning);
            document.getElementById('spinBtn').textContent = isSpinning ? 'Stop Spin' : 'Start Spin';
            log("Spin: " + isSpinning);
        }
        
        // Toggle protein visibility
        function toggleProtein() {
            if (!proteinComponent) return;
            
            const isVisible = proteinComponent.visible;
            proteinComponent.setVisibility(!isVisible);
            document.getElementById('proteinBtn').textContent = !isVisible ? 'Hide Protein' : 'Show Protein';
            log("Protein visibility: " + !isVisible);
        }

        // Update button states
        function updateButtons() {
            const buttons = document.querySelectorAll('.ligand-button');
            buttons.forEach(button => {
                const ligandName = button.dataset.ligand;
                if (activeLigands.has(ligandName)) {
                    button.classList.add('active');
                } else {
                    button.classList.remove('active');
                }
            });
        }

        // Create buttons for each ligand
        function createLigandButtons() {
            const container = document.getElementById('ligand-buttons');
            container.innerHTML = ''; // Clear any existing buttons
            
            ligandData.forEach(ligandInfo => {
                const button = document.createElement('button');
                button.className = 'ligand-button';
                button.dataset.ligand = ligandInfo.name;
                // Use normal string concatenation for JS template
                button.innerHTML = 'Rank ' + ligandInfo.rank + '<span class="confidence"> (Confidence: ' + ligandInfo.confidence.toFixed(2) + ')</span>';
                
                button.onclick = function() {
                    toggleLigand(ligandInfo.name);
                };
                
                container.appendChild(button);
            });
            
            log("Created " + ligandData.length + " ligand buttons");
        }

        // Initialize the buttons
        createLigandButtons();
        
        // Load the top-ranked ligand by default
        if (ligandData.length > 0) {
            log("Loading top-ranked ligand by default");
            setTimeout(() => {
                loadLigand(ligandData[0]).catch(error => {
                    log("Error loading default ligand: " + error);
                });
            }, 500);
        }
        
        // Add console output display for debugging
        window.onerror = function(message, source, lineno, colno, error) {
            log("ERROR: " + message + " at " + source + ":" + lineno);
            return false;
        };
        
        log("Visualization initialized");
    </script>
</body>
</html>
""" 