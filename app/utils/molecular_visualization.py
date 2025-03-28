import os
import logging
import tempfile
from typing import Dict, Any, Optional, List, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import base64
from io import BytesIO

logger = logging.getLogger(__name__)

def smiles_to_2d_image(smile: str, width: int = 300, height: int = 200, 
                      highlight_atoms: Optional[List[int]] = None) -> Optional[str]:
    """
    Convert a SMILES string to a 2D image in base64 format.
    
    Args:
        smile: SMILES string representation of a molecule
        width: Width of the image in pixels
        height: Height of the image in pixels
        highlight_atoms: Optional list of atom indices to highlight
        
    Returns:
        Base64 encoded PNG image as a string, or None if conversion fails
    """
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            logger.error(f"Failed to convert SMILES to molecule: {smile}")
            return None
            
        # Compute 2D coordinates
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        mol = Chem.RemoveHs(mol)
        
        # Create molecular drawing
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().additionalAtomLabelPadding = 0.3
        
        # Add highlight if provided
        if highlight_atoms:
            highlight_bond_indices = []
            colors = {}
            for idx in highlight_atoms:
                colors[idx] = (1, 0, 0)  # Red color for highlights
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, 
                               highlightBonds=highlight_bond_indices,
                               highlightAtomColors=colors)
        else:
            drawer.DrawMolecule(mol)
        
        drawer.FinishDrawing()
        
        # Get PNG data and convert to base64
        png_data = drawer.GetDrawingText()
        
        # Return just the base64 encoded string without the prefix
        # The frontend will add the prefix as needed
        return base64.b64encode(png_data).decode('utf-8')
    
    except Exception as e:
        logger.error(f"Error creating 2D visualization for SMILES {smile}: {e}")
        return None

def smiles_to_3d_coordinates(smile: str) -> Optional[Dict[str, Any]]:
    """
    Convert a SMILES string to 3D coordinates for visualization.
    
    Args:
        smile: SMILES string representation of the molecule
        
    Returns:
        Dictionary with atom positions, bonds, and atom types, or None if conversion fails
    """
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            logger.error(f"Failed to convert SMILES to molecule: {smile}")
            return None
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Optimize the structure
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Get atom positions
        conf = mol.GetConformer()
        atom_positions = []
        atom_types = []
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atom_positions.append({
                "x": pos.x,
                "y": pos.y,
                "z": pos.z
            })
            atom_types.append(atom.GetSymbol())
        
        # Get bonds
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                "source": bond.GetBeginAtomIdx(),
                "target": bond.GetEndAtomIdx(),
                "order": int(bond.GetBondTypeAsDouble())
            })
        
        return {
            "atoms": atom_positions,
            "atom_types": atom_types,
            "bonds": bonds
        }
        
    except Exception as e:
        logger.error(f"Error creating 3D coordinates for SMILES {smile}: {e}")
        return None

def smiles_to_molecule_properties(smile: str) -> Dict[str, Any]:
    """
    Get basic molecular properties from a SMILES string.
    
    Args:
        smile: SMILES string representation of the molecule
        
    Returns:
        Dictionary with molecular properties
    """
    properties = {
        "formula": "",
        "molecular_weight": 0,
        "num_atoms": 0,
        "num_bonds": 0,
        "num_rings": 0,
        "is_valid": False
    }
    
    try:
        # Convert SMILES to molecule
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            return properties
            
        properties["is_valid"] = True
        properties["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
        properties["molecular_weight"] = Chem.Descriptors.MolWt(mol)
        properties["num_atoms"] = mol.GetNumAtoms()
        properties["num_bonds"] = mol.GetNumBonds()
        properties["num_rings"] = Chem.rdMolDescriptors.CalcNumRings(mol)
        
        # Additional descriptors
        properties["logp"] = Chem.Descriptors.MolLogP(mol)
        properties["tpsa"] = Chem.Descriptors.TPSA(mol)
        properties["hba"] = Chem.Descriptors.NumHAcceptors(mol)
        properties["hbd"] = Chem.Descriptors.NumHDonors(mol)
        properties["rotatable_bonds"] = Chem.Descriptors.NumRotatableBonds(mol)
        
        return properties
        
    except Exception as e:
        logger.error(f"Error calculating properties for SMILES {smile}: {e}")
        return properties

def compare_molecules_similarity(smile1: str, smile2: str) -> float:
    """
    Calculate similarity between two molecules based on their SMILES representations.
    
    Args:
        smile1: First SMILES string
        smile2: Second SMILES string
        
    Returns:
        Tanimoto similarity score (0-1)
    """
    try:
        mol1 = Chem.MolFromSmiles(smile1)
        mol2 = Chem.MolFromSmiles(smile2)
        
        if mol1 is None or mol2 is None:
            return 0.0
            
        # Calculate Morgan fingerprints
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        
        # Calculate Tanimoto similarity
        similarity = Chem.DataStructs.TanimotoSimilarity(fp1, fp2)
        
        return similarity
        
    except Exception as e:
        logger.error(f"Error calculating similarity between {smile1} and {smile2}: {e}")
        return 0.0

def create_molecule_grid_image(smiles_list: List[str], 
                              legends: Optional[List[str]] = None,
                              mols_per_row: int = 4,
                              sub_img_size: Tuple[int, int] = (200, 200)) -> Optional[str]:
    """
    Create a grid image of multiple molecules from their SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        legends: Optional list of legends for each molecule
        mols_per_row: Number of molecules per row in the grid
        sub_img_size: Size of each molecular image in the grid
        
    Returns:
        Base64 encoded PNG image as a string, or None if conversion fails
    """
    try:
        # Convert SMILES to molecules
        mols = []
        for smile in smiles_list:
            mol = Chem.MolFromSmiles(smile)
            if mol is not None:
                # Add 2D coordinates
                AllChem.Compute2DCoords(mol)
                mols.append(mol)
        
        if not mols:
            logger.error("No valid molecules to create grid image")
            return None
            
        # Create legends if not provided
        if legends is None:
            legends = [f"Molecule {i+1}" for i in range(len(mols))]
        
        # Ensure legends match the number of valid molecules
        legends = legends[:len(mols)]
        
        # Create grid image
        img = Draw.MolsToGridImage(mols, molsPerRow=mols_per_row, subImgSize=sub_img_size,
                                 legends=legends, useSVG=False)
        
        # Convert to base64
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        
        # Return just the base64 encoded string without the prefix
        # The frontend will add the prefix as needed
        return base64.b64encode(buffered.getvalue()).decode('utf-8')
        
    except Exception as e:
        logger.error(f"Error creating grid image: {e}")
        return None 