#!/usr/bin/env python3
"""
Script 1: Protein Target Preparation for Molecular Docking

Purpose: Download, prepare, and validate protein structures for docking.
Converts PDB to PDBQT format, defines binding sites, generates docking grids.

Usage:
    python 01_prepare_protein.py --pdb_id 3ODU --binding_residues 110,114,121,125
"""

import argparse
import numpy as np
from prody import *
import os

def prepare_protein(pdb_id, binding_site_residues=None, ph=7.4, output_dir='./'):
    """
    Download and prepare protein structure for docking.
    
    Args:
        pdb_id: PDB identifier (e.g., '3ODU')
        binding_site_residues: List of residue numbers defining pocket (None for blind docking)
        ph: Physiological pH (default 7.4)
        output_dir: Output directory for prepared structures
    
    Returns:
        Prepared PDB filename, box_center, box_size
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Downloading PDB structure {pdb_id}...")
    pdb = fetchPDB(pdb_id)
    ag = parsePDB(pdb)
    
    print(f"Loaded {ag.numAtoms()} atoms, {ag.numResidues()} residues")
    
    # Remove heteroatoms (water, ligands) except essential cofactors
    ag_cleaned = ag.select('not heteroatm or (resname ZN or resname MG)')
    
    # Save prepared structure
    output_pdb = os.path.join(output_dir, f'{pdb_id}_prepared.pdb')
    writePDB(output_pdb, ag_cleaned)
    print(f"Saved prepared structure: {output_pdb}")
    
    # Define docking grid
    if binding_site_residues:
        # Targeted docking: define pocket
        residue_str = ' or '.join([f'residue {r}' for r in binding_site_residues])
        pocket_atoms = ag_cleaned.select(residue_str)
        if pocket_atoms is None:
            print(f"Warning: Residues {binding_site_residues} not found")
            pocket_atoms = ag_cleaned
    else:
        # Blind docking: use entire protein
        pocket_atoms = ag_cleaned
    
    # Calculate box dimensions
    coords = pocket_atoms.getCoords()
    box_center = coords.mean(axis=0)
    box_size = np.ptp(coords, axis=0) + np.array([6, 6, 6])  # Add padding
    
    print(f"\nDocking Grid Parameters:")
    print(f"Box center: ({box_center[0]:.2f}, {box_center[1]:.2f}, {box_center[2]:.2f})")
    print(f"Box size: ({box_size[0]:.2f}, {box_size[1]:.2f}, {box_size[2]:.2f})")
    
    return output_pdb, box_center, box_size

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare protein structure for docking')
    parser.add_argument('--pdb_id', required=True, help='PDB identifier (e.g., 3ODU)')
    parser.add_argument('--binding_residues', help='Comma-separated list of binding-site residues')
    parser.add_argument('--output_dir', default='./', help='Output directory')
    
    args = parser.parse_args()
    
    binding_residues = None
    if args.binding_residues:
        binding_residues = [int(r.strip()) for r in args.binding_residues.split(',')]
    
    pdb_file, center, size = prepare_protein(args.pdb_id, binding_residues, output_dir=args.output_dir)
    print(f"\n✓ Protein preparation complete: {pdb_file}")
