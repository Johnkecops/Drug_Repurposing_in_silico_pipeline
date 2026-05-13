#!/usr/bin/env python3
"""
Script 2: Chemical Library Preparation and Drug-Likeness Filtering

Purpose: Standardise SMILES, add hydrogens, apply Lipinski's Rule of Five,
generate 3D structures for docking.

Usage:
    python 02_prepare_library.py --smiles_file compounds.smi --output_sdf library_prepared.sdf
"""

import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, AllChem
import os

def prepare_ligand_library(smiles_file, output_sdf, apply_ro5=True, verbose=True):
    """
    Prepare ligand library: standardise SMILES, add hydrogens, apply filters.
    
    Args:
        smiles_file: Text file with SMILES strings (one per line)
        output_sdf: Output SDF file with prepared 3D structures
        apply_ro5: Apply Lipinski's Rule of Five filtering
        verbose: Print progress messages
    
    Returns:
        Output SDF filename, number of retained compounds, statistics dict
    """
    mols = []
    rejected_counts = {'invalid_smiles': 0, 'ro5_violation': 0, 'other': 0}
    
    if verbose:
        print(f"Reading SMILES from {smiles_file}...")
    
    with open(smiles_file) as f:
        total = sum(1 for _ in f)
    
    with open(smiles_file) as f:
        for i, line in enumerate(f):
            smiles = line.strip()
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                rejected_counts['invalid_smiles'] += 1
                continue
            
            # Apply Lipinski's Rule of Five
            if apply_ro5:
                mw = Descriptors.MolWt(mol)
                logp = Crippen.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                rotatable = Descriptors.NumRotatableBonds(mol)
                
                # Check RO5 criteria
                if mw > 500 or logp > 5 or hbd > 5 or hba > 10 or rotatable > 10:
                    rejected_counts['ro5_violation'] += 1
                    continue
            
            try:
                # Add hydrogens
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.UMMFFOptimizeMolecule(mol)
                
                mols.append(mol)
            except:
                rejected_counts['other'] += 1
                continue
            
            if verbose and (i + 1) % 10000 == 0:
                print(f"  Processed {i+1}/{total} compounds...")
    
    # Write SDF
    if verbose:
        print(f"Writing {len(mols)} compounds to {output_sdf}...")
    
    writer = Chem.SDWriter(output_sdf)
    for mol in mols:
        writer.write(mol)
    writer.close()
    
    # Statistics
    stats = {
        'input_compounds': total,
        'output_compounds': len(mols),
        'invalid_smiles': rejected_counts['invalid_smiles'],
        'ro5_violations': rejected_counts['ro5_violation'],
        'other_errors': rejected_counts['other'],
        'retention_rate': 100 * len(mols) / total if total > 0 else 0
    }
    
    if verbose:
        print(f"\n✓ Library preparation complete:")
        print(f"  Input compounds: {stats['input_compounds']}")
        print(f"  Output compounds: {stats['output_compounds']}")
        print(f"  Retention rate: {stats['retention_rate']:.1f}%")
        print(f"  Invalid SMILES: {stats['invalid_smiles']}")
        print(f"  RO5 violations: {stats['ro5_violations']}")
    
    return output_sdf, len(mols), stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare chemical library for docking')
    parser.add_argument('--smiles_file', required=True, help='Input SMILES file')
    parser.add_argument('--output_sdf', required=True, help='Output SDF file')
    parser.add_argument('--no_ro5', action='store_true', help='Skip Lipinski RO5 filtering')
    
    args = parser.parse_args()
    
    sdf_file, n_mols, stats = prepare_ligand_library(
        args.smiles_file, 
        args.output_sdf,
        apply_ro5=not args.no_ro5
    )
