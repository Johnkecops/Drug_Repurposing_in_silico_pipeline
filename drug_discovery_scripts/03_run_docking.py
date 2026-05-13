#!/usr/bin/env python3
"""
Script 3: Virtual Screening with AutoDock Vina

Purpose: Perform high-throughput docking to identify hit compounds.
Requires: AutoDock Vina, prepare_receptor4.py, prepare_ligand4.py (from ADFR)

Usage:
    python 03_run_docking.py --receptor protein.pdbqt --ligands library.sdf \\
                             --center 0,0,0 --size 30,30,30
"""

import argparse
import subprocess
import pandas as pd
from rdkit import Chem
import os
import re

def run_vina_docking(receptor_pdbqt, ligand_sdf, grid_center, grid_size,
                     output_dir='docking_results', exhaustiveness=8, verbose=True):
    """
    Run AutoDock Vina for ligand library screening.
    
    Args:
        receptor_pdbqt: Prepared receptor in PDBQT format
        ligand_sdf: Ligand library SDF file
        grid_center: Tuple (x, y, z) for docking box center
        grid_size: Tuple (sx, sy, sz) for docking box dimensions
        exhaustiveness: Search intensity (1-32; default 8)
        verbose: Print progress
    
    Returns:
        DataFrame with docking results (compound_id, smiles, binding_energy)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if verbose:
        print(f"Starting virtual screening with AutoDock Vina...")
        print(f"Receptor: {receptor_pdbqt}")
        print(f"Grid center: {grid_center}")
        print(f"Grid size: {grid_size}")
    
    results = []
    suppl = Chem.SDMolSupplier(ligand_sdf)
    total = len(suppl)
    
    if verbose:
        print(f"Screening {total} compounds...\n")
    
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        
        # Report progress
        if verbose and (i + 1) % 1000 == 0:
            print(f"  Processed {i+1}/{total} compounds...")
        
        # Prepare ligand PDBQT
        ligand_pdb = f'{output_dir}/lig_{i}.pdb'
        ligand_pdbqt = f'{output_dir}/lig_{i}.pdbqt'
        
        try:
            Chem.MolToPDBFile(mol, ligand_pdb)
        except:
            continue
        
        # Convert to PDBQT (requires prepare_ligand4.py from ADFR)
        try:
            subprocess.run(
                f'prepare_ligand4.py -l {ligand_pdb} -o {ligand_pdbqt}',
                shell=True, capture_output=True, timeout=5
            )
        except:
            continue
        
        # Run Vina docking
        output_pdbqt = f'{output_dir}/lig_{i}_out.pdbqt'
        cmd = (
            f'vina --receptor {receptor_pdbqt} '
            f'--ligand {ligand_pdbqt} '
            f'--center_x {grid_center[0]} --center_y {grid_center[1]} '
            f'--center_z {grid_center[2]} '
            f'--size_x {grid_size[0]} --size_y {grid_size[1]} '
            f'--size_z {grid_size[2]} '
            f'--out {output_pdbqt} --exhaustiveness {exhaustiveness}'
        )
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                   text=True, timeout=60)
            
            # Extract binding energy from output
            if 'kcal/mol' in result.stdout:
                lines = result.stdout.split('\n')
                for line in lines:
                    if 'kcal/mol' in line and line.strip()[0].isdigit():
                        try:
                            energy = float(line.split()[0])
                            smiles = Chem.MolToSmiles(mol)
                            results.append({
                                'compound_id': i,
                                'smiles': smiles,
                                'binding_energy_kcal_mol': energy
                            })
                            break
                        except:
                            pass
        except:
            continue
    
    # Create results dataframe and sort
    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values('binding_energy_kcal_mol')
    
    if verbose:
        print(f"\n✓ Docking complete!")
        print(f"  Successful poses: {len(df)}/{total}")
        if len(df) > 0:
            print(f"  Best binding energy: {df['binding_energy_kcal_mol'].min():.2f} kcal/mol")
            print(f"  Worst binding energy: {df['binding_energy_kcal_mol'].max():.2f} kcal/mol")
            print(f"\nTop 10 hits:")
            print(df.head(10).to_string(index=False))
    
    # Save results
    csv_output = os.path.join(output_dir, 'docking_results.csv')
    df.to_csv(csv_output, index=False)
    if verbose:
        print(f"\nResults saved: {csv_output}")
    
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run AutoDock Vina screening')
    parser.add_argument('--receptor', required=True, help='Receptor PDBQT file')
    parser.add_argument('--ligands', required=True, help='Ligand library SDF file')
    parser.add_argument('--center', required=True, help='Grid center (x,y,z)')
    parser.add_argument('--size', required=True, help='Grid size (sx,sy,sz)')
    parser.add_argument('--output_dir', default='docking_results', help='Output directory')
    parser.add_argument('--exhaustiveness', type=int, default=8, help='Search exhaustiveness')
    
    args = parser.parse_args()
    
    center = tuple(float(x) for x in args.center.split(','))
    size = tuple(float(x) for x in args.size.split(','))
    
    results = run_vina_docking(args.receptor, args.ligands, center, size,
                              args.output_dir, args.exhaustiveness)
