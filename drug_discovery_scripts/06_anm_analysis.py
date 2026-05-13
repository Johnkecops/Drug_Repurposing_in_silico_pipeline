#!/usr/bin/env python3
"""
Script 6: Anisotropic Network Model (ANM) Analysis with ProDy

Purpose: Analyse binding-pocket flexibility and conformational dynamics
using coarse-grained elastic network model. Calculate RMSF, identify
collective motions, assess pocket stability.

Usage:
    python 06_anm_analysis.py --pdb protein.pdb --ligand_xyz 0,0,0 \\
                              --output_dir anm_results
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from prody import *
import os

def analyze_binding_pocket_anm(protein_pdb, ligand_center=None, 
                               pocket_radius=5.0, n_modes=20,
                               output_dir='anm_results', verbose=True):
    """
    Perform Anisotropic Network Model (ANM) analysis of binding pocket.
    
    ANM is a coarse-grained elastic network model using Cα atoms connected
    by harmonic springs. Diagonalisation of the Hessian yields normal modes
    capturing collective thermal fluctuations at physiological temperature
    without explicit solvent or force-field parametrisation.
    
    Args:
        protein_pdb: Protein structure (PDB file)
        ligand_center: Ligand center coordinates (x, y, z) or None
        pocket_radius: Radius (Å) around ligand for pocket definition
        n_modes: Number of slowest modes to analyse
        output_dir: Directory for output files and plots
        verbose: Print progress messages
    
    Returns:
        ANM object, per-residue RMSF, pocket flexibility metrics
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if verbose:
        print(f"Loading protein structure from {protein_pdb}...")
    
    # Load protein
    protein = parsePDB(protein_pdb)
    ca = protein.select('ca')
    
    if verbose:
        print(f"Found {ca.numAtoms()} Cα atoms")
    
    # Define binding pocket
    if ligand_center is not None:
        ligand_center = np.array(ligand_center)
        # Select Cα atoms within pocket_radius of ligand
        distances = np.linalg.norm(ca.getCoords() - ligand_center, axis=1)
        pocket_indices = np.where(distances < pocket_radius)[0]
        pocket_atoms = ca[pocket_indices]
        
        if verbose:
            print(f"Binding pocket: {len(pocket_atoms)} Cα atoms within {pocket_radius} Å")
    else:
        pocket_atoms = ca
        if verbose:
            print("Analysing entire protein (no specific pocket defined)")
    
    # Build ANM
    if verbose:
        print(f"Building Anisotropic Network Model (cutoff=15.0 Å)...")
    
    anm = ANM('protein_pocket')
    anm.buildHessian(pocket_atoms, cutoff=15.0, gamma=1.0)
    
    if verbose:
        print(f"Computing {n_modes} slowest normal modes...")
    
    anm.calcModes(n_modes=n_modes)
    
    # Calculate RMSF from slowest modes
    rmsf = calcSqFluct(anm[:n_modes])
    
    # Per-residue flexibility
    pocket_rmsf = rmsf[pocket_atoms.getIndices()]
    
    # Statistics
    rmsf_mean = pocket_rmsf.mean()
    rmsf_std = pocket_rmsf.std()
    rmsf_min = pocket_rmsf.min()
    rmsf_max = pocket_rmsf.max()
    
    # Classify flexibility
    if rmsf_mean < 0.8:
        flexibility_class = "RIGID (favourable for stable binding)"
    elif rmsf_mean < 1.5:
        flexibility_class = "MODERATE (ideal for most drug targets)"
    else:
        flexibility_class = "FLEXIBLE (entropy-driven binding possible)"
    
    if verbose:
        print(f"\n✓ ANM Analysis Complete!")
        print(f"\nBinding Pocket Flexibility (RMSF from {n_modes} slowest modes):")
        print(f"  Mean:   {rmsf_mean:.2f} Å")
        print(f"  Std:    {rmsf_std:.2f} Å")
        print(f"  Min:    {rmsf_min:.2f} Å")
        print(f"  Max:    {rmsf_max:.2f} Å")
        print(f"  Class:  {flexibility_class}")
    
    # Generate visualisation
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Per-residue RMSF
    residue_indices = pocket_atoms.getResindices()
    axes[0].plot(residue_indices, pocket_rmsf, 'o-', linewidth=2, markersize=4)
    axes[0].axhline(rmsf_mean, color='r', linestyle='--', label=f'Mean: {rmsf_mean:.2f} Å')
    axes[0].set_xlabel('Residue Index')
    axes[0].set_ylabel('RMSF (Å)')
    axes[0].set_title('Binding Pocket Flexibility (ANM Normal Modes)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Mode 1 projection
    mode_1_energy = anm[0].getEigval()
    axes[1].bar(range(1, n_modes+1), 
               [anm[i].getEigval() for i in range(n_modes)],
               color='steelblue')
    axes[1].set_xlabel('Mode Number')
    axes[1].set_ylabel('Eigenvalue (λ)')
    axes[1].set_title('ANM Mode Eigenvalues (Collective Motion Frequencies)')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, 'anm_analysis.png')
    plt.savefig(plot_file, dpi=150)
    if verbose:
        print(f"\nPlot saved: {plot_file}")
    
    # Save results to CSV
    import pandas as pd
    results_df = pd.DataFrame({
        'residue_index': residue_indices,
        'rmsf_angstrom': pocket_rmsf
    })
    csv_file = os.path.join(output_dir, 'anm_rmsf.csv')
    results_df.to_csv(csv_file, index=False)
    if verbose:
        print(f"RMSF values saved: {csv_file}")
    
    return anm, rmsf, {
        'mean_rmsf': rmsf_mean,
        'std_rmsf': rmsf_std,
        'flexibility_class': flexibility_class,
        'n_modes': n_modes,
        'pocket_radius': pocket_radius
    }

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform ANM analysis')
    parser.add_argument('--pdb', required=True, help='Protein PDB file')
    parser.add_argument('--ligand_xyz', help='Ligand center (x,y,z)')
    parser.add_argument('--pocket_radius', type=float, default=5.0,
                       help='Pocket radius (Å)')
    parser.add_argument('--n_modes', type=int, default=20,
                       help='Number of slowest modes to analyse')
    parser.add_argument('--output_dir', default='anm_results',
                       help='Output directory')
    
    args = parser.parse_args()
    
    ligand_center = None
    if args.ligand_xyz:
        ligand_center = tuple(float(x) for x in args.ligand_xyz.split(','))
    
    anm, rmsf, metrics = analyze_binding_pocket_anm(
        args.pdb, ligand_center, args.pocket_radius, args.n_modes, 
        args.output_dir
    )
