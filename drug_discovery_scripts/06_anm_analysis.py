#!/usr/bin/env python3
"""
Script 6: Anisotropic Network Model (ANM) Analysis with ProDy

Purpose: Analyse binding-pocket flexibility and collective motions with a
coarse-grained elastic network model on Cα atoms.

Output units: per-residue mean-square fluctuation (MSF) in ANM's relative
units, not Å. ANM is not calibrated to absolute B-factors here, so compare
values within one pocket and one protein, not against fixed thresholds.

Usage:
    python 06_anm_analysis.py --pdb protein.pdb --ligand_xyz 12.3,-4.5,20.1 \\
                              --pocket_radius 8.0 --output_dir anm_results
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from prody import ANM, calcSqFluct, parsePDB


def analyze_binding_pocket_anm(protein_pdb, ligand_center=None,
                               pocket_radius=8.0, n_modes=20,
                               output_dir='anm_results', verbose=True):
    """
    Perform Anisotropic Network Model (ANM) analysis of a binding pocket.

    ANM is a coarse-grained elastic network model using Cα atoms connected
    by harmonic springs. Diagonalising the Hessian yields normal modes that
    describe collective fluctuations without explicit solvent or force-field
    parametrisation.

    Args:
        protein_pdb: Protein structure (PDB file)
        ligand_center: Ligand centre coordinates (x, y, z), or None for the whole protein
        pocket_radius: Radius (Å) around the ligand centre that defines the pocket
        n_modes: Number of slowest modes to use
        output_dir: Directory for output files and plots
        verbose: Print progress messages

    Returns:
        ANM object, per-residue DataFrame, summary metrics dict
    """
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print(f"Loading protein structure from {protein_pdb}...")

    protein = parsePDB(protein_pdb)
    ca = protein.select('ca')
    if verbose:
        print(f"Found {ca.numAtoms()} Cα atoms")

    if ligand_center is not None:
        distances = np.linalg.norm(ca.getCoords() - np.asarray(ligand_center), axis=1)
        pocket_atoms = ca[np.where(distances < pocket_radius)[0]]
        if verbose:
            print(f"Binding pocket: {pocket_atoms.numAtoms()} Cα atoms within {pocket_radius} Å")
    else:
        pocket_atoms = ca
        if verbose:
            print("Analysing entire protein (no pocket defined)")

    if pocket_atoms.numAtoms() < 10:
        raise ValueError(f"Only {pocket_atoms.numAtoms()} Cα atoms in the pocket; "
                         "increase --pocket_radius (8 to 10 Å is typical)")

    if verbose:
        print("Building Anisotropic Network Model (cutoff=15.0 Å, gamma=1.0)...")
    anm = ANM('pocket')
    anm.buildHessian(pocket_atoms, cutoff=15.0, gamma=1.0)
    anm.calcModes(n_modes=n_modes)

    # calcSqFluct on the pocket ANM returns one value per pocket Cα, in pocket order.
    # Use pocket-local positions; residue numbers come from the same atoms.
    msf = calcSqFluct(anm[:n_modes])
    resnums = pocket_atoms.getResnums()

    metrics = {
        'mean_msf_relative': float(msf.mean()),
        'std_msf_relative': float(msf.std()),
        'min_msf_relative': float(msf.min()),
        'max_msf_relative': float(msf.max()),
        'n_pocket_calpha': int(pocket_atoms.numAtoms()),
        'n_modes': n_modes,
        'pocket_radius': pocket_radius,
    }

    if verbose:
        print(f"\n✓ ANM analysis complete ({n_modes} modes, relative units):")
        for key, value in metrics.items():
            print(f"  {key}: {value:.4g}" if isinstance(value, float) else f"  {key}: {value}")

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    axes[0].plot(resnums, msf, 'o-', linewidth=2, markersize=4)
    axes[0].axhline(msf.mean(), color='r', linestyle='--',
                    label=f"Mean: {msf.mean():.3g} (relative)")
    axes[0].set_xlabel('Residue number')
    axes[0].set_ylabel('MSF (ANM relative units)')
    axes[0].set_title('Pocket Cα fluctuations (ANM)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    eigvals = [anm[i].getEigval() for i in range(n_modes)]
    axes[1].bar(range(1, n_modes + 1), eigvals, color='steelblue')
    axes[1].set_xlabel('Mode number')
    axes[1].set_ylabel('Eigenvalue')
    axes[1].set_title('ANM mode eigenvalues (increasing = stiffer modes)')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(output_dir, 'anm_analysis.png')
    plt.savefig(plot_file, dpi=150)
    plt.close(fig)

    results_df = pd.DataFrame({'resnum': resnums, 'msf_relative': msf})
    csv_file = os.path.join(output_dir, 'anm_rmsf.csv')
    results_df.to_csv(csv_file, index=False)
    if verbose:
        print(f"\nPlot saved: {plot_file}")
        print(f"Per-residue values saved: {csv_file}")

    return anm, results_df, metrics


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform ANM analysis')
    parser.add_argument('--pdb', required=True, help='Protein PDB file')
    parser.add_argument('--ligand_xyz', help='Ligand centre (x,y,z)')
    parser.add_argument('--pocket_radius', type=float, default=8.0,
                        help='Pocket radius (Å); 8 to 10 is typical')
    parser.add_argument('--n_modes', type=int, default=20,
                        help='Number of slowest modes to use')
    parser.add_argument('--output_dir', default='anm_results',
                        help='Output directory')

    args = parser.parse_args()

    ligand_center = None
    if args.ligand_xyz:
        ligand_center = tuple(float(x) for x in args.ligand_xyz.split(','))

    anm, results_df, metrics = analyze_binding_pocket_anm(
        args.pdb, ligand_center, args.pocket_radius, args.n_modes,
        args.output_dir
    )
