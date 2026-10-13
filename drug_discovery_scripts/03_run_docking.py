#!/usr/bin/env python3
"""
Script 3: Virtual Screening with AutoDock Vina

Purpose: Perform high-throughput docking to identify hit compounds.
Requires: AutoDock Vina and prepare_ligand4.py (from ADFR) on PATH.

Usage:
    python 03_run_docking.py --receptor protein.pdbqt --ligands library.sdf \\
                             --center 0,0,0 --size 30,30,30
"""

import argparse
import os
import shutil
import subprocess
from collections import Counter

import pandas as pd
from rdkit import Chem


def read_vina_affinity(out_pdbqt):
    """Return the best (mode 1) affinity in kcal/mol from a Vina output PDBQT, or None."""
    with open(out_pdbqt) as f:
        for line in f:
            # Format: "REMARK VINA RESULT:  -7.2  0.000  0.000"
            if line.startswith('REMARK VINA RESULT:'):
                return float(line.split()[3])
    return None


def run_vina_docking(receptor_pdbqt, ligand_sdf, grid_center, grid_size,
                     output_dir='docking_results', exhaustiveness=8,
                     prep_timeout=60, vina_timeout=600, verbose=True):
    """
    Run AutoDock Vina for ligand library screening.

    Args:
        receptor_pdbqt: Prepared receptor in PDBQT format
        ligand_sdf: Ligand library SDF file (from 02_prepare_library.py)
        grid_center: Tuple (x, y, z) for docking box center
        grid_size: Tuple (sx, sy, sz) for docking box dimensions
        exhaustiveness: Search intensity (1-32; default 8)
        prep_timeout: Seconds allowed for prepare_ligand4.py per ligand
        vina_timeout: Seconds allowed for vina per ligand
        verbose: Print progress

    Returns:
        DataFrame with docking results (compound_id, smiles, binding_energy_kcal_mol)
    """
    for tool in ('vina', 'prepare_ligand4.py'):
        if shutil.which(tool) is None:
            raise FileNotFoundError(f"'{tool}' not found on PATH; install AutoDock Vina / ADFR")

    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        print(f"Starting virtual screening with AutoDock Vina...")
        print(f"Receptor: {receptor_pdbqt}")
        print(f"Grid center: {grid_center}")
        print(f"Grid size: {grid_size}")

    results = []
    skipped = Counter()
    suppl = Chem.SDMolSupplier(ligand_sdf)
    total = len(suppl)

    if verbose:
        print(f"Screening {total} compounds...\n")

    for i, mol in enumerate(suppl):
        if mol is None:
            skipped['unreadable_sdf_record'] += 1
            continue

        # Stable ID from stage 2; fall back to record index for hand-made SDFs
        cid = mol.GetIntProp('compound_id') if mol.HasProp('compound_id') else i

        if verbose and (i + 1) % 1000 == 0:
            print(f"  Processed {i+1}/{total} compounds...")

        ligand_pdb = os.path.join(output_dir, f'lig_{cid}.pdb')
        ligand_pdbqt = os.path.join(output_dir, f'lig_{cid}.pdbqt')
        output_pdbqt = os.path.join(output_dir, f'lig_{cid}_out.pdbqt')

        Chem.MolToPDBFile(mol, ligand_pdb)

        # List arguments, not shell=True: the project path contains spaces.
        try:
            prep = subprocess.run(
                ['prepare_ligand4.py', '-l', ligand_pdb, '-o', ligand_pdbqt],
                capture_output=True, text=True, timeout=prep_timeout)
        except subprocess.TimeoutExpired:
            skipped['prepare_ligand_timeout'] += 1
            continue
        if prep.returncode != 0 or not os.path.exists(ligand_pdbqt):
            skipped['prepare_ligand_failed'] += 1
            continue

        try:
            subprocess.run(
                ['vina', '--receptor', receptor_pdbqt, '--ligand', ligand_pdbqt,
                 '--center_x', str(grid_center[0]), '--center_y', str(grid_center[1]),
                 '--center_z', str(grid_center[2]),
                 '--size_x', str(grid_size[0]), '--size_y', str(grid_size[1]),
                 '--size_z', str(grid_size[2]),
                 '--out', output_pdbqt, '--exhaustiveness', str(exhaustiveness)],
                capture_output=True, text=True, timeout=vina_timeout, check=True)
        except subprocess.TimeoutExpired:
            skipped['vina_timeout'] += 1
            continue
        except subprocess.CalledProcessError:
            skipped['vina_failed'] += 1
            continue

        energy = read_vina_affinity(output_pdbqt) if os.path.exists(output_pdbqt) else None
        if energy is None:
            skipped['no_vina_result'] += 1
            continue

        results.append({
            'compound_id': cid,
            'smiles': Chem.MolToSmiles(Chem.RemoveHs(mol)),
            'binding_energy_kcal_mol': energy
        })

    df = pd.DataFrame(results, columns=['compound_id', 'smiles', 'binding_energy_kcal_mol'])
    if len(df) > 0:
        df = df.sort_values('binding_energy_kcal_mol')

    if verbose:
        print(f"\n✓ Docking complete!")
        print(f"  Successful poses: {len(df)}/{total}")
        if skipped:
            print(f"  Skipped: {dict(skipped)}")
        if len(df) > 0:
            print(f"  Best binding energy: {df['binding_energy_kcal_mol'].min():.2f} kcal/mol")
            print(f"  Worst binding energy: {df['binding_energy_kcal_mol'].max():.2f} kcal/mol")
            print(f"\nTop 10 hits:")
            print(df.head(10).to_string(index=False))

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
    parser.add_argument('--vina_timeout', type=int, default=600,
                        help='Seconds allowed per ligand for Vina')

    args = parser.parse_args()

    center = tuple(float(x) for x in args.center.split(','))
    size = tuple(float(x) for x in args.size.split(','))

    results = run_vina_docking(args.receptor, args.ligands, center, size,
                               args.output_dir, args.exhaustiveness,
                               vina_timeout=args.vina_timeout)
