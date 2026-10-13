#!/usr/bin/env python3
"""
Script 2: Chemical Library Preparation and Drug-Likeness Filtering

Purpose: Standardise SMILES, add hydrogens, apply Lipinski's Rule of Five,
generate 3D structures for docking. Each output record carries the stable
'compound_id' property (see smiles_io.py) used by later stages.

Usage:
    python 02_prepare_library.py --smiles_file compounds.smi --output_sdf library_prepared.sdf
"""

import argparse
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, AllChem

from smiles_io import read_smiles


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
    records = read_smiles(smiles_file)
    total = len(records)
    mols = []
    rejected_counts = {'invalid_smiles': 0, 'ro5_violation': 0, 'other': 0}

    if verbose:
        print(f"Reading {total} SMILES from {smiles_file}...")

    for cid, smiles in records:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            rejected_counts['invalid_smiles'] += 1
            continue

        # Lipinski's Rule of Five. Rotatable bonds are not part of RO5 (Veber's rule).
        if apply_ro5:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            if mw > 500 or logp > 5 or hbd > 5 or hba > 10:
                rejected_counts['ro5_violation'] += 1
                continue

        mol = Chem.AddHs(mol)
        # EmbedMolecule returns -1 when no 3D conformer can be generated.
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            rejected_counts['other'] += 1
            continue
        AllChem.MMFFOptimizeMolecule(mol)

        mol.SetIntProp('compound_id', cid)
        mols.append(mol)

    if verbose:
        print(f"Writing {len(mols)} compounds to {output_sdf}...")

    writer = Chem.SDWriter(output_sdf)
    for mol in mols:
        writer.write(mol)
    writer.close()

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
        print(f"  Embedding failures: {stats['other_errors']}")

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
