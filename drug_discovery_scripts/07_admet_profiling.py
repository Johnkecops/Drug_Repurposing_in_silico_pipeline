#!/usr/bin/env python3
"""
Script 7: ADMET Profiling and Drug-Likeness Assessment

Purpose: Compute drug-likeness descriptors (Lipinski's Rule of Five, TPSA,
rule-of-thumb BBB and solubility heuristics). These are not measured ADMET:
no hERG, metabolism, toxicity, or PAINS screen is performed.

Usage:
    python 07_admet_profiling.py --smiles_file compounds.smi \\
                                 --output_csv admet_results.csv
"""

import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import os

from smiles_io import read_smiles

def assess_admet_properties(smiles_list, output_csv='admet_results.csv', 
                           verbose=True):
    """
    Evaluate Lipinski's Rule of Five and ADMET properties for compounds.
    
    Args:
        smiles_list: List of SMILES strings (or filename)
        output_csv: Output CSV file with ADMET results
        verbose: Print progress messages
    
    Returns:
        DataFrame with ADMET properties and drug-likeness assessment
    """
    results = []

    # Stable compound IDs shared with stages 2 and 3 (see smiles_io.py)
    if isinstance(smiles_list, str) and os.path.isfile(smiles_list):
        if verbose:
            print(f"Reading SMILES from {smiles_list}...")
        records = read_smiles(smiles_list)
    else:
        records = list(enumerate(smiles_list))

    if verbose:
        print(f"Assessing {len(records)} compounds...")

    for i, smiles in records:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        
        # Lipinski's Rule of Five
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Count RO5 violations
        ro5_violations = 0
        if mw > 500:
            ro5_violations += 1
        if logp > 5:
            ro5_violations += 1
        if hbd > 5:
            ro5_violations += 1
        if hba > 10:
            ro5_violations += 1
        
        # Drug-likeness classification
        drug_like = 'YES' if ro5_violations <= 1 else 'NO'
        
        # Additional ADMET properties
        tpsa = Descriptors.TPSA(mol)  # Polar surface area (BBB penetration)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        sp3_fraction = Descriptors.FractionCSP3(mol)
        
        # BBB penetrance prediction (TPSA < 60 and LogP 1-5)
        bbb_penetrant = 'YES' if (tpsa < 60 and 1 < logp < 5) else 'NO'
        
        # Aqueous solubility class (simplified)
        if logp > 5:
            solubility = 'LOW'
        elif logp < 0:
            solubility = 'HIGH'
        else:
            solubility = 'MODERATE'
        
        results.append({
            'compound_id': i,
            'SMILES': smiles,
            'MW_Da': mw,
            'LogP': logp,
            'HBD': hbd,
            'HBA': hba,
            'RotBonds': rotatable_bonds,
            'TPSA_A2': tpsa,
            'AromaRings': aromatic_rings,
            'Sp3Frac': sp3_fraction,
            'RO5_Violations': ro5_violations,
            'Drug_Like': drug_like,
            'BBB_Penetrant': bbb_penetrant,
            'Solubility': solubility
        })
    
    df = pd.DataFrame(results)
    
    # Save results
    df.to_csv(output_csv, index=False)
    
    if verbose:
        print(f"\n✓ ADMET Assessment Complete!")
        print(f"\nSummary Statistics:")
        print(f"  Total compounds evaluated: {len(df)}")
        print(f"  Drug-like (RO5 ≤ 1): {(df['Drug_Like']=='YES').sum()} "
              f"({100*(df['Drug_Like']=='YES').sum()/len(df):.1f}%)")
        print(f"  BBB penetrant: {(df['BBB_Penetrant']=='YES').sum()} "
              f"({100*(df['BBB_Penetrant']=='YES').sum()/len(df):.1f}%)")
        print(f"\nRO5 Violation Distribution:")
        print(df['RO5_Violations'].value_counts().sort_index())
        print(f"\nResults saved: {output_csv}")
    
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assess ADMET properties')
    parser.add_argument('--smiles_file', required=True, 
                       help='Input SMILES file or list')
    parser.add_argument('--output_csv', default='admet_results.csv',
                       help='Output CSV file')
    
    args = parser.parse_args()
    
    df = assess_admet_properties(args.smiles_file, args.output_csv)
