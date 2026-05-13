#!/usr/bin/env python3
"""
Script 8: Multi-Criteria Decision Analysis (MCDA) for Candidate Ranking

Purpose: Integrate docking, molecular dynamics, and ADMET metrics into
a unified ranking system to identify optimal lead candidates.

Ranking criteria (weights):
- Binding energy (25%): Primary pharmacophore requirement
- Conformational stability (20%): RMSD plateau & persistence
- Pocket flexibility (15%): ANM RMSF metrics
- Drug-likeness (20%): Lipinski RO5 compliance
- Toxicity/Safety (10%): hERG, PAINS alerts
- Synthesis tractability (10%): SAscore

Usage:
    python 08_rank_candidates.py --docking docking_results.csv \\
                                 --admet admet_results.csv \\
                                 --output ranked_candidates.csv
"""

import argparse
import pandas as pd
import numpy as np
import os

def rank_candidates_mcda(docking_df, admet_df, dynamics_dict=None,
                        weights=None, output_csv='ranked_candidates.csv',
                        verbose=True):
    """
    Multi-Criteria Decision Analysis: rank candidates by composite score.
    
    Args:
        docking_df: DataFrame from virtual screening (binding_energy)
        admet_df: DataFrame from ADMET assessment (RO5_Violations, etc.)
        dynamics_dict: Dict mapping compound_id -> dynamics metrics (optional)
        weights: Dict of criterion weights (uses defaults if None)
        output_csv: Output CSV file
        verbose: Print progress
    
    Returns:
        Ranked DataFrame of top candidates
    """
    if dynamics_dict is None:
        dynamics_dict = {}
    
    if weights is None:
        weights = {
            'binding_energy': 0.25,
            'stability': 0.20,
            'flexibility': 0.15,
            'druggability': 0.20,
            'safety': 0.10,
            'synthesis': 0.10
        }
    
    if verbose:
        print(f"Merging datasets...")
        print(f"  Docking: {len(docking_df)} compounds")
        print(f"  ADMET: {len(admet_df)} compounds")
        print(f"  Dynamics: {len(dynamics_dict)} compounds")
    
    candidates = []
    
    # Iterate over docking results
    for idx, dock_row in docking_df.iterrows():
        cid = dock_row['compound_id']
        
        # Match with ADMET data
        admet_matches = admet_df[admet_df['compound_id'] == cid]
        if admet_matches.empty:
            continue
        
        admet_row = admet_matches.iloc[0]
        
        # Normalize binding energy (lower = better)
        # Assume range -12 to 0 kcal/mol (typical for good binders)
        be = dock_row['binding_energy_kcal_mol']
        be_norm = max(0, min(1.0, (be + 12) / 12))  # Clamp to [0,1]
        
        # Normalize stability (RMSD; lower = better)
        rmsd = dynamics_dict.get(cid, {}).get('rmsd', 2.0)
        stability_norm = max(0, 1.0 - (rmsd / 3.0))  # Assume max RMSD ~3 Å
        
        # Normalize flexibility (RMSF; moderate is best)
        rmsf = dynamics_dict.get(cid, {}).get('pocket_rmsf', 1.0)
        # Optimal ~0.8 Å; penalise if too rigid or too flexible
        flexibility_norm = 1.0 - abs(rmsf - 0.8) / 1.5
        
        # Drug-likeness (RO5)
        ro5_viol = admet_row['RO5_Violations']
        druggability = 1.0 if ro5_viol <= 1 else max(0, 1.0 - ro5_viol/3)
        
        # Safety (hERG block risk, PAINS alerts; simplified)
        # Proxy: LogP > 5 increases hERG risk
        logp = admet_row['LogP']
        safety_score = 1.0 if logp <= 5 else max(0.3, 1.0 - (logp - 5)/5)
        
        # Synthesis tractability (simplified; use SAscore if available)
        # Proxy: Number of aromatic rings and rotatable bonds
        aroma_rings = admet_row['AromaRings']
        rot_bonds = admet_row['RotBonds']
        synthesis_score = 1.0 - (aroma_rings + rot_bonds) / 20  # Rough estimate
        synthesis_score = max(0.2, min(1.0, synthesis_score))
        
        # Composite score
        composite_score = (
            weights['binding_energy'] * be_norm +
            weights['stability'] * stability_norm +
            weights['flexibility'] * flexibility_norm +
            weights['druggability'] * druggability +
            weights['safety'] * safety_score +
            weights['synthesis'] * synthesis_score
        )
        
        candidates.append({
            'rank_order': 0,  # Will set after sorting
            'compound_id': cid,
            'smiles': admet_row.get('SMILES', ''),
            'binding_energy_kcal_mol': be,
            'rmsd': rmsd,
            'pocket_rmsf': rmsf,
            'logp': logp,
            'ro5_violations': ro5_viol,
            'drug_like': 'YES' if ro5_viol <= 1 else 'NO',
            'be_score': be_norm,
            'stability_score': stability_norm,
            'flexibility_score': flexibility_norm,
            'druggability_score': druggability,
            'safety_score': safety_score,
            'synthesis_score': synthesis_score,
            'composite_score': composite_score
        })
    
    # Create DataFrame and rank
    df_ranked = pd.DataFrame(candidates).sort_values('composite_score', ascending=False)
    df_ranked['rank_order'] = range(1, len(df_ranked) + 1)
    
    # Reorder columns
    cols = ['rank_order', 'compound_id', 'smiles', 'binding_energy_kcal_mol',
            'composite_score', 'be_score', 'stability_score', 'flexibility_score',
            'druggability_score', 'safety_score', 'synthesis_score']
    df_ranked = df_ranked[cols + [c for c in df_ranked.columns if c not in cols]]
    
    # Save results
    df_ranked.to_csv(output_csv, index=False)
    
    if verbose:
        print(f"\n✓ Ranking Complete!")
        print(f"\nTop 10 Lead Candidates (MCDA Ranking):")
        print(df_ranked[['rank_order', 'compound_id', 'binding_energy_kcal_mol',
                         'composite_score']].head(10).to_string(index=False))
        print(f"\nDetailed scores for top 5:")
        for idx in range(min(5, len(df_ranked))):
            row = df_ranked.iloc[idx]
            print(f"\n  Rank #{int(row['rank_order'])}: Compound {int(row['compound_id'])}")
            print(f"    Binding Energy: {row['binding_energy_kcal_mol']:.2f} kcal/mol")
            print(f"    Composite Score: {row['composite_score']:.3f}")
            print(f"    Components: BE={row['be_score']:.2f}, Stab={row['stability_score']:.2f}, "
                  f"Flex={row['flexibility_score']:.2f}")
        
        print(f"\nResults saved: {output_csv}")
    
    return df_ranked

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rank candidates via MCDA')
    parser.add_argument('--docking', required=True,
                       help='Docking results CSV')
    parser.add_argument('--admet', required=True,
                       help='ADMET results CSV')
    parser.add_argument('--dynamics', help='Dynamics metrics JSON or CSV')
    parser.add_argument('--output', default='ranked_candidates.csv',
                       help='Output CSV file')
    
    args = parser.parse_args()
    
    # Load dataframes
    docking_df = pd.read_csv(args.docking)
    admet_df = pd.read_csv(args.admet)
    
    dynamics_dict = {}
    if args.dynamics and os.path.isfile(args.dynamics):
        if args.dynamics.endswith('.csv'):
            dyn_df = pd.read_csv(args.dynamics)
            for idx, row in dyn_df.iterrows():
                dynamics_dict[int(row['compound_id'])] = row.to_dict()
    
    ranked = rank_candidates_mcda(docking_df, admet_df, dynamics_dict, 
                                 output_csv=args.output)
