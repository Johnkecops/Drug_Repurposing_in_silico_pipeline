#!/usr/bin/env python3
"""
Script 8: Multi-Criteria Decision Analysis (MCDA) for Candidate Ranking

Purpose: Combine docking energy, drug-likeness descriptors, and (optionally)
per-compound dynamics into one weighted score.

Terms and default weights:
- Binding energy (0.25): docking score, linearly scaled from -12 to 0 kcal/mol
- Druggability (0.20): Lipinski RO5 violations
- Safety (0.10): LogP proxy only (not a hERG or PAINS screen)
- Synthesis (0.10): ring and rotor count proxy (not SAscore)
- Stability (0.20): RMSD per compound. Needs --dynamics.
- Flexibility (0.15): pocket MSF per compound. Needs --dynamics.

Without --dynamics, stability and flexibility are dropped and the remaining
weights are renormalised. With --dynamics, every docked compound needs a row.

Dynamics CSV columns: compound_id, rmsd (Å), pocket_msf (ANM relative units).

Usage:
    python 08_rank_candidates.py --docking docking_results/docking_results.csv \\
                                 --admet admet_results.csv \\
                                 --dynamics dynamics.csv \\
                                 --output ranked_candidates.csv
"""

import argparse
import os

import pandas as pd

DEFAULT_WEIGHTS = {
    'binding_energy': 0.25,
    'stability': 0.20,
    'flexibility': 0.15,
    'druggability': 0.20,
    'safety': 0.10,
    'synthesis': 0.10,
}
DYNAMICS_TERMS = ('stability', 'flexibility')


def rank_candidates_mcda(docking_df, admet_df, dynamics_dict=None,
                         weights=None, output_csv='ranked_candidates.csv',
                         verbose=True):
    """
    Multi-Criteria Decision Analysis: rank candidates by composite score.

    Args:
        docking_df: DataFrame with compound_id, binding_energy_kcal_mol
        admet_df: DataFrame with compound_id, SMILES, LogP, RO5_Violations, AromaRings, RotBonds
        dynamics_dict: Dict compound_id -> {'rmsd': float, 'pocket_msf': float}, or None
        weights: Dict of criterion weights (defaults above if None)
        output_csv: Output CSV file
        verbose: Print progress

    Returns:
        Ranked DataFrame
    """
    if docking_df.empty:
        raise ValueError("Docking results are empty; nothing to rank. "
                         "Check the skip counts printed by 03_run_docking.py.")

    dynamics_dict = dynamics_dict or {}
    weights = dict(DEFAULT_WEIGHTS if weights is None else weights)

    if dynamics_dict:
        missing = [int(c) for c in docking_df['compound_id'] if int(c) not in dynamics_dict]
        if missing:
            raise ValueError(f"Dynamics missing for {len(missing)} compounds "
                             f"(first few: {missing[:5]}). Provide every compound or omit --dynamics.")
    else:
        for term in DYNAMICS_TERMS:
            weights.pop(term)

    total_weight = sum(weights.values())
    weights = {k: v / total_weight for k, v in weights.items()}

    if verbose:
        print(f"Merging datasets...")
        print(f"  Docking: {len(docking_df)} compounds")
        print(f"  ADMET: {len(admet_df)} compounds")
        print(f"  Dynamics: {len(dynamics_dict)} compounds")
        print(f"  Active terms: {', '.join(weights)}")

    admet_by_id = admet_df.drop_duplicates('compound_id').set_index('compound_id')

    candidates = []
    for _, dock_row in docking_df.iterrows():
        cid = int(dock_row['compound_id'])
        if cid not in admet_by_id.index:
            continue
        admet_row = admet_by_id.loc[cid]

        be = float(dock_row['binding_energy_kcal_mol'])
        be_norm = max(0.0, min(1.0, (be + 12) / 12))

        ro5_viol = int(admet_row['RO5_Violations'])
        druggability = 1.0 if ro5_viol <= 1 else max(0.0, 1.0 - ro5_viol / 3)

        # LogP above 5 is penalised as a crude proxy for hERG risk
        logp = float(admet_row['LogP'])
        safety = 1.0 if logp <= 5 else max(0.3, 1.0 - (logp - 5) / 5)

        # Ring and rotor count as a rough synthesis-effort proxy
        synthesis = max(0.2, min(1.0, 1.0 - (int(admet_row['AromaRings']) + int(admet_row['RotBonds'])) / 20))

        terms = {
            'binding_energy': be_norm,
            'druggability': druggability,
            'safety': safety,
            'synthesis': synthesis,
        }
        row = {
            'compound_id': cid,
            'smiles': admet_row.get('SMILES', ''),
            'binding_energy_kcal_mol': be,
            'logp': logp,
            'ro5_violations': ro5_viol,
        }

        if dynamics_dict:
            dyn = dynamics_dict[cid]
            rmsd = float(dyn['rmsd'])
            msf = float(dyn['pocket_msf'])
            terms['stability'] = max(0.0, 1.0 - rmsd / 3.0)  # assumes RMSD of 3 Å is the worst case
            # Centre 0.8 is unvalidated. Calibrate on your own pocket before relying on this term.
            terms['flexibility'] = 1.0 - abs(msf - 0.8) / 1.5
            row.update({'rmsd': rmsd, 'pocket_msf': msf})

        row.update({f'{k}_score': v for k, v in terms.items()})
        row['composite_score'] = sum(weights[k] * terms[k] for k in weights)
        candidates.append(row)

    if not candidates:
        raise ValueError("No docked compound has a matching ADMET row. Check that compound IDs "
                         "come from the same SMILES file (see smiles_io.py).")

    df_ranked = pd.DataFrame(candidates).sort_values('composite_score', ascending=False)
    df_ranked.insert(0, 'rank_order', range(1, len(df_ranked) + 1))
    df_ranked.to_csv(output_csv, index=False)

    if verbose:
        print(f"\n✓ Ranking complete ({len(df_ranked)} compounds):")
        print(df_ranked[['rank_order', 'compound_id', 'binding_energy_kcal_mol',
                         'composite_score']].head(10).to_string(index=False))
        print(f"\nResults saved: {output_csv}")

    return df_ranked


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rank candidates via MCDA')
    parser.add_argument('--docking', required=True, help='Docking results CSV')
    parser.add_argument('--admet', required=True, help='ADMET results CSV')
    parser.add_argument('--dynamics', help='Per-compound dynamics CSV (compound_id, rmsd, pocket_msf)')
    parser.add_argument('--output', default='ranked_candidates.csv', help='Output CSV file')

    args = parser.parse_args()

    docking_df = pd.read_csv(args.docking)
    admet_df = pd.read_csv(args.admet)

    dynamics_dict = None
    if args.dynamics:
        if not os.path.isfile(args.dynamics) or not args.dynamics.endswith('.csv'):
            raise ValueError(f"--dynamics must be an existing .csv file: {args.dynamics}")
        dyn_df = pd.read_csv(args.dynamics)
        dynamics_dict = {int(row['compound_id']): row.to_dict() for _, row in dyn_df.iterrows()}

    rank_candidates_mcda(docking_df, admet_df, dynamics_dict, output_csv=args.output)
