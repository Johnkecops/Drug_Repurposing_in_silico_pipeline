#!/usr/bin/env python3
"""
Self-check for the pipeline fixes. Run from this folder:
    python test_pipeline_fixes.py

Covers logic that does not need ProDy, OpenMM, or Vina: compound IDs,
library prep, Vina output parsing, ADMET IDs, and MCDA weighting.
"""

import importlib
import os
import tempfile

import pandas as pd
from rdkit import Chem

from smiles_io import read_smiles

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = os.path.join(HERE, 'sample_compounds.smi')

lib = importlib.import_module('02_prepare_library')
admet = importlib.import_module('07_admet_profiling')
docking = importlib.import_module('03_run_docking')
ranking = importlib.import_module('08_rank_candidates')


def test_read_smiles_skips_comments_and_numbers_stably():
    records = read_smiles(SAMPLE)
    assert records, "sample file produced no records"
    assert [cid for cid, _ in records] == list(range(len(records)))
    assert all(not s.startswith('#') for _, s in records)


def test_library_prep_keeps_stable_ids():
    ids = {cid for cid, _ in read_smiles(SAMPLE)}
    with tempfile.TemporaryDirectory() as tmp:
        sdf = os.path.join(tmp, 'lib.sdf')
        _, n_out, stats = lib.prepare_ligand_library(SAMPLE, sdf, verbose=False)
        assert n_out > 0, f"nothing retained: {stats}"
        seen = [int(m.GetIntProp('compound_id')) for m in Chem.SDMolSupplier(sdf) if m]
    assert len(seen) == n_out
    assert set(seen) <= ids


def test_vina_affinity_parsed_from_remark_line():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, 'out.pdbqt')
        with open(out, 'w') as f:
            f.write("MODEL 1\nREMARK VINA RESULT:    -7.2      0.000      0.000\nENDMDL\n")
        assert docking.read_vina_affinity(out) == -7.2
        empty = os.path.join(tmp, 'empty.pdbqt')
        open(empty, 'w').close()
        assert docking.read_vina_affinity(empty) is None


def test_docking_refuses_to_run_without_vina():
    saved = os.environ.get('PATH', '')
    os.environ['PATH'] = ''
    try:
        raised = False
        try:
            docking.run_vina_docking('r.pdbqt', SAMPLE, (0, 0, 0), (20, 20, 20),
                                     output_dir=tempfile.mkdtemp(), verbose=False)
        except FileNotFoundError:
            raised = True
        assert raised, "expected FileNotFoundError when vina is missing"
    finally:
        os.environ['PATH'] = saved


def test_admet_ids_match_library_ids():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, 'admet.csv')
        df = admet.assess_admet_properties(SAMPLE, out, verbose=False)
    valid = [cid for cid, s in read_smiles(SAMPLE) if Chem.MolFromSmiles(s) is not None]
    assert list(df['compound_id']) == valid


def _toy_inputs():
    docking_df = pd.DataFrame({'compound_id': [0, 1],
                               'binding_energy_kcal_mol': [-9.0, -7.0]})
    admet_df = pd.DataFrame({'compound_id': [0, 1], 'SMILES': ['C', 'CC'],
                             'LogP': [1.0, 2.0], 'RO5_Violations': [0, 0],
                             'AromaRings': [1, 0], 'RotBonds': [0, 1]})
    return docking_df, admet_df


def test_ranking_drops_dynamics_terms_without_dynamics():
    docking_df, admet_df = _toy_inputs()
    with tempfile.TemporaryDirectory() as tmp:
        df = ranking.rank_candidates_mcda(docking_df, admet_df, None,
                                          output_csv=os.path.join(tmp, 'r.csv'), verbose=False)
    assert 'stability_score' not in df.columns
    assert 'flexibility_score' not in df.columns
    assert df['composite_score'].between(0, 1).all()


def test_ranking_requires_full_dynamics_coverage():
    docking_df, admet_df = _toy_inputs()
    dyn = {0: {'rmsd': 1.0, 'pocket_msf': 0.8}}  # compound 1 missing
    try:
        ranking.rank_candidates_mcda(docking_df, admet_df, dyn,
                                     output_csv=os.devnull, verbose=False)
    except ValueError as err:
        assert 'missing' in str(err)
    else:
        raise AssertionError("expected ValueError for incomplete dynamics")


def test_ranking_rejects_empty_docking():
    _, admet_df = _toy_inputs()
    try:
        ranking.rank_candidates_mcda(pd.DataFrame(columns=['compound_id', 'binding_energy_kcal_mol']),
                                     admet_df, output_csv=os.devnull, verbose=False)
    except ValueError:
        return
    raise AssertionError("expected ValueError for empty docking results")


if __name__ == '__main__':
    tests = [v for k, v in dict(globals()).items() if k.startswith('test_')]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\n{len(tests)} checks passed.")
