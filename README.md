# Drug Discovery Pipeline - Python Scripts

The scripts were developed for this following journal manuscript:

Dhea Priskila , Angelo Christiano Aouad Tomodok, Beatrice Valerie Basuki, Carissa Putri Pratama, Carolyn Nathaniel, Aishah Muhsin, Cheryl Yang, Daisy Jemima Bappedyanto, **Arli Aditya Parikesit**. 2026. Molecular Simulation-Based Drug Repurposing of Thioridazine and Losartan for Targeting CXCR4 and CCR5 Chemokine Receptors in Breast Cancer. _(in press/accepted)_. Journal of Science and Technology. UTHM Publisher, Malaysia.

Complete computational pipeline for structure-based virtual screening, lead characterisation, and optimisation.


## Overview

This package contains 8 modular Python scripts implementing a full drug discovery workflow:

1. **01_prepare_protein.py** — Protein target preparation for docking
2. **02_prepare_library.py** — Chemical library curation and filtering
3. **03_run_docking.py** — Virtual screening with AutoDock Vina
4. **04_run_md_simulation.py** — Molecular dynamics simulations (OpenMM)
5. **05_binding_energy.py** — Binding free energy calculation (MM-PBSA)
6. **06_anm_analysis.py** — Normal mode analysis (ProDy/ANM)
7. **07_admet_profiling.py** — ADMET and drug-likeness assessment
8. **08_rank_candidates.py** — Multi-criteria candidate ranking (MCDA)

## Installation Requirements

### Core Dependencies

```bash
# Docking and structure preparation
pip install autodock-vina openbabel-wheel rdkit

# Molecular dynamics
conda install -c conda-forge openmm mdanalysis

# Protein dynamics
pip install prody

# ADMET and descriptors
pip install pandas numpy matplotlib scipy scikit-learn
```

### System Requirements

- CPU: 4+ cores (8-16 recommended)
- RAM: 8-16 GB (for MD simulations)
- GPU: Optional (10-100x speedup for MD with CUDA/OpenCL)
- Storage: 10-50 GB (for trajectories from 100+ compounds)

### External Tools

- **AutoDock Vina** (http://vina.scripps.edu/)
- **ADFR Tools** (prepare_receptor4.py, prepare_ligand4.py)
- **AmberTools** (optional, for MMPBSA.py)

## Workflow

### Phase 0: Target and Library Preparation

```bash
# Prepare protein structure
python 01_prepare_protein.py --pdb_id 3ODU \
                             --binding_residues 110,114,121,125 \
                             --output_dir protein_prep

# Prepare ligand library (RO5-filtered)
python 02_prepare_library.py --smiles_file compounds.smi \
                            --output_sdf library_prepared.sdf
```

### Phase 1: Virtual Screening

```bash
# Run molecular docking (AutoDock Vina)
python 03_run_docking.py --receptor protein.pdbqt \
                         --ligands library_prepared.sdf \
                         --center 0,0,0 \
                         --size 30,30,30 \
                         --output_dir docking_results
```

### Phase 2: Molecular Dynamics Refinement

```bash
# Run MD simulations for top hits (e.g., top 50)
python 04_run_md_simulation.py --pdb complex.pdb \
                               --output_dir md_output \
                               --production_time_ns 200
```

### Phase 3: Binding Free Energy

```bash
# Estimate binding free energy from MD
python 05_binding_energy.py --protein protein.pdb \
                            --ligand_smiles "CC(C)Cc1ccc(cc1)C(C)C(O)=O" \
                            --output_csv binding_energy.csv

# Print MMPBSA.py template for full calculation
python 05_binding_energy.py --mmpbsa_template
```

### Phase 3b: Pocket Dynamics Analysis (ANM)

```bash
# Analyse binding-pocket flexibility
python 06_anm_analysis.py --pdb protein.pdb \
                          --ligand_xyz 0,0,0 \
                          --pocket_radius 5.0 \
                          --n_modes 20 \
                          --output_dir anm_results
```

### Phase 4: ADMET Assessment

```bash
# Evaluate ADMET properties
python 07_admet_profiling.py --smiles_file top_hits.smi \
                            --output_csv admet_results.csv
```

### Phase 5: Candidate Ranking

```bash
# Rank candidates via multi-criteria decision analysis
python 08_rank_candidates.py --docking docking_results/docking_results.csv \
                            --admet admet_results.csv \
                            --output ranked_candidates.csv
```

## Example Complete Workflow

```bash
# 1. Prepare protein (CXCR4)
python 01_prepare_protein.py --pdb_id 3ODU \
                             --binding_residues 110,114,121,125

# 2. Prepare library
python 02_prepare_library.py --smiles_file nci_compounds.smi \
                            --output_sdf library.sdf

# 3. Dock (requires prepare_receptor4.py and Vina)
python 03_run_docking.py --receptor 3ODU_prepared.pdb \
                         --ligands library.sdf \
                         --center 0,0,0 --size 30,30,30

# 4. MD simulations for top 50
for i in {0..49}; do
  python 04_run_md_simulation.py --pdb complex_${i}.pdb \
                                 --production_time_ns 200
done

# 5. ANM analysis
python 06_anm_analysis.py --pdb 3ODU_prepared.pdb \
                          --ligand_xyz 0,0,0

# 6. ADMET filtering
python 07_admet_profiling.py --smiles_file top_50_hits.smi

# 7. Rank candidates
python 08_rank_candidates.py --docking docking_results/docking_results.csv \
                            --admet admet_results.csv
```

## Input File Formats

### SMILES File (for 02, 07)
```
CC(C)Cc1ccc(cc1)C(C)C(O)=O
CC(=O)Oc1ccccc1C(=O)O
...
```

### PDB/PDBQT Files
Standard PDB format (protein structures)

### CSV Output Format

**docking_results.csv:**
```
compound_id,smiles,binding_energy_kcal_mol
0,SMILES_STRING,-9.5
1,SMILES_STRING,-8.8
...
```

**admet_results.csv:**
```
compound_id,MW_Da,LogP,HBD,HBA,RO5_Violations,Drug_Like
0,123.4,3.2,1,5,0,YES
1,234.5,5.1,2,8,1,NO
...
```

**ranked_candidates.csv:**
```
rank_order,compound_id,binding_energy_kcal_mol,composite_score
1,15,-9.2,0.78
2,32,-8.9,0.75
3,8,-8.7,0.72
...
```

## Key Parameters to Customize

### Docking Grid (03_run_docking.py)
- `--center`: Binding pocket center (x, y, z)
- `--size`: Grid box dimensions (default 30×30×30 Å)
- `--exhaustiveness`: Search quality (1-32; default 8)

### MD Simulation (04_run_md_simulation.py)
- `--production_time_ns`: Simulation length (default 200 ns)
- Force field: amber14-all.xml + tip3pfb.xml (changeable)
- Temperature: 300 K (physiological)

### ANM Analysis (06_anm_analysis.py)
- `--pocket_radius`: Distance from ligand center (default 5 Å)
- `--n_modes`: Number of slowest modes (default 20)
- `--cutoff`: ANM contact cutoff (default 15 Å)

### ADMET Thresholds (07_admet_profiling.py)
- Lipinski MW ≤ 500 Da
- Lipinski LogP ≤ 5
- Lipinski HBD ≤ 5, HBA ≤ 10
- TPSA < 60 Å² for BBB penetrance

### MCDA Weights (08_rank_candidates.py)
- Binding energy: 25%
- Stability: 20%
- Flexibility: 15%
- Drug-likeness: 20%
- Safety: 10%
- Synthesis: 10%

(Customise weights in script or via future command-line arguments)

## Output Files

Each script generates output files in specified directories:

- `protein_prep/`: Prepared protein structures
- `docking_results/`: Docking poses (PDBQT) and energy CSV
- `md_output/`: MD trajectories (DCD) and energy logs
- `anm_results/`: ANM RMSF plots and CSV
- `admet_results.csv`: ADMET properties
- `ranked_candidates.csv`: Final ranked list

## Performance Tips

1. **Parallelise docking**: Use batch submission on HPC for multiple compounds
2. **GPU acceleration**: Enable CUDA for OpenMM (10-100x faster MD)
3. **Screen reduction**: Apply ADMET filters early to reduce downstream calculations
4. **Check convergence**: Verify MD RMSD plateaus before proceeding
5. **Validate positive controls**: Include known actives to benchmark docking accuracy

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Docking fails | Missing ADFR tools | Install: `conda install -c bioconda mgltools-prepare_ligand4` |
| MD diverges | Poor equilibration | Increase NVT/NPT equilibration (100 ps → 1 ns) |
| Low docking success | Misaligned binding site | Validate pocket geometry; try blind docking |
| ANM-MD mismatch | Contact cutoff too short | Increase ANM cutoff (15 → 20 Å); analyse more modes |

## Citation

If you use this pipeline in your research, please cite:

- AutoDock Vina: Trott & Olson (2010) J Comput Chem 31:455-461
- ProDy: Bakan et al. (2011) Bioinformatics 27:1575-1577
- OpenMM: Eastman et al. (2017) PLoS Comput Biol 13:e1005659
- RDKit: Landrum & Co. (Zenodo DOI:10.5281/zenodo.591637)

## License

MIT License — See individual scripts for details.

## Contact

For questions or suggestions, contact:
Dr. Arli Aditya Parikesit
i3L University Jakarta, Indonesia
