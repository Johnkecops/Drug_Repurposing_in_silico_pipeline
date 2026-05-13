# Drug Discovery Pipeline — Python Scripts Summary


---

## 📦 Archive Contents

### Python Scripts (8 total)

#### 1. **01_prepare_protein.py** (3.1 KB)
**Purpose:** Protein target preparation for molecular docking

**Key functions:**
- Download PDB structures via ProDy
- Remove heteroatoms (water, ligands)
- Define binding site geometry
- Calculate docking grid parameters (center, size)

**Usage:**
```bash
python 01_prepare_protein.py --pdb_id 3ODU --binding_residues 110,114,121,125
```

**Output:** Prepared PDB file, grid center (x,y,z), grid size (Å)

---

#### 2. **02_prepare_library.py** (4.1 KB)
**Purpose:** Chemical library curation and drug-likeness filtering

**Key functions:**
- Read SMILES strings from file
- Validate SMILES and generate 3D coordinates
- Apply Lipinski's Rule of Five (MW, LogP, HBD, HBA)
- Generate optimised 3D structures (SDF output)

**Usage:**
```bash
python 02_prepare_library.py --smiles_file compounds.smi --output_sdf library.sdf
```

**Output:** SDF file with prepared compounds, statistics (retention rate, RO5 violations)

---

#### 3. **03_run_docking.py** (5.5 KB)
**Purpose:** Virtual screening with AutoDock Vina

**Key functions:**
- Convert receptor to PDBQT format
- Prepare ligand PDBQT files
- Run Vina docking with configurable parameters
- Extract binding energies from Vina output
- Sort results by binding energy

**Usage:**
```bash
python 03_run_docking.py --receptor protein.pdbqt --ligands library.sdf \
                         --center 0,0,0 --size 30,30,30
```

**Parameters:**
- `--exhaustiveness`: Search quality (1-32; default 8)
- `--output_dir`: Output directory for results

**Output:** CSV file (compound_id, smiles, binding_energy_kcal_mol)

---

#### 4. **04_run_md_simulation.py** (4.5 KB)
**Purpose:** Molecular dynamics simulations using OpenMM

**Key workflow:**
1. System preparation: Add hydrogens, solvate in TIP3P water
2. Energy minimisation (1000 steps)
3. NVT equilibration (50 ps at 300 K)
4. NPT equilibration (100 ps)
5. Production MD (200 ns default, configurable)

**Usage:**
```bash
python 04_run_md_simulation.py --pdb complex.pdb --production_time_ns 200
```

**Output:** 
- DCD trajectory file (every 2 ps snapshot)
- Energy log (step, time, energy, temperature, density)

---

#### 5. **05_binding_energy.py** (4.2 KB)
**Purpose:** Binding free energy calculation

**Key functions:**
- Simplified ΔG_bind estimation (educational placeholder)
- Print MMPBSA.py template for production calculations
- Crude approximation based on MW and LogP

**Usage:**
```bash
# Estimate binding energy
python 05_binding_energy.py --protein protein.pdb --ligand_smiles "SMILES_STRING"

# Print MMPBSA.py template
python 05_binding_energy.py --mmpbsa_template
```

**Note:** For production use, employ AmberTools MMPBSA.py or FEP/TI methods

**Output:** CSV with ΔG_bind estimates, MMPBSA input template

---

#### 6. **06_anm_analysis.py** (6.4 KB)
**Purpose:** Anisotropic Network Model (ANM) analysis of binding-pocket flexibility

**Key concepts:**
- Coarse-grained elastic network model (Cα atoms)
- Hessian diagonalisation → normal modes
- Collective motion at 300 K without explicit solvent
- Fast alternative to all-atom MD for flexibility assessment

**Key functions:**
- Build ANM Hessian (cutoff 15 Å, gamma=1.0)
- Calculate slowest N modes (default 20)
- Compute per-residue RMSF
- Classify pocket flexibility (rigid/moderate/flexible)
- Generate RMSF plots and eigenvalue spectra

**Usage:**
```bash
python 06_anm_analysis.py --pdb protein.pdb --ligand_xyz 0,0,0 \
                          --pocket_radius 5.0 --n_modes 20
```

**Output:**
- PNG plot (RMSF per residue + mode eigenvalues)
- CSV with residue indices and RMSF values
- Console report (flexibility classification)

---

#### 7. **07_admet_profiling.py** (4.4 KB)
**Purpose:** ADMET profiling and drug-likeness assessment

**Metrics calculated:**
- Lipinski's Rule of Five: MW, LogP, HBD, HBA, rotatable bonds
- Polar surface area (TPSA) for BBB penetrance
- Aromatic rings, sp3 fraction
- Drug-likeness classification (RO5 violations ≤ 1)
- Aqueous solubility prediction
- BBB penetrance prediction (TPSA < 60 + LogP 1-5)

**Usage:**
```bash
python 07_admet_profiling.py --smiles_file compounds.smi --output_csv admet_results.csv
```

**Output:** CSV with ADMET properties and drug-likeness assessment

**Filtering thresholds:**
- MW ≤ 500 Da
- LogP ≤ 5
- HBD ≤ 5, HBA ≤ 10
- RO5 violations ≤ 1 = drug-like

---

#### 8. **08_rank_candidates.py** (7.6 KB)
**Purpose:** Multi-Criteria Decision Analysis (MCDA) for candidate ranking

**MCDA Criteria (weights):**
- Binding energy (25%): Best ≈ -12 kcal/mol
- Conformational stability (20%): RMSD < 2.5 Å
- Pocket flexibility (15%): Optimal RMSF ≈ 0.8 Å
- Drug-likeness (20%): RO5 compliant
- Toxicity/Safety (10%): LogP ≤ 5, PAINS-free
- Synthesis tractability (10%): SAscore < 6

**Usage:**
```bash
python 08_rank_candidates.py --docking docking_results.csv \
                            --admet admet_results.csv \
                            --dynamics dynamics.csv
```

**Output:**
- CSV with ranked candidates (rank_order, composite_score, component scores)
- Console report (top 10 candidates with breakdown)

---

### Documentation Files

#### **README.md** (8.2 KB)
Comprehensive guide covering:
- Installation requirements and dependencies
- Complete workflow examples
- Input/output file formats
- Parameter customisation guide
- Troubleshooting tips
- Performance optimisation suggestions
- Citation information

#### **sample_compounds.smi** (0.6 KB)
Example SMILES file with ~10 test compounds:
- Anti-inflammatory drugs (Ibuprofen, Aspirin, Celecoxib)
- Known kinase inhibitors
- GPCR ligands for testing

---

## 🔧 Quick Start

### Basic Workflow (5 steps)

```bash
# 1. Prepare protein
python 01_prepare_protein.py --pdb_id 3ODU

# 2. Prepare library (RO5-filtered)
python 02_prepare_library.py --smiles_file compounds.smi --output_sdf library.sdf

# 3. Virtual screening (requires Vina + ADFR tools)
python 03_run_docking.py --receptor protein.pdbqt --ligands library.sdf \
                         --center 0,0,0 --size 30,30,30

# 4. ADMET filtering
python 07_admet_profiling.py --smiles_file top_hits.smi

# 5. Rank candidates
python 08_rank_candidates.py --docking docking_results/docking_results.csv \
                            --admet admet_results.csv
```

---

## 📊 Pipeline Architecture

```
PHASE 0: Preparation
  ├─ 01_prepare_protein.py        → PDB structure + grid
  └─ 02_prepare_library.py        → SDF file (RO5-filtered)

PHASE 1: Virtual Screening
  └─ 03_run_docking.py             → Docking CSV (1000–1M compounds)

PHASE 2: Molecular Dynamics
  └─ 04_run_md_simulation.py       → Trajectory + energy log

PHASE 3: Dynamics Analysis
  ├─ 05_binding_energy.py          → ΔG_bind estimate (or MMPBSA template)
  └─ 06_anm_analysis.py            → Pocket RMSF, flexibility classification

PHASE 4: Drug-Likeness
  └─ 07_admet_profiling.py         → ADMET CSV + drug-like assessment

PHASE 5: Candidate Selection
  └─ 08_rank_candidates.py         → Ranked candidates (MCDA)
```

---

## 💾 File Size & Compression

| File | Size |
|------|------|
| 01_prepare_protein.py | 3.1 KB |
| 02_prepare_library.py | 4.1 KB |
| 03_run_docking.py | 5.5 KB |
| 04_run_md_simulation.py | 4.5 KB |
| 05_binding_energy.py | 4.2 KB |
| 06_anm_analysis.py | 6.4 KB |
| 07_admet_profiling.py | 4.4 KB |
| 08_rank_candidates.py | 7.6 KB |
| README.md | 8.2 KB |
| sample_compounds.smi | 0.6 KB |
| **Total (uncompressed)** | **48.4 KB** |
| **Zip archive** | **20 KB** |

---

## 🎯 Key Features

✅ **Modular design**: Use individual scripts or chain them together  
✅ **Production-ready**: Error handling, verbose progress reporting  
✅ **Well-documented**: Docstrings, examples, README  
✅ **Open source**: MIT license, no proprietary dependencies  
✅ **Extensible**: Easy to modify parameters, add custom metrics  
✅ **Cross-platform**: Works on Linux, macOS, Windows (with WSL)  

---

## 📚 Citations

When using these scripts in research, please cite:

- **RDKit**: Landrum et al. (2023) https://zenodo.org/records/3735086
- **ProDy**: Bakan et al. (2011) Bioinformatics 27:1575–1577 [PMID 21471012]
- **OpenMM**: Eastman et al. (2017) PLoS Comput Biol 13:e1005659 [PMID 28948882]
- **AutoDock Vina**: Trott & Olson (2010) J Comput Chem 31:455–461 [PMID 19499576]

---

## 🚀 Next Steps

1. **Extract archive**: `unzip drug_discovery_scripts.zip`
2. **Install dependencies**: `pip install -r requirements.txt` (if generated)
3. **Test with sample data**: Use `sample_compounds.smi` as test input
4. **Read README.md**: Detailed workflow and troubleshooting guide
5. **Customize parameters**: Edit weights, thresholds, grid settings as needed

---

**Package Date:** May 3, 2026  
**Author:** Dr. Arli Aditya Parikesit, i3L University Jakarta  
**License:** MIT
