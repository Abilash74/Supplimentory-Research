# Docking-DNAGyrase
# Automated Molecular Docking Pipeline

A fully automated virtual screening pipeline using **GNINA 1.0** with CNN-based rescoring.  
Reads ligand SMILES from an Excel file, prepares 3D structures, runs docking, and extracts CNN_VS scores.

---

## Scripts

| Script | Purpose |
|---|---|
| `docking_pipeline_excel.py` | End-to-end ligand prep + GNINA docking for all compounds in Excel |
| `cnn_vs_extractor.py` | Post-processing: extracts CNN scores from docked SDF files, computes CNN_VS |

---

## Requirements

### Python packages
```bash
pip install rdkit pandas openpyxl
```

### System dependencies
```bash
# Ubuntu / Debian
sudo apt install openbabel

# Conda
conda install -c conda-forge openbabel
```

### GNINA binary
Download the latest release from:
```
https://github.com/gnina/gnina/releases
```
Place the binary in your working directory:
```bash
chmod +x gnina
```

---

## Input Files

Place all of the following in the same working directory before running:

```
your_project/
├── ligands.xlsx          ← Excel file with ligand SMILES
├── proteingyr.pdbqt      ← Prepared receptor (PDBQT format)
├── autoboxlig.pdbqt      ← Reference ligand to define the binding box
├── gnina                 ← GNINA binary (executable)
├── docking_pipeline_excel.py
└── cnn_vs_extractor.py
```

### Excel file format (`ligands.xlsx`)

The Excel file must have at minimum a **SMILES column**. A **Name column** is recommended:

| Name | SMILES |
|---|---|
| Ciprofloxacin | OB(O)c1ccc(cc1)-c1ccc2c(F)... |
| Compound_2 | CC(=O)Oc1ccccc1C(=O)O |
| Compound_3 | CN1CCC[C@H]1c2cccnc2 |

> Column names are configurable at the top of `docking_pipeline_excel.py`.

### Receptor preparation

Your receptor must be in **PDBQT format** with hydrogens and charges. Prepare using:
```bash
# Using AutoDockTools
prepare_receptor4.py -r protein.pdb -o proteingyr.pdbqt -A hydrogens

# Or using Open Babel
obabel protein.pdb -O proteingyr.pdbqt -xh --partialcharge gasteiger
```

### Autobox ligand

A reference ligand in PDBQT format that defines the binding site search space.  
This can be a known co-crystal ligand or any small molecule placed at the binding site:
```bash
obabel reference_ligand.sdf -O autoboxlig.pdbqt --partialcharge gasteiger
```

---

## Step 1 — Run Docking

### Configuration

Open `docking_pipeline_excel.py` and edit the **CONFIGURATION** block at the top:

```python
INPUT_EXCEL     = "ligands.xlsx"       # Your Excel file
SMILES_COL      = "SMILES"            # Column name containing SMILES strings
NAME_COL        = "Name"              # Column name for compound names (or None)
SHEET_NAME      = 0                   # Sheet index or name

RECEPTOR_PDBQT  = "proteingyr.pdbqt"
AUTOBOX_LIGAND  = "autoboxlig.pdbqt"
GNINA_BINARY    = "./gnina"

OUTPUT_DIR      = "docking_results"   # Where docked SDF files are saved
LIGAND_DIR      = "ligand_pdbqts"     # Intermediate PDB/PDBQT files

# Docking parameters
EXHAUSTIVENESS  = 32
NUM_MODES       = 10
FLEXDIST        = 3
SEED            = 0
```

### Run

```bash
python docking_pipeline_excel.py
```

### What it does (per ligand)

1. Reads SMILES from Excel
2. Generates a 3D conformer using RDKit (ETKDGv3 + MMFF optimization)
3. Converts PDB → PDBQT with Gasteiger partial charges via Open Babel
4. Runs GNINA flexible docking (flexible side chains within `FLEXDIST` Å of the ligand)
5. Decompresses the `.sdf.gz` output to a plain `.sdf` file
6. Writes a `docking_summary.xlsx` with best affinity, CNNscore, and CNNaffinity per compound

### Output structure after docking

```
your_project/
├── docking_results/
│   ├── Ciprofloxacin_docked.sdf
│   ├── Compound_2_docked.sdf
│   └── ...
├── ligand_pdbqts/
│   ├── Ciprofloxacin.pdb
│   ├── Ciprofloxacin.pdbqt
│   └── ...
└── docking_summary.xlsx
```

---

## Step 2 — CNN_VS Score Extraction

### Configuration

Open `cnn_vs_extractor.py` and edit the top:

```python
SDF_DIR        = "docking_results"          # Folder with *_docked.sdf files
OUTPUT_LOG_DIR = "log_files_with_CNN_VS"    # Per-compound .log files
SUMMARY_EXCEL  = "CNN_VS_summary.xlsx"      # Combined results Excel
```

### Run

```bash
python cnn_vs_extractor.py
```

### What it does

For every `.sdf` (or `.sdf.gz`) file in `docking_results/`:

1. Reads all docked poses
2. Extracts `minimizedAffinity`, `CNNscore`, and `CNNaffinity` per pose
3. Calculates **CNN_VS** = `CNNscore × CNNaffinity`
4. Writes a formatted `.log` file per compound
5. Compiles everything into `CNN_VS_summary.xlsx` with two sheets:
   - **All Poses** — every pose from every compound
   - **Best Per Compound** — highest CNN_VS pose per compound, sorted

### Output structure after extraction

```
your_project/
├── log_files_with_CNN_VS/
│   ├── Ciprofloxacin/
│   │   └── Ciprofloxacin.log
│   ├── Compound_2/
│   │   └── Compound_2.log
│   └── ...
└── CNN_VS_summary.xlsx
```

### Example `.log` file

```
Compound: Ciprofloxacin
Pose  MinimizedAffinity   CNNscore       CNNaffinity    CNN_VS
===========================================================================
1     -8.432              0.8821         7.341          6.4742
2     -7.981              0.7654         6.892          5.2779
3     -7.123              0.6201         6.105          3.7861
```

---

## Understanding the Scores

GNINA outputs three key scoring metrics per pose (McNutt et al., 2021):

| Score | Description | Range | Interpretation |
|---|---|---|---|
| **minimizedAffinity** | Vina-based binding affinity (kcal/mol) | Negative values | More negative = stronger binding |
| **CNNscore** | CNN pose quality score | 0–1 | Closer to 1 = better pose quality |
| **CNNaffinity** | CNN-predicted binding affinity (pK units) | Positive values | Higher = stronger predicted binding |
| **CNN_VS** | Virtual screening score = CNNscore × CNNaffinity | Derived | Higher = more likely to be an active |

> **CNN_VS** combines pose confidence with predicted affinity and is recommended for ranking compounds in virtual screening campaigns.

According to McNutt et al. (2021), a CNNscore > 0.8 indicates at least a **56% probability** (cross-docking) or **79% probability** (redocking) that the pose is within 2 Å RMSD of the true binding pose.

---

## Docking Parameters Explained

Based on the GNINA 1.0 paper (McNutt et al., 2021):

| Parameter | Default Used | Recommended Range | Notes |
|---|---|---|---|
| `exhaustiveness` | 32 | 8–64 | Higher = more thorough sampling. Use 32+ for virtual screening |
| `num_modes` | 10 | 9–100 | Number of output poses per ligand |
| `flexdist` | 3 | 3–3.5 Å | Distance threshold for flexible side chain selection |
| `seed` | 0 | Any integer | For reproducibility |

> For **whole-protein docking** (unknown binding site), use `exhaustiveness ≥ 32` as recommended by the paper.

---

## Troubleshooting

### `could not open receptor.pdbqt`
Ensure your receptor file is in the same directory and the filename in the config matches exactly.

### `obabel: command not found`
Install Open Babel: `sudo apt install openbabel` or `conda install -c conda-forge openbabel`

### `3D embedding failed`
Some SMILES strings may contain unusual valences or unsupported chemistries. Check the SMILES with RDKit:
```python
from rdkit import Chem
mol = Chem.MolFromSmiles("your_smiles_here")
print(mol)  # None = invalid SMILES
```

### `BadGzipFile` errors
The `cnn_vs_extractor.py` script automatically detects whether files are genuinely gzipped using magic number checking, and falls back to reading them as plain SDF. This handles GNINA versions that output uncompressed SDF even with a `.sdf.gz` extension.

### Skipped ligands in docking
Check the terminal output and `docking_summary.xlsx` — the `Error` column records the failure reason per compound.

---

## Full Pipeline — Quick Reference

```bash
# 1. Install dependencies
pip install rdkit pandas openpyxl
sudo apt install openbabel

# 2. Place files in working directory
#    ligands.xlsx, proteingyr.pdbqt, autoboxlig.pdbqt, ./gnina

# 3. Edit configuration blocks at the top of each script

# 4. Run docking
python docking_pipeline_excel.py

# 5. Extract CNN scores
python cnn_vs_extractor.py

# 6. Review results
#    docking_summary.xlsx     ← best scores per compound
#    CNN_VS_summary.xlsx      ← all poses + best per compound
#    log_files_with_CNN_VS/   ← detailed per-compound logs
```

---

## Citation

If you use this pipeline in your research, please cite:

**GNINA 1.0:**
> McNutt, A.T., Francoeur, P., Aggarwal, R., Masuda, T., Meli, R., Ragoza, M., Sunseri, J., & Koes, D.R. (2021).
> GNINA 1.0: molecular docking with deep learning.
> *Journal of Cheminformatics*, 13, 43.
> https://doi.org/10.1186/s13321-021-00522-2

**RDKit (ligand 3D generation):**
> RDKit: Open-source cheminformatics. http://www.rdkit.org

**Open Babel (format conversion & Gasteiger charges):**
> O'Boyle, N.M., Banck, M., James, C.A., Morley, C., Vandermeersch, T., & Hutchison, G.R. (2011).
> Open Babel: An open chemical toolbox.
> *Journal of Cheminformatics*, 3, 33.
> https://doi.org/10.1186/1758-2946-3-33

---

## License

These scripts are provided for academic and research use.
GNINA is released under GPL2/Apache License — see https://github.com/gnina/gnina for details.
