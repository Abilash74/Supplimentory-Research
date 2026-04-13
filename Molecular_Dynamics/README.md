# MD Simulation & MM/GBSA Analysis Pipeline

A complete pipeline for **Molecular Dynamics (MD) post-processing** and **MM/GBSA binding free energy** calculations using AMBER/AmberTools.

---

## Scripts Overview

| Script | Purpose |
|---|---|
| `run_mmpbsa.sh` | Prepares stripped topology and trajectory files for MM/GBSA (Steps 1–6) |
| `run_mpi.sh` | Launches the actual MM/GBSA calculation with MPI parallelization |
| `run_analysis.sh` | Runs MD trajectory analyses: RMSD, RMSF, H-Bonds, Radius of Gyration, SASA |

---

## Requirements

- **AMBER / AmberTools** (v20 or later recommended)
  - `cpptraj` — trajectory processing and analysis
  - `parmed` — topology manipulation
  - `ante-MMPBSA.py` — topology splitting
  - `MMPBSA.py.MPI` — parallel MM/GBSA calculation
- **OpenMPI** — for MPI parallelization
- AMBER environment sourced: `source /usr/local/amber.sh`

---

## Input Files Required

Place all files in your working directory before running:

```
project/
├── solvated.parm7       ← Solvated system topology (AMBER parm7 format)
├── prod.nc              ← Production MD trajectory (NetCDF format)
├── run_mmpbsa.sh        ← Step 1: topology/trajectory preparation
├── run_analysis.sh      ← MD trajectory analysis
└── MMPBSA/
    ├── mmpbsa.in        ← MMPBSA input parameters (you create this)
    └── run_mpi.sh       ← Step 2: launch MPI calculation
```

---

## Workflow

### Step A — MD Trajectory Analysis

Run this at any point after production MD is complete.

**1. Edit settings in `run_analysis.sh`:**

```bash
PARM="solvated.parm7"         # Your topology file
TRAJ="prod.nc"                 # Your trajectory file
PROTEIN_RESIDUES="1-300"       # Residue range of your protein
ACTIVE_SITE_CA="50,80,100"     # Active site residue numbers (CA atoms, comma-separated)
LIGAND="LIG"                   # Ligand residue name in topology
```

**2. Run:**

```bash
chmod +x run_analysis.sh
./run_analysis.sh
```

**3. Outputs:**

```
analysis/
├── rmsd/
│   ├── rmsd_to_first_complex.dat    ← Whole complex RMSD vs frame 1
│   ├── rmsd_to_first_ligand.dat     ← Ligand-only RMSD vs frame 1
│   └── rmsd_to_first_AS_CA.dat      ← Active site CA RMSD vs frame 1
├── rmsf/
│   └── rmsf.dat                     ← Per-residue fluctuation
├── hbonds/
│   ├── hbonds.dat                   ← H-bond occupancy per frame
│   └── hbonds_avg.dat               ← Averaged H-bond statistics
├── radgyr/
│   └── radgyr.dat                   ← Radius of gyration over time
└── sasa/
    └── sasa.dat                     ← Solvent accessible surface area over time
```

---

### Step B — MM/GBSA Preparation

Strips solvent from topology and trajectory, then prepares mbondi2 topology files.

**1. Edit settings in `run_mmpbsa.sh`:**

```bash
SOLVATED_PARM="solvated.parm7"       # Your solvated topology
TRAJ="prod.nc"                         # Your production trajectory
LIGAND_MASK=":LIG"                     # Ligand residue mask
SOLVENT_MASK=":WAT,HOH,Na+,Cl-,CLA,POT"   # Solvent/ions to strip
```

**2. Run from your main project directory:**

```bash
chmod +x run_mmpbsa.sh
./run_mmpbsa.sh
```

**3. Outputs created:**

```
project/
├── complex_mb2.prmtop     ← Stripped complex with mbondi2 radii
├── rec_mb2.prmtop         ← Receptor-only topology with mbondi2
├── ligand_mb2.prmtop      ← Ligand-only topology with mbondi2
└── prod_stripped.nc       ← Solvent-stripped trajectory
```

---

### Step C — MM/GBSA Calculation (MPI)

**1. Create your `mmpbsa.in` file** inside the `MMPBSA/` subdirectory.

Example `mmpbsa.in` for GB calculation:

```
&general
  startframe=1, endframe=10000, interval=5,
  verbose=2,
/
&gb
  igb=2, saltcon=0.150,
/
&decomp
  idecomp=2, dec_verbose=3,
  print_res="within 6",
/
```

**2. Edit settings in `run_mpi.sh`:**

```bash
NCORES=10                          # Number of MPI cores (match your HPC allocation)
export AMBERHOME=/usr/local        # Path to your AMBER installation
export PATH="/home/user/openmpi/bin:$AMBERHOME/bin:$PATH"
```

**3. Run from inside the `MMPBSA/` directory:**

```bash
cd MMPBSA/
chmod +x run_mpi.sh
./run_mpi.sh
```

**4. Monitor progress:**

```bash
tail -f MMPBSA/mmpbsa_mpi.log
ps aux | grep MMPBSA
```

**5. Output files (when complete):**

```
MMPBSA/
├── FINAL_RESULTS_MMPBSA.dat     ← Summary: delta G, component energies
├── FINAL_RESULTS_MMPBSA.csv     ← Per-frame energy breakdown
└── DECOMP_MMPBSA.dat            ← Per-residue energy decomposition
```

---

## Understanding the Output

### RMSD (run_analysis.sh)
- **Complex RMSD** — overall structural stability of the protein-ligand complex
- **Ligand RMSD** — stability of ligand position within the binding site
- **Active site CA RMSD** — local structural changes at the binding pocket

A stable simulation typically shows RMSD convergence within 2–3 Å for proteins.

### RMSF (run_analysis.sh)
- Per-residue flexibility; high RMSF indicates flexible/disordered regions
- Useful for identifying loop regions and rigid binding site residues

### H-Bonds (run_analysis.sh)
- Protein-ligand hydrogen bonds per frame
- `hbonds_avg.dat` shows occupancy (%) — bonds with >30% occupancy are considered significant

### Radius of Gyration (run_analysis.sh)
- Measures protein compactness over the simulation
- Stable Rg indicates the protein maintains its folded structure

### SASA (run_analysis.sh)
- Solvent Accessible Surface Area — changes indicate conformational shifts or burial of residues

### MM/GBSA Results (run_mpi.sh)
- **DeltaG_bind** = MM_energy + GB_solvation + SA_nonpolar − TΔS (optional)
- **Per-residue decomposition** identifies key contributing residues to binding

---

## Troubleshooting

### `AMBERHOME is not set`
```bash
source /usr/local/amber.sh
# or wherever your AMBER is installed
```

### `cpptraj not found`
Ensure AMBER bin is in your PATH:
```bash
export PATH="$AMBERHOME/bin:$PATH"
```

### `Topology and trajectory mismatch`
The atom counts must match exactly. Ensure both were created from the same system. Re-run `run_mmpbsa.sh` to regenerate consistent files.

### `mmpbsa.in not found`
Create your `mmpbsa.in` file inside the `MMPBSA/` folder before running `run_mpi.sh`.

### MPI calculation fails immediately
Check the log:
```bash
tail -50 MMPBSA/mmpbsa_mpi.log
```
Common causes: wrong number of cores, topology/trajectory mismatch, or missing input parameters in `mmpbsa.in`.

### `ante-MMPBSA.py` does not write `com.prmtop`
This is expected when the input topology is already solvent-stripped. The script automatically copies `complex_mb2.prmtop` → `com.prmtop`.

---

## Full Pipeline Quick Reference

```bash
# 1. Source AMBER environment
source /usr/local/amber.sh

# 2. Edit USER SETTINGS in each script

# 3. MD trajectory analysis (can run independently)
chmod +x run_analysis.sh
./run_analysis.sh

# 4. Prepare MM/GBSA topologies
chmod +x run_mmpbsa.sh
./run_mmpbsa.sh

# 5. Create MMPBSA/ subdirectory and add mmpbsa.in
mkdir -p MMPBSA
cp run_mpi.sh MMPBSA/
# (create mmpbsa.in inside MMPBSA/)

# 6. Launch MPI calculation
cd MMPBSA/
chmod +x run_mpi.sh
./run_mpi.sh

# 7. Monitor
tail -f mmpbsa_mpi.log
```

---

## Citations

If you use this pipeline in your work, please cite the following:

---

### AMBER (Molecular Dynamics Engine)

> D.A. Case, H.M. Aktulga, K. Belfon, I.Y. Ben-Shalom, J.T. Berryman, S.R. Brozell,
> D.S. Cerutti, T.E. Cheatham III, G.A. Cisneros, V.W.D. Cruzeiro, T.A. Darden,
> R.E. Duke, N.H. Giese, H. Gohlke, A.W. Goetz, J. Harris, S. Izadi, S.A. Izmailov,
> K. Kasavajhala, M.C. Kaymak, E. King, A. Kovalenko, T. Kurtzman, T.S. Lee,
> S. LeGrand, P. Li, C. Lin, J. Liu, T. Luchko, R. Luo, M. Machado, V. Man,
> M. Manathunga, K.M. Merz, Y. Miao, O. Mikhailovskii, G. Monard, H. Nguyen,
> K.A. O'Hearn, A. Onufriev, F. Pan, S. Pantano, R. Qi, A. Rahnamoun, D.R. Roe,
> A. Roitberg, C. Sagui, S. Schott-Verdugo, A. Shajan, J. Shen, C.L. Simmerling,
> N.R. Skrynnikov, J. Smith, J. Swails, R.C. Walker, J. Wang, J. Wang, H. Wei,
> R.M. Wolf, X. Wu, Y. Xiong, Y. Xue, D.M. York, S. Zhao, and P.A. Kollman (2023).
> *AMBER 2023*. University of California, San Francisco.
> https://ambermd.org

---

### AmberTools (cpptraj, parmed, ante-MMPBSA.py)

> D.A. Case, H.M. Aktulga, K. Belfon, I.Y. Ben-Shalom, J.T. Berryman, S.R. Brozell,
> D.S. Cerutti, T.E. Cheatham III, G.A. Cisneros, V.W.D. Cruzeiro, T.A. Darden,
> R.E. Duke, N.H. Giese, H. Gohlke, A.W. Goetz, J. Harris, S. Izadi, S.A. Izmailov,
> K. Kasavajhala, M.C. Kaymak, E. King, A. Kovalenko, T. Kurtzman, T.S. Lee,
> S. LeGrand, P. Li, C. Lin, J. Liu, T. Luchko, R. Luo, M. Machado, V. Man,
> M. Manathunga, K.M. Merz, Y. Miao, O. Mikhailovskii, G. Monard, H. Nguyen,
> K.A. O'Hearn, A. Onufriev, F. Pan, S. Pantano, R. Qi, A. Rahnamoun, D.R. Roe,
> A. Roitberg, C. Sagui, S. Schott-Verdugo, A. Shajan, J. Shen, C.L. Simmerling,
> N.R. Skrynnikov, J. Smith, J. Swails, R.C. Walker, J. Wang, J. Wang, H. Wei,
> R.M. Wolf, X. Wu, Y. Xiong, Y. Xue, D.M. York, S. Zhao, and P.A. Kollman (2023).
> *AmberTools 2023*. University of California, San Francisco.
> https://ambermd.org

---

### cpptraj (Trajectory Analysis)

> Roe, D.R. & Cheatham, T.E. (2013).
> PTRAJ and CPPTRAJ: Software for Processing and Analysis of Molecular Dynamics Trajectory Data.
> *Journal of Chemical Theory and Computation*, 9(7), 3084–3095.
> https://doi.org/10.1021/ct400341p

---

### MMPBSA.py (MM/GBSA Free Energy Calculations)

> Miller, B.R., McGee, T.D., Swails, J.M., Homeyer, N., Gohlke, H., & Roitberg, A.E. (2012).
> MMPBSA.py: An Efficient Program for End-State Free Energy Calculations.
> *Journal of Chemical Theory and Computation*, 8(9), 3314–3321.
> https://doi.org/10.1021/ct300418h

---

### MM/PBSA Method (Theoretical Basis)

> Kollman, P.A., Massova, I., Reyes, C., Kuhn, B., Huo, S., Chong, L., Lee, M.,
> Lee, T., Duan, Y., Wang, W., Donini, O., Cieplak, P., Srinivasan, J., Case, D.A.,
> & Cheatham, T.E. (2000).
> Calculating Structures and Free Energies of Complex Molecules: Combining Molecular
> Mechanics and Continuum Models.
> *Accounts of Chemical Research*, 33(12), 889–897.
> https://doi.org/10.1021/ar000033j

> Srinivasan, J., Cheatham, T.E., Cieplak, P., Kollman, P.A., & Case, D.A. (1998).
> Continuum Solvent Studies of the Stability of DNA, RNA, and Phosphoramidate-DNA Helices.
> *Journal of the American Chemical Society*, 120(37), 9401–9409.
> https://doi.org/10.1021/ja981844+

---

### GB Solvation Model (igb=2, mbondi2 radii)

> Onufriev, A., Bashford, D., & Case, D.A. (2004).
> Exploring protein native states and large-scale conformational changes with a modified
> generalized Born model.
> *Proteins: Structure, Function, and Bioinformatics*, 55(2), 383–394.
> https://doi.org/10.1002/prot.20033

> Tsui, V. & Case, D.A. (2000).
> Theory and Applications of the Generalized Born Solvation Model in Macromolecular Simulations.
> *Biopolymers*, 56(4), 275–291.
> https://doi.org/10.1002/1097-0282(2000)56:4<275::AID-BIP10024>3.0.CO;2-E

---

### ParmEd (Topology Manipulation)

> Swails, J., Hernandez, C., Mobley, D.L., Nguyen, H., Wang, L.P., & Janowski, P. (2016).
> ParmEd: Cross-program parameter and topology file editor and molecular mechanical simulator engine.
> https://github.com/ParmEd/ParmEd

---

### RMSD / RMSF / SASA Analysis Methods

> Humphrey, W., Dalke, A., & Schulten, K. (1996).
> VMD: Visual Molecular Dynamics.
> *Journal of Molecular Graphics*, 14(1), 33–38.
> https://doi.org/10.1016/0263-7855(96)00018-5
> *(Reference for analysis methodology context)*

---

### GNINA 1.0 (Docking — for input complex structures)

> McNutt, A.T., Francoeur, P., Aggarwal, R., Masuda, T., Meli, R., Ragoza, M.,
> Sunseri, J., & Koes, D.R. (2021).
> GNINA 1.0: molecular docking with deep learning.
> *Journal of Cheminformatics*, 13, 43.
> https://doi.org/10.1186/s13321-021-00522-2

---

### RDKit (Ligand 3D Preparation)

> RDKit: Open-source cheminformatics.
> https://www.rdkit.org

---

### Open Babel (Format Conversion & Gasteiger Charges)

> O'Boyle, N.M., Banck, M., James, C.A., Morley, C., Vandermeersch, T., & Hutchison, G.R. (2011).
> Open Babel: An open chemical toolbox.
> *Journal of Cheminformatics*, 3, 33.
> https://doi.org/10.1186/1758-2946-3-33

---

## License

These scripts are provided for academic and research use.
AMBER is subject to its own licensing terms — see https://ambermd.org/GetAmber.php
