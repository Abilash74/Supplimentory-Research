# MD Simulation & MM/GBSA Analysis Pipeline
## MD Simulation → MM/GBSA → Analysis & Plotting

A complete pipeline for molecular dynamics post-processing, binding free energy calculations, and publication-quality plotting.

---

## Scripts Overview

| Script | Language | Purpose |
|---|---|---|
| `run_mmpbsa.sh` | Bash | Strips trajectory/topology, applies mbondi2 radii (MM/GBSA prep) |
| `run_mpi.sh` | Bash | Launches MM/GBSA calculation with MPI parallelization |
| `run_analysis.sh` | Bash | Runs MD trajectory analyses: RMSD, RMSF, H-Bonds, Rg, SASA |
| `plot_md_full.py` | Python | Generates all MD + decomposition plots and text analysis report |

---

## Requirements

### System tools
- **AMBER / AmberTools** (v20 or later)
  - `cpptraj` — trajectory processing
  - `parmed` — topology manipulation
  - `ante-MMPBSA.py` — topology splitting
  - `MMPBSA.py.MPI` — parallel MM/GBSA
- **OpenMPI** — MPI parallelization
### Python packages
```bash
pip install pandas plotly numpy matplotlib kaleido
```

---

## Full Project Structure

```
project/
├── solvated.parm7                <- Solvated MD topology
├── prod.nc                       <- Production MD trajectory
│
├── run_analysis.sh               <- Step 1: MD analysis
├── run_mmpbsa.sh                 <- Step 2: MM/GBSA prep
│
├── MMPBSA/
│   ├── mmpbsa.in                 <- MM/GBSA input parameters
│   └── run_mpi.sh                <- Step 3: MPI calculation
│
└── plot_md_full.py               <- Step 4: All plots + decomp report
```

---

## Pipeline Workflow

---

### Step 1 — MD Trajectory Analysis (`run_analysis.sh`)

Runs cpptraj analyses on the production MD trajectory.

**Edit settings at the top of the script:**
```bash
PARM="solvated.parm7"
TRAJ="prod.nc"
PROTEIN_RESIDUES="1-300"          # Residue range for SASA and Rg
ACTIVE_SITE_CA="50,80,100,120"    # Active site residues for RMSD (CA atoms)
LIGAND="LIG"                       # Ligand residue name in topology
```

**Run:**
```bash
source /usr/local/amber.sh
chmod +x run_analysis.sh
./run_analysis.sh
```

**Outputs:**
```
analysis/
├── rmsd/
│   ├── rmsd_to_first_complex.dat    <- Whole complex RMSD vs frame 1
│   ├── rmsd_to_first_ligand.dat     <- Ligand RMSD vs frame 1
│   └── rmsd_to_first_AS_CA.dat      <- Active site CA RMSD vs frame 1
├── rmsf/
│   └── rmsf.dat                     <- Per-residue RMSF
├── hbonds/
│   ├── hbonds.dat                   <- H-bonds per frame
│   └── hbonds_avg.dat               <- Averaged H-bond statistics
├── radgyr/
│   └── radgyr.dat                   <- Radius of gyration over time
└── sasa/
    └── sasa.dat                     <- SASA over time
```

---

### Step 2 — MM/GBSA Preparation (`run_mmpbsa.sh`)

Strips solvent from topology and trajectory, then prepares mbondi2 topology files required for MM/GBSA.

**Edit settings:**
```bash
SOLVATED_PARM="solvated.parm7"
TRAJ="prod.nc"
LIGAND_MASK=":LIG"
SOLVENT_MASK=":WAT,HOH,Na+,Cl-,CLA,POT"
```

**Run from your main project directory:**
```bash
chmod +x run_mmpbsa.sh
./run_mmpbsa.sh
```

**Outputs:**
```
complex_mb2.prmtop    <- Stripped complex with mbondi2 radii
rec_mb2.prmtop        <- Receptor-only topology
ligand_mb2.prmtop     <- Ligand-only topology
prod_stripped.nc      <- Solvent-stripped trajectory
```

**Steps performed internally:**
1. Strip solvent from trajectory to `prod_stripped.nc`
2. Strip solvent from topology to `complex.prmtop`
3. Verify topology-trajectory compatibility
4. Apply mbondi2 radii to `complex_mb2.prmtop`
5. Split with `ante-MMPBSA.py` into receptor and ligand topologies
6. Enforce mbondi2 on all topology files

---

### Step 3 — MM/GBSA MPI Calculation (`run_mpi.sh`)

Launches the actual MM/GBSA binding free energy calculation using MPI parallelization.

**1. Create `mmpbsa.in` inside the `MMPBSA/` folder:**

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
NCORES=10
export AMBERHOME=/usr/local
export PATH="/home/user/openmpi/bin:$AMBERHOME/bin:$PATH"
```

**3. Run from inside `MMPBSA/`:**
```bash
cd MMPBSA/
chmod +x run_mpi.sh
./run_mpi.sh
```

**4. Monitor:**
```bash
tail -f MMPBSA/mmpbsa_mpi.log
ps aux | grep MMPBSA
```

**Outputs:**
```
MMPBSA/
├── FINAL_RESULTS_MMPBSA.dat     <- Summary delta G and energy components
├── FINAL_RESULTS_MMPBSA.csv     <- Per-frame energy breakdown
└── DECOMP_MMPBSA.dat            <- Per-residue energy decomposition
```

---

### Step 4 — Plotting & Analysis (`plot_md_full.py`)

Generates all MD trajectory plots and MM/GBSA decomposition figures in one run.

**Edit the configuration block:**
```python
# MD settings
MD_BASE_PATH = "."
COMPOUNDS    = ["Compound_A", "Compound_B", ...]
TOTAL_NS     = 100.0
RMSF_RES_MIN = 0
RMSF_RES_MAX = 300

# Decomposition settings
DECOMP_FILES = {
    "Compound_A": "Compound_A_FINAL_DECOMP_MMPBSA.dat",
    ...
}
DECOMP_SHORT = {"Compound_A": "CmpA", ...}   # Short labels for plots
DATA_START   = 8       # Line where residue data starts in .dat
DATA_END     = 203     # Line where residue data ends
TOP_N        = 15      # Top residues per compound for union set
DECOMP_OUT   = "DECOMP_PLOTS"
```

**Expected folder structure for MD data:**
```
MD_BASE_PATH/
├── Compound_A/
│   ├── rmsd_to_first_ligand.dat
│   ├── rmsd_to_first_complex.dat
│   ├── rmsd_to_first_AS_CA.dat
│   ├── sasa.dat
│   ├── radgyr.dat
│   ├── hbonds.dat
│   └── rmsf.dat
└── Compound_B/
    └── ...
```

**Run:**
```bash
python plot_md_full.py
```

**Outputs — MD trajectory plots (HTML + PNG + SVG):**
```
RMSD_PLOTS/ligand.*          <- Ligand RMSD vs time
RMSD_PLOTS/complex.*         <- Complex RMSD vs time
RMSD_PLOTS/AS.*              <- Active site RMSD vs time
SASA_PLOTS/SASA.*            <- SASA vs time
RG_PLOTS/Rg.*                <- Radius of gyration vs time
HBOND_PLOTS/HBonds.*         <- H-bonds vs time
RMSF_PLOTS/RMSF.*            <- RMSF vs residue number
```

**Outputs — MM/GBSA decomposition figures:**
```
DECOMP_PLOTS/
├── fig1_heatmaps.png            <- 4-panel heatmap (Total/vdW/NetElec/NP)
├── fig2_summary.png             <- Grand total + components + radar
├── fig3_grouped_bar.png         <- Top-15 grouped bar with error bars
├── interactive_heatmaps.html    <- Interactive 4-panel heatmap
├── interactive_grouped_bar.html <- Interactive grouped bar chart
├── interactive_radar.html       <- Interactive binding profile radar
├── interactive_summary_bar.html <- Interactive grand total bar
└── decomp_analysis_report.txt   <- Full text analysis report
```

**Decomposition report sections:**
1. Ligand self-energy contributions (vdW, Elec, Polar, NP per compound)
2. Grand total delta G ranking (most to least favorable)
3. Top 15 protein residues by |TOTAL| per compound with interaction type
4. Cross-compound comparison table of union top residues
5. Dominant interaction type per compound (vdW/hydrophobic vs electrostatic)

---

## Understanding the Output

### RMSD
- **Complex RMSD** — overall structural stability; convergence within 2-3 Angstrom indicates a stable simulation
- **Ligand RMSD** — stability of ligand position within the binding pocket
- **Active site CA RMSD** — local flexibility at the binding pocket

### RMSF
- High RMSF = flexible/disordered regions (loops, termini)
- Low RMSF at active site residues = rigid, well-defined binding pocket

### H-Bonds
- Occupancy >30% is considered significant
- `hbonds_avg.dat` reports average occupancy across the full trajectory

### Radius of Gyration
- Stable Rg indicates the protein maintains a compact, folded structure throughout the simulation

### SASA
- Decreasing SASA may indicate burial of hydrophobic residues
- Sudden changes can indicate conformational transitions

### MM/GBSA Energy Terms

| Term | Description |
|---|---|
| `VDWAALS` | van der Waals contribution |
| `EEL` | Electrostatic contribution |
| `EGB` | Generalized Born solvation |
| `ESURF` | Non-polar surface area solvation |
| `DELTA G binding` | Total estimated binding free energy |

### Decomposition Energy Terms

| Term | Meaning |
|---|---|
| `vdW` | van der Waals interaction |
| `Elec` | Electrostatic interaction |
| `Polar` | Polar solvation penalty |
| `NP` | Non-polar solvation |
| `TOTAL` | Sum of all components — use for ranking |

Residues with TOTAL < -0.05 kcal/mol are favorable; those > +0.05 kcal/mol are unfavorable.

---

## Full Pipeline Quick Reference

```bash
# Source AMBER environment
source /usr/local/amber.sh

# Step 1: MD trajectory analysis
chmod +x run_analysis.sh
./run_analysis.sh

# Step 2: MM/GBSA topology prep
chmod +x run_mmpbsa.sh
./run_mmpbsa.sh

# Step 3: MM/GBSA MPI calculation
mkdir -p MMPBSA
cp run_mpi.sh MMPBSA/
# Create mmpbsa.in inside MMPBSA/
cd MMPBSA/
chmod +x run_mpi.sh
./run_mpi.sh
tail -f mmpbsa_mpi.log
cd ..

# Step 4: All plots + decomp analysis
python plot_md_full.py
```

---

## Troubleshooting

**`AMBERHOME is not set`**
```bash
source /usr/local/amber.sh
```

**`cpptraj not found`**
```bash
export PATH="$AMBERHOME/bin:$PATH"
```

**Topology and trajectory mismatch** — Re-run `run_mmpbsa.sh` to regenerate consistent files.

**`mmpbsa.in not found`** — Create it inside the `MMPBSA/` folder before running `run_mpi.sh`.

**`ante-MMPBSA.py` does not write `com.prmtop`** — Expected for already-stripped topologies; the script copies `complex_mb2.prmtop` automatically.

**Decomposition plots show no data** — Check that `DATA_START` and `DATA_END` match the actual line range of the residue data block in your `.dat` file.

**PNG/SVG export fails** — Install kaleido: `pip install kaleido`. HTML plots still work without it.

---

## Citations

If you use this pipeline in your work, please cite the following:

---

### AMBER 2023 (Molecular Dynamics Engine)

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

### AmberTools 2023 (cpptraj, parmed, ante-MMPBSA.py)

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
> *Journal of Chemical Theory and Computation*, 9(7), 3084-3095.
> https://doi.org/10.1021/ct400341p

---

### MMPBSA.py (MM/GBSA Free Energy Calculations)

> Miller, B.R., McGee, T.D., Swails, J.M., Homeyer, N., Gohlke, H., & Roitberg, A.E. (2012).
> MMPBSA.py: An Efficient Program for End-State Free Energy Calculations.
> *Journal of Chemical Theory and Computation*, 8(9), 3314-3321.
> https://doi.org/10.1021/ct300418h

---

### MM/PBSA Method (Theoretical Basis)

> Kollman, P.A., Massova, I., Reyes, C., Kuhn, B., Huo, S., Chong, L., Lee, M.,
> Lee, T., Duan, Y., Wang, W., Donini, O., Cieplak, P., Srinivasan, J., Case, D.A.,
> & Cheatham, T.E. (2000).
> Calculating Structures and Free Energies of Complex Molecules: Combining Molecular
> Mechanics and Continuum Models.
> *Accounts of Chemical Research*, 33(12), 889-897.
> https://doi.org/10.1021/ar000033j

> Srinivasan, J., Cheatham, T.E., Cieplak, P., Kollman, P.A., & Case, D.A. (1998).
> Continuum Solvent Studies of the Stability of DNA, RNA, and Phosphoramidate-DNA Helices.
> *Journal of the American Chemical Society*, 120(37), 9401-9409.
> https://doi.org/10.1021/ja981844+

---

### GB Solvation Model (igb=2, mbondi2 radii)

> Onufriev, A., Bashford, D., & Case, D.A. (2004).
> Exploring protein native states and large-scale conformational changes with a modified
> generalized Born model.
> *Proteins: Structure, Function, and Bioinformatics*, 55(2), 383-394.
> https://doi.org/10.1002/prot.20033

> Tsui, V. & Case, D.A. (2000).
> Theory and Applications of the Generalized Born Solvation Model in Macromolecular Simulations.
> *Biopolymers*, 56(4), 275-291.
> https://doi.org/10.1002/1097-0282(2000)56:4<275::AID-BIP10024>3.0.CO;2-E

---

### ParmEd (Topology Manipulation)

> Swails, J., Hernandez, C., Mobley, D.L., Nguyen, H., Wang, L.P., & Janowski, P. (2016).
> ParmEd: Cross-program parameter and topology file editor and molecular mechanical simulator engine.
> https://github.com/ParmEd/ParmEd

---

### Plotly (Interactive Visualisation)

> Plotly Technologies Inc. (2015).
> Collaborative data science. Plotly Technologies Inc., Montreal, QC.
> https://plot.ly

---

### Matplotlib (Static Figures)

> Hunter, J.D. (2007).
> Matplotlib: A 2D graphics environment.
> *Computing in Science and Engineering*, 9(3), 90-95.
> https://doi.org/10.1109/MCSE.2007.55

---

### NumPy

> Harris, C.R., et al. (2020).
> Array programming with NumPy.
> *Nature*, 585, 357-362.
> https://doi.org/10.1038/s41586-020-2649-2

---

### pandas

> McKinney, W. (2010).
> Data Structures for Statistical Computing in Python.
> *Proceedings of the 9th Python in Science Conference*, 56-61.
> https://doi.org/10.25080/Majora-92bf1922-00a

---

## License

These scripts are provided for academic and research use.
AMBER is subject to its own licensing terms — see https://ambermd.org/GetAmber.php
