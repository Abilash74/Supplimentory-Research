"""
Automated Multi-Ligand Docking Pipeline — Excel SMILES Input
=============================================================
1. Reads ligand names + SMILES from an Excel file
2. Generates 3D conformers + MMFF optimization for each (RDKit)
3. Adds Gasteiger partial charges, converts to PDBQT (Open Babel)
4. Runs GNINA docking for every ligand
5. Saves docked poses as individual .sdf files
6. Writes a summary results Excel file

Usage:
    python docking_pipeline_excel.py

Requirements:
    pip install rdkit pandas openpyxl
    sudo apt install openbabel        # or: conda install -c conda-forge openbabel
    # GNINA binary: https://github.com/gnina/gnina/releases
"""

import os
import sys
import gzip
import shutil
import subprocess
from pathlib import Path

# ──────────────────────────────────────────────────────────────
# CONFIGURATION — edit these to match your setup
# ──────────────────────────────────────────────────────────────

# --- Input Excel file ---
INPUT_EXCEL     = "ligands.xlsx"        # Path to your Excel file
SMILES_COL      = "SMILES"             # Column header containing SMILES
NAME_COL        = "Name"               # Column header for ligand names
                                        # (set to None to auto-generate names)
SHEET_NAME      = 0                    # Sheet index (0) or name e.g. "Sheet1"

# --- Docking inputs ---
RECEPTOR_PDBQT  = "proteingyr.pdbqt"
AUTOBOX_LIGAND  = "autoboxlig.pdbqt"
GNINA_BINARY    = "./gnina"

# --- Output ---
OUTPUT_DIR      = "docking_results"    # Docked SDF files go here
LIGAND_DIR      = "ligand_pdbqts"      # Intermediate PDB/PDBQT files
RESULTS_EXCEL   = "docking_summary.xlsx"

# --- Docking parameters ---
EXHAUSTIVENESS  = 32
NUM_MODES       = 10
FLEXDIST        = 3
SEED            = 0

# ──────────────────────────────────────────────────────────────
# STEP 0 — Install missing packages
# ──────────────────────────────────────────────────────────────
def ensure_packages():
    for pkg, import_name in [("rdkit", "rdkit"), ("pandas", "pandas"), ("openpyxl", "openpyxl")]:
        try:
            __import__(import_name)
        except ImportError:
            print(f"[*] Installing {pkg}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])

# ──────────────────────────────────────────────────────────────
# STEP 1 — Read SMILES from Excel
# ──────────────────────────────────────────────────────────────
def load_smiles_from_excel(excel_path: str):
    import pandas as pd

    df = pd.read_excel(excel_path, sheet_name=SHEET_NAME, dtype=str)
    df.columns = df.columns.str.strip()

    if SMILES_COL not in df.columns:
        raise ValueError(
            f"Column '{SMILES_COL}' not found in Excel.\n"
            f"Available columns: {list(df.columns)}\n"
            f"Update SMILES_COL in the config section."
        )

    ligands = []
    for i, row in df.iterrows():
        smiles = str(row[SMILES_COL]).strip()
        if not smiles or smiles.lower() == "nan":
            continue

        if NAME_COL and NAME_COL in df.columns:
            name = str(row[NAME_COL]).strip()
            if not name or name.lower() == "nan":
                name = f"ligand_{i+1:04d}"
        else:
            name = f"ligand_{i+1:04d}"

        # Sanitize name for filesystem
        name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
        ligands.append((name, smiles))

    print(f"[✓] Loaded {len(ligands)} ligand(s) from {excel_path}")
    return ligands

# ──────────────────────────────────────────────────────────────
# STEP 2 — SMILES → 3D PDB (RDKit)
# ──────────────────────────────────────────────────────────────
def smiles_to_pdb(smiles: str, name: str, work_dir: str) -> str:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # ETKDGv3 → ETKDGv2 fallback
    for params in [AllChem.ETKDGv3(), AllChem.ETKDGv2()]:
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == 0:
            break
    else:
        raise RuntimeError("3D embedding failed for all methods")

    # MMFF → UFF fallback
    if AllChem.MMFFOptimizeMolecule(mol, maxIters=2000) == -1:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)

    out_pdb = os.path.join(work_dir, f"{name}.pdb")
    writer = Chem.PDBWriter(out_pdb)
    writer.write(mol)
    writer.close()
    return out_pdb

# ──────────────────────────────────────────────────────────────
# STEP 3 — PDB → PDBQT with Gasteiger charges (Open Babel)
# ──────────────────────────────────────────────────────────────
def pdb_to_pdbqt(pdb_path: str, pdbqt_path: str):
    cmd = ["obabel", pdb_path, "-h", "--partialcharge", "gasteiger", "-O", pdbqt_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.isfile(pdbqt_path):
        raise RuntimeError(f"obabel failed:\n{result.stderr.strip()}")

# ──────────────────────────────────────────────────────────────
# STEP 4 — GNINA docking
# ──────────────────────────────────────────────────────────────
def run_gnina(ligand_pdbqt: str, output_gz: str) -> str:
    cmd = [
        GNINA_BINARY,
        "-r",  RECEPTOR_PDBQT,
        "-l",  ligand_pdbqt,
        "--autobox_ligand",  AUTOBOX_LIGAND,
        "--exhaustiveness",  str(EXHAUSTIVENESS),
        "--flexdist_ligand", ligand_pdbqt,
        "--flexdist",        str(FLEXDIST),
        "--num_modes",       str(NUM_MODES),
        "--seed",            str(SEED),
        "--out",             output_gz,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"GNINA error:\n{result.stderr.strip()}")
    return result.stdout

# ──────────────────────────────────────────────────────────────
# STEP 5 — Decompress .sdf.gz → .sdf
# ──────────────────────────────────────────────────────────────
def decompress_gz(gz_path: str, out_path: str):
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(gz_path)

# ──────────────────────────────────────────────────────────────
# Parse best score from GNINA stdout
# ──────────────────────────────────────────────────────────────
def parse_best_score(gnina_log: str):
    affinity, cnn_score, cnn_affinity = "N/A", "N/A", "N/A"
    for line in gnina_log.splitlines():
        stripped = line.strip()
        if stripped.startswith("1 "):
            parts = stripped.split()
            if len(parts) >= 4:
                affinity     = parts[1]
                cnn_score    = parts[2]
                cnn_affinity = parts[3]
            break
    return affinity, cnn_score, cnn_affinity

# ──────────────────────────────────────────────────────────────
# Pre-flight checks
# ──────────────────────────────────────────────────────────────
def preflight_checks():
    errors = []
    if not os.path.isfile(INPUT_EXCEL):
        errors.append(f"Input Excel not found: {INPUT_EXCEL}")
    if not os.path.isfile(RECEPTOR_PDBQT):
        errors.append(f"Receptor not found: {RECEPTOR_PDBQT}")
    if not os.path.isfile(AUTOBOX_LIGAND):
        errors.append(f"Autobox ligand not found: {AUTOBOX_LIGAND}")
    if not os.path.isfile(GNINA_BINARY):
        errors.append(
            f"GNINA binary not found: '{GNINA_BINARY}'\n"
            "  → Download: https://github.com/gnina/gnina/releases"
        )
    if errors:
        print("\n[✗] Pre-flight failed:")
        for e in errors:
            print(f"    • {e}")
        sys.exit(1)
    print("[✓] Pre-flight checks passed.\n")

# ──────────────────────────────────────────────────────────────
# Write summary Excel
# ──────────────────────────────────────────────────────────────
def write_summary(results: list):
    import pandas as pd

    df = pd.DataFrame(results, columns=[
        "Name", "SMILES", "Status",
        "Affinity (kcal/mol)", "CNN Score", "CNN Affinity",
        "Output SDF", "Error"
    ])
    df.to_excel(RESULTS_EXCEL, index=False)
    print(f"[✓] Summary saved to: {RESULTS_EXCEL}")

# ──────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────
def main():
    ensure_packages()
    preflight_checks()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LIGAND_DIR, exist_ok=True)

    ligands = load_smiles_from_excel(INPUT_EXCEL)
    total   = len(ligands)
    results = []

    print(f"{'='*60}")
    print(f"  Starting pipeline for {total} ligand(s)")
    print(f"{'='*60}\n")

    for idx, (name, smiles) in enumerate(ligands, 1):
        print(f"[{idx}/{total}] {name}")
        print(f"        SMILES: {smiles[:60]}{'...' if len(smiles)>60 else ''}")
        row = [name, smiles, "FAILED", "N/A", "N/A", "N/A", "", ""]

        try:
            # 3D structure
            pdb_path = smiles_to_pdb(smiles, name, LIGAND_DIR)
            print(f"        3D conformer   → {pdb_path}")

            # Gasteiger PDBQT
            pdbqt_path = os.path.join(LIGAND_DIR, f"{name}.pdbqt")
            pdb_to_pdbqt(pdb_path, pdbqt_path)
            print(f"        PDBQT (Gast.)  → {pdbqt_path}")

            # Docking
            out_gz  = os.path.join(OUTPUT_DIR, f"{name}_docked.sdf.gz")
            out_sdf = os.path.join(OUTPUT_DIR, f"{name}_docked.sdf")
            log     = run_gnina(pdbqt_path, out_gz)

            # Decompress
            decompress_gz(out_gz, out_sdf)
            print(f"        Docked poses   → {out_sdf}")

            # Parse scores
            affinity, cnn_score, cnn_aff = parse_best_score(log)
            print(f"        Best affinity  = {affinity} kcal/mol  |  CNNscore = {cnn_score}  |  CNNaffinity = {cnn_aff}")

            row = [name, smiles, "SUCCESS", affinity, cnn_score, cnn_aff, out_sdf, ""]

        except Exception as e:
            row[-1] = str(e)
            print(f"        [✗] FAILED — {e}")

        results.append(row)
        print()

    # Summary
    succeeded = sum(1 for r in results if r[2] == "SUCCESS")
    print(f"{'='*60}")
    print(f"  Done — {succeeded}/{total} ligands docked successfully")
    print(f"{'='*60}\n")

    write_summary(results)
    print(f"\n  Docked SDF files : ./{OUTPUT_DIR}/")
    print(f"  Ligand PDBQTs    : ./{LIGAND_DIR}/")
    print(f"  Summary Excel    : ./{RESULTS_EXCEL}\n")


if __name__ == "__main__":
    main()