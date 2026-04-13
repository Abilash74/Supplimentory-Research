"""
MM-GBSA Per-Residue Decomposition Analysis — Multi-Ligand
==========================================================
Parses AMBER MMPBSA.py decomposition .dat files for multiple ligands,
extracts per-residue energy contributions, and prints:
  - Ligand self-energy summary
  - Grand total ΔG (protein + ligand contributions)
  - Top contributing residues per ligand
  - Shared hot-spot residue comparison table

Usage:
    python decomp_analysis.py

Dependencies: numpy, pandas
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ─── CONFIGURATION — edit these paths and names ──────────────────────────────
LIGAND_FILES = {
    "Ciprofloxacin":      "Ciprofloxacin_FINAL_DECOMP_MMPBSA.dat",
    "Diallyl Trisulfide": "Diallyl_trisulfide_FINAL_DECOMP_MMPBSA.dat",
    "Cinnamaldehyde":     "Cinnamaldehyde_DECOMP_MMPBSA.dat",
    "Eugenol":            "Eugenol_DECOMP_MMPBSA.dat",
    "Zingiberene":        "Zingiberene_DECOMP_MMPBSA.dat",
}
DATA_START = 8    # line index where residue data begins (after 2 header rows)
DATA_END   = 203  # line index where residue data ends
TOP_N      = 15   # number of top residues to report per ligand
SIG_THRESH = 0.05 # |TOTAL| threshold (kcal/mol) for "significant" residues

COLS = [
    "Residue", "Location",
    "Int_Avg",   "Int_Std",   "Int_SEM",
    "vdW_Avg",   "vdW_Std",   "vdW_SEM",
    "Elec_Avg",  "Elec_Std",  "Elec_SEM",
    "Polar_Avg", "Polar_Std", "Polar_SEM",
    "NP_Avg",    "NP_Std",    "NP_SEM",
    "TOTAL_Avg", "TOTAL_Std", "TOTAL_SEM",
]


# ─── PARSER ──────────────────────────────────────────────────────────────────
def parse_decomp(filepath: str, start: int, end: int) -> pd.DataFrame:
    """Parse GB Total Energy Decomposition block from a .dat file."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath) as fh:
        lines = fh.readlines()
    rows = []
    for line in lines[start:end]:
        line = line.strip()
        if not line or (not line[0].isalpha() and not line[0].isdigit()):
            break
        vals = line.split(",")
        if len(vals) >= 20:
            rows.append(vals[:20])
    df = pd.DataFrame(rows, columns=COLS)
    df["Residue"] = df["Residue"].str.strip()
    for col in df.columns[2:]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# ─── LOAD ALL LIGANDS ────────────────────────────────────────────────────────
all_data = {}
for name, fname in LIGAND_FILES.items():
    all_data[name] = parse_decomp(fname, DATA_START, DATA_END)

LIGANDS = list(all_data.keys())


# ─── HELPER ──────────────────────────────────────────────────────────────────
def get_ligand_row(df: pd.DataFrame) -> pd.Series:
    lig = df[df["Residue"].str.startswith("LIG")]
    if lig.empty:
        raise ValueError("No LIG residue found in decomp data.")
    return lig.iloc[0]


def get_protein(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df["Residue"].str.startswith("LIG")].copy()


# ─── 1. LIGAND SELF-ENERGY SUMMARY ──────────────────────────────────────────
sep  = "=" * 72
sep2 = "-" * 72

print(sep)
print("  MM-GBSA DECOMPOSITION ANALYSIS — MULTI-LIGAND")
print(sep)

print("\n[ 1. Ligand self-energy contributions (LIG 195) ]")
print(f"\n  {'Ligand':<22} {'TOTAL':>9} {'±Std':>7}  {'vdW':>9} {'Elec':>9} {'Polar':>9} {'NP':>9}")
print(f"  {'-'*72}")
for name in LIGANDS:
    lig = get_ligand_row(all_data[name])
    net_elec = lig["Elec_Avg"] + lig["Polar_Avg"]
    print(f"  {name:<22} {lig['TOTAL_Avg']:>9.3f} {lig['TOTAL_Std']:>7.3f}  "
          f"{lig['vdW_Avg']:>9.3f} {net_elec:>9.3f} "
          f"{lig['Polar_Avg']:>9.3f} {lig['NP_Avg']:>9.3f}")


# ─── 2. GRAND TOTAL & PROTEIN SUM ────────────────────────────────────────────
print(f"\n[ 2. Grand total ΔG (protein + ligand contributions) ]")
print(f"\n  {'Ligand':<22} {'Ligand ΔG':>10} {'Protein sum':>12} {'Grand sum':>10}  Rank")
print(f"  {'-'*60}")
grand_vals = {}
for name in LIGANDS:
    lig   = get_ligand_row(all_data[name])
    prot  = get_protein(all_data[name])
    ps    = prot["TOTAL_Avg"].sum()
    gs    = all_data[name]["TOTAL_Avg"].sum()
    grand_vals[name] = gs
    print(f"  {name:<22} {lig['TOTAL_Avg']:>10.3f} {ps:>12.3f} {gs:>10.3f}")

ranked = sorted(grand_vals.items(), key=lambda x: x[1])
print(f"\n  Ranking (most favorable → least):")
for i, (name, val) in enumerate(ranked, 1):
    print(f"    {i}. {name:<22}  ΔG = {val:.3f} kcal/mol")


# ─── 3. TOP CONTRIBUTING RESIDUES PER LIGAND ────────────────────────────────
print(f"\n[ 3. Top {TOP_N} protein residues by |TOTAL| per ligand ]")
for name in LIGANDS:
    prot = get_protein(all_data[name]).copy()
    prot["abs_total"] = prot["TOTAL_Avg"].abs()
    top = prot.nlargest(TOP_N, "abs_total")
    print(f"\n  {name}:")
    print(f"  {'Residue':<12} {'TOTAL':>8} {'±Std':>7}  {'vdW':>8} {'Elec':>8} {'Polar':>8} {'NP':>8}  Type")
    print(f"  {'-'*78}")
    for _, r in top.iterrows():
        net = r["Elec_Avg"] + r["Polar_Avg"]
        if r["TOTAL_Avg"] > SIG_THRESH:
            rtype = "UNFAV"
        elif abs(r["vdW_Avg"]) > abs(net) * 1.5:
            rtype = "hydrophobic"
        elif abs(net) > abs(r["vdW_Avg"]):
            rtype = "electrostatic"
        else:
            rtype = "mixed"
        flag = " *" if r["TOTAL_Avg"] < -SIG_THRESH else ("  !" if r["TOTAL_Avg"] > SIG_THRESH else "")
        print(f"  {r['Residue']:<12} {r['TOTAL_Avg']:>+8.4f} {r['TOTAL_Std']:>7.4f}  "
              f"{r['vdW_Avg']:>8.4f} {r['Elec_Avg']:>8.4f} "
              f"{r['Polar_Avg']:>8.4f} {r['NP_Avg']:>8.4f}  {rtype}{flag}")


# ─── 4. SHARED HOT-SPOT COMPARISON TABLE ─────────────────────────────────────
print(f"\n[ 4. Cross-ligand comparison — union of top-{TOP_N} residues ]")
union_top = set()
for name in LIGANDS:
    prot = get_protein(all_data[name]).copy()
    prot["abs_total"] = prot["TOTAL_Avg"].abs()
    union_top |= set(prot.nlargest(TOP_N, "abs_total")["Residue"].tolist())

all_res_ord = all_data["Ciprofloxacin"]["Residue"].tolist()
sorted_top  = [r for r in all_res_ord if r in union_top]

col_w = 10
header = f"  {'Residue':<12}" + "".join(f"{n[:col_w]:>{col_w}}" for n in LIGANDS)
print(f"\n{header}")
print(f"  {'-'*(12 + col_w*len(LIGANDS))}")
for res in sorted_top:
    row_str = f"  {res:<12}"
    for name in LIGANDS:
        df = all_data[name]
        v  = df[df["Residue"] == res]["TOTAL_Avg"].values
        val = v[0] if len(v) > 0 else 0.0
        marker = "*" if abs(val) > SIG_THRESH else " "
        row_str += f"{val:>+9.4f}{marker}"
    print(row_str)
print(f"\n  * = |TOTAL| > {SIG_THRESH} kcal/mol (significant)")


# ─── 5. ENERGY TYPE SUMMARY ──────────────────────────────────────────────────
print(f"\n[ 5. Dominant interaction type per ligand ]")
for name in LIGANDS:
    prot = get_protein(all_data[name]).copy()
    prot["net_elec"] = prot["Elec_Avg"] + prot["Polar_Avg"]
    prot["abs_total"] = prot["TOTAL_Avg"].abs()
    sig = prot[prot["abs_total"] > SIG_THRESH]

    sum_vdw  = prot["vdW_Avg"].sum()
    sum_elec = prot["net_elec"].sum()
    sum_np   = prot["NP_Avg"].sum()
    n_fav    = (sig["TOTAL_Avg"] < 0).sum()
    n_unfav  = (sig["TOTAL_Avg"] > 0).sum()

    dom = "vdW/hydrophobic" if abs(sum_vdw) > abs(sum_elec) else "electrostatic"
    print(f"\n  {name}:")
    print(f"    Protein vdW sum:         {sum_vdw:+.4f} kcal/mol")
    print(f"    Protein net-elec sum:    {sum_elec:+.4f} kcal/mol")
    print(f"    Protein non-polar sum:   {sum_np:+.4f} kcal/mol")
    print(f"    Significant residues:    {n_fav} favorable, {n_unfav} unfavorable")
    print(f"    Dominant mechanism:      {dom}")

print(f"\n{sep}")
print("  END OF REPORT")
print(sep)
