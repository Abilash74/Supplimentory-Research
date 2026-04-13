"""
CNN_VS Score Extractor
======================
Reads all docked .sdf or .sdf.gz files from a directory,
extracts minimizedAffinity, CNNscore, CNNaffinity per pose,
calculates CNN_VS = CNNscore × CNNaffinity,
writes a .log file per compound AND a combined Excel summary.

Usage:
    python cnn_vs_extractor.py

Requirements:
    pip install rdkit pandas openpyxl
"""

import os
import gzip
import sys
import subprocess

# ──────────────────────────────────────────────────────────────
# CONFIGURATION — edit these
# ──────────────────────────────────────────────────────────────
SDF_DIR         = "docking_results"              # Folder with *_docked.sdf or *_docked.sdf.gz
OUTPUT_LOG_DIR  = "log_files_with_CNN_VS"        # Per-compound .log files go here
SUMMARY_EXCEL   = "CNN_VS_summary.xlsx"          # Combined results Excel
# ──────────────────────────────────────────────────────────────


def ensure_packages():
    for pkg in ["rdkit", "pandas", "openpyxl"]:
        try:
            __import__(pkg)
        except ImportError:
            print(f"[*] Installing {pkg}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])


def is_gzipped(path: str) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def open_sdf(sdf_path: str):
    """Return an RDKit supplier for .sdf or .sdf.gz (real or fake gz)."""
    from rdkit import Chem
    if sdf_path.endswith(".sdf.gz") and is_gzipped(sdf_path):
        fh = gzip.open(sdf_path, "rb")
        return Chem.ForwardSDMolSupplier(fh), fh
    else:
        # Plain SDF (including misnamed .sdf.gz)
        return Chem.SDMolSupplier(sdf_path), None


def process_file(sdf_path: str, compound_name: str, log_dir: str) -> list:
    """
    Parse one docked SDF, write a .log file, return list of pose dicts.
    """
    compound_folder = os.path.join(log_dir, compound_name)
    os.makedirs(compound_folder, exist_ok=True)
    log_path = os.path.join(compound_folder, f"{compound_name}.log")

    suppl, fh = open_sdf(sdf_path)
    pose_count = 0
    pose_rows  = []

    try:
        with open(log_path, "w") as log:
            log.write(f"Compound: {compound_name}\n")
            log.write(
                f"{'Pose':<6}{'MinimizedAffinity':<20}"
                f"{'CNNscore':<15}{'CNNaffinity':<15}{'CNN_VS':<15}\n"
            )
            log.write("=" * 75 + "\n")

            for mol in suppl:
                if mol is None:
                    continue
                pose_count += 1

                min_aff   = mol.GetProp("minimizedAffinity") if mol.HasProp("minimizedAffinity") else "N/A"
                cnn_score = mol.GetProp("CNNscore")          if mol.HasProp("CNNscore")          else "N/A"
                cnn_aff   = mol.GetProp("CNNaffinity")       if mol.HasProp("CNNaffinity")       else "N/A"

                try:
                    cnn_vs = round(float(cnn_score) * float(cnn_aff), 4)
                except Exception:
                    cnn_vs = "N/A"

                log.write(
                    f"{pose_count:<6}{min_aff:<20}"
                    f"{cnn_score:<15}{cnn_aff:<15}{str(cnn_vs):<15}\n"
                )

                pose_rows.append({
                    "Compound":           compound_name,
                    "Pose":               pose_count,
                    "MinimizedAffinity":  min_aff,
                    "CNNscore":           cnn_score,
                    "CNNaffinity":        cnn_aff,
                    "CNN_VS":             cnn_vs,
                })

            if pose_count == 0:
                log.write("No valid docking poses found.\n")

    finally:
        if fh:
            fh.close()

    return pose_rows


def write_summary_excel(all_rows: list, out_path: str):
    import pandas as pd

    df = pd.DataFrame(all_rows, columns=[
        "Compound", "Pose",
        "MinimizedAffinity", "CNNscore", "CNNaffinity", "CNN_VS"
    ])

    # Numeric conversion for sorting
    for col in ["MinimizedAffinity", "CNNscore", "CNNaffinity", "CNN_VS"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Best-pose summary sheet (best CNN_VS per compound)
    best = (
        df.sort_values("CNN_VS", ascending=False)
          .groupby("Compound", as_index=False)
          .first()
          .sort_values("CNN_VS", ascending=False)
    )

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df.to_excel(writer,   sheet_name="All Poses",    index=False)
        best.to_excel(writer, sheet_name="Best Per Compound", index=False)

    print(f"[✓] Summary Excel saved: {out_path}")
    print(f"    • All Poses sheet      : {len(df)} rows")
    print(f"    • Best Per Compound    : {len(best)} compounds")


def main():
    ensure_packages()

    if not os.path.isdir(SDF_DIR):
        print(f"[✗] SDF directory not found: {SDF_DIR}")
        sys.exit(1)

    os.makedirs(OUTPUT_LOG_DIR, exist_ok=True)

    # Collect .sdf and .sdf.gz files
    files = sorted([
        f for f in os.listdir(SDF_DIR)
        if f.endswith(".sdf") or f.endswith(".sdf.gz")
    ])

    if not files:
        print(f"[!] No .sdf or .sdf.gz files found in: {SDF_DIR}")
        sys.exit(0)

    print(f"[✓] Found {len(files)} SDF file(s) in '{SDF_DIR}'\n")

    all_rows  = []
    skipped   = []

    for i, fname in enumerate(files, 1):
        compound_name = fname.replace("_docked.sdf.gz", "").replace("_docked.sdf", "")
        sdf_path      = os.path.join(SDF_DIR, fname)
        print(f"[{i}/{len(files)}] {compound_name}", end=" ... ", flush=True)

        try:
            rows = process_file(sdf_path, compound_name, OUTPUT_LOG_DIR)
            all_rows.extend(rows)
            print(f"{len(rows)} pose(s)")
        except Exception as e:
            print(f"SKIPPED — {e}")
            skipped.append((fname, str(e)))

    print(f"\n{'='*55}")
    print(f"  Processed : {len(files) - len(skipped)}/{len(files)} files")
    if skipped:
        print(f"  Skipped   : {len(skipped)}")
        for fname, err in skipped:
            print(f"    • {fname}: {err}")
    print(f"{'='*55}\n")

    if all_rows:
        write_summary_excel(all_rows, SUMMARY_EXCEL)

    print(f"\n  Log files : ./{OUTPUT_LOG_DIR}/")
    print(f"  Summary   : ./{SUMMARY_EXCEL}\n")


if __name__ == "__main__":
    main()