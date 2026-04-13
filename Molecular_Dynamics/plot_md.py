#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MD Trajectory Plotting Script
==============================
Generates interactive Plotly plots for all MD analyses:
  1. Ligand RMSD      → RMSD_PLOTS/ligand.html / .png / .svg
  2. Complex RMSD     → RMSD_PLOTS/complex.html / .png / .svg
  3. Active Site RMSD → RMSD_PLOTS/AS.html / .png / .svg
  4. SASA             → SASA_PLOTS/SASA.html / .png / .svg
  5. Radius of Gyration → RG_PLOTS/Rg.html / .png / .svg
  6. H-Bonds          → HBOND_PLOTS/HBonds.html / .png / .svg
  7. RMSF             → RMSF_PLOTS/RMSF.html / .png / .svg

Usage:
    python plot_md.py

Requirements:
    pip install plotly numpy kaleido
"""

from __future__ import division
import os
import sys
import subprocess
import numpy as np

# Auto-install missing packages
for pkg, import_name in [("plotly", "plotly"), ("numpy", "numpy")]:
    try:
        __import__(import_name)
    except ImportError:
        print(f"[*] Installing {pkg}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])

try:
    import plotly.graph_objs as go
    import plotly.io as pio
except ImportError:
    from plotly import graph_objs as go
    pio = None

# ============================================================
# USER SETTINGS — EDIT HERE
# ============================================================

BASE_PATH = "."          # Root folder containing one subfolder per compound
                         # e.g. ./Compound_A/rmsd_to_first_ligand.dat

COMPOUNDS = [            # Subfolder names inside BASE_PATH
    "Compound_A",
    "Compound_B",
    "Compound_C",
    "Compound_D",
    "Compound_E",
]

TOTAL_NS = 100.0         # Total simulation time in nanoseconds

# RMSF residue range filter
RMSF_RES_MIN = 0         # First residue to include in RMSF plot
RMSF_RES_MAX = 300       # Last residue to include in RMSF plot

# Plot colours — one per compound (add more if needed)
COLORS = [
    "#1f77b4",   # blue
    "#ff7f0e",   # orange
    "#2ca02c",   # green
    "#d62728",   # red
    "#9467bd",   # purple
    "#8c564b",   # brown
    "#e377c2",   # pink
    "#7f7f7f",   # grey
    "#bcbd22",   # yellow-green
    "#17becf",   # cyan
]

# ============================================================
# PLOT CONFIGURATIONS
# Each entry: (filename, title, y_axis_label, tag, output_subdir,
#              smooth_window, target_frames, is_rmsf)
# ============================================================
PLOTS = [
    {
        "filename"      : "rmsd_to_first_ligand.dat",
        "title"         : "<b>Ligand RMSD (Å) vs Time (ns)</b>",
        "y_label"       : "<b>RMSD (Å)</b>",
        "tag"           : "ligand",
        "output_subdir" : "RMSD_PLOTS",
        "smooth_window" : 1,
        "target_frames" : 8000,
        "is_rmsf"       : False,
    },
    {
        "filename"      : "rmsd_to_first_complex.dat",
        "title"         : "<b>Complex RMSD (Å) vs Time (ns)</b>",
        "y_label"       : "<b>RMSD (Å)</b>",
        "tag"           : "complex",
        "output_subdir" : "RMSD_PLOTS",
        "smooth_window" : 1,
        "target_frames" : 8000,
        "is_rmsf"       : False,
    },
    {
        "filename"      : "rmsd_to_first_AS_CA.dat",
        "title"         : "<b>Active Site RMSD (Å) vs Time (ns)</b>",
        "y_label"       : "<b>RMSD (Å)</b>",
        "tag"           : "AS",
        "output_subdir" : "RMSD_PLOTS",
        "smooth_window" : 1,
        "target_frames" : 1000,
        "is_rmsf"       : False,
    },
    {
        "filename"      : "sasa.dat",
        "title"         : "<b>Solvent Accessible Surface Area (SASA) vs Time (ns)</b>",
        "y_label"       : "<b>SASA (Å²)</b>",
        "tag"           : "SASA",
        "output_subdir" : "SASA_PLOTS",
        "smooth_window" : 1,
        "target_frames" : 1000,
        "is_rmsf"       : False,
    },
    {
        "filename"      : "radgyr.dat",
        "title"         : "<b>Radius of Gyration (Rg) vs Time (ns)</b>",
        "y_label"       : "<b>Radius of Gyration (Å)</b>",
        "tag"           : "Rg",
        "output_subdir" : "RG_PLOTS",
        "smooth_window" : 3,
        "target_frames" : 1000,
        "is_rmsf"       : False,
    },
    {
        "filename"      : "hbonds.dat",
        "title"         : "<b>Hydrogen Bonds vs Time (ns)</b>",
        "y_label"       : "<b>Number of Hydrogen Bonds</b>",
        "tag"           : "HBonds",
        "output_subdir" : "HBOND_PLOTS",
        "smooth_window" : 5,
        "target_frames" : 1000,
        "is_rmsf"       : False,
    },
    {
        "filename"      : "rmsf.dat",
        "title"         : "<b>RMSF vs Residue Number</b>",
        "y_label"       : "<b>RMSF (Å)</b>",
        "tag"           : "RMSF",
        "output_subdir" : "RMSF_PLOTS",
        "smooth_window" : 5,
        "target_frames" : None,   # no downsampling for RMSF (residue-based)
        "is_rmsf"       : True,
    },
]

# ============================================================
# HELPERS
# ============================================================

def read_dat(path, max_frame=None):
    """Read a cpptraj .dat file → (frame/residue array, value array)."""
    t, y = [], []
    if not os.path.isfile(path):
        return np.array(t), np.array(y)
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("@"):
                continue
            p = s.split()
            if len(p) < 2:
                continue
            try:
                col1 = int(float(p[0]))
                col2 = float(p[1])
            except ValueError:
                continue
            if max_frame and col1 > max_frame:
                break
            t.append(col1)
            y.append(col2)
    return np.array(t), np.array(y)


def read_rmsf(path):
    """Read RMSF .dat, filter to RMSF_RES_MIN–RMSF_RES_MAX."""
    residues, values = [], []
    if not os.path.isfile(path):
        return np.array(residues), np.array(values)
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("@"):
                continue
            p = s.split()
            if len(p) < 2:
                continue
            try:
                res = int(float(p[0]))
                val = float(p[1])
            except ValueError:
                continue
            if RMSF_RES_MIN <= res <= RMSF_RES_MAX:
                residues.append(res)
                values.append(val)
    return np.array(residues), np.array(values)


def moving_average(y, window):
    if window <= 1 or len(y) < window:
        return y
    return np.convolve(y, np.ones(window) / float(window), mode="valid")


def downsample(t, y, target):
    if target is None or len(y) <= target:
        return t, y
    idx = np.linspace(0, len(y) - 1, target).astype(int)
    return t[idx], y[idx]


def build_figure(title, x_label, y_label):
    axis_style = dict(
        ticks="inside",
        ticklen=14,
        tickwidth=3,
        showline=True,
        mirror="all",
        linewidth=2,
        linecolor="black",
        showgrid=False,
    )
    fig = go.Figure()
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=40)),
        font=dict(family="Helvetica", size=20, color="black"),
        xaxis=dict(title=x_label, **axis_style),
        yaxis=dict(title=y_label, **axis_style),
        width=1000,
        height=650,
        template="plotly_white",
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-0.18,
            font=dict(size=20),
        ),
        margin=dict(l=70, r=40, t=80, b=80),
    )
    return fig


def save_figure(fig, out_dir, tag):
    os.makedirs(out_dir, exist_ok=True)

    html_path = os.path.join(out_dir, f"{tag}.html")
    fig.write_html(html_path)
    print(f"  [ok] HTML  → {html_path}")

    if pio:
        try:
            pio.write_image(fig, os.path.join(out_dir, f"{tag}.png"),
                            width=1920, height=1080, scale=5)
            pio.write_image(fig, os.path.join(out_dir, f"{tag}.svg"))
            print(f"  [ok] PNG + SVG → {out_dir}/")
        except Exception as e:
            print(f"  [info] Static export skipped (install kaleido): {e}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("  MD Trajectory Plotting Pipeline")
    print(f"  Base path : {os.path.abspath(BASE_PATH)}")
    print(f"  Compounds : {COMPOUNDS}")
    print(f"  Plots     : {len(PLOTS)}")
    print("=" * 60)
    print()

    for plot in PLOTS:
        tag        = plot["tag"]
        title      = plot["title"]
        y_label    = plot["y_label"]
        filename   = plot["filename"]
        out_dir    = os.path.join(BASE_PATH, plot["output_subdir"])
        smooth_win = plot["smooth_window"]
        target_fr  = plot["target_frames"]
        is_rmsf    = plot["is_rmsf"]

        x_label = "<b>Residue Number</b>" if is_rmsf else "<b>Time (ns)</b>"
        fig = build_figure(title, x_label, y_label)

        print(f"[{tag}] {title.replace('<b>','').replace('</b>','').strip()}")

        traces_added = 0
        for i, compound in enumerate(COMPOUNDS):
            color = COLORS[i % len(COLORS)]
            path  = os.path.join(BASE_PATH, compound, filename)

            if is_rmsf:
                x, y = read_rmsf(path)
            else:
                x, y = read_dat(path)

            if x.size == 0:
                print(f"  [!] Missing or empty: {path}")
                continue

            # Smooth
            y_sm = moving_average(y, smooth_win)
            x_sm = x[:len(y_sm)]

            if is_rmsf:
                x_plot = x_sm
            else:
                # Downsample
                x_ds, y_sm = downsample(x_sm, y_sm, target_fr)
                # Normalise frames → nanoseconds
                x_plot = (x_ds / float(x_ds.max())) * TOTAL_NS

            fig.add_trace(go.Scatter(
                x=x_plot,
                y=y_sm,
                mode="lines",
                name=compound,
                line=dict(color=color, width=2.5),
            ))
            traces_added += 1

        if traces_added == 0:
            print(f"  [!] No data found for any compound — skipping {tag}")
            print()
            continue

        save_figure(fig, out_dir, tag)
        print()

    print("=" * 60)
    print("  All plots complete.")
    print(f"  Output folders inside: {os.path.abspath(BASE_PATH)}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
