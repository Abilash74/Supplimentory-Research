"""
MM-GBSA Decomposition — Publication Plots (Plotly + Matplotlib)
================================================================
Generates 6 publication-quality figures:
  Fig 1  — 4-panel heatmap (Total/vdW/NetElec/NP)     [matplotlib PNG]
  Fig 2  — Summary panel: grand total + components + radar  [matplotlib PNG]
  Fig 3  — Top-15 grouped bar with error bars           [matplotlib PNG]
  HTML-A — Interactive heatmaps                         [Plotly HTML]
  HTML-B — Interactive grouped bar                      [Plotly HTML]
  HTML-C — Interactive radar chart                      [Plotly HTML]
  HTML-D — Interactive grand-total summary bar          [Plotly HTML]

Usage:
    python decomp_plots.py

Dependencies: numpy, pandas, matplotlib, plotly
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ─── CONFIGURATION ───────────────────────────────────────────────────────────
LIGAND_FILES = {
    "Ciprofloxacin":      "Ciprofloxacin_FINAL_DECOMP_MMPBSA.dat",
    "Diallyl Trisulfide": "Diallyl_trisulfide_FINAL_DECOMP_MMPBSA.dat",
    "Cinnamaldehyde":     "Cinnamaldehyde_DECOMP_MMPBSA.dat",
    "Eugenol":            "Eugenol_DECOMP_MMPBSA.dat",
    "Zingiberene":        "Zingiberene_DECOMP_MMPBSA.dat",
}
DATA_START  = 8
DATA_END    = 203
TOP_PER_LIG = 15    # residues per ligand for union set
OUTPUT_DIR  = "."   # change to your output directory

LIGANDS = list(LIGAND_FILES.keys())
COLORS  = ["#1565C0", "#C62828", "#6A1B9A", "#2E7D32", "#E65100"]
SHORT   = {
    "Ciprofloxacin":      "Cipro",
    "Diallyl Trisulfide": "DTS",
    "Cinnamaldehyde":     "Cinn",
    "Eugenol":            "Eug",
    "Zingiberene":        "Zing",
}

COLS = [
    "Residue", "Location",
    "Int_Avg",   "Int_Std",   "Int_SEM",
    "vdW_Avg",   "vdW_Std",   "vdW_SEM",
    "Elec_Avg",  "Elec_Std",  "Elec_SEM",
    "Polar_Avg", "Polar_Std", "Polar_SEM",
    "NP_Avg",    "NP_Std",    "NP_SEM",
    "TOTAL_Avg", "TOTAL_Std", "TOTAL_SEM",
]

# Light theme constants
BG      = "#FFFFFF"
FG      = "#1a1a2e"
SUB_FG  = "#555566"
GRID_C  = "#E8E8F0"
CARD_BG = "#F7F8FC"

DIV_COLORS = [
    "#2166ac", "#4393c3", "#74add1", "#abd9e9", "#ffffff",
    "#fee090", "#fdae61", "#d73027", "#a50026",
]
DIV_SCALE_PLOTLY = [
    [0.000, "#2166ac"], [0.125, "#4393c3"], [0.250, "#74add1"],
    [0.375, "#abd9e9"], [0.500, "#ffffff"],  [0.625, "#fee090"],
    [0.750, "#fdae61"], [0.875, "#d73027"],  [1.000, "#a50026"],
]


# ─── DATA LOADING ────────────────────────────────────────────────────────────
def parse_decomp(filepath: str) -> pd.DataFrame:
    with open(filepath) as fh:
        lines = fh.readlines()
    rows = []
    for line in lines[DATA_START:DATA_END]:
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


all_data = {name: parse_decomp(fname) for name, fname in LIGAND_FILES.items()}


# ─── DERIVE UNION TOP RESIDUES ───────────────────────────────────────────────
top_per = {}
for name, df in all_data.items():
    prot = df[~df["Residue"].str.startswith("LIG")].copy()
    prot["abs_total"] = prot["TOTAL_Avg"].abs()
    top_per[name] = set(prot.nlargest(TOP_PER_LIG, "abs_total")["Residue"].tolist())

union_top = set()
for s in top_per.values():
    union_top |= s

all_res_ord = all_data["Ciprofloxacin"]["Residue"].tolist()
top_res  = [r for r in all_res_ord if r in union_top and not r.startswith("LIG")]
short_res = [r.strip() for r in top_res]


# ─── MATRIX BUILDER ──────────────────────────────────────────────────────────
def get_matrix(col: str) -> np.ndarray:
    """Return (n_ligands × n_residues) array for a given column."""
    mat = []
    for name in LIGANDS:
        df = all_data[name]
        row = []
        for r in top_res:
            v = df[df["Residue"] == r][col].values
            row.append(float(v[0]) if len(v) > 0 else 0.0)
        mat.append(row)
    return np.array(mat)


total_mat = get_matrix("TOTAL_Avg")
total_std = get_matrix("TOTAL_Std")
vdw_mat   = get_matrix("vdW_Avg")
net_elec  = get_matrix("Elec_Avg") + get_matrix("Polar_Avg")
np_mat    = get_matrix("NP_Avg")
max_abs   = np.max(np.abs(total_mat), axis=0)


# ─── LIGAND SUMMARY DICT ─────────────────────────────────────────────────────
lig_summary = {}
for name, df in all_data.items():
    lig  = df[df["Residue"].str.startswith("LIG")].iloc[0]
    prot = df[~df["Residue"].str.startswith("LIG")]
    lig_summary[name] = {
        "lig_total": lig["TOTAL_Avg"], "lig_std":  lig["TOTAL_Std"],
        "lig_vdw":   lig["vdW_Avg"],  "lig_elec": lig["Elec_Avg"],
        "lig_polar": lig["Polar_Avg"],"lig_np":   lig["NP_Avg"],
        "prot_sum":  prot["TOTAL_Avg"].sum(),
        "grand_sum": df["TOTAL_Avg"].sum(),
    }


# ══════════════════════════════════════════════════════════════════════════════
#  MATPLOTLIB FIGURES
# ══════════════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec

matplotlib.rcParams['font.family'] = ['Helvetica', 'Helvetica Neue', 'Arial', 'sans-serif']

DIV_CMAP = LinearSegmentedColormap.from_list("div", DIV_COLORS)


# ── Figure 1: 4-panel heatmap ─────────────────────────────────────────────────
def plot_heatmaps(outpath: str):
    fig, axes = plt.subplots(4, 1, figsize=(22, 14), facecolor=BG)
    fig.subplots_adjust(hspace=0.55, left=0.08, right=0.94, top=0.92, bottom=0.12)

    panels = [
        (total_mat, "Total ΔG contribution"),
        (vdw_mat,   "van der Waals"),
        (net_elec,  "Net electrostatic (EEL + Polar)"),
        (np_mat,    "Non-polar solvation"),
    ]
    for ax, (mat, title) in zip(axes, panels):
        ax.set_facecolor(BG)
        lim  = max(abs(mat.min()), abs(mat.max()), 0.01)
        norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
        im   = ax.imshow(mat, aspect="auto", interpolation="nearest",
                         cmap=DIV_CMAP, norm=norm)
        ax.set_yticks(range(5))
        ax.set_yticklabels(LIGANDS, color=FG, fontsize=10, fontweight="bold")
        ax.set_xticks(range(len(short_res)))
        ax.set_xticklabels(short_res, rotation=45, ha="right",
                           fontsize=7.5, color=SUB_FG)
        ax.tick_params(left=False, bottom=False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.01)
        cb.ax.tick_params(labelsize=7.5, colors=SUB_FG)
        cb.outline.set_visible(False)
        cb.set_label("kcal/mol", color=SUB_FG, fontsize=8)
        ax.set_title(title, color=FG, fontsize=11, fontweight="bold",
                     loc="left", pad=5)
        # Annotate significant cells
        for i in range(5):
            for j in range(len(short_res)):
                v = mat[i, j]
                if abs(v) > 0.05:
                    tc = "white" if abs(v) > 0.5 else FG
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            fontsize=5.5, color=tc, fontweight="bold")

    fig.suptitle(
        "MM-GBSA Per-Residue Energy Decomposition — 5 Ligands vs DNA Gyrase",
        color=FG, fontsize=15, fontweight="bold", y=0.975,
    )
    plt.savefig(outpath, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


# ── Figure 2: Summary (grand total + components + radar) ─────────────────────
def plot_summary(outpath: str):
    fig = plt.figure(figsize=(22, 7), facecolor=BG)
    gs  = GridSpec(1, 3, figure=fig, wspace=0.38,
                   left=0.06, right=0.97, top=0.86, bottom=0.22)
    ax_grand = fig.add_subplot(gs[0, 0])
    ax_comp  = fig.add_subplot(gs[0, 1])
    ax_radar = fig.add_subplot(gs[0, 2], projection="polar")

    for ax in [ax_grand, ax_comp]:
        ax.set_facecolor(CARD_BG)
    ax_radar.set_facecolor(CARD_BG)

    # Grand total horizontal bars
    names_s  = sorted(LIGANDS, key=lambda x: lig_summary[x]["grand_sum"])
    grand_v  = [lig_summary[n]["grand_sum"] for n in names_s]
    cols_s   = [COLORS[LIGANDS.index(n)] for n in names_s]
    bars = ax_grand.barh(names_s, grand_v, color=cols_s, alpha=0.88,
                         height=0.6, zorder=3)
    ax_grand.axvline(0, color=SUB_FG, linewidth=0.8, linestyle="--", zorder=2)
    for bar, v in zip(bars, grand_v):
        ax_grand.text(
            v + (0.3 if v > 0 else -0.3),
            bar.get_y() + bar.get_height() / 2,
            f"{v:.2f}", va="center",
            ha="left" if v > 0 else "right",
            color=FG, fontsize=9.5, fontweight="bold",
        )
    ax_grand.set_xlabel("ΔG (kcal/mol)", color=SUB_FG, fontsize=10)
    ax_grand.set_title("Grand total ΔG", color=FG, fontsize=11,
                       fontweight="bold", pad=8)
    ax_grand.tick_params(colors=FG, labelsize=9)
    for sp in ax_grand.spines.values():
        sp.set_edgecolor(GRID_C)
    ax_grand.grid(axis="x", color=GRID_C, linewidth=0.5, zorder=0)

    # Ligand self-energy components
    x = np.arange(5)
    w = 0.18
    comp_map = [
        ("vdW",   "lig_vdw",   "#00897B"),
        ("Elec",  "lig_elec",  "#E53935"),
        ("Polar", "lig_polar", "#1E88E5"),
        ("NP",    "lig_np",    "#FB8C00"),
    ]
    for k, (label, key, color) in enumerate(comp_map):
        vals = [lig_summary[n][key] for n in LIGANDS]
        ax_comp.bar(x + k * w, vals, w, label=label,
                    color=color, alpha=0.85, zorder=3)
    ax_comp.axhline(0, color=SUB_FG, linewidth=0.8, linestyle="--", zorder=2)
    ax_comp.set_xticks(x + 1.5 * w)
    ax_comp.set_xticklabels([SHORT[n] for n in LIGANDS],
                             color=FG, fontsize=9.5, rotation=20)
    ax_comp.set_ylabel("ΔG (kcal/mol)", color=SUB_FG, fontsize=10)
    ax_comp.set_title("Ligand self-energy components", color=FG, fontsize=11,
                      fontweight="bold", pad=8)
    ax_comp.legend(fontsize=8, framealpha=0.9, facecolor=BG,
                   edgecolor=GRID_C, labelcolor=FG, loc="upper right")
    ax_comp.tick_params(colors=FG, labelsize=9)
    for sp in ax_comp.spines.values():
        sp.set_edgecolor(GRID_C)
    ax_comp.grid(axis="y", color=GRID_C, linewidth=0.5, zorder=0)

    # Radar chart
    top10_idx    = np.argsort(max_abs)[-10:][::-1]
    radar_res    = [short_res[i] for i in top10_idx]
    angles       = np.linspace(0, 2 * np.pi, len(radar_res), endpoint=False).tolist()
    angles_c     = angles + [angles[0]]
    radar_res_c  = radar_res + [radar_res[0]]
    for i, (name, color) in enumerate(zip(LIGANDS, COLORS)):
        vals   = [abs(total_mat[i, idx]) for idx in top10_idx]
        vals_c = vals + [vals[0]]
        ax_radar.plot(angles_c, vals_c, color=color, linewidth=2,
                      label=SHORT[name])
        ax_radar.fill(angles_c, vals_c, color=color, alpha=0.1)
    ax_radar.set_xticks(angles)
    ax_radar.set_xticklabels(radar_res, color=FG, fontsize=7.5)
    ax_radar.set_yticklabels([])
    ax_radar.grid(color=GRID_C, linewidth=0.5)
    ax_radar.spines["polar"].set_color(GRID_C)
    ax_radar.set_title("Binding profile\n(top 10 residues)", color=FG,
                       fontsize=10, fontweight="bold", pad=18)
    ax_radar.legend(fontsize=7.5, framealpha=0.9, facecolor=BG,
                    edgecolor=GRID_C, labelcolor=FG,
                    loc="upper right", bbox_to_anchor=(1.35, 1.12))

    fig.suptitle("Binding Free Energy Summary — 5 Ligands vs DNA Gyrase",
                 color=FG, fontsize=14, fontweight="bold", y=0.98)
    plt.savefig(outpath, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


# ── Figure 3: Grouped bar (top 15) ───────────────────────────────────────────
def plot_grouped_bar(outpath: str):
    top15_idx = np.argsort(max_abs)[-15:][::-1]
    top15_res = [short_res[i] for i in top15_idx]
    fig, ax   = plt.subplots(figsize=(22, 6), facecolor=BG)
    ax.set_facecolor(CARD_BG)
    x = np.arange(len(top15_res))
    w = 0.14
    for k, (name, color) in enumerate(zip(LIGANDS, COLORS)):
        vals = [total_mat[k, i] for i in top15_idx]
        stds = [total_std[k, i] for i in top15_idx]
        ax.bar(x + k * w, vals, w, label=SHORT[name],
               color=color, alpha=0.85, zorder=3)
        ax.errorbar(x + k * w, vals, yerr=stds, fmt="none",
                    ecolor="#333344", elinewidth=0.8,
                    capsize=2, capthick=0.8, zorder=4, alpha=0.7)
    ax.axhline(0, color=SUB_FG, linewidth=0.8, linestyle="--", zorder=2)
    ax.set_xticks(x + 2 * w)
    ax.set_xticklabels(top15_res, rotation=40, ha="right",
                       color=FG, fontsize=9.5)
    ax.set_ylabel("ΔG contribution (kcal/mol)", color=SUB_FG, fontsize=11)
    ax.set_title(
        "Top 15 key residues — per-ligand ΔG contribution ± std dev",
        color=FG, fontsize=13, fontweight="bold", pad=10,
    )
    ax.legend(fontsize=9, framealpha=0.9, facecolor=BG,
              edgecolor=GRID_C, labelcolor=FG, loc="lower right", ncol=5)
    ax.tick_params(colors=FG, labelsize=9)
    for sp in ax.spines.values():
        sp.set_edgecolor(GRID_C)
    ax.grid(axis="y", color=GRID_C, linewidth=0.5, zorder=0)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180, bbox_inches="tight", facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════════════
#  PLOTLY INTERACTIVE HTML FIGURES
# ══════════════════════════════════════════════════════════════════════════════
import plotly.graph_objects as go
from plotly.subplots import make_subplots

AXIS_STYLE = dict(
    ticks="inside",
    ticklen=10,
    tickwidth=2,
    showline=True,
    mirror="all",
    linewidth=2,
    linecolor="black",
    showgrid=False,
    tickfont=dict(family="Helvetica", size=14, color="black"),
    title_font=dict(family="Helvetica", size=16, color="black"),
)

AXIS_STYLE_BASE = dict(
    ticks="inside",
    ticklen=14,
    tickwidth=3,
    showline=True,
    mirror="all",
    linewidth=2,
    linecolor="black",
    showgrid=False,
    tickfont=dict(family="Helvetica", size=14, color="black"),
    title_font=dict(family="Helvetica", size=16, color="black"),
)

def build_figure(title: str) -> go.Figure:
    """
    Return a Figure with shared publication style — exactly as requested.
    Title, font, bgcolor, and default axis style are set here.
    Callers override xaxis/yaxis titles and any extra layout via
    a separate fig.update_layout() call after build_figure().
    """
    fig = go.Figure()
    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b>",
            x=0.5,
            font=dict(family="Helvetica", size=40, color="black"),
        ),
        font=dict(family="Helvetica", size=20, color="black"),
        plot_bgcolor=BG,
        paper_bgcolor=BG,
        xaxis=dict(**AXIS_STYLE_BASE),
        yaxis=dict(**AXIS_STYLE_BASE),
    )
    return fig


def _style_axes(fig, rows=1):
    """Apply shared axis style to all subplot rows (for make_subplots figures)."""
    for i in range(1, rows + 1):
        fig.update_xaxes(AXIS_STYLE_BASE, row=i, col=1)
        fig.update_yaxes(AXIS_STYLE_BASE, row=i, col=1)


def plot_interactive_heatmaps(outpath: str):
    fig = make_subplots(
        rows=4, cols=1,
        subplot_titles=["Total ΔG", "van der Waals",
                        "Net electrostatic", "Non-polar solvation"],
        vertical_spacing=0.08,
    )
    for ri, (mat, _) in enumerate(
        [(total_mat, ""), (vdw_mat, ""), (net_elec, ""), (np_mat, "")], 1
    ):
        lim   = max(abs(mat.min()), abs(mat.max()), 0.01)
        hover = [
            [f"<b>{short_res[j]}</b><br>{LIGANDS[i]}<br>ΔG = {mat[i,j]:.4f} kcal/mol"
             for j in range(len(short_res))]
            for i in range(5)
        ]
        fig.add_trace(
            go.Heatmap(
                z=mat, x=short_res, y=LIGANDS,
                colorscale=DIV_SCALE_PLOTLY, zmid=0, zmin=-lim, zmax=lim,
                colorbar=dict(
                    len=0.22, y=1 - (ri - 1) * 0.265, thickness=10,
                    tickfont=dict(size=8, color=SUB_FG),
                    title=dict(text="kcal/mol", font=dict(size=9, color=SUB_FG)),
                ),
                hovertext=hover,
                hovertemplate="%{hovertext}<extra></extra>",
            ),
            row=ri, col=1,
        )
    fig.update_layout(
        height=900, width=1400,
        title=dict(text="<b>MM-GBSA Per-Residue Decomposition — 5 Ligands</b>",
                   font=dict(family="Helvetica", size=22, color="black"), x=0.5),
        font=dict(family="Helvetica", size=14, color="black"),
        plot_bgcolor=BG, paper_bgcolor=BG,
        margin=dict(l=150, r=80, t=80, b=50),
    )
    _style_axes(fig, rows=4)
    for i in range(1, 5):
        fig.update_xaxes(tickangle=-45, tickfont=dict(family="Helvetica", size=10, color="black"),
                         row=i, col=1)
        fig.update_yaxes(tickfont=dict(family="Helvetica", size=12, color="black"), row=i, col=1)
    fig.write_html(outpath)
    print(f"  Saved: {outpath}")


def plot_interactive_grouped_bar(outpath: str):
    top15_idx = np.argsort(max_abs)[-15:][::-1]
    top15_res = [short_res[i] for i in top15_idx]
    fig = build_figure("Top 15 residues — per-ligand ΔG ± std dev")
    fig.update_layout(
        barmode="group",
        height=500, width=1300,
        margin=dict(l=70, r=40, t=130, b=130),
        xaxis=dict(**AXIS_STYLE_BASE, title="<b>Residue</b>", tickangle=-40),
        yaxis=dict(**AXIS_STYLE_BASE, title="<b>ΔG contribution (kcal/mol)</b>"),
        legend=dict(
            orientation="h", y=1.08, x=0.5, xanchor="center",
            font=dict(family="Helvetica", size=14, color="black"),
            bgcolor="rgba(255,255,255,0.9)", bordercolor="#CCCCDD", borderwidth=1,
        ),
    )
    for i, (name, color) in enumerate(zip(LIGANDS, COLORS)):
        vals = [total_mat[i, j] for j in top15_idx]
        stds = [total_std[i, j] for j in top15_idx]
        fig.add_trace(go.Bar(
            name=name, x=top15_res, y=vals,
            marker_color=color, opacity=0.85,
            error_y=dict(type="data", array=stds, visible=True,
                         thickness=1.2, width=3,
                         color="rgba(60,60,60,0.6)"),
            hovertemplate=f"<b>%{{x}}</b><br>{name}<br>ΔG = %{{y:.4f}} kcal/mol<extra></extra>",
        ))
    fig.add_hline(y=0, line=dict(color="black", width=1, dash="dot"))
    fig.write_html(outpath)
    print(f"  Saved: {outpath}")


def plot_interactive_radar(outpath: str):
    top10_idx   = np.argsort(max_abs)[-10:][::-1]
    radar_res   = [short_res[i] for i in top10_idx]
    radar_c     = radar_res + [radar_res[0]]
    fig = go.Figure()
    for i, (name, color) in enumerate(zip(LIGANDS, COLORS)):
        vals = [abs(total_mat[i, j]) for j in top10_idx]
        raw  = [total_mat[i, j] for j in top10_idx]
        hover = [f"{radar_res[k]}: {raw[k]:.4f} kcal/mol"
                 for k in range(len(radar_res))] + [""]
        fig.add_trace(go.Scatterpolar(
            r=vals + [vals[0]], theta=radar_c,
            fill="toself", name=name,
            line=dict(color=color, width=2),
            fillcolor=color, opacity=0.15,
            hovertemplate=f"<b>{name}</b><br>%{{text}}<extra></extra>",
            text=hover,
        ))
    fig.update_layout(
        polar=dict(
            bgcolor=CARD_BG,
            radialaxis=dict(
                visible=True,
                tickfont=dict(family="Helvetica", size=11, color="black"),
                gridcolor=GRID_C, linecolor="black",
            ),
            angularaxis=dict(
                tickfont=dict(family="Helvetica", size=12, color="black"),
                linecolor="black",
            ),
        ),
        title=dict(
            text="<b>Binding profile radar — top 10 key residues</b>",
            font=dict(family="Helvetica", size=40, color="black"), x=0.5,
        ),
        height=620, width=720,
        paper_bgcolor=BG,
        font=dict(family="Helvetica", size=20, color="black"),
        legend=dict(
            font=dict(family="Helvetica", size=13, color="black"),
            bgcolor="rgba(255,255,255,0.9)", bordercolor="#CCCCDD", borderwidth=1,
        ),
        margin=dict(l=80, r=180, t=80, b=80),
    )
    fig.write_html(outpath)
    print(f"  Saved: {outpath}")


def plot_interactive_summary_bar(outpath: str):
    names_s = sorted(LIGANDS, key=lambda x: lig_summary[x]["grand_sum"])
    grand_v = [lig_summary[n]["grand_sum"] for n in names_s]
    bar_c   = [COLORS[LIGANDS.index(n)] for n in names_s]
    fig = build_figure("Grand total ΔG per ligand")
    fig.update_layout(
        height=400, width=820,
        margin=dict(l=170, r=130, t=120, b=70),
        xaxis=dict(**AXIS_STYLE_BASE, title="<b>ΔG (kcal/mol)</b>"),
        yaxis=dict(**AXIS_STYLE_BASE),
    )
    fig.add_trace(go.Bar(
        x=grand_v, y=names_s, orientation="h",
        marker_color=bar_c, opacity=0.85, name="Grand total",
        text=[f"{v:.3f}" for v in grand_v], textposition="outside",
        textfont=dict(family="Helvetica", size=13, color="black"),
        hovertemplate="<b>%{y}</b><br>ΔG = %{x:.3f} kcal/mol<extra></extra>",
    ))
    fig.add_vline(x=0, line=dict(color="black", width=1, dash="dash"))
    fig.write_html(outpath)
    print(f"  Saved: {outpath}")


# ─── MAIN ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    out = Path(OUTPUT_DIR)
    out.mkdir(parents=True, exist_ok=True)

    print("\n[ Matplotlib static figures ]")
    plot_heatmaps(str(out / "fig1_heatmaps.png"))
    plot_summary(str(out / "fig2_summary.png"))
    plot_grouped_bar(str(out / "fig3_grouped_bar.png"))

    print("\n[ Plotly interactive HTML figures ]")
    plot_interactive_heatmaps(str(out / "interactive_heatmaps.html"))
    plot_interactive_grouped_bar(str(out / "interactive_grouped_bar.html"))
    plot_interactive_radar(str(out / "interactive_radar.html"))
    plot_interactive_summary_bar(str(out / "interactive_summary_bar.html"))

    print("\nDone — all figures saved.")
