#!/bin/bash
# =============================================================================
#  MD Analysis Automation Script
#  Analyses: RMSD, RMSF, H-Bonds, Radius of Gyration, SASA
#  All outputs → analysis/<analysis_name>/
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                     USER SETTINGS — EDIT HERE                          ║
# ╠══════════════════════════════════════════════════════════════════════════╣

PARM="solvated.parm7"        # Topology file
TRAJ="prod.nc"               # Trajectory file
CPPTRAJ="cpptraj"            # cpptraj executable (full path if not in PATH)

# ── Residue Selections ───────────────────────────────────────────────────────
PROTEIN_RESIDUES="1-300"     # Residue range for SASA and Radius of Gyration
                             # e.g. "1-300" = residues 1 to 300 of your protein

ACTIVE_SITE_CA="50,80,100,120,150"
                             # Active site residue numbers for RMSD (CA atoms)
                             # Comma-separated, no spaces

LIGAND="LIG"                 # Ligand residue name as it appears in the topology

# ╚══════════════════════════════════════════════════════════════════════════╝

# ── Colour helpers ────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; CYAN='\033[0;36m'; RED='\033[0;31m'; NC='\033[0m'
info()  { echo -e "${CYAN}[INFO]${NC}  $*"; }
ok()    { echo -e "${GREEN}[ OK ]${NC}  $*"; }
fail()  { echo -e "${RED}[FAIL]${NC}  $*"; }

# ── Validate inputs ───────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  MD Analysis Pipeline"
echo "  Topology  : $PARM"
echo "  Trajectory: $TRAJ"
echo "  Protein residues (SASA/RadGyr): :${PROTEIN_RESIDUES}"
echo "  Active site CA residues (RMSD) : :${ACTIVE_SITE_CA}@CA"
echo "  Ligand (RMSD)                  : :${LIGAND}"
echo "============================================================"
echo ""

[[ ! -f "$PARM" ]] && { fail "Topology not found: $PARM"; exit 1; }
[[ ! -f "$TRAJ" ]] && { fail "Trajectory not found: $TRAJ"; exit 1; }
command -v "$CPPTRAJ" &>/dev/null || { fail "cpptraj not found: $CPPTRAJ"; exit 1; }

ANALYSIS_DIR="analysis"
mkdir -p "$ANALYSIS_DIR"
info "Output root -> ${ANALYSIS_DIR}/"
echo ""

# =============================================================================
#  HELPER FUNCTION
# =============================================================================
run_analysis() {
    local label="$1"
    local script_content="$2"
    local out_subdir="$3"

    local out_path="${ANALYSIS_DIR}/${out_subdir}"
    mkdir -p "$out_path"

    local tmp_script
    tmp_script=$(mktemp /tmp/cpptraj_XXXXXX.in)
    printf "%s\n" "$script_content" > "$tmp_script"

    info "Running ${label} ..."
    if "$CPPTRAJ" -i "$tmp_script" > "${out_path}/${out_subdir}.log" 2>&1; then
        ok "${label} complete  ->  ${out_path}/"
    else
        fail "${label} FAILED. Check log: ${out_path}/${out_subdir}.log"
    fi

    rm -f "$tmp_script"
}

# =============================================================================
#  1. HBONDS
# =============================================================================
run_analysis "H-Bonds" "
parm ${PARM}
trajin ${TRAJ}
hbond out ${ANALYSIS_DIR}/hbonds/hbonds.dat avgout ${ANALYSIS_DIR}/hbonds/hbonds_avg.dat
run
" "hbonds"

# =============================================================================
#  2. RADIUS OF GYRATION
# =============================================================================
run_analysis "Radius of Gyration" "
parm ${PARM}
trajin ${TRAJ}
radgyr out ${ANALYSIS_DIR}/radgyr/radgyr.dat :${PROTEIN_RESIDUES}
run
" "radgyr"

# =============================================================================
#  3. RMSD
# =============================================================================
run_analysis "RMSD" "
parm ${PARM}
trajin ${TRAJ}
autoimage
strip :WAT,HOH,Na+,Cl-
resinfo :${LIGAND}
rms tofirst_complex   !@H=                  first mass out ${ANALYSIS_DIR}/rmsd/rmsd_to_first_complex.dat
rms tofirst_ligand    :${LIGAND}&!@H=       first mass out ${ANALYSIS_DIR}/rmsd/rmsd_to_first_ligand.dat
rms tofirst_activeCA  :${ACTIVE_SITE_CA}@CA first mass out ${ANALYSIS_DIR}/rmsd/rmsd_to_first_AS_CA.dat
run
" "rmsd"

# =============================================================================
#  4. RMSF
# =============================================================================
run_analysis "RMSF" "
parm ${PARM}
trajin ${TRAJ}
atomicfluct out ${ANALYSIS_DIR}/rmsf/rmsf.dat byres
run
" "rmsf"

# =============================================================================
#  5. SASA
# =============================================================================
run_analysis "SASA" "
parm ${PARM}
trajin ${TRAJ}
autoimage
strip :WAT,HOH,Na+,Cl-
surf :${PROTEIN_RESIDUES} out ${ANALYSIS_DIR}/sasa/sasa.dat
run
" "sasa"

# =============================================================================
#  SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo "  All analyses finished. Results tree:"
echo "============================================================"
find "$ANALYSIS_DIR" -type f | sort | sed 's/^/  /'
echo ""
ok "Done! All outputs saved inside -> ${ANALYSIS_DIR}/"
