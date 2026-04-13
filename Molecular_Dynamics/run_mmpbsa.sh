#!/bin/bash

# ==========================================
# MM/GBSA PREPARATION PIPELINE
# Strips trajectory, prepares topologies
# Stops after topology preparation (STEP 6)
# ==========================================

set -e

echo "=========================================="
echo "   MM/GBSA PREPARATION (STOP AT STEP 6)"
echo "=========================================="

# ==================================================
# USER SETTINGS — EDIT HERE
# ==================================================

SOLVATED_PARM="solvated.parm7"       # Solvated topology file
TRAJ="prod.nc"                        # Production trajectory
LIGAND_MASK=":LIG"                    # Ligand residue mask (e.g. :LIG, :MOL)
SOLVENT_MASK=":WAT,HOH,Na+,Cl-,CLA,POT"  # Solvent/ion mask to strip

# ==================================================
# 0. CHECK AMBER ENVIRONMENT
# ==================================================

if [ -z "$AMBERHOME" ]; then
    echo "ERROR: AMBERHOME is not set."
    echo "  Run: source /usr/local/amber.sh"
    exit 1
fi

for cmd in cpptraj parmed ante-MMPBSA.py; do
    if ! command -v $cmd &> /dev/null; then
        echo "ERROR: $cmd not found in PATH."
        exit 1
    fi
done

echo "AMBER environment OK"
echo "AMBERHOME = $AMBERHOME"
echo ""

# ==================================================
# 1. CLEAN OLD FILES
# ==================================================

echo "Cleaning old MM/GBSA files..."
rm -f *.prmtop prod_stripped.nc *.dat *.csv *.log

# ==================================================
# 2. STRIP TRAJECTORY
# ==================================================

echo "Stripping trajectory..."

cpptraj <<EOF
parm $SOLVATED_PARM
trajin $TRAJ
strip $SOLVENT_MASK
trajout prod_stripped.nc
run
EOF

echo "Stripped trajectory saved: prod_stripped.nc"
echo ""

# ==================================================
# 3. STRIP TOPOLOGY
# ==================================================

echo "Creating stripped topology..."

cpptraj <<EOF
parm $SOLVATED_PARM
parmstrip $SOLVENT_MASK
parmwrite out complex.prmtop
run
EOF

echo "Stripped topology saved: complex.prmtop"
echo ""

# ==================================================
# 4. VERIFY MATCH
# ==================================================

echo "Verifying topology-trajectory compatibility..."

cpptraj -p complex.prmtop -y prod_stripped.nc <<EOF
quit
EOF

echo "Topology and trajectory match."
echo ""

# ==================================================
# 5. APPLY MBONDI2 RADII TO COMPLEX
# ==================================================

echo "Applying mbondi2 radii to complex..."

parmed complex.prmtop <<EOF
changeRadii mbondi2
outparm complex_mb2.prmtop
quit
EOF

echo "mbondi2 topology saved: complex_mb2.prmtop"
echo ""

# ==================================================
# 6. SPLIT TOPOLOGY + ENFORCE MBONDI2 FOR ALL
# ==================================================

echo "Splitting topology into receptor and ligand..."

ante-MMPBSA.py \
 -p complex_mb2.prmtop \
 -c com.prmtop \
 -r rec.prmtop \
 -l ligand.prmtop \
 -n $LIGAND_MASK

# ante-MMPBSA.py may skip com.prmtop for already-stripped input
if [ ! -f com.prmtop ]; then
    echo "  (com.prmtop not written — using complex_mb2.prmtop directly)"
    cp complex_mb2.prmtop com.prmtop
fi

echo "Enforcing mbondi2 radii on all topology files..."

for TOP in com rec ligand; do
    echo "  -> Converting ${TOP}.prmtop"
    parmed ${TOP}.prmtop <<EOF
changeRadii mbondi2
outparm ${TOP}_mb2.prmtop
quit
EOF
done

echo ""
echo "=========================================="
echo "  TOPOLOGY PREPARATION COMPLETE"
echo "  Files created:"
echo "    complex_mb2.prmtop  (complex with mbondi2)"
echo "    rec_mb2.prmtop      (receptor with mbondi2)"
echo "    ligand_mb2.prmtop   (ligand with mbondi2)"
echo "    prod_stripped.nc    (solvent-stripped trajectory)"
echo ""
echo "  Next step: run run_mpi.sh in the MMPBSA/ folder"
echo "  (Stopped at STEP 6 intentionally)"
echo "=========================================="
