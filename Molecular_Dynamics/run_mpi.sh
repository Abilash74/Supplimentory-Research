#!/bin/bash
# =============================================================================
#  MMPBSA MPI Calculation Launcher
#  Run this script from inside your MMPBSA/ subdirectory.
#  Input topologies and trajectory live in the parent directory.
# =============================================================================

set -e

echo "========================================"
echo "MMPBSA MPI Calculation Launcher"
echo "Directory : $(pwd)"
echo "Date      : $(date)"
echo "========================================"
echo ""

# =============================================================================
#  USER SETTINGS — EDIT HERE
# =============================================================================

NCORES=10       # Number of MPI ranks (adjust to your system)

# AMBER + OpenMPI paths — update if your installation differs
export AMBERHOME=/usr/local
export PATH="/home/user/openmpi/bin:$AMBERHOME/bin:$PATH"
export LD_LIBRARY_PATH="/home/user/openmpi/lib:$AMBERHOME/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$AMBERHOME/local/lib/python3.12/dist-packages:$PYTHONPATH"

# =============================================================================
#  ENVIRONMENT CHECK
# =============================================================================

echo "Environment:"
echo "  AMBERHOME : $AMBERHOME"
echo "  mpirun    : $(which mpirun)"
echo "  MMPBSA.py : $(which MMPBSA.py.MPI)"
echo ""

if ! command -v mpirun &> /dev/null; then
    echo "ERROR: mpirun not found in PATH."
    exit 1
fi

if ! command -v MMPBSA.py.MPI &> /dev/null; then
    echo "ERROR: MMPBSA.py.MPI not found in PATH."
    exit 1
fi

# =============================================================================
#  FILE CHECKS
# =============================================================================

# Topologies and trajectory are in the parent directory
DATA_DIR="$(cd .. && pwd)"

echo "Checking required files in: $DATA_DIR"
for FILE in complex_mb2.prmtop rec_mb2.prmtop ligand_mb2.prmtop prod_stripped.nc; do
    if [ ! -f "$DATA_DIR/$FILE" ]; then
        echo "ERROR: Required file not found: $DATA_DIR/$FILE"
        echo "  Run run_mmpbsa.sh first to prepare these files."
        exit 1
    fi
    echo "  OK  $FILE"
done

if [ ! -f mmpbsa.in ]; then
    echo "ERROR: mmpbsa.in not found in $(pwd)"
    echo "  Create your mmpbsa.in input file first."
    exit 1
fi
echo "  OK  mmpbsa.in"
echo ""

# =============================================================================
#  SHOW CONFIGURATION
# =============================================================================

echo "Configuration:"
echo "  MPI ranks   : $NCORES"
echo "  Input file  : mmpbsa.in"
echo "  Complex     : $DATA_DIR/complex_mb2.prmtop"
echo "  Receptor    : $DATA_DIR/rec_mb2.prmtop"
echo "  Ligand      : $DATA_DIR/ligand_mb2.prmtop"
echo "  Trajectory  : $DATA_DIR/prod_stripped.nc"
echo "  Output dir  : $(pwd)"
echo ""

echo "MMPBSA Input Parameters:"
echo "----------------------------------------"
cat mmpbsa.in
echo "----------------------------------------"
echo ""

# =============================================================================
#  LAUNCH CALCULATION
# =============================================================================

echo "Launching MMPBSA.py.MPI with $NCORES cores..."
echo ""

nohup mpirun -np $NCORES MMPBSA.py.MPI -O \
  -i mmpbsa.in \
  -cp "$DATA_DIR/complex_mb2.prmtop" \
  -rp "$DATA_DIR/rec_mb2.prmtop" \
  -lp "$DATA_DIR/ligand_mb2.prmtop" \
  -y  "$DATA_DIR/prod_stripped.nc" \
  -o  FINAL_RESULTS_MMPBSA.dat \
  -eo FINAL_RESULTS_MMPBSA.csv \
  -do DECOMP_MMPBSA.dat \
  > mmpbsa_mpi.log 2>&1 &

MPI_PID=$!
sleep 3

# =============================================================================
#  STATUS CHECK
# =============================================================================

if ps -p $MPI_PID > /dev/null 2>&1; then
    echo "MPI calculation started successfully!"
    echo "  Process ID : $MPI_PID"
    echo ""
    echo "Monitor progress:"
    echo "  tail -f $(pwd)/mmpbsa_mpi.log"
    echo "  ps aux | grep MMPBSA"
    echo ""
    echo "Output files (when complete):"
    echo "  FINAL_RESULTS_MMPBSA.dat   <- Summary binding free energy"
    echo "  FINAL_RESULTS_MMPBSA.csv   <- Per-frame energies"
    echo "  DECOMP_MMPBSA.dat          <- Per-residue decomposition"
    echo ""
    echo "========================================"
    echo "MPI calculation running in background"
    echo "========================================"
else
    echo "ERROR: Failed to start MPI calculation."
    echo "Last 30 lines of mmpbsa_mpi.log:"
    tail -30 mmpbsa_mpi.log 2>/dev/null || echo "(log file empty)"
    exit 1
fi
