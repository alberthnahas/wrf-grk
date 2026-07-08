#!/bin/bash
# =============================================================================
# Run real.exe to generate wrfinput_d01 and wrfbdy_d01
# Must run AFTER 00_setup_simulation.sh
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
SIM_DIR=$WRF_GRK/simulations/IDN_BB_2015
RUN_DIR=$SIM_DIR/run
LOG_DIR=$SIM_DIR/logs

mkdir -p $LOG_DIR

echo "=========================================="
echo " Running real.exe"
echo "=========================================="
echo "Run dir: $RUN_DIR"

cd $RUN_DIR

# real.exe needs namelist.input in the current directory
if [ ! -f namelist.input ]; then
    echo "ERROR: namelist.input not found. Run 00_setup_simulation.sh first."
    exit 1
fi
if [ ! -f real.exe ]; then
    echo "ERROR: real.exe not found. Run 00_setup_simulation.sh first."
    exit 1
fi
if [ ! -f met_em.d01.2015-07-25_00:00:00.nc ]; then
    echo "ERROR: met_em files not found. Run WPS first."
    exit 1
fi

echo "Running real.exe with 16 MPI tasks..."
mpirun --use-hwthread-cpus -np 16 ./real.exe > $LOG_DIR/real.log 2>&1

# Check output
echo ""
echo "Checking output..."
for f in wrfinput_d01 wrfbdy_d01; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  [OK] $f ($SIZE)"
        # Copy to input dir for preprocessing scripts
        cp "$f" "$SIM_DIR/$f"
    else
        echo "  [MISSING] $f"
        echo "  Check: $LOG_DIR/real.log"
        tail -30 $LOG_DIR/real.log
        exit 1
    fi
done

echo ""
grep -i "SUCCESS\|success\|error\|Error" $LOG_DIR/real.log | tail -5
echo ""
echo "real.exe complete. wrfinput_d01 and wrfbdy_d01 ready."
echo "Next: run preprocessing scripts, then 02_run_wrf.sh"
