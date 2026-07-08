#!/bin/bash
# =============================================================================
# Run all preprocessing scripts in parallel (after real.exe completes)
# Usage: bash scripts/run/03_run_preprocessing.sh [--wait]
#   --wait : block until all screens finish (polls every 30s)
# Requires: wrfinput_d01 and wrfbdy_d01 in simulations/IDN_BB_2015/run/
# =============================================================================

WRF_GRK=/home/igrk/WRF-GRK
VENV=$WRF_GRK/.venv/bin/activate
LOG_DIR=$WRF_GRK/simulations/IDN_BB_2015/logs
RUN_DIR=$WRF_GRK/simulations/IDN_BB_2015/run

mkdir -p $LOG_DIR

# Check prerequisites
if [ ! -f "$RUN_DIR/wrfinput_d01" ] && [ ! -f "$WRF_GRK/simulations/IDN_BB_2015/wrfinput_d01" ]; then
    echo "ERROR: wrfinput_d01 not found. Run 01_run_real.sh first."
    exit 1
fi

echo "================================================"
echo " Launching preprocessing scripts in parallel"
echo " Logs: $LOG_DIR/"
echo "================================================"

# Kill any stale screens with same names
for name in edgar fire ocean carbontracker; do
    screen -S $name -X quit 2>/dev/null || true
done

# 1. EDGAR anthropogenic emissions
screen -dmS edgar bash -c "
    source $VENV
    cd $WRF_GRK
    echo '[EDGAR] Started at \$(date)' > $LOG_DIR/edgar.log
    python preprocessing/anthropogenic/01_process_edgar.py >> $LOG_DIR/edgar.log 2>&1
    echo '[EDGAR] Exit code: \$?' >> $LOG_DIR/edgar.log
    echo '[EDGAR] Done at \$(date)' >> $LOG_DIR/edgar.log
"
echo "[LAUNCHED] edgar  → $LOG_DIR/edgar.log"

# 2. Fire emissions (FINN + hotspot)
screen -dmS fire bash -c "
    source $VENV
    cd $WRF_GRK
    echo '[FIRE] Started at \$(date)' > $LOG_DIR/fire.log
    python preprocessing/fire/01_process_finn_hotspot.py >> $LOG_DIR/fire.log 2>&1
    echo '[FIRE] Exit code: \$?' >> $LOG_DIR/fire.log
    echo '[FIRE] Done at \$(date)' >> $LOG_DIR/fire.log
"
echo "[LAUNCHED] fire   → $LOG_DIR/fire.log"

# 3. Ocean CO2 flux (SOMFFN)
screen -dmS ocean bash -c "
    source $VENV
    cd $WRF_GRK
    echo '[OCEAN] Started at \$(date)' > $LOG_DIR/ocean.log
    python preprocessing/ocean/01_process_somffn.py >> $LOG_DIR/ocean.log 2>&1
    echo '[OCEAN] Exit code: \$?' >> $LOG_DIR/ocean.log
    echo '[OCEAN] Done at \$(date)' >> $LOG_DIR/ocean.log
"
echo "[LAUNCHED] ocean  → $LOG_DIR/ocean.log"

# 4. CarbonTracker BCs/ICs
screen -dmS carbontracker bash -c "
    source $VENV
    cd $WRF_GRK
    echo '[CT] Started at \$(date)' > $LOG_DIR/carbontracker.log
    python preprocessing/boundary_conditions/01_process_carbontracker.py >> $LOG_DIR/carbontracker.log 2>&1
    echo '[CT] Exit code: \$?' >> $LOG_DIR/carbontracker.log
    echo '[CT] Done at \$(date)' >> $LOG_DIR/carbontracker.log
"
echo "[LAUNCHED] carbontracker → $LOG_DIR/carbontracker.log"

# Note: VPRM launched separately after MODIS completes
echo ""
echo "NOTE: VPRM script must be launched manually after MODIS download completes:"
echo "  screen -dmS vprm bash -c 'source $VENV && cd $WRF_GRK && python preprocessing/vprm/01_process_vprm.py >> $LOG_DIR/vprm.log 2>&1'"
echo ""

if [ "$1" == "--wait" ]; then
    echo "Waiting for all preprocessing to complete..."
    while true; do
        RUNNING=$(screen -ls | grep -cE "edgar|fire|ocean|carbontracker" || true)
        if [ "$RUNNING" -eq 0 ]; then
            echo "All preprocessing screens finished."
            break
        fi
        echo "  $RUNNING screens still running... ($(date +%H:%M:%S))"
        # Show last line of each log
        for name in edgar fire ocean carbontracker; do
            last=$(tail -1 $LOG_DIR/${name}.log 2>/dev/null)
            [ -n "$last" ] && echo "    [$name] $last"
        done
        sleep 30
    done
    echo ""
    echo "--- Final status ---"
    for name in edgar fire ocean carbontracker; do
        echo "[$name] $(tail -2 $LOG_DIR/${name}.log 2>/dev/null | tr '\n' ' ')"
    done
fi
