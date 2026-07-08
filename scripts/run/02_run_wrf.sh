#!/bin/bash
# =============================================================================
# Run wrf.exe with restart cycling for IDN_BB_2015
# Period: 2015-07-25 to 2015-12-01 (5 segments of ~30 days)
# Uses 28 MPI tasks (leave 4 for OS/IO)
# =============================================================================

WRF_GRK=/home/igrk/WRF-GRK
SIM_DIR=$WRF_GRK/simulations/IDN_BB_2015
RUN_DIR=$SIM_DIR/run
LOG_DIR=$SIM_DIR/logs
OUT_DIR=$SIM_DIR/output
RST_DIR=$SIM_DIR/restart

mkdir -p $LOG_DIR $OUT_DIR $RST_DIR

NP=28

# Simulation segments (restart cycling)
# Format: START END SEGMENT_NAME
SEGMENTS=(
    "2015-07-25_00:00:00 2015-08-24_00:00:00 seg1"
    "2015-08-24_00:00:00 2015-09-23_00:00:00 seg2"
    "2015-09-23_00:00:00 2015-10-23_00:00:00 seg3"
    "2015-10-23_00:00:00 2015-11-22_00:00:00 seg4"
    "2015-11-22_00:00:00 2015-11-30_18:00:00 seg5"
)

echo "=========================================="
echo " WRF-GHG (chem_opt=17) IDN_BB_2015"
echo " Period: 2015-07-25 to 2015-12-01"
echo " MPI tasks: $NP"
echo "=========================================="

cd $RUN_DIR

if [ ! -f wrf.exe ]; then
    echo "ERROR: wrf.exe not found. Run 00_setup_simulation.sh first."
    exit 1
fi

for SEG in "${SEGMENTS[@]}"; do
    START=$(echo $SEG | awk '{print $1}')
    END=$(echo $SEG | awk '{print $2}')
    NAME=$(echo $SEG | awk '{print $3}')

    START_Y=$(echo $START | cut -c1-4)
    START_M=$(echo $START | cut -c6-7)
    START_D=$(echo $START | cut -c9-10)
    START_H=$(echo $START | cut -c12-13)
    END_Y=$(echo $END | cut -c1-4)
    END_M=$(echo $END | cut -c6-7)
    END_D=$(echo $END | cut -c9-10)
    END_H=$(echo $END | cut -c12-13)

    # Calculate run hours
    START_EPOCH=$(date -d "${START_Y}-${START_M}-${START_D} ${START_H}:00:00" +%s)
    END_EPOCH=$(date -d "${END_Y}-${END_M}-${END_D} ${END_H}:00:00" +%s)
    RUN_HOURS=$(( (END_EPOCH - START_EPOCH) / 3600 ))

    echo ""
    echo "=========================================="
    echo " Segment $NAME: $START → $END ($RUN_HOURS h)"
    echo "=========================================="

    # Update namelist.input
    sed -i "s/run_hours\s*=.*/run_hours = $RUN_HOURS,/" namelist.input
    sed -i "s/start_year\s*=.*/start_year = $START_Y,/" namelist.input
    sed -i "s/start_month\s*=.*/start_month = $START_M,/" namelist.input
    sed -i "s/start_day\s*=.*/start_day = $START_D,/" namelist.input
    sed -i "s/start_hour\s*=.*/start_hour = $START_H,/" namelist.input
    sed -i "s/end_year\s*=.*/end_year = $END_Y,/" namelist.input
    sed -i "s/end_month\s*=.*/end_month = $END_M,/" namelist.input
    sed -i "s/end_day\s*=.*/end_day = $END_D,/" namelist.input
    sed -i "s/end_hour\s*=.*/end_hour = $END_H,/" namelist.input

    # Set restart flag (false for first segment, true for subsequent)
    if [ "$NAME" = "seg1" ]; then
        sed -i "s/^ *restart  *= .*/restart                             = .false.,/" namelist.input
    else
        sed -i "s/^ *restart  *= .*/restart                             = .true.,/" namelist.input
        # Link restart file
        RST_FILE="wrfrst_d01_${START}"
        if [ ! -f "$RST_FILE" ]; then
            # Check restart dir
            if [ -f "$RST_DIR/$RST_FILE" ]; then
                ln -sf "$RST_DIR/$RST_FILE" "$RUN_DIR/$RST_FILE"
            else
                echo "ERROR: Restart file not found: $RST_FILE"
                exit 1
            fi
        fi
    fi

    # Run wrf.exe
    echo "Starting wrf.exe at $(date)..."
    mpirun --use-hwthread-cpus -np $NP ./wrf.exe > "$LOG_DIR/wrf_${NAME}.log" 2>&1
    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: wrf.exe failed for $NAME (exit $EXIT_CODE)"
        tail -30 "$LOG_DIR/wrf_${NAME}.log"
        exit 1
    fi

    # Check for SUCCESS
    if ! grep -q "wrf: SUCCESS COMPLETE WRF" "$LOG_DIR/wrf_${NAME}.log" 2>/dev/null; then
        if ! tail -5 "$LOG_DIR/wrf_${NAME}.log" | grep -qi "success\|complete"; then
            echo "WARNING: Could not confirm WRF success for $NAME"
            tail -10 "$LOG_DIR/wrf_${NAME}.log"
        fi
    fi
    echo "  Segment $NAME complete at $(date)"

    # Move output files
    echo "  Moving output files..."
    for f in wrfout_d01_*; do
        [ -f "$f" ] && mv "$f" "$OUT_DIR/"
    done

    # Move restart files
    for f in wrfrst_d01_*; do
        [ -f "$f" ] && mv "$f" "$RST_DIR/"
    done

    echo "  Output in: $OUT_DIR"
    ls "$OUT_DIR" | tail -5
done

echo ""
echo "=========================================="
echo " WRF-GHG simulation COMPLETE"
echo " Output in: $OUT_DIR"
echo "=========================================="
ls -lh "$OUT_DIR"
