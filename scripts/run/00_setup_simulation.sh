#!/bin/bash
# =============================================================================
# Setup simulation run directory
# Links WRF executables, data tables, met_em files, namelist, emissions
# Must run AFTER WRF compiled and preprocessing complete
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
WRF_SRC=$WRF_GRK/WRF
SIM_DIR=$WRF_GRK/simulations/IDN_BB_2015
INPUT_DIR=$SIM_DIR/input
MET_EM_DIR=$SIM_DIR/met_em
RUN_DIR=$SIM_DIR/run

mkdir -p $RUN_DIR $SIM_DIR/output $SIM_DIR/restart $SIM_DIR/logs

echo "Setting up run directory: $RUN_DIR"

# ---- 1. Link WRF executables ----
echo "[1] Linking executables..."
for exe in wrf.exe real.exe ndown.exe tc.exe; do
    if [ -f "$WRF_SRC/main/$exe" ]; then
        ln -sf "$WRF_SRC/main/$exe" "$RUN_DIR/$exe"
        echo "  [OK] $exe"
    fi
done

# ---- 2. Link WRF data tables ----
echo "[2] Linking data tables..."
for f in "$WRF_SRC/run/"*.TBL "$WRF_SRC/run/"*.tbl \
         "$WRF_SRC/run/"*.formatted "$WRF_SRC/run/"ETAMPNEW_DATA* \
         "$WRF_SRC/run/RRTM"* "$WRF_SRC/run/CAM"* \
         "$WRF_SRC/run/ozone"* "$WRF_SRC/run/aerosol"*; do
    [ -f "$f" ] && ln -sf "$f" "$RUN_DIR/" 2>/dev/null || true
done
# Also link chemistry tables
for f in "$WRF_SRC/run/"*.dat "$WRF_SRC/run/"*.inp \
         "$WRF_SRC/run/tr"* "$WRF_SRC/run/SOILPARM"*; do
    [ -f "$f" ] && ln -sf "$f" "$RUN_DIR/" 2>/dev/null || true
done
echo "  [OK] Data tables linked"

# ---- 3. Copy namelist.input ----
echo "[3] Copying namelist.input..."
cp "$INPUT_DIR/namelist.input" "$RUN_DIR/namelist.input"
echo "  [OK] namelist.input"

# ---- 4. Link met_em files ----
echo "[4] Linking met_em files..."
N_MET=$(ls "$MET_EM_DIR"/met_em.d01.* 2>/dev/null | wc -l)
if [ "$N_MET" -eq 0 ]; then
    echo "  WARNING: No met_em files found in $MET_EM_DIR"
    echo "  Run preprocessing/domain/run_wps.sh first"
else
    for f in "$MET_EM_DIR"/met_em.d01.*; do
        ln -sf "$f" "$RUN_DIR/"
    done
    echo "  [OK] $N_MET met_em files linked"
fi

# ---- 5. Link emission files ----
echo "[5] Linking emission files..."
# wrfchemi (anthropogenic)
N_CHEMI=$(ls "$INPUT_DIR"/wrfchemi_d01_* 2>/dev/null | wc -l)
if [ "$N_CHEMI" -gt 0 ]; then
    for f in "$INPUT_DIR"/wrfchemi_d01_*; do ln -sf "$f" "$RUN_DIR/"; done
    echo "  [OK] $N_CHEMI wrfchemi files"
else
    echo "  WARNING: No wrfchemi files (run preprocessing/anthropogenic/01_process_edgar.py)"
fi

# wrffirechemi (fire)
N_FIRE=$(ls "$INPUT_DIR"/wrffirechemi_d01_* 2>/dev/null | wc -l)
if [ "$N_FIRE" -gt 0 ]; then
    for f in "$INPUT_DIR"/wrffirechemi_d01_*; do ln -sf "$f" "$RUN_DIR/"; done
    echo "  [OK] $N_FIRE wrffirechemi files"
else
    echo "  WARNING: No wrffirechemi files (run preprocessing/fire/01_process_finn_hotspot.py)"
fi

# wrfoce (ocean)
N_OCE=$(ls "$INPUT_DIR"/wrfoce_d01_* 2>/dev/null | wc -l)
if [ "$N_OCE" -gt 0 ]; then
    for f in "$INPUT_DIR"/wrfoce_d01_*; do ln -sf "$f" "$RUN_DIR/"; done
    echo "  [OK] $N_OCE wrfoce files"
fi

# vprm_input
N_VPRM=$(ls "$INPUT_DIR"/vprm_input_d01_* 2>/dev/null | wc -l)
if [ "$N_VPRM" -gt 0 ]; then
    for f in "$INPUT_DIR"/vprm_input_d01_*; do ln -sf "$f" "$RUN_DIR/"; done
    echo "  [OK] $N_VPRM vprm_input files"
fi

echo ""
echo "Run directory ready: $RUN_DIR"
ls "$RUN_DIR" | head -20
