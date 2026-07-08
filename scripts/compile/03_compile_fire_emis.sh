#!/bin/bash
# =============================================================================
# Compile NCAR fire_emis tool for WRF-GHG fire emissions
# NOTE: fire_emis source must be downloaded manually from:
#   https://www2.acom.ucar.edu/wrf-chem/wrf-chem-tools-community
# Download fire_emis.tar.gz and place in /home/igrk/WRF-GRK/
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
FIRE_EMIS_SRC=$WRF_GRK/fire_emis

echo "=========================================="
echo " fire_emis Tool Compilation"
echo "=========================================="

# ---- 1. Check source ---------------------------------------------------------
if [ ! -s "$WRF_GRK/fire_emis.tar.gz" ]; then
    echo "ERROR: fire_emis.tar.gz not found at $WRF_GRK/fire_emis.tar.gz"
    echo "Please download from:"
    echo "  https://www2.acom.ucar.edu/wrf-chem/wrf-chem-tools-community"
    exit 1
fi

# ---- 2. Extract --------------------------------------------------------------
echo "[1/3] Extracting fire_emis..."
TMPDIR_EXTRACT=$(mktemp -d /tmp/fire_emis_XXXXXX)
tar -xzf "$WRF_GRK/fire_emis.tar.gz" -C "$TMPDIR_EXTRACT/"
EXTRACTED=$(find "$TMPDIR_EXTRACT" -maxdepth 1 -mindepth 1 -type d | head -1)
mkdir -p "$FIRE_EMIS_SRC"
cp -a "$EXTRACTED/." "$FIRE_EMIS_SRC/"
rm -rf "$TMPDIR_EXTRACT"

# ---- 3. Set environment and compile ------------------------------------------
echo "[2/3] Compiling fire_emis..."
export NETCDF=/usr
export FC=gfortran
export FFLAGS="-O2"

cd "$FIRE_EMIS_SRC"

if [ -f Makefile ]; then
    make clean > /dev/null 2>&1 || true
    make -j $(nproc) > "$FIRE_EMIS_SRC/compile.log" 2>&1
else
    echo "ERROR: No Makefile found in $FIRE_EMIS_SRC"
    exit 1
fi

# ---- 4. Verify ---------------------------------------------------------------
echo "[3/3] Verifying..."
if [ -f "$FIRE_EMIS_SRC/fire_emis" ]; then
    echo ""
    echo "=============================="
    echo " fire_emis compilation SUCCESS"
    echo " Executable: $FIRE_EMIS_SRC/fire_emis"
    echo "=============================="
else
    echo "ERROR: fire_emis executable not found. Check compile.log"
    tail -20 "$FIRE_EMIS_SRC/compile.log"
    exit 1
fi
