#!/bin/bash
# =============================================================================
# Compile WPS 4.1 with Jasper support (GRIB2)
# Must run AFTER 00_build_jasper.sh and 01_compile_wrf.sh
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
WPS_SRC=$WRF_GRK/WPS
WRF_SRC=$WRF_GRK/WRF
JASPER_DIR=$WRF_GRK/libs/jasper
WPS_VERSION=4.1

echo "=========================================="
echo " WPS ${WPS_VERSION} Compilation Script"
echo "=========================================="

# ---- 1. Check prerequisites --------------------------------------------------
if [ ! -f "$WRF_SRC/main/real.exe" ]; then
    echo "ERROR: WRF not compiled yet. Run 01_compile_wrf.sh first."
    exit 1
fi
if [ ! -d "$JASPER_DIR/lib" ]; then
    echo "ERROR: Jasper not built. Run 00_build_jasper.sh first."
    exit 1
fi

# ---- 2. Download WPS ---------------------------------------------------------
if [ ! -s "$WRF_GRK/WPS-${WPS_VERSION}.tar.gz" ]; then
    echo "[1/4] Downloading WPS ${WPS_VERSION}..."
    wget -O "$WRF_GRK/WPS-${WPS_VERSION}.tar.gz" \
        "https://github.com/wrf-model/WPS/archive/refs/tags/v${WPS_VERSION}.tar.gz"
else
    echo "[1/4] WPS tarball already present."
fi

# ---- 3. Extract --------------------------------------------------------------
echo "[2/4] Extracting WPS..."
TMPDIR_EXTRACT=$(mktemp -d /tmp/wps_extract_XXXXXX)
tar -xzf "$WRF_GRK/WPS-${WPS_VERSION}.tar.gz" -C "$TMPDIR_EXTRACT/"
EXTRACTED=$(find "$TMPDIR_EXTRACT" -maxdepth 1 -mindepth 1 -type d | head -1)
mkdir -p "$WPS_SRC"
cp -a "$EXTRACTED/." "$WPS_SRC/"
rm -rf "$TMPDIR_EXTRACT"

# ---- 4. Set environment ------------------------------------------------------
echo "[3/4] Setting environment..."
export WRF_DIR=$WRF_SRC
export NETCDF=/usr
export JASPERLIB=$JASPER_DIR/lib
export JASPERINC=$JASPER_DIR/include
export PATH=/usr/lib/x86_64-linux-gnu/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$JASPER_DIR/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:$LD_LIBRARY_PATH

# ---- 5. Configure and compile ------------------------------------------------
echo "[4/4] Configuring and compiling WPS..."
cd "$WPS_SRC"

if [ -f Makefile ]; then
    make clean > /dev/null 2>&1 || true
fi

# Option 3 = Linux x86_64, gfortran (serial) – WPS doesn't use MPI
echo "3\n" | ./configure > /tmp/wps_configure.log 2>&1

./compile > "$WPS_SRC/compile.log" 2>&1

# Verify
echo ""
MISSING=0
for exe in geogrid.exe ungrib.exe metgrid.exe; do
    if [ -f "$WPS_SRC/${exe}" ]; then
        echo "  [OK] $exe"
    else
        echo "  [MISSING] $exe"
        MISSING=1
    fi
done

if [ $MISSING -eq 0 ]; then
    echo ""
    echo "=============================="
    echo " WPS compilation SUCCESS"
    echo "=============================="
else
    echo ""
    echo "=============================="
    echo " WPS compilation FAILED"
    echo " Check: $WPS_SRC/compile.log"
    echo "=============================="
    tail -30 "$WPS_SRC/compile.log"
    exit 1
fi
