#!/bin/bash
# =============================================================================
# Compile WRF-Chem 4.1.5 for WRF-GRK (WRF-GHG) simulation
# Target: Indonesia biomass burning episode Aug-Nov 2015
# System: Ubuntu 22.04, Intel Xeon E5-2650 v2, OpenMPI + gfortran
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
WRF_SRC=$WRF_GRK/WRF
WRF_VERSION=4.1.5

echo "=========================================="
echo " WRF-Chem $WRF_VERSION Compilation Script"
echo "=========================================="

# ---- 1. Download WRF source --------------------------------------------------
if [ ! -s "$WRF_GRK/WRFV${WRF_VERSION}.tar.gz" ]; then
    echo "[1/5] Downloading WRF-Chem ${WRF_VERSION} source..."
    wget -O "$WRF_GRK/WRFV${WRF_VERSION}.tar.gz" \
        "https://github.com/wrf-model/WRF/archive/refs/tags/v${WRF_VERSION}.tar.gz"
else
    echo "[1/5] WRF source tarball already present."
fi

# ---- 2. Extract --------------------------------------------------------------
echo "[2/5] Extracting WRF source..."
# Extract tarball to a temporary location, then move to $WRF_SRC
TMPDIR_EXTRACT=$(mktemp -d /tmp/wrf_extract_XXXXXX)
tar -xzf $WRF_GRK/WRFV${WRF_VERSION}.tar.gz -C $TMPDIR_EXTRACT/
# The tarball extracts as WRF-4.1.5/ or WRFV4.1.5/
EXTRACTED=$(find $TMPDIR_EXTRACT -maxdepth 1 -mindepth 1 -type d | head -1)
if [ -z "$EXTRACTED" ]; then
    echo "ERROR: Could not find extracted WRF directory in $TMPDIR_EXTRACT"
    exit 1
fi
echo "  Extracted: $EXTRACTED → $WRF_SRC"
mkdir -p $WRF_SRC
# Move extracted contents into WRF_SRC
cp -a "$EXTRACTED/." "$WRF_SRC/"
rm -rf "$TMPDIR_EXTRACT"

# ---- 3. Set environment variables --------------------------------------------
echo "[3/5] Setting environment..."
export NETCDF=/usr
export NETCDF_classic=1
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

# WRF-Chem specific
export WRF_CHEM=1
export WRF_KPP=0
export WRF_EM_CORE=1
export WRF_NMM_CORE=0

# Compiler flags for stability
export FCFLAGS="-m64 -fallow-argument-mismatch -fallow-invalid-boz"
export FFLAGS="-m64 -fallow-argument-mismatch -fallow-invalid-boz"

# OpenMPI path
export PATH=/usr/lib/x86_64-linux-gnu/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/openmpi/lib:$LD_LIBRARY_PATH

# NetCDF Fortran (Ubuntu 22.04 netcdff)
export NETCDFF=/usr

# ---- 4. Configure WRF --------------------------------------------------------
echo "[4/5] Configuring WRF-Chem..."
cd $WRF_SRC

# Clean any previous build
if [ -f configure.wrf ]; then
    echo "  Cleaning previous build..."
    ./clean -a > /dev/null 2>&1 || true
fi

# Configure: option 34 = GNU/gfortran + OpenMPI (dmpar) on Linux
# Option 35 = dm+sm (shared+distributed memory); use 34 for pure MPI
echo "  Running ./configure (selecting option 34: GNU/OpenMPI dmpar)..."
printf "34\n1\n" | ./configure > /tmp/wrf_configure.log 2>&1

# Verify configure succeeded
if [ ! -f configure.wrf ]; then
    echo "ERROR: configure.wrf not created. See /tmp/wrf_configure.log"
    exit 1
fi

# Patch configure.wrf for Ubuntu 22.04 NetCDF paths
# NetCDF-Fortran library is libnetcdff (separate in Ubuntu)
sed -i 's/-lnetcdf$/-lnetcdf -lnetcdff/' configure.wrf
echo "  configure.wrf patched for Ubuntu 22.04 (added -lnetcdff)."

# ---- 5. Compile WRF-Chem ----------------------------------------------------
echo "[5/5] Compiling WRF-Chem (em_real with chem)..."
echo "  This will take 30-90 minutes on 32 cores..."
echo "  Log: $WRF_SRC/compile.log"

./compile -j 16 em_real > $WRF_SRC/compile.log 2>&1

# Verify executables
echo ""
echo "Checking compiled executables..."
MISSING=0
for exe in wrf.exe real.exe ndown.exe; do
    if [ -f "$WRF_SRC/main/${exe}" ]; then
        echo "  [OK] main/${exe}"
    else
        echo "  [MISSING] main/${exe}"
        MISSING=1
    fi
done

if [ $MISSING -eq 0 ]; then
    echo ""
    echo "=============================="
    echo " WRF-Chem compilation SUCCESS"
    echo "=============================="
    echo "Executables in: $WRF_SRC/main/"
else
    echo ""
    echo "=============================="
    echo " WRF-Chem compilation FAILED"
    echo " Check: $WRF_SRC/compile.log"
    echo "=============================="
    tail -30 $WRF_SRC/compile.log
    exit 1
fi
