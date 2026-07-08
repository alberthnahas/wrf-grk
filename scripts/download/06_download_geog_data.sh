#!/bin/bash
# =============================================================================
# Download WPS geographical static data
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
GEOG_DIR=$WRF_GRK/rawdata/geog
mkdir -p $GEOG_DIR

WPS_BASE="https://www2.mmm.ucar.edu/wrf/src/wps_files"

echo "Downloading WPS geographical data..."

# High-resolution datasets needed for the simulation
FILES=(
    "geog_high_res_mandatory.tar.gz"
    "modis_landuse_20class_30s_with_lakes.tar.bz2"
    "modis_landuse_20class_15s.tar.bz2"
    "topo_gmted2010_30s.tar.bz2"
    "soiltype_top_30s.tar.bz2"
    "soiltype_bot_30s.tar.bz2"
    "albedo_modis.tar.bz2"
    "maxsnowalb_modis.tar.bz2"
    "greenfrac_fpar_modis.tar.bz2"
    "lai_modis_30s.tar.bz2"
    "lai_modis_10m.tar.bz2"
    "orogwd_2deg.tar.bz2"
    "orogwd_1deg.tar.bz2"
    "soiltemp_1deg.tar.bz2"
    "vegparm.tar.bz2"
    "soilparm.tar.bz2"
)

for f in "${FILES[@]}"; do
    DEST=$GEOG_DIR/$f
    if [ -s "$DEST" ]; then
        echo "[SKIP] $f"
        continue
    fi
    echo "Downloading $f..."
    wget -q --show-progress -P $GEOG_DIR ${WPS_BASE}/$f || echo "  WARNING: $f not found (may not exist)"
done

# Extract all
echo ""
echo "Extracting geographical data..."
cd $GEOG_DIR
for f in *.tar.gz *.tar.bz2; do
    [ -s "$f" ] || continue
    echo "  Extracting $f..."
    tar -xf "$f" --no-same-owner 2>/dev/null || true
done

echo ""
echo "Geog data ready in: $GEOG_DIR"
ls -la $GEOG_DIR/ | head -20
