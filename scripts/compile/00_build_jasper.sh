#!/bin/bash
# =============================================================================
# Build Jasper 2.0.33 from source
# Required for WPS GRIB2 support (libjasper-dev not in Ubuntu 22.04)
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
JASPER_VERSION=2.0.33
JASPER_DIR=$WRF_GRK/libs/jasper
JASPER_SRC=/tmp/jasper-${JASPER_VERSION}

echo "Building Jasper ${JASPER_VERSION}..."

# Build dependencies (cmake libjpeg-dev libpng-dev zlib1g-dev) assumed installed

# Download
if [ ! -f /tmp/jasper-${JASPER_VERSION}.tar.gz ]; then
    wget -O /tmp/jasper-${JASPER_VERSION}.tar.gz \
        "https://github.com/jasper-software/jasper/releases/download/version-${JASPER_VERSION}/jasper-${JASPER_VERSION}.tar.gz"
fi

# Extract
rm -rf $JASPER_SRC
tar -xzf /tmp/jasper-${JASPER_VERSION}.tar.gz -C /tmp/
# jasper archives as jasper-version-NUMBER, find it
EXTRACTED=$(find /tmp -maxdepth 1 -name "jasper-${JASPER_VERSION}*" -type d | head -1)
if [ "$EXTRACTED" != "$JASPER_SRC" ]; then
    mv "$EXTRACTED" "$JASPER_SRC"
fi

# Build
mkdir -p $JASPER_SRC/build
cd $JASPER_SRC/build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$JASPER_DIR \
    -DJAS_ENABLE_OPENGL=false \
    -DJAS_ENABLE_AUTOMATIC_DEPENDENCIES=false \
    -DCMAKE_BUILD_TYPE=Release

make -j $(nproc)
make install

echo ""
echo "Jasper installed to: $JASPER_DIR"
echo "  Headers: $JASPER_DIR/include"
echo "  Libs:    $JASPER_DIR/lib"
ls $JASPER_DIR/lib/libjasper* 2>/dev/null && echo "[OK] libjasper found" || echo "[WARNING] libjasper not found"
