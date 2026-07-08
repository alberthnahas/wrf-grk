#!/bin/bash
# =============================================================================
# Run WPS (geogrid, ungrib, metgrid) for IDN_BB_2015
# Must run AFTER WPS compiled and ERA5 downloaded
# =============================================================================
set -e

WRF_GRK=/home/igrk/WRF-GRK
WPS_SRC=$WRF_GRK/WPS
ERA5_DIR=$WRF_GRK/rawdata/era5
MET_EM_DIR=$WRF_GRK/simulations/IDN_BB_2015/met_em
WPS_WORK=$WRF_GRK/preprocessing/domain

START_DATE="2015-07-25"
END_DATE="2015-12-01"

echo "=========================================="
echo " Running WPS for IDN_BB_2015"
echo " Period: $START_DATE to $END_DATE"
echo "=========================================="

mkdir -p $MET_EM_DIR
cd $WPS_WORK

# Link WPS executables if not already linked
for exe in geogrid.exe ungrib.exe metgrid.exe; do
    [ -f $exe ] || ln -sf $WPS_SRC/$exe .
done

# Link Vtable for ERA5 (ECMWF GRIB)
[ -f Vtable ] || ln -sf $WPS_SRC/ungrib/Variable_Tables/Vtable.ECMWF Vtable

# ---- 1. geogrid ----
echo "[1/3] Running geogrid..."
./geogrid.exe > $WRF_GRK/simulations/IDN_BB_2015/logs/geogrid.log 2>&1
if grep -q "Successful completion" $WRF_GRK/simulations/IDN_BB_2015/logs/geogrid.log; then
    echo "  [OK] geogrid complete"
else
    echo "  ERROR: geogrid failed. Check logs/geogrid.log"
    tail -20 $WRF_GRK/simulations/IDN_BB_2015/logs/geogrid.log
    exit 1
fi

# ---- 2. ungrib ----
echo "[2/3] Running ungrib..."
# Link ERA5 GRIB files
rm -f GRIBFILE.* FILE:*
$WPS_SRC/link_grib.csh $ERA5_DIR/era5_pl_2015*.grib $ERA5_DIR/era5_sl_2015*.grib

./ungrib.exe > $WRF_GRK/simulations/IDN_BB_2015/logs/ungrib.log 2>&1
if grep -q "Successful completion" $WRF_GRK/simulations/IDN_BB_2015/logs/ungrib.log; then
    echo "  [OK] ungrib complete"
else
    echo "  ERROR: ungrib failed. Check logs/ungrib.log"
    tail -20 $WRF_GRK/simulations/IDN_BB_2015/logs/ungrib.log
    exit 1
fi

# ---- 3. metgrid ----
echo "[3/3] Running metgrid..."
ln -sf $WPS_SRC/metgrid/METGRID.TBL.ARW METGRID.TBL
mpirun -np 4 ./metgrid.exe > $WRF_GRK/simulations/IDN_BB_2015/logs/metgrid.log 2>&1
if grep -q "Successful completion" $WRF_GRK/simulations/IDN_BB_2015/logs/metgrid.log; then
    echo "  [OK] metgrid complete"
else
    echo "  ERROR: metgrid failed. Check logs/metgrid.log"
    tail -20 $WRF_GRK/simulations/IDN_BB_2015/logs/metgrid.log
    exit 1
fi

echo ""
echo "WPS complete. met_em files:"
ls $MET_EM_DIR/met_em.d01.* | wc -l
echo "files in $MET_EM_DIR"
