#!/bin/bash
# Convert WRF output files to HYSPLIT ARL format for STILT
# Uses hard-linked .nc files (colon-free) from stilt/met/wrfout_links/
#
# Usage: bash stilt/convert_wrf2arl.sh
# Output: stilt/met/2015-10-22.arl ... stilt/met/2015-10-29.arl

EXE_DIR=/home/igrk/WRF-GRK/stilt/exe
MET_DIR=/home/igrk/WRF-GRK/stilt/met
LINK_DIR=$MET_DIR/wrfout_links

echo "Starting WRF->ARL conversion"
echo "EXE_DIR : $EXE_DIR"
echo "MET_DIR : $MET_DIR"
echo "Links   : $LINK_DIR"
echo ""

for NCFILE in "$LINK_DIR"/wrfout_d01_2015-10-*.nc; do
  DATE=$(basename "$NCFILE" | sed 's/wrfout_d01_\([0-9-]*\)\.nc/\1/')
  ARLFILE="$MET_DIR/${DATE}.arl"

  if [[ -f "$ARLFILE" ]]; then
    echo "[SKIP] $DATE -- already converted ($(du -h "$ARLFILE" | cut -f1))"
    continue
  fi

  echo "[RUN ] $DATE -- converting $NCFILE ..."
  START=$(date +%s)

  # arw2arl must be run from EXE_DIR so it finds WRFDATA.CFG and ARLDATA.CFG
  # It writes output to ARLDATA.BIN in the current directory by default
  cd "$EXE_DIR"
  ./arw2arl "$NCFILE" 2>&1 | tee "$MET_DIR/${DATE}_arw2arl.log" | grep -E "Completed|ERROR|WARNING" | tail -5

  if [[ -f "$EXE_DIR/ARLDATA.BIN" ]]; then
    mv "$EXE_DIR/ARLDATA.BIN" "$ARLFILE"
    END=$(date +%s)
    echo "[DONE] $DATE -- $(du -h "$ARLFILE" | cut -f1) in $((END-START))s"
  else
    echo "[FAIL] $DATE -- ARLDATA.BIN not produced, see $MET_DIR/${DATE}_arw2arl.log"
  fi

  cd /home/igrk/WRF-GRK
done

echo ""
echo "=== Conversion complete. ARL files: ==="
ls -lh "$MET_DIR"/*.arl 2>/dev/null || echo "No ARL files found."
