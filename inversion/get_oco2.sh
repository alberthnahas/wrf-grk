#!/bin/bash
# Download OCO-2 Lite FP v11.1r granules for Oct 23-26 2015
# over Indonesia (local overpass times ~02-06 UTC)
#
# REQUIRES: a valid NASA Earthdata login stored in ~/.netrc
#   (register/manage at https://urs.earthdata.nasa.gov/).
#   Credentials are read from ~/.netrc and are never stored in this repo.
#
# Run: bash inversion/get_oco2.sh

set -e
BASE="https://data.gesdisc.earthdata.nasa.gov/data/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2015"
OUTDIR="/home/igrk/WRF-GRK/inversion/data"

FILES=(
  "oco2_LtCO2_151023_B11100Ar_230525020827s.nc4"
  "oco2_LtCO2_151024_B11100Ar_230525020939s.nc4"
  "oco2_LtCO2_151025_B11100Ar_230525021120s.nc4"
  "oco2_LtCO2_151026_B11100Ar_230525021134s.nc4"
)

mkdir -p "$OUTDIR"

for f in "${FILES[@]}"; do
    echo "Downloading $f ..."
    curl -n -L \
         --netrc-file ~/.netrc \
         -c ~/.urs_cookies -b ~/.urs_cookies \
         -o "$OUTDIR/$f" \
         "$BASE/$f"
    sz=$(stat -c %s "$OUTDIR/$f" 2>/dev/null || echo 0)
    if [[ "$sz" -lt 1000 ]]; then
        echo "  WARNING: $f is only $sz bytes — download may have failed"
    else
        echo "  OK: $f  ($sz bytes)"
    fi
done

echo ""
echo "After successful download, run the inversion:"
echo "  source /home/igrk/WRF-GRK/.venv/bin/activate"
echo "  python3 /home/igrk/WRF-GRK/inversion/run_inversion.py"
