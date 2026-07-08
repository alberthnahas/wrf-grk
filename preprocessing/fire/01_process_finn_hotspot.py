#!/usr/bin/env python3
"""
Process FINN v2.5 + BMKG hotspot data → wrffirechemi_d01
Hybrid approach:
  - Where FINN has signal: use FINN * SCALE_ALPHA scaled by hotspot density
  - Where FINN=0 but hotspots exist: inject baseline from hotspot density
Output: daily files with 24 hourly time frames
Units: mol km-2 hr-1
Period: 2015-07-25 to 2015-11-30
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
from datetime import datetime, timedelta

WRF_GRK = "/home/igrk/WRF-GRK"
FINN_DIR = os.path.join(WRF_GRK, "rawdata/finn_fire")
HOTSPOT_FILE = os.path.join(WRF_GRK, "rawdata/hotspot/archived_hotspot_idn.csv")
OUT_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/input")
WRFINPUT = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/wrfinput_d01")

os.makedirs(OUT_DIR, exist_ok=True)

SCALE_ALPHA = 0.5     # hotspot density scaling factor
# BMKG `archived_hotspot_idn.csv` is already pre-filtered to high-confidence
# detections; the `confidence` column is a constant marker, not a filter key.
# Use all rows.
BASELINE_CO2 = 500.0  # mol km-2 hr-1 baseline CO2 fire rate when injecting
BASELINE_CO  = 50.0   # mol km-2 hr-1 baseline CO
BASELINE_CH4 = 2.0    # mol km-2 hr-1 baseline CH4

# Diurnal fire profile (peak 13-15 LT ≈ 05-07 UTC for Indonesia ~UTC+7)
# UTC: fire peaks 05:00-09:00 UTC (noon-14:00 local)
FIRE_DIURNAL = np.array([
    0.2, 0.2, 0.2, 0.3, 0.5, 1.5,
    2.5, 2.8, 2.5, 2.0, 1.5, 1.2,
    1.0, 0.8, 0.6, 0.5, 0.4, 0.3,
    0.3, 0.3, 0.3, 0.2, 0.2, 0.2,
], dtype=np.float32)
FIRE_DIURNAL = FIRE_DIURNAL / FIRE_DIURNAL.mean()

# FINN unit: molecules cm-2 s-1 → mol km-2 hr-1
# 1 mol = 6.022e23 molecules, 1 km² = 1e10 cm², 1 hr = 3600 s
def finn_to_wrf(data_molecules_cm2_s, mw=None):
    return data_molecules_cm2_s * 3600.0 * 1e10 / 6.022e23

# --- Load WRF grid ---
print("Loading WRF grid...")
ds_wrf = xr.open_dataset(WRFINPUT)
wrf_lat = ds_wrf["XLAT"].values[0]
wrf_lon = ds_wrf["XLONG"].values[0]
NY, NX = wrf_lat.shape
ds_wrf.close()
print(f"  Grid: {NY} x {NX}")

# --- Load hotspot data ---
print("Loading hotspot data...")
hs = pd.read_csv(HOTSPOT_FILE)
# Normalize column names
hs.columns = [c.strip().lower() for c in hs.columns]
# Expected columns: lat, lon, month, day, year, confidence (or similar)
# CSV is already pre-filtered to high-confidence detections — use all 2015 rows
hs = hs[hs["year"] == 2015].copy()
print(f"  {len(hs)} hotspot detections (2015)")

# --- Load FINN data ---
print("Loading FINN data...")
finn_files = {
    "CO2": sorted([f for f in os.listdir(FINN_DIR) if "CO2" in f and f.endswith(".nc")]),
    "CO":  sorted([f for f in os.listdir(FINN_DIR) if "_CO_" in f and f.endswith(".nc")]),
    "CH4": sorted([f for f in os.listdir(FINN_DIR) if "CH4" in f and f.endswith(".nc")]),
}

finn_ds = {}
for sp, files in finn_files.items():
    if not files:
        print(f"  WARNING: No FINN file for {sp}")
        finn_ds[sp] = None
    else:
        fpath = os.path.join(FINN_DIR, files[0])
        print(f"  {sp}: {os.path.basename(fpath)}")
        finn_ds[sp] = xr.open_dataset(fpath)

# --- Regrid FINN to WRF ---
def regrid_nearest(data2d, src_lat, src_lon, dst_lat, dst_lon):
    """Nearest-neighbor regrid from regular lat/lon to WRF grid."""
    from scipy.interpolate import RegularGridInterpolator
    if src_lat[0] > src_lat[-1]:
        src_lat = src_lat[::-1]
        data2d = data2d[::-1, :]
    interp = RegularGridInterpolator(
        (src_lat, src_lon), data2d.astype(np.float32),
        method="nearest", bounds_error=False, fill_value=0.0
    )
    pts = np.column_stack([dst_lat.ravel(), dst_lon.ravel()])
    return interp(pts).reshape(dst_lat.shape).astype(np.float32)

# --- Build hotspot density grid for a given date ---
def hotspot_density(date, wrf_lat, wrf_lon, dx_km=27.0):
    """Count hotspots per WRF cell for a given date."""
    day_hs = hs[(hs["month"] == date.month) & (hs["day"] == date.day)]
    grid = np.zeros(wrf_lat.shape, dtype=np.float32)
    if len(day_hs) == 0:
        return grid
    # Bin hotspots into WRF cells using nearest neighbor
    for _, row in day_hs.iterrows():
        # Find nearest grid cell
        dist = (wrf_lat - row["lat"])**2 + (wrf_lon - row["lon"])**2
        iy, ix = np.unravel_index(np.argmin(dist), dist.shape)
        grid[iy, ix] += 1
    return grid

# --- Process each day ---
START = datetime(2015, 7, 25)
END   = datetime(2015, 11, 30)

current = START
count = 0

while current <= END:
    fname = f"wrffirechemi_d01_{current.strftime('%Y-%m-%d')}_00:00:00"
    fpath = os.path.join(OUT_DIR, fname)

    if os.path.exists(fpath):
        current += timedelta(days=1)
        continue

    # Get FINN day index (days since Jan 1 2015)
    day_of_year = (current - datetime(2015, 1, 1)).days  # 0-based

    # Load FINN fields for this day
    finn_fields = {}
    for sp in ["CO2", "CO", "CH4"]:
        ds = finn_ds[sp]
        if ds is None:
            finn_fields[sp] = np.zeros((NY, NX), dtype=np.float32)
            continue
        # Find time index
        times = ds["time"].values
        # Try to find matching day
        tidx = min(day_of_year, len(times) - 1)
        var_name = [v for v in ds.data_vars if sp in v.upper() or "emi" in v.lower()]
        if not var_name:
            var_name = list(ds.data_vars)
        vname = var_name[0]
        data = ds[vname].values[tidx].squeeze()
        # Handle fill values
        if hasattr(ds[vname], "_FillValue"):
            data = np.where(data == ds[vname]._FillValue, 0.0, data)
        data = np.nan_to_num(data, nan=0.0)
        data = np.maximum(data, 0.0)
        # Regrid
        f_lat = ds["lat"].values if "lat" in ds else ds["latitude"].values
        f_lon = ds["lon"].values if "lon" in ds else ds["longitude"].values
        finn_fields[sp] = regrid_nearest(data, f_lat, f_lon, wrf_lat, wrf_lon)
        # Convert units: molecules cm-2 s-1 → mol km-2 hr-1
        finn_fields[sp] = finn_to_wrf(finn_fields[sp])

    # Hotspot density for this day
    hs_grid = hotspot_density(current, wrf_lat, wrf_lon)
    hs_norm = hs_grid / (hs_grid.max() + 1e-10)  # 0–1

    # Hybrid: scale FINN where hotspots exist; inject where FINN=0 but hotspots
    final_fields = {}
    for sp in ["CO2", "CO", "CH4"]:
        f = finn_fields[sp].copy()
        # Scale existing FINN signal
        has_finn = f > 0
        has_hs   = hs_grid > 0
        # Where both: scale
        f = np.where(has_finn, f * (1.0 + SCALE_ALPHA * hs_norm), f)
        # Where only hotspot (FINN=0): inject baseline
        baselines = {"CO2": BASELINE_CO2, "CO": BASELINE_CO, "CH4": BASELINE_CH4}
        inject = (~has_finn) & has_hs
        f = np.where(inject, baselines[sp] * hs_norm, f)
        final_fields[sp] = f.astype(np.float32)

    # Write NetCDF with 24 hourly frames
    ds_out = nc.Dataset(fpath, "w")
    ds_out.TITLE = "WRF-Chem EMISSIONS"
    ds_out.createDimension("Time", 24)
    ds_out.createDimension("south_north", NY)
    ds_out.createDimension("west_east", NX)
    ds_out.createDimension("emissions_zdim", 1)
    ds_out.createDimension("DateStrLen", 19)

    times_var = ds_out.createVariable("Times", "S1", ("Time", "DateStrLen"))
    for h in range(24):
        t = current + timedelta(hours=h)
        times_var[h, :] = list(t.strftime("%Y-%m-%d_%H:%M:%S"))

    vmap = {"CO2": "ebu_in_co2", "CO": "ebu_in_co", "CH4": "ebu_in_ch4"}
    for sp, vname in vmap.items():
        v = ds_out.createVariable(vname, "f4", ("Time", "emissions_zdim", "south_north", "west_east"))
        v.units = "mol km-2 hr-1"
        v.description = f"Fire emissions {sp} (FINN+hotspot hybrid)"
        v.FieldType = 104
        v.MemoryOrder = "XYZ"
        v.stagger = "Z"
        daily = final_fields[sp]
        for h in range(24):
            v[h, 0, :, :] = daily * FIRE_DIURNAL[h]

    ds_out.close()
    count += 1
    if count % 10 == 0:
        print(f"  {current.strftime('%Y-%m-%d')} done ({count} files)")

    current += timedelta(days=1)

print(f"\nDone. {count} wrffirechemi files written to {OUT_DIR}")
