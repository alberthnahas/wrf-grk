#!/usr/bin/env python3
"""
Process SOMFFN ocean CO2 flux → wrfoce_d01
Monthly files, variable: ebio_co2oce
Units: mol km-2 hr-1
Period: 2015-07 to 2015-12
"""
import os
import numpy as np
import xarray as xr
import netCDF4 as nc
from datetime import datetime

WRF_GRK = "/home/igrk/WRF-GRK"
OCEAN_DIR = os.path.join(WRF_GRK, "rawdata/ocean")
OUT_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/input")
WRFINPUT = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/wrfinput_d01")

os.makedirs(OUT_DIR, exist_ok=True)

# SOMFFN flux unit: mol m-2 yr-1 (or mol m-2 day-1 depending on version)
# Convert mol m-2 yr-1 → mol km-2 hr-1: * 1e6 / 8760
# Convert mol m-2 day-1 → mol km-2 hr-1: * 1e6 / 24
def flux_to_wrf(flux_mol_m2_yr):
    """mol m-2 yr-1 → mol km-2 hr-1"""
    return flux_mol_m2_yr * 1e6 / 8760.0

def flux_daily_to_wrf(flux_mol_m2_day):
    """mol m-2 day-1 → mol km-2 hr-1"""
    return flux_mol_m2_day * 1e6 / 24.0

print("Loading WRF grid...")
ds_wrf = xr.open_dataset(WRFINPUT)
wrf_lat = ds_wrf["XLAT"].values[0]
wrf_lon = ds_wrf["XLONG"].values[0]
NY, NX = wrf_lat.shape
ds_wrf.close()
print(f"  Grid: {NY} x {NX}")

# Find SOMFFN file
somffn_files = [f for f in os.listdir(OCEAN_DIR) if f.endswith(".nc")]
if not somffn_files:
    print("ERROR: No SOMFFN NetCDF file found in", OCEAN_DIR)
    print("Run scripts/download/04_download_somffn.py first.")
    exit(1)

somffn_path = os.path.join(OCEAN_DIR, sorted(somffn_files)[-1])
print(f"Using: {os.path.basename(somffn_path)}")

ds = xr.open_dataset(somffn_path)
print("Variables:", list(ds.data_vars))

# Find flux variable (look for fgco2 or similar)
flux_var = None
for v in ds.data_vars:
    if "fgco2" in v.lower() or "flux" in v.lower() or "co2_flux" in v.lower() or "sfco2" in v.lower():
        flux_var = v
        break
if flux_var is None:
    flux_var = list(ds.data_vars)[0]
print(f"  Flux variable: {flux_var}")

# The MPI-SOM-FFN file encodes time as decimal years stored as "seconds since 2000-01-01"
# (a broken encoding where 1982.0 means Jan 1982, 1982.0833 means Feb 1982, etc.)
# We must use raw netCDF4 time values and compute indices manually.
import netCDF4 as _nc4
_nc4ds = _nc4.Dataset(somffn_path)
_time_raw = _nc4ds.variables["time"][:]   # decimal years, e.g. 1982.0, 1982.0833, ...
_nc4ds.close()

def get_month_index(time_raw, year, month):
    """Return index of closest time step matching year/month in decimal-year time axis."""
    target = year + (month - 1) / 12.0
    diffs = np.abs(time_raw - target)
    return int(np.argmin(diffs))

def regrid_to_wrf(data2d, src_lat, src_lon, dst_lat, dst_lon):
    from scipy.interpolate import RegularGridInterpolator
    if src_lat[0] > src_lat[-1]:
        src_lat = src_lat[::-1]
        data2d = data2d[::-1, :]
    # Handle 0-360 longitude
    if src_lon.max() > 180:
        src_lon = np.where(src_lon > 180, src_lon - 360, src_lon)
        # Re-sort
        idx = np.argsort(src_lon)
        src_lon = src_lon[idx]
        data2d = data2d[:, idx]
    data2d = np.nan_to_num(data2d, nan=0.0).astype(np.float32)
    interp = RegularGridInterpolator(
        (src_lat, src_lon), data2d,
        method="linear", bounds_error=False, fill_value=0.0
    )
    pts = np.column_stack([dst_lat.ravel(), dst_lon.ravel()])
    return interp(pts).reshape(dst_lat.shape).astype(np.float32)

# Get lat/lon from dataset
src_lat = ds["lat"].values if "lat" in ds else ds["latitude"].values
src_lon = ds["lon"].values if "lon" in ds else ds["longitude"].values

MONTHS = [7, 8, 9, 10, 11, 12]
YEAR = 2015

for mo in MONTHS:
    fname = f"wrfoce_d01_{YEAR}-{mo:02d}-01_00:00:00"
    fpath = os.path.join(OUT_DIR, fname)
    if os.path.exists(fpath):
        print(f"[SKIP] {fname}")
        continue

    # Extract monthly flux using correct decimal-year time indexing
    tidx = get_month_index(_time_raw, YEAR, mo)
    flux_mo = ds[flux_var].isel(time=tidx)
    print(f"  {YEAR}-{mo:02d}: time_raw={_time_raw[tidx]:.4f} (index {tidx})")

    data = flux_mo.values.squeeze()
    # SOMFFN stores data as (lon, lat, time); after isel we have (lon, lat).
    # RegularGridInterpolator expects (lat, lon), so transpose.
    if data.shape == (len(src_lon), len(src_lat)):
        data = data.T  # now shape (lat, lon)

    # Mask fill/missing values (SOMFFN uses 1e20 as fill)
    data = np.where(np.abs(data) > 1e10, np.nan, data)

    # Detect units and convert
    units = ds[flux_var].attrs.get("units", "")
    if "yr" in units or "year" in units:
        flux_wrf = flux_to_wrf(data)
    elif "day" in units:
        flux_wrf = flux_daily_to_wrf(data)
    else:
        # Assume mol m-2 yr-1
        flux_wrf = flux_to_wrf(data)

    # Regrid
    flux_grid = regrid_to_wrf(flux_wrf, src_lat.copy(), src_lon.copy(), wrf_lat, wrf_lon)

    # Land cells should be zero (ocean flux only)
    # wrfinput LANDMASK: 1=land, 0=water
    ds_wrf2 = xr.open_dataset(WRFINPUT)
    if "LANDMASK" in ds_wrf2:
        landmask = ds_wrf2["LANDMASK"].values[0]
        flux_grid = np.where(landmask > 0.5, 0.0, flux_grid)
    ds_wrf2.close()

    # Write NetCDF
    ds_out = nc.Dataset(fpath, "w")
    ds_out.TITLE = "WRF-Chem EMISSIONS"
    ds_out.createDimension("west_east", NX)
    ds_out.createDimension("south_north", NY)
    ds_out.createDimension("Time", 1)
    ds_out.createDimension("DateStrLen", 19)

    times_var = ds_out.createVariable("Times", "S1", ("Time", "DateStrLen"))
    t = datetime(YEAR, mo, 1)
    times_var[0, :] = list(t.strftime("%Y-%m-%d_%H:%M:%S"))

    # Variable name MUST be UPPERCASE: WRF-Chem's registry preprocessor uppercases
    # the DataName for scalar state vars (allocs.inc emits 'EBIO_CO2OCE'), and
    # NetCDF nc_inq_varid is case-sensitive. Lowercase here means CO2_OCE = 0
    # in wrfout (silent failure). See Bug 16 in docs/simulation.md.
    v = ds_out.createVariable("EBIO_CO2OCE", "f4", ("Time", "south_north", "west_east"))
    v.units = "mol km-2 hr-1"
    v.description = "Ocean CO2 air-sea flux (SOMFFN)"
    v.FieldType = 104
    v.MemoryOrder = "XY"
    v.stagger = ""
    v[0, :, :] = flux_grid

    ds_out.close()
    print(f"  Wrote {fname}")

ds.close()
print(f"\nDone. wrfoce files written to {OUT_DIR}")
