#!/usr/bin/env python3
"""
Process CarbonTracker → wrfinput_d01 (ICs) + wrfbdy_d01 (BCs)
Variables: CO2_BCK, CH4_BCK (+ CO_BCK proxy)
Must run AFTER real.exe (needs wrfinput_d01 and wrfbdy_d01)
"""
import os
import numpy as np
import xarray as xr
import netCDF4 as nc
from datetime import datetime, timedelta
from scipy.interpolate import RegularGridInterpolator

WRF_GRK = "/home/igrk/WRF-GRK"
CT_DIR = os.path.join(WRF_GRK, "rawdata/carbontracker")
RUN_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/run")
WRFINPUT = os.path.join(RUN_DIR, "wrfinput_d01")
WRFBDY   = os.path.join(RUN_DIR, "wrfbdy_d01")

SIM_START = datetime(2015, 7, 28, 0)
SIM_END   = datetime(2015, 12, 1, 0)

# CO background (ppm, approximate tropospheric)
CO_BACKGROUND_PPM = 0.08  # 80 ppb

def find_ct_file(ct_dir, species, year, month):
    """Find CarbonTracker file for given species/date."""
    prefix = "CT2025" if species == "CO2" else "CTCH4_2024"
    pattern = f"{prefix}.molefrac_glb3x2_{year}-{month:02d}"
    files = sorted([f for f in os.listdir(ct_dir) if f.startswith(pattern)])
    return [os.path.join(ct_dir, f) for f in files]

def pick_ct_file(files, target_dt):
    """Pick the daily CarbonTracker file closest to target_dt."""
    target = target_dt.strftime("%Y-%m-%d")
    for f in files:
        if target in os.path.basename(f):
            return f

    dated = []
    for f in files:
        stem = os.path.basename(f).rsplit(".", 1)[0]
        date_s = stem[-10:]
        try:
            dated.append((abs((datetime.strptime(date_s, "%Y-%m-%d") - target_dt).total_seconds()), f))
        except ValueError:
            pass
    if dated:
        return min(dated)[1]
    return files[len(files) // 2]

def load_ct_field(ct_dir, species, target_dt):
    """Load CarbonTracker 3D mole fraction field nearest to target datetime."""
    files = find_ct_file(ct_dir, species, target_dt.year, target_dt.month)
    # Also check adjacent month
    if target_dt.day <= 7:
        prev = target_dt.replace(day=1) - timedelta(days=1)
        files = find_ct_file(ct_dir, species, prev.year, prev.month) + files
    if not files:
        print(f"  WARNING: No CarbonTracker file for {species} {target_dt.strftime('%Y-%m')}")
        return None, None, None, None

    # Load closest daily file. Older versions accidentally used the middle
    # file of the month, which could shift IC/BC fields by ~2 weeks.
    ds = xr.open_dataset(pick_ct_file(files, target_dt))
    var = "co2" if species == "CO2" else "ch4"
    var_names = [v for v in ds.data_vars if var in v.lower() and "total" in v.lower()]
    if not var_names:
        var_names = [v for v in ds.data_vars if var in v.lower()]
    if not var_names:
        var_names = list(ds.data_vars)
    vname = var_names[0]
    data = ds[vname].values  # (time, lev, lat, lon) or similar

    # CT uses pressure levels, get lat/lon/lev
    lat = ds["latitude"].values if "latitude" in ds else ds["lat"].values
    lon = ds["longitude"].values if "longitude" in ds else ds["lon"].values
    # Derive mid-layer pressure profile (hPa) for vertical interpolation
    if "pressure" in ds:
        # CT2025 CO2: pressure is (time, boundary, lat, lon) in Pa
        p_edge = ds["pressure"].values
        if p_edge.ndim == 4:
            p_edge_t = p_edge[0]  # (boundary, lat, lon)
        else:
            p_edge_t = p_edge
        p_mid = 0.5 * (p_edge_t[:-1] + p_edge_t[1:])  # (level, lat, lon)
        lev = p_mid.mean(axis=(1, 2)) / 100.0  # Pa → hPa, shape (level,)
    elif "at" in ds and "bt" in ds and "surf_pressure" in ds:
        # CTCH4: hybrid sigma levels. p = at + bt * psfc  (Pa)
        at = ds["at"].values       # (boundary,)
        bt = ds["bt"].values       # (boundary,)
        psfc = ds["surf_pressure"].values[0]  # (lat, lon), first time step
        p_edge = at[:, None, None] + bt[:, None, None] * psfc[None, :, :]
        p_mid = 0.5 * (p_edge[:-1] + p_edge[1:])  # (level, lat, lon)
        lev = p_mid.mean(axis=(1, 2)) / 100.0  # Pa → hPa
    elif "level" in ds:
        lev = ds["level"].values.astype(float)
    else:
        lev = None
    # Get nearest time
    if "time" in ds.dims or "time" in ds.coords:
        tidx = 0  # first time step
        data = data[tidx]
    data = data.squeeze()
    ds.close()
    return data, lat, lon, lev

def regrid_3d_to_wrf(ct_data, ct_lat, ct_lon, ct_lev,
                     wrf_lat, wrf_lon, wrf_pres):
    """
    Regrid CarbonTracker (lat,lon,lev) to WRF (south_north, west_east, bottom_top).
    ct_lev: pressure levels in hPa
    wrf_pres: WRF pressure array (bottom_top, south_north, west_east)
    """
    NY, NX = wrf_lat.shape
    NZ = wrf_pres.shape[0]

    # Ensure lat ascending
    if ct_lat[0] > ct_lat[-1]:
        ct_lat = ct_lat[::-1]
        ct_data = ct_data[::-1, :, :] if ct_data.ndim == 3 else ct_data

    result = np.zeros((NZ, NY, NX), dtype=np.float32)

    # For each WRF level, interpolate horizontally from CT
    # CT has shape (lev, lat, lon) or (lat, lon) for surface
    if ct_data.ndim == 2:
        # Surface only - replicate to all levels
        interp = RegularGridInterpolator(
            (ct_lat, ct_lon), ct_data.astype(np.float32),
            method="linear", bounds_error=False, fill_value=None
        )
        pts = np.column_stack([wrf_lat.ravel(), wrf_lon.ravel()])
        surface_val = interp(pts).reshape(NY, NX)
        for k in range(NZ):
            result[k] = surface_val
        return result

    # 3D case: (lev, lat, lon)
    n_ct_lev = ct_data.shape[0]
    # Horizontal regrid for each CT level
    ct_on_wrf = np.zeros((n_ct_lev, NY, NX), dtype=np.float32)
    for k in range(n_ct_lev):
        interp = RegularGridInterpolator(
            (ct_lat, ct_lon), ct_data[k].astype(np.float32),
            method="linear", bounds_error=False, fill_value=None
        )
        pts = np.column_stack([wrf_lat.ravel(), wrf_lon.ravel()])
        ct_on_wrf[k] = interp(pts).reshape(NY, NX)

    # Vertical interpolation to WRF levels
    if ct_lev is None:
        for k in range(NZ):
            result[k] = ct_on_wrf[0]
        return result

    # Ensure pressure decreasing (top of atm = small pressure)
    if float(ct_lev[0]) < float(ct_lev[-1]):
        ct_lev = ct_lev[::-1]
        ct_on_wrf = ct_on_wrf[::-1]

    for j in range(NY):
        for i in range(NX):
            wrf_col = wrf_pres[:, j, i]
            ct_col  = ct_on_wrf[:, j, i]
            result[:, j, i] = np.interp(wrf_col, ct_lev[::-1], ct_col[::-1])

    return result

print("Loading wrfinput and wrfbdy...")
if not os.path.exists(WRFINPUT):
    print(f"ERROR: {WRFINPUT} not found. Run real.exe first.")
    exit(1)

ds_in = xr.open_dataset(WRFINPUT)
wrf_lat  = ds_in["XLAT"].values[0]
wrf_lon  = ds_in["XLONG"].values[0]
NY, NX   = wrf_lat.shape
# WRF pressure (approximate from P + PB)
P  = ds_in["P"].values[0]   # perturbation pressure (Pa)
PB = ds_in["PB"].values[0]  # base pressure (Pa)
wrf_pres = (P + PB) / 100.0  # hPa, (bottom_top, south_north, west_east)
NZ = wrf_pres.shape[0]
ds_in.close()
print(f"  Grid: {NZ} x {NY} x {NX}")

# ---- 1. Patch wrfinput (initial conditions) ----
print(f"\nPatching wrfinput with CarbonTracker ICs for {SIM_START}...")

for species in ["CO2", "CH4"]:
    print(f"  {species}...")
    ct_data, ct_lat, ct_lon, ct_lev = load_ct_field(CT_DIR, species, SIM_START)
    if ct_data is None:
        print(f"  Skipping {species} (no data)")
        continue
    field_3d = regrid_3d_to_wrf(ct_data, ct_lat, ct_lon, ct_lev, wrf_lat, wrf_lon, wrf_pres)

    with nc.Dataset(WRFINPUT, "r+") as ds_nc:
        vname = f"{species}_BCK"
        if vname not in ds_nc.variables:
            v = ds_nc.createVariable(vname, "f4", ("Time", "bottom_top", "south_north", "west_east"))
            v.units = "ppm" if species == "CO2" else "ppb"
            v.description = f"{species} background mole fraction (CarbonTracker)"
            v.FieldType = 104
            v.MemoryOrder = "XYZ"
            v.stagger = ""
        # Convert CH4 from mol/mol to ppb if needed
        if species == "CH4" and field_3d.max() < 0.01:
            field_3d = field_3d * 1e9  # mol/mol → ppb
        elif species == "CO2" and field_3d.max() < 0.01:
            field_3d = field_3d * 1e6  # mol/mol → ppm
        ds_nc.variables[vname][0] = field_3d
    print(f"    {vname} range: {field_3d.min():.2f}–{field_3d.max():.2f}")

# CO background (simple constant)
with nc.Dataset(WRFINPUT, "r+") as ds_nc:
    if "CO_BCK" not in ds_nc.variables:
        v = ds_nc.createVariable("CO_BCK", "f4", ("Time", "bottom_top", "south_north", "west_east"))
        v.units = "ppb"
        v.description = "CO background mole fraction"
        v.FieldType = 104
        v.MemoryOrder = "XYZ"
        v.stagger = ""
    ds_nc.variables["CO_BCK"][0] = np.full((NZ, NY, NX), CO_BACKGROUND_PPM * 1000, np.float32)
print("  CO_BCK set to constant background")

# ---- 2. Patch wrfbdy (boundary conditions) ----
print(f"\nPatching wrfbdy with CarbonTracker BCs...")
if not os.path.exists(WRFBDY):
    print(f"ERROR: {WRFBDY} not found.")
    exit(1)

def make_bdy_slices(field_3d, bdy_width):
    """
    Extract boundary slices from 3D field (NZ, NY, NX).
    WRF wrfbdy convention:
      BXS = west  boundary: dims (bdy_width, bottom_top, south_north)
      BXE = east  boundary: dims (bdy_width, bottom_top, south_north)
      BYS = south boundary: dims (bdy_width, bottom_top, west_east)
      BYE = north boundary: dims (bdy_width, bottom_top, west_east)
    """
    # field_3d: (NZ, NY, NX)
    # BXS/BXE: west/east → slice along NX axis, keep NZ and NY
    bxs = np.transpose(field_3d[:, :, :bdy_width], (2, 0, 1))   # (bdy_width, NZ, NY)
    bxe = np.transpose(field_3d[:, :, -bdy_width:], (2, 0, 1))  # (bdy_width, NZ, NY)
    # BYS/BYE: south/north → slice along NY axis, keep NZ and NX
    bys = np.transpose(field_3d[:, :bdy_width, :], (1, 0, 2))   # (bdy_width, NZ, NX)
    bye = np.transpose(field_3d[:, -bdy_width:, :], (1, 0, 2))  # (bdy_width, NZ, NX)
    return bxs, bxe, bys, bye

with nc.Dataset(WRFBDY, "r+") as ds_bdy:
    n_times = ds_bdy.dimensions["Time"].size
    bdy_width = ds_bdy.dimensions["bdy_width"].size

    # Load one representative CT field per species (SIM_START).
    # CO2_BCK and CH4_BCK change slowly (~10-15 ppm / ~50 ppb over 5 months),
    # so using a constant boundary field (BT tendency = 0) is a good approximation.
    bdy_fields = {}
    for species, tracers in [("CO2", ["BCK"]), ("CH4", ["BCK"]), ("CO", ["BCK"])]:
        print(f"  {species} boundaries...")
        for tracer in tracers:
            vname_base = f"{species}_{tracer}"
            if species == "CO":
                field_3d = np.full((NZ, NY, NX), CO_BACKGROUND_PPM * 1000, np.float32)
            else:
                ct_data, ct_lat, ct_lon, ct_lev = load_ct_field(CT_DIR, species, SIM_START)
                if ct_data is None:
                    print(f"    WARNING: no CT data for {species}, skipping")
                    continue
                field_3d = regrid_3d_to_wrf(ct_data, ct_lat, ct_lon, ct_lev,
                                            wrf_lat, wrf_lon, wrf_pres)
                if species == "CH4" and field_3d.max() < 0.01:
                    field_3d = field_3d * 1e9   # mol/mol → ppb
                elif species == "CO2" and field_3d.max() < 0.01:
                    field_3d = field_3d * 1e6   # mol/mol → ppm

            bxs, bxe, bys, bye = make_bdy_slices(field_3d, bdy_width)
            print(f"    {vname_base}: range {field_3d.min():.2f}–{field_3d.max():.2f}")

            for suffix, data in [("BXS", bxs), ("BXE", bxe), ("BYS", bys), ("BYE", bye)]:
                vname = f"{vname_base}_{suffix}"
                if vname in ds_bdy.variables:
                    # Fill all time steps with the same constant field
                    for t in range(n_times):
                        ds_bdy.variables[vname][t] = data

            # Set tendency (BT) to zero — constant BCs
            for suffix in ["BTXS", "BTXE", "BTYS", "BTYE"]:
                vname = f"{vname_base}_{suffix}"
                if vname in ds_bdy.variables:
                    ds_bdy.variables[vname][:] = 0.0

print("\nCarbonTracker BC processing complete.")
