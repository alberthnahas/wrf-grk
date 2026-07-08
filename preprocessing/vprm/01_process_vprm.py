#!/usr/bin/env python3
"""
Process MODIS MOD09A1 + MCD12Q1 → vprm_input_d01
VPRM requires: EVI, LSWI, EVI_MIN, EVI_MAX, LSWI_MIN, LSWI_MAX,
               VEGFRA_VPRM (8 VPRM classes), WETMAP, TANN, CPOOL
Output: 8-day files with single time frame
Period: 2015 (full year for seasonal min/max)
"""
import os
import numpy as np
import glob
import xarray as xr
import netCDF4 as nc
from datetime import datetime, timedelta

WRF_GRK = "/home/igrk/WRF-GRK"
MODIS_DIR = os.path.join(WRF_GRK, "rawdata/modis")
OUT_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/input")
WRFINPUT = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/wrfinput_d01")

os.makedirs(OUT_DIR, exist_ok=True)

# VPRM vegetation classes (8 classes for vprm_opt=1 US table, applied globally)
# MCD12Q1 IGBP → VPRM class mapping
IGBP_TO_VPRM = {
    1:  1,  # Evergreen Needleleaf → EverNeedleleaf
    2:  2,  # Evergreen Broadleaf → EverBroadleaf (tropical forests)
    3:  3,  # Deciduous Needleleaf → DecidNeedleleaf
    4:  4,  # Deciduous Broadleaf → DecidBroadleaf
    5:  5,  # Mixed Forest → MixedForest
    6:  6,  # Closed Shrublands → Shrubland
    7:  6,  # Open Shrublands → Shrubland
    8:  7,  # Woody Savannas → Savanna
    9:  7,  # Savannas → Savanna
    10: 8,  # Grasslands → Grassland
    11: 0,  # Permanent Wetlands → Water/NA
    12: 8,  # Croplands → Grassland (proxy)
    13: 0,  # Urban → no biosphere
    14: 8,  # Cropland/Natural → Grassland
    15: 0,  # Snow/Ice
    16: 0,  # Barren
    17: 0,  # Water
}

def get_wrf_grid():
    ds = xr.open_dataset(WRFINPUT)
    lat = ds["XLAT"].values[0]
    lon = ds["XLONG"].values[0]
    ds.close()
    return lat, lon

print("Loading WRF grid...")
wrf_lat, wrf_lon = get_wrf_grid()
NY, NX = wrf_lat.shape
print(f"  Grid: {NY} x {NX}")

# Check for pyhdf
try:
    from pyhdf.SD import SD, SDC
    HDF4_AVAILABLE = True
    print("  pyhdf available: will read actual MODIS HDF4 files")
except ImportError:
    HDF4_AVAILABLE = False
    print("  WARNING: pyhdf not installed. Using synthetic VPRM fields.")
    print("  Install with: pip install pyhdf")

def read_modis_hdf(filepath, sds_name):
    """Read a Scientific Dataset from MODIS HDF4 file."""
    from pyhdf.SD import SD, SDC
    hdf = SD(filepath, SDC.READ)
    sds = hdf.select(sds_name)
    data = sds.get().astype(np.float32)
    attrs = sds.attributes()
    # Apply scale factor and fill value only if scale_factor is present in file
    fill = attrs.get("_FillValue", -28672)
    scale = attrs.get("scale_factor", None)  # None means no scale in file
    offset = attrs.get("add_offset", 0.0)
    if scale is not None:
        data = np.where(data == fill, np.nan, data * scale + offset)
    else:
        # No scale factor: integer/categorical data (e.g. LC_Type1), just mask fill
        data = np.where(data == fill, np.nan, data)
    hdf.end()
    return data

def mosaic_tiles(tile_files, sds_name, tile_res=500):
    """Read and mosaic MODIS tiles into a global-ish array."""
    # This is a simplified mosaic - in production use GDAL or pymodis
    data_list = []
    lat_list = []
    lon_list = []
    for f in tile_files:
        try:
            data = read_modis_hdf(f, sds_name)
            # Extract tile info from filename e.g. h27v08
            import re
            m = re.search(r'h(\d+)v(\d+)', os.path.basename(f))
            if not m:
                continue
            h, v = int(m.group(1)), int(m.group(2))
            # MODIS sinusoidal → approximate lat/lon for tile center
            # Tile size: 10 degrees
            tile_lon_min = h * 10 - 180
            tile_lat_max = 90 - v * 10
            tile_lat_min = tile_lat_max - 10
            ny, nx = data.shape
            lats = np.linspace(tile_lat_max, tile_lat_min, ny)
            lons = np.linspace(tile_lon_min, tile_lon_min + 10, nx)
            data_list.append(data)
            lat_list.append(lats)
            lon_list.append(lons)
        except Exception as e:
            print(f"  WARNING: Could not read {os.path.basename(f)}: {e}")
    return data_list, lat_list, lon_list

def regrid_to_wrf(data_list, lat_list, lon_list, wrf_lat, wrf_lon):
    """Regrid mosaicked tiles to WRF grid."""
    from scipy.interpolate import RegularGridInterpolator
    result = np.full(wrf_lat.shape, np.nan, dtype=np.float32)
    for data, lats, lons in zip(data_list, lat_list, lon_list):
        # Only fill NaN cells with data from this tile
        # Sort lats ascending
        if lats[0] > lats[-1]:
            lats = lats[::-1]
            data = data[::-1, :]
        # Check if any WRF points fall in this tile
        mask = ((wrf_lat >= lats.min()) & (wrf_lat <= lats.max()) &
                (wrf_lon >= lons.min()) & (wrf_lon <= lons.max()) &
                np.isnan(result))
        if not np.any(mask):
            continue
        interp = RegularGridInterpolator(
            (lats, lons), np.nan_to_num(data, nan=0.0),
            method="nearest", bounds_error=False, fill_value=0.0
        )
        pts = np.column_stack([wrf_lat[mask], wrf_lon[mask]])
        result[mask] = interp(pts)
    result = np.nan_to_num(result, nan=0.0)
    return result

# 8-day periods for 2015
def get_8day_dates(year):
    d = datetime(year, 1, 1)
    dates = []
    while d.year == year:
        dates.append(d)
        d += timedelta(days=8)
    return dates

dates_2015 = get_8day_dates(2015)
# Filter to simulation range
sim_dates = [d for d in dates_2015 if d >= datetime(2015, 7, 17) and d <= datetime(2015, 12, 1)]

if HDF4_AVAILABLE:
    # --- Real MODIS processing ---
    print(f"\nProcessing {len(dates_2015)} 8-day composites for seasonal min/max...")

    all_evi   = []
    all_lswi  = []
    all_dates = []

    for d in dates_2015:
        doy = d.timetuple().tm_yday
        doy_str = f"{doy:03d}"
        year_str = str(d.year)

        # Find MOD09A1 files for this date
        pattern = os.path.join(MODIS_DIR, f"MOD09A1.A{year_str}{doy_str}.*.hdf")
        tile_files = glob.glob(pattern)
        if not tile_files:
            all_evi.append(None)
            all_lswi.append(None)
            all_dates.append(d)
            continue

        # Read bands: B1=red, B2=NIR, B5=SWIR1, B6=SWIR2
        b2_list, lat_list, lon_list = mosaic_tiles(tile_files, "sur_refl_b02")  # NIR
        b1_list, _, _ = mosaic_tiles(tile_files, "sur_refl_b01")                # Red
        b6_list, _, _ = mosaic_tiles(tile_files, "sur_refl_b06")                # SWIR

        b2_wrf = regrid_to_wrf(b2_list, lat_list, lon_list, wrf_lat, wrf_lon)
        b1_wrf = regrid_to_wrf(b1_list, lat_list, lon_list, wrf_lat, wrf_lon)
        b6_wrf = regrid_to_wrf(b6_list, lat_list, lon_list, wrf_lat, wrf_lon)

        # EVI = 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)  [using Blue=0]
        evi = 2.5 * (b2_wrf - b1_wrf) / (b2_wrf + 6*b1_wrf - 7.5*0.03 + 1 + 1e-6)
        evi = np.clip(evi, -1.0, 1.0)

        # LSWI = (NIR - SWIR) / (NIR + SWIR)
        lswi = (b2_wrf - b6_wrf) / (b2_wrf + b6_wrf + 1e-6)
        lswi = np.clip(lswi, -1.0, 1.0)

        all_evi.append(evi)
        all_lswi.append(lswi)
        all_dates.append(d)
        print(f"  {d.strftime('%Y-%m-%d')}: EVI range {evi.min():.3f}–{evi.max():.3f}")

    # Compute annual min/max
    valid_evi  = [x for x in all_evi  if x is not None]
    valid_lswi = [x for x in all_lswi if x is not None]
    evi_min  = np.minimum.reduce(valid_evi)  if valid_evi  else np.zeros((NY, NX), np.float32)
    evi_max  = np.maximum.reduce(valid_evi)  if valid_evi  else np.ones((NY, NX), np.float32) * 0.7
    lswi_min = np.minimum.reduce(valid_lswi) if valid_lswi else np.zeros((NY, NX), np.float32)
    lswi_max = np.maximum.reduce(valid_lswi) if valid_lswi else np.ones((NY, NX), np.float32) * 0.3

    # Land cover (MCD12Q1) for VEGFRA_VPRM
    lc_files = glob.glob(os.path.join(MODIS_DIR, "MCD12Q1.A2015*.hdf"))
    if lc_files:
        lc_list, lc_lat, lc_lon = mosaic_tiles(lc_files, "LC_Type1")
        lc_wrf = regrid_to_wrf(lc_list, lc_lat, lc_lon, wrf_lat, wrf_lon).astype(int)
    else:
        print("  WARNING: No MCD12Q1 files. Using broadleaf (class 2) as default.")
        lc_wrf = np.full((NY, NX), 2, dtype=int)

    # Build VEGFRA_VPRM: fraction of each class per cell (binary nearest-neighbor)
    vprm_class = np.zeros((NY, NX), dtype=int)
    for igbp, vprm in IGBP_TO_VPRM.items():
        mask = lc_wrf == igbp
        vprm_class[mask] = vprm

else:
    # --- Synthetic fields (no pyhdf) ---
    print("\nGenerating synthetic VPRM fields (pyhdf not available)...")
    # Tropical Indonesia: mostly broadleaf forest (high EVI) with seasonal signal
    # This is a placeholder - install pyhdf for real data
    evi_background = np.full((NY, NX), 0.55, dtype=np.float32)
    evi_min   = np.full((NY, NX), 0.35, dtype=np.float32)
    evi_max   = np.full((NY, NX), 0.75, dtype=np.float32)
    lswi_min  = np.full((NY, NX), 0.10, dtype=np.float32)
    lswi_max  = np.full((NY, NX), 0.45, dtype=np.float32)
    vprm_class = np.full((NY, NX), 2, dtype=int)  # EverBroadleaf

# TANN: annual mean surface temperature (K), spatially varying.
# Read from wrfinput_d01's TMN field — Noah's annual climatological deep-soil
# temperature, set by real.exe from ERA5 climatology. Realistic land range
# 280–301 K (cold mountains → warm lowlands). A constant placeholder removes
# the spatial gradient from VPRM's Q10 respiration scaling — see Bug 17 in
# docs/simulation.md.
WRFINPUT = os.path.join(os.path.dirname(OUT_DIR), "run", "wrfinput_d01")
if os.path.exists(WRFINPUT):
    with nc.Dataset(WRFINPUT) as _ds:
        tann = _ds.variables["TMN"][0].astype(np.float32)
    print(f"  T_ANN ← wrfinput TMN: min={tann.min():.2f} max={tann.max():.2f} mean={tann.mean():.2f} K")
else:
    print("  WARN: wrfinput_d01 not found; using fallback constant T_ANN = 300.15 K")
    tann = np.full((NY, NX), 300.15, dtype=np.float32)

# WETMAP: wetland fraction (small over most of Indonesia)
wetmap = np.zeros((NY, NX), dtype=np.float32)

# CPOOL: carbon pool (g C m-2), tropical forest ~20 kg C m-2
cpool = np.full((NY, NX), 20000.0, dtype=np.float32)

# --- Write vprm_input files ---
print(f"\nWriting vprm_input files for {len(sim_dates)} 8-day periods...")
for d in sim_dates:
    fname = f"vprm_input_d01_{d.strftime('%Y-%m-%d')}_00:00:00"
    fpath = os.path.join(OUT_DIR, fname)
    if os.path.exists(fpath):
        continue

    if HDF4_AVAILABLE:
        # Find matching EVI/LSWI composite
        idx = None
        for i, ad in enumerate(all_dates):
            if ad == d and all_evi[i] is not None:
                idx = i
                break
        if idx is not None:
            evi_day  = all_evi[idx]
            lswi_day = all_lswi[idx]
        else:
            evi_day  = evi_background if 'evi_background' in dir() else np.full((NY,NX), 0.55, np.float32)
            lswi_day = np.full((NY, NX), 0.3, dtype=np.float32)
    else:
        # Synthetic seasonal signal (small variation)
        doy_frac = (d.timetuple().tm_yday - 1) / 365.0
        seasonal = 0.05 * np.sin(2 * np.pi * doy_frac)
        evi_day  = evi_background + seasonal
        lswi_day = np.full((NY, NX), 0.28 + seasonal * 0.5, dtype=np.float32)

    ds_out = nc.Dataset(fpath, "w")
    ds_out.TITLE = "WRF-Chem EMISSIONS"
    ds_out.createDimension("west_east", NX)
    ds_out.createDimension("south_north", NY)
    ds_out.createDimension("num_vprm_classes", 8)
    ds_out.createDimension("Time", 1)
    ds_out.createDimension("DateStrLen", 19)

    times_var = ds_out.createVariable("Times", "S1", ("Time", "DateStrLen"))
    times_var[0, :] = list(d.strftime("%Y-%m-%d_%H:%M:%S"))

    def add_var_2d(name, data, desc, units=""):
        """Write 2D (Time, south_north, west_east) variable."""
        v = ds_out.createVariable(name, "f4", ("Time", "south_north", "west_east"))
        v.description = desc
        v.units = units
        v.FieldType = 104
        v.MemoryOrder = "XY"
        v.stagger = ""
        v[0, :, :] = data.astype(np.float32)

    def add_var_vprm(name, data, desc, units=""):
        """Write 4D (Time, num_vprm_classes, south_north, west_east) variable.
        'data' is 2D (NY, NX); the same values are broadcast across all 8 classes."""
        v = ds_out.createVariable(name, "f4", ("Time", "num_vprm_classes", "south_north", "west_east"))
        v.description = desc
        v.units = units
        v.FieldType = 104
        v.MemoryOrder = "XYZ"
        v.stagger = ""
        data_3d = np.broadcast_to(data.astype(np.float32), (8, NY, NX)).copy()
        v[0, :, :, :] = data_3d

    add_var_vprm("EVI",      evi_day,   "Enhanced Vegetation Index")
    add_var_vprm("EVI_MIN",  evi_min,   "Annual minimum EVI")
    add_var_vprm("EVI_MAX",  evi_max,   "Annual maximum EVI")
    add_var_vprm("LSWI",     lswi_day,  "Land Surface Water Index")
    add_var_vprm("LSWI_MIN", lswi_min,  "Annual minimum LSWI")
    add_var_vprm("LSWI_MAX", lswi_max,  "Annual maximum LSWI")
    add_var_2d("T_ANN",    tann,      "Annual mean temperature", "K")
    add_var_2d("WETMAP",   wetmap,    "Wetland fraction")
    add_var_2d("CPOOL",    cpool,     "Carbon pool", "g C m-2")

    # VEGFRA_VPRM: 8 vegetation class fractions (binary)
    vf = ds_out.createVariable("VEGFRA_VPRM", "f4", ("Time", "num_vprm_classes", "south_north", "west_east"))
    vf.description = "VPRM vegetation class fractions"
    vf.FieldType = 104
    vf.MemoryOrder = "XYZ"
    vf.stagger = ""
    vegfra = np.zeros((8, NY, NX), dtype=np.float32)
    for cls in range(1, 9):
        vegfra[cls - 1] = (vprm_class == cls).astype(np.float32)
    vf[0, :, :, :] = vegfra

    ds_out.close()

print(f"\nDone. vprm_input files written to {OUT_DIR}")
