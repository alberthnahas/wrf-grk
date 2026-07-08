#!/usr/bin/env python3
"""
Process EDGAR v8.0 → wrfchemi_d01_YYYY-MM-DD_HH:00:00
Anthropogenic emissions: CO2, CH4, CO
Output: hourly files covering 2015-07-25 to 2015-12-01
Units: mol km-2 hr-1
Must run AFTER real.exe (needs wrfinput_d01 for grid)

CO2/CH4: EDGAR v8.0 per-sector annual files → sector-specific temporal profiles
CO:      EDGAR HTAP v3 monthly sector file → annual mean + generic profile
AWB sector excluded from CH4 (covered by FINN fire emissions).
"""
import os
import glob
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import netCDF4 as nc

WRF_GRK = "/home/igrk/WRF-GRK"
EDGAR_DIR = os.path.join(WRF_GRK, "rawdata/edgar")
OUT_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015/input")
RUN_DIR = os.path.join(WRF_GRK, "simulations/IDN_BB_2015")
WRFINPUT = os.path.join(RUN_DIR, "wrfinput_d01")

os.makedirs(OUT_DIR, exist_ok=True)

# Molecular weights (g/mol)
MW = {"CO2": 44.01, "CH4": 16.04, "CO": 28.01}

# EDGAR v8 unit: Tonnes/year/0.1°-cell → mol/km²/hr
# Cell area at 0.1° depends on latitude: A = (0.1 * π/180 * R)² * cos(lat)
# R_earth = 6371 km, Δlat = Δlon = 0.1°
def edgar_to_wrf(data_tonnes_yr, elat, mw):
    """
    Convert EDGAR v8 Tonnes/yr/0.1°-cell → mol/km²/hr.
    data_tonnes_yr: 2D array (nlat, nlon)
    elat:           1D latitude array (nlat,) in degrees
    mw:             molecular weight in g/mol
    """
    R = 6371.0  # km
    deg2rad = np.pi / 180.0
    cell_area_km2 = (0.1 * deg2rad * R) ** 2 * np.cos(np.radians(elat))  # (nlat,)
    cell_area_km2 = cell_area_km2[:, np.newaxis]  # broadcast over lon
    # 1 tonne = 1e6 g; 8760 hr/yr
    return (data_tonnes_yr * 1e6 / mw / 8760.0 / cell_area_km2).astype(np.float32)

# =============================================================================
# Temporal profiles — calibrated for Indonesia (UTC+7)
# Local 07:00 = 00 UTC; local 17:00 = 10 UTC; local 20:00 = 13 UTC
# =============================================================================

def _norm(arr):
    a = np.array(arr, dtype=np.float32)
    return a / a.mean()

MONTHLY_PROFILES = {
    "energy":      _norm([1.08,1.06,1.03,0.99,0.96,0.92,0.90,0.91,0.96,1.01,1.06,1.10]),
    "transport":   _norm([0.96,0.96,1.00,1.02,1.04,1.02,0.98,0.97,1.04,1.04,0.99,0.96]),
    "residential": _norm([1.10,1.08,1.05,1.00,0.96,0.93,0.92,0.93,0.97,1.00,1.03,1.06]),
    "industry":    _norm([1.02,1.01,1.01,1.00,1.00,0.99,0.98,0.99,1.00,1.01,1.00,1.01]),
    "agriculture": _norm([1.00]*12),
    "waste":       _norm([1.01,1.01,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.01]),
    "shipping":    _norm([1.00]*12),
    "aviation":    _norm([0.97,0.96,0.97,0.99,1.02,1.04,1.05,1.04,1.02,1.00,0.98,0.97]),
}

DIURNAL_PROFILES = {
    # Power plants: peaks at local 07:00 (00 UTC) and 18:00 (11 UTC)
    "energy":      _norm([1.3,1.2,1.1,1.0,0.9,0.9,0.8,0.8,0.9,1.0,1.1,1.3,
                          1.2,1.1,0.9,0.8,0.7,0.7,0.7,0.7,0.8,0.9,1.0,1.2]),
    # Road transport: commute peaks at local 08:00 (01 UTC) and 17:00 (10 UTC)
    "transport":   _norm([0.8,1.5,1.8,1.5,1.0,0.8,0.8,0.8,0.8,0.9,1.6,1.8,
                          1.4,0.9,0.7,0.6,0.6,0.5,0.5,0.5,0.5,0.5,0.5,0.6]),
    # Residential: cooking peaks at local 07:00 (00 UTC) and 19:00 (12 UTC)
    "residential": _norm([1.5,1.6,1.4,0.9,0.7,0.6,0.6,0.6,0.6,0.7,0.8,1.0,
                          1.6,1.7,1.4,1.0,0.7,0.6,0.5,0.5,0.5,0.5,0.6,0.8]),
    # Industry: daytime working hours local 07-17 (00-10 UTC)
    "industry":    _norm([1.1,1.3,1.3,1.3,1.2,1.2,1.1,1.1,1.0,1.0,0.8,0.7,
                          0.7,0.6,0.6,0.6,0.7,0.7,0.7,0.8,0.8,0.8,0.9,1.0]),
    "flat":        _norm([1.0]*24),
    "shipping":    _norm([0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1.1,1.1,1.1,1.1,
                          1.1,1.1,1.1,1.1,1.0,1.0,1.0,1.0,0.9,0.9,0.9,0.9]),
    "aviation":    _norm([0.5,0.4,0.5,0.7,1.0,1.3,1.5,1.6,1.5,1.4,1.3,1.2,
                          1.1,1.1,1.2,1.2,1.1,0.9,0.8,0.7,0.7,0.7,0.6,0.5]),
}

# Sector → (monthly_profile, diurnal_profile)
CO2_SECTOR_PROFILES = {
    "ENE":              ("energy",      "energy"),
    "IND":              ("industry",    "industry"),
    "CHE":              ("industry",    "industry"),
    "NMM":              ("industry",    "industry"),
    "IRO":              ("industry",    "industry"),
    "NFE":              ("industry",    "industry"),
    "NEU":              ("industry",    "flat"),
    "TRO":              ("transport",   "transport"),
    "RCO":              ("residential", "residential"),
    "REF_TRF":          ("industry",    "industry"),
    "PRO_FFF":          ("energy",      "flat"),
    "PRU_SOL":          ("industry",    "flat"),
    "AGS":              ("agriculture", "flat"),
    "SWD_INC":          ("waste",       "industry"),
    "TNR_Ship":         ("shipping",    "shipping"),
    "TNR_Aviation_LTO": ("aviation",    "aviation"),
    "TNR_Aviation_CDS": ("aviation",    "aviation"),
    "TNR_Aviation_CRS": ("aviation",    "aviation"),
    # TNR_Aviation_SPS: not available in EDGAR v8.0
    "TNR_Other":        ("transport",   "transport"),
}

CH4_SECTOR_PROFILES = {
    "ENE":              ("energy",      "energy"),
    "PRO_FFF":          ("energy",      "flat"),
    "PRO_COAL":         ("energy",      "flat"),
    "PRO_GAS":          ("energy",      "flat"),
    "PRO_OIL":          ("energy",      "flat"),
    "IND":              ("industry",    "industry"),
    "CHE":              ("industry",    "industry"),
    "IRO":              ("industry",    "industry"),
    "REF_TRF":          ("industry",    "industry"),
    "RCO":              ("residential", "residential"),
    "TRO":              ("transport",   "transport"),
    "ENF":              ("agriculture", "flat"),   # enteric fermentation
    "MNM":              ("agriculture", "flat"),   # manure management
    "AGS":              ("agriculture", "flat"),
    "SWD_LDF":          ("waste",       "flat"),   # landfill
    "SWD_INC":          ("waste",       "industry"),
    "WWT":              ("waste",       "flat"),
    # AWB excluded — covered by FINN fire emissions
    "TNR_Ship":         ("shipping",    "shipping"),
    "TNR_Aviation_LTO": ("aviation",    "aviation"),
    "TNR_Aviation_CDS": ("aviation",    "aviation"),
    "TNR_Aviation_CRS": ("aviation",    "aviation"),
    # TNR_Aviation_SPS: not available in EDGAR v8.0
    "TNR_Other":        ("transport",   "transport"),
}

# Generic fallback profiles (CO HTAP + missing sectors)
MONTHLY_GENERIC = {
    "CO2": MONTHLY_PROFILES["energy"],
    "CH4": MONTHLY_PROFILES["agriculture"],
    "CO":  _norm([1.10,1.08,1.05,0.98,0.93,0.88,0.87,0.88,0.95,1.01,1.06,1.11]),
}
DIURNAL_GENERIC = _norm([
    0.35,0.30,0.28,0.28,0.30,0.50,0.90,1.30,1.45,1.40,1.35,1.35,
    1.40,1.45,1.45,1.40,1.35,1.30,1.20,1.10,0.95,0.80,0.65,0.50,
])

# =============================================================================
# Grid utilities
# =============================================================================

def get_grid_info(wrfinput_path):
    ds = xr.open_dataset(wrfinput_path)
    lat = ds["XLAT"].values[0]
    lon = ds["XLONG"].values[0]
    ds.close()
    return lat, lon

def regrid_edgar(edgar_var, edgar_lat, edgar_lon, wrf_lat, wrf_lon):
    from scipy.interpolate import RegularGridInterpolator
    if edgar_lat[0] > edgar_lat[-1]:
        edgar_lat = edgar_lat[::-1]
        edgar_var = edgar_var[::-1, :]
    interp = RegularGridInterpolator(
        (edgar_lat, edgar_lon), np.nan_to_num(edgar_var, nan=0.0).astype(np.float32),
        method="nearest", bounds_error=False, fill_value=0.0
    )
    pts = np.column_stack([wrf_lat.ravel(), wrf_lon.ravel()])
    return interp(pts).reshape(wrf_lat.shape).astype(np.float32)

print("Loading grid from wrfinput_d01...")
wrf_lat, wrf_lon = get_grid_info(WRFINPUT)
NY, NX = wrf_lat.shape
print(f"  Grid: {NY} x {NX}")

# =============================================================================
# Load sector-resolved fields for CO2 and CH4
# =============================================================================

def load_edgar_v8_sector(species, sector, wrf_lat, wrf_lon):
    pattern = os.path.join(EDGAR_DIR, f"v8.0_FT2022_GHG_{species}_2015_{sector}_emi.nc")
    matches = glob.glob(pattern)
    if not matches:
        return None
    ds = xr.open_dataset(matches[0])
    var_names = [v for v in ds.data_vars if "emi" in v.lower() or species.lower() in v.lower()]
    if not var_names:
        var_names = list(ds.data_vars)
    data = ds[var_names[0]].values.squeeze()
    if data.ndim > 2:
        data = data.mean(axis=0)
    elat = ds["lat"].values if "lat" in ds else ds["latitude"].values
    elon = ds["lon"].values if "lon" in ds else ds["longitude"].values
    ds.close()
    return regrid_edgar(edgar_to_wrf(data, elat, MW[species]), elat, elon, wrf_lat, wrf_lon)

print("\nLoading EDGAR v8.0 sector files for CO2 and CH4...")
edgar_sectors = {}
for species, sector_profiles in [("CO2", CO2_SECTOR_PROFILES), ("CH4", CH4_SECTOR_PROFILES)]:
    edgar_sectors[species] = {}
    for sector in sector_profiles:
        field = load_edgar_v8_sector(species, sector, wrf_lat, wrf_lon)
        if field is not None:
            edgar_sectors[species][sector] = field
    n = len(edgar_sectors[species])
    print(f"  {species}: {n}/{len(sector_profiles)} sectors loaded", end="")
    if n == 0:
        print(" → will fall back to TOTALS")
    else:
        total = sum(edgar_sectors[species].values())
        print(f"  (total range: {total.min():.2e} – {total.max():.2e} mol km-2 hr-1)")

# Load CO (HTAP v3: sum all sectors, annual mean, then apply generic profile)
print("\nLoading EDGAR HTAP v3 CO...")
edgar_co_field = None
# Prefer dedicated HTAP CO file; fall back to any file with standalone "CO" in name
htap_co = os.path.join(EDGAR_DIR, "edgar_HTAPv3_2015_CO.nc")
if os.path.exists(htap_co):
    co_files = [os.path.basename(htap_co)]
else:
    co_files = [f for f in os.listdir(EDGAR_DIR)
                if f.endswith(".nc") and ("HTAPv3" in f or "HTAP" in f) and "_CO." in f]
if co_files:
    fpath = os.path.join(EDGAR_DIR, sorted(co_files)[-1])
    print(f"  CO: {os.path.basename(fpath)}")
    ds = xr.open_dataset(fpath)
    total = None
    for v in ds.data_vars:
        arr = ds[v].values
        if arr.ndim == 3:
            arr = arr.mean(axis=0)
        total = arr if total is None else total + arr
    elat = ds["lat"].values if "lat" in ds else ds["latitude"].values
    elon = ds["lon"].values if "lon" in ds else ds["longitude"].values
    ds.close()
    # HTAP CO is in ton/month; annual mean (over 12 months) × 12 → Tonnes/year equivalent
    total_tonnes_yr = total * 12.0
    edgar_co_field = regrid_edgar(edgar_to_wrf(total_tonnes_yr, elat, MW["CO"]), elat, elon, wrf_lat, wrf_lon)
    print(f"    CO range: {edgar_co_field.min():.2e} – {edgar_co_field.max():.2e} mol km-2 hr-1")
else:
    print("  WARNING: No CO HTAP file found")

# Fallback TOTALS for any species with no sector files loaded
edgar_totals = {}
for species in ["CO2", "CH4"]:
    if edgar_sectors.get(species):
        continue
    tot_files = [f for f in os.listdir(EDGAR_DIR)
                 if species in f and "TOTALS" in f and f.endswith(".nc")]
    if tot_files:
        fpath = os.path.join(EDGAR_DIR, sorted(tot_files)[-1])
        print(f"  Fallback TOTALS {species}: {os.path.basename(fpath)}")
        ds = xr.open_dataset(fpath)
        data = ds["emissions"].values.squeeze()
        elat = ds["lat"].values if "lat" in ds else ds["latitude"].values
        elon = ds["lon"].values if "lon" in ds else ds["longitude"].values
        ds.close()
        edgar_totals[species] = regrid_edgar(edgar_to_wrf(data, elat, MW[species]), elat, elon, wrf_lat, wrf_lon)

# =============================================================================
# Write hourly wrfchemi files
# =============================================================================

START = datetime(2015, 7, 25, 0)
END   = datetime(2015, 12, 1, 23)

print(f"\nWriting wrfchemi files: {START} to {END}")
current = START
count = 0

while current <= END:
    mo = current.month
    hr = current.hour
    fname = f"wrfchemi_d01_{current.strftime('%Y-%m-%d_%H:%M:%S')}"
    fpath = os.path.join(OUT_DIR, fname)

    if os.path.exists(fpath):
        current += timedelta(hours=1)
        continue

    ds = nc.Dataset(fpath, "w")
    ds.TITLE = "WRF-Chem EMISSIONS"
    ds.createDimension("west_east", NX)
    ds.createDimension("south_north", NY)
    ds.createDimension("emissions_zdim", 1)
    ds.createDimension("Time", 1)
    ds.createDimension("DateStrLen", 19)

    times_var = ds.createVariable("Times", "S1", ("Time", "DateStrLen"))
    times_var[0, :] = list(current.strftime("%Y-%m-%d_%H:%M:%S"))

    for species in ["CO2", "CH4", "CO"]:
        if species == "CO":
            if edgar_co_field is not None:
                scale = MONTHLY_GENERIC["CO"][mo - 1] * DIURNAL_GENERIC[hr]
                field = edgar_co_field * scale
            else:
                field = np.zeros((NY, NX), dtype=np.float32)

        elif edgar_sectors.get(species):
            sector_map = CO2_SECTOR_PROFILES if species == "CO2" else CH4_SECTOR_PROFILES
            field = np.zeros((NY, NX), dtype=np.float32)
            for sector, (mo_prof, di_prof) in sector_map.items():
                if sector not in edgar_sectors[species]:
                    continue
                field += (edgar_sectors[species][sector]
                          * MONTHLY_PROFILES[mo_prof][mo - 1]
                          * DIURNAL_PROFILES[di_prof][hr])

        elif species in edgar_totals:
            scale = MONTHLY_GENERIC[species][mo - 1] * DIURNAL_GENERIC[hr]
            field = edgar_totals[species] * scale

        else:
            field = np.zeros((NY, NX), dtype=np.float32)

        field = np.maximum(field, 0.0).astype(np.float32)

        for vname in [f"E_{species}", f"E_{species}TST"]:
            v = ds.createVariable(vname, "f4", ("Time", "emissions_zdim", "south_north", "west_east"))
            v.units = "mol km-2 hr-1"
            v.description = f"EDGAR {species} anthropogenic emissions"
            v.FieldType = 104
            v.MemoryOrder = "XYZ"
            v.stagger = ""
            v[0, 0, :, :] = field

    ds.close()
    count += 1
    if count % 100 == 0:
        print(f"  Written {count} files... ({current})")

    current += timedelta(hours=1)

print(f"\nDone. {count} wrfchemi files written to {OUT_DIR}")

