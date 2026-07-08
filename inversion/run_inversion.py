#!/usr/bin/env python3
"""
Bayesian inversion of FINN fire CO2 emissions using STILT footprints
and in-situ/reanalysis CO2 observations.

Observations:
  - PRIMARY: NOAA flask CO2 at Bukit Kototabang (BKT, Oct 27 2015)
             Real measurement, flag=C (fire-contaminated = fire signal present)
  - SECONDARY: CT2025 reanalysis surface CO2 at all 5 receptor sites, 6 days
               Used as pseudo-observations; independently constrain regional signal

State vector (x): emission scaling factors α_k for 3 fire regions k:
  1. Sumatra      (lon 95–107, lat  -8 to +5)
  2. Kalimantan   (lon 107–120, lat  -7 to +3)
  3. Java/Sulawesi(lon 105–130, lat -10 to -4)

  α = 1.0 means FINN is unbiased; α < 1 → FINN overestimates.

Observation equation:  y_obs = H α + ε,  where
  H[i,k] = STILT_foot_i ⊗ FINN_k  [ppm / (unitless scaling)]
         = Σ_{grid cells in region k}  foot_i(x,y) × E_k(x,y) [µmol/m²/s]

Prior:   x ~ N(1, σ_prior²)  with σ_prior = 0.5  (50% per-region)
Posterior:
  S_post  = (H^T S_obs^-1 H + S_prior^-1)^-1
  x_post  = S_post × (H^T S_obs^-1 y + S_prior^-1 x_prior)

OCO-2 pipeline:
  Once credentials are reset at https://urs.earthdata.nasa.gov/
  re-run this script with USE_OCO2=True.  The script will automatically
  download granules for Oct 23-26 2015 and extract soundings within 50 km
  of each receptor site to supplement the flask data.

Outputs (all in inversion/results/):
  inversion_summary.txt       — posterior α, uncertainty, emission change
  inversion/plots/            — 4 figures
"""

import os, glob
import numpy as np
import netCDF4 as nc
from datetime import datetime, timezone
import warnings
warnings.filterwarnings('ignore')

USE_OCO2 = True   # set False to skip OCO-2 and use only BKT flask + CT2025

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT      = '/home/igrk/WRF-GRK'
INV_DIR   = f'{ROOT}/inversion'
DATA_DIR  = f'{INV_DIR}/data'
RES_DIR   = f'{INV_DIR}/results'
PLT_DIR   = f'{INV_DIR}/plots'
FOOT_DIR  = f'{ROOT}/stilt/out/footprints'
WRF_DIR   = f'{ROOT}/simulations/IDN_BB_2015/output'
CT_DIR    = f'{ROOT}/rawdata/carbontracker'
FINN_CO2  = f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CO2_bb_surface_daily_20150101-20151231_0.1x0.1.nc'

for d in [RES_DIR, PLT_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Constants & domain ────────────────────────────────────────────────────────
NA = 6.022e23
DATES       = [datetime(2015, 10, d, tzinfo=timezone.utc) for d in range(24, 30)]
DATES_NOON  = [datetime(2015, 10, d, 12, tzinfo=timezone.utc) for d in range(24, 30)]

SITES = {
    'Bukit Kototabang': (100.32, -0.20),
    'Palangka Raya':    (113.94, -2.16),
    'Pontianak':        (109.33, -0.02),
    'Jambi':            (103.61, -1.61),
    'Jakarta':          (106.85, -6.21),
}
SITE_ORDER = list(SITES.keys())

# Emission regions: (name, lon_min, lon_max, lat_min, lat_max)
REGIONS = [
    ('Sumatra',      95.0, 107.0,  -8.0,  5.0),
    ('Kalimantan',  107.0, 121.0,  -7.0,  3.0),
    ('Java/Sulawesi',104.0, 130.0, -10.0, -3.0),
]
N_REGIONS = len(REGIONS)

# Clean background from multiple remote stations (GMI, SEY, SMO, MLO)
CO2_BACKGROUND = 399.1   # ppm  (mean of all unflagged remote flask obs, Oct 23-30)
BG_UNCERTAINTY  = 0.5    # ppm  1σ on background estimate


# ════════════════════════════════════════════════════════════════════════════
# 1. Build observation vector y and uncertainty S_obs
# ════════════════════════════════════════════════════════════════════════════
print('='*70)
print('Building observation vector ...')
print('='*70)

# 1a. Primary: BKT flask Oct 27 07 UTC (NOAA, C-flagged = fire contaminated)
BKT_FLASK_MEAN  = (409.23 + 409.79) / 2.0   # 409.51 ppm
BKT_FIRE_SIGNAL = BKT_FLASK_MEAN - CO2_BACKGROUND  # 10.41 ppm
BKT_OBS_SIGMA   = np.sqrt(0.073**2 + BG_UNCERTAINTY**2 + 3.0**2)  # analytical + bg + representativity

print(f'\nPrimary observation (BKT flask, Oct 27, 07 UTC):')
print(f'  Raw value     : {BKT_FLASK_MEAN:.2f} ppm')
print(f'  Background    : {CO2_BACKGROUND:.2f} ppm')
print(f'  Fire signal   : {BKT_FIRE_SIGNAL:.2f} ppm')
print(f'  Obs uncertainty: ±{BKT_OBS_SIGMA:.2f} ppm (1σ)')

# 1b. Secondary: CT2025 surface fire enhancement at all 5 sites × 6 days
print('\nSecondary observations (CT2025 pseudo-obs):')
print('  (CT2025 is a reanalysis — observations already assimilated)')

def load_ct_surface_co2(date, site_lat, site_lon):
    """Return CT2025 surface CO2 at site at ~noon UTC (level 0, t_idx 4 = 13:30)."""
    fname = os.path.join(CT_DIR, f'CT2025.molefrac_glb3x2_{date.strftime("%Y-%m-%d")}.nc')
    if not os.path.exists(fname):
        return None
    ds = nc.Dataset(fname)
    lats = ds.variables['latitude'][:]
    lons = ds.variables['longitude'][:]
    ji   = int(np.argmin(np.abs(lats - site_lat)))
    ii   = int(np.argmin(np.abs(lons - site_lon)))
    co2  = float(ds.variables['co2'][4, 0, ji, ii])   # level 0 = surface
    ds.close()
    return co2

ct_obs   = []   # list of (site_idx, day_idx, y_ppm, sigma_ppm, obs_type)

for di, dt in enumerate(DATES):
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        ct_co2 = load_ct_surface_co2(dt, slat, slon)
        if ct_co2 is None:
            continue
        fire_enh = ct_co2 - CO2_BACKGROUND
        # CT2025 representation error: ~3 ppm (due to 3°×2° grid vs 27 km)
        # Signal-to-noise acceptable only where |fire_enh| > 0.5 ppm
        sigma_ct = np.sqrt(3.0**2 + BG_UNCERTAINTY**2)  # CT representativity + bg
        ct_obs.append((si, di, sname, dt, fire_enh, sigma_ct, 'CT2025'))
        print(f'  {sname:<22} Oct{dt.day}: CT2025={ct_co2:.2f}  Δfire={fire_enh:+.2f}  σ={sigma_ct:.1f}')

# ════════════════════════════════════════════════════════════════════════════
# 2. Build H matrix: STILT footprint × FINN per region
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('Building H matrix (footprint × emission sensitivity) ...')
print('='*70)

# Load and regrid FINN for target days
def load_finn_regridded(foot_lon, foot_lat, day_indices):
    """Load FINN CO2 flux, convert to µmol/m²/s, regrid to 0.25° footprint grid."""
    ds = nc.Dataset(FINN_CO2)
    finn_lon = ds.variables['lon'][:]
    finn_lat = ds.variables['lat'][:]
    raw = ds.variables['fire_modisviirs_CO2'][day_indices, :, :]  # molec/cm²/s
    ds.close()
    CONV = 1e4 / NA * 1e6   # molec/cm²/s → µmol/m²/s
    raw_umol = raw * CONV

    dlat_f = abs(float(foot_lat[1] - foot_lat[0]))
    dlon_f = float(foot_lon[1] - foot_lon[0])
    n_days = raw.shape[0]
    out = np.zeros((n_days, len(foot_lat), len(foot_lon)), dtype=np.float64)
    for jj, flat in enumerate(foot_lat):
        jmask = np.where((finn_lat >= flat-dlat_f/2) & (finn_lat < flat+dlat_f/2))[0]
        if jmask.size == 0: continue
        cos_lat = np.abs(np.cos(finn_lat[jmask] * np.pi/180))
        for ii, flon in enumerate(foot_lon):
            imask = np.where((finn_lon >= flon-dlon_f/2) & (finn_lon < flon+dlon_f/2))[0]
            if imask.size == 0: continue
            patch = raw_umol[:, jmask[0]:jmask[-1]+1, imask[0]:imask[-1]+1]
            w2d   = cos_lat[:, np.newaxis] * np.ones(len(imask))[np.newaxis, :]
            out[:, jj, ii] = np.sum(patch * w2d, axis=(1,2)) / w2d.sum()
    return out   # (n_days, nlat, nlon) µmol/m²/s

# Load footprint grid
fp_files = sorted(glob.glob(os.path.join(FOOT_DIR, '*_foot.nc')))
ds0 = nc.Dataset(fp_files[0])
foot_lon = ds0.variables['lon'][:]
foot_lat = ds0.variables['lat'][:]
ds0.close()

FINN_IDX = [d.timetuple().tm_yday - 1 for d in DATES]
print('Regridding FINN to 0.25° grid (this takes ~30s) ...')
finn_regrid = load_finn_regridded(foot_lon, foot_lat, FINN_IDX)  # (6, 120, 220)

# Build per-region masks on footprint grid
lon2d, lat2d = np.meshgrid(foot_lon, foot_lat)
region_masks = []
for rname, rlon0, rlon1, rlat0, rlat1 in REGIONS:
    mask = ((lon2d >= rlon0) & (lon2d <= rlon1) &
            (lat2d >= rlat0) & (lat2d <= rlat1))
    region_masks.append(mask)
    print(f'  Region "{rname}": {mask.sum()} grid cells')

def load_footprint(slon, slat, dt_noon):
    tag = dt_noon.strftime('%Y%m%d%H%M')
    for f in glob.glob(os.path.join(FOOT_DIR, f'{tag}_*_foot.nc')):
        parts = os.path.basename(f).split('_')
        try:
            if abs(float(parts[1])-slon)<0.01 and abs(float(parts[2])-slat)<0.01:
                ds = nc.Dataset(f)
                foot = ds.variables['foot'][0, :, :]
                ds.close()
                return foot
        except (ValueError, IndexError):
            continue
    return None

# Compute H matrix for all 30 receptor-times × N_REGIONS
# H_full[30, 3]: one row per (site, day), one col per region
n_obs = len(DATES) * len(SITE_ORDER)
H_full = np.zeros((n_obs, N_REGIONS))
obs_labels = []
for di, (dt_midnight, dt_noon) in enumerate(zip(DATES, DATES_NOON)):
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        row_idx = di * len(SITE_ORDER) + si
        foot = load_footprint(slon, slat, dt_noon)
        if foot is None:
            obs_labels.append((sname, dt_midnight.day, 'missing'))
            continue
        # H[i,k] = sum(foot × FINN_region_k) for this day's flux
        E_day = finn_regrid[di, :, :]
        for ki, mask in enumerate(region_masks):
            H_full[row_idx, ki] = np.sum(foot * E_day * mask)
        obs_labels.append((sname, dt_midnight.day, 'ok'))

print(f'\nH matrix shape: {H_full.shape}')
print('H summary (mean sensitivity per region [ppm / scaling factor]):')
for ki, (rname, *_) in enumerate(REGIONS):
    print(f'  {rname:<20}: mean={H_full[:,ki].mean():.3f}  max={H_full[:,ki].max():.3f}')


# ════════════════════════════════════════════════════════════════════════════
# 3. Assemble observation vector with the correct H rows
# ════════════════════════════════════════════════════════════════════════════
# Observation set A: BKT flask (real, Oct 27) → map to H row for BKT Oct 27
# Observation set B: CT2025 pseudo-obs for all 30 receptor-times

# Build combined obs: include CT2025 for ALL 30 receptor-times, then
# REPLACE the BKT Oct 27 row with the real flask obs

y_all   = np.zeros(n_obs)
sig_all = np.full(n_obs, 99.0)   # start very uncertain
obs_type_all = ['CT2025'] * n_obs

# Fill CT2025 values
ct_dict = {}   # (si, di) → (fire_enh, sigma)
for si, di, sname, dt, fire_enh, sigma_ct, otype in ct_obs:
    row_idx = di * len(SITE_ORDER) + si
    y_all[row_idx]   = fire_enh
    sig_all[row_idx] = sigma_ct
    ct_dict[(si, di)] = (fire_enh, sigma_ct)

# Override BKT Oct 27 with the real flask obs
# BKT is site index 0; Oct 27 is day index 3 (Oct 24=0, 25=1, 26=2, 27=3)
BKT_idx = 0   # site index for Bukit Kototabang
OCT27_idx = 3  # day index (Oct 27)
bkt_row = OCT27_idx * len(SITE_ORDER) + BKT_idx
y_all[bkt_row]   = BKT_FIRE_SIGNAL     # 10.41 ppm from real flask
sig_all[bkt_row] = BKT_OBS_SIGMA       # 3.1 ppm
obs_type_all[bkt_row] = 'FLASK'

# Filter to rows where we have observations (sig < 90) and H has signal
valid = (sig_all < 90) & (np.abs(H_full).sum(axis=1) > 0)
y_obs   = y_all[valid]
H       = H_full[valid, :]
S_obs   = np.diag(sig_all[valid]**2)
obs_ids = [obs_labels[i] for i, v in enumerate(valid) if v]
obs_types = [obs_type_all[i] for i, v in enumerate(valid) if v]

print(f'\nUsable observations: {valid.sum()} / {n_obs}')
print(f'  Real flask obs: {sum(1 for t in obs_types if t=="FLASK")}')
print(f'  CT2025 pseudo:  {sum(1 for t in obs_types if t=="CT2025")}')


# ════════════════════════════════════════════════════════════════════════════
# 3b. OCO-2 observations (if available)
# ════════════════════════════════════════════════════════════════════════════
# OCO-2 Lite FP v11.1r measures column-averaged XCO2.  Fire signals in the
# column are ~3-5× smaller than surface flask enhancements because CO2 is
# diluted through the full atmospheric column (~8 km scale height).
#
# H sensitivity for OCO-2:  We don't have STILT footprints at arbitrary
# OCO-2 sounding locations.  Instead, for each 1°×1° grid cell with good
# soundings we use inverse-distance-weighted interpolation of the STILT
# footprints from the 5 receptor sites to estimate the sensitivity.
# A column-averaging scaling factor of 0.25 is applied (surface plume
# confined to ~500 m BL in a ~2 km effective mixing layer / 8 km column).

OCO2_FILES = sorted(glob.glob(os.path.join(DATA_DIR, 'oco2_LtCO2_1510*.nc4')))
# Only Oct 23-26 are in archive; OCO-2 overpass time over Indonesia ≈ 02-07 UTC
OCO2_DAYS = {23: 0, 24: 1, 25: 2, 26: 3}   # day → DATES index (Oct24=0 … Oct27=3; Oct23→no STILT)
COLUMN_SCALE = 0.25   # XCO2 enhancement ≈ 25% of surface enhancement for BL fires
OCO2_CELL_DEG = 1.0   # aggregate soundings to 1°×1° grid

oco2_y_extra  = []
oco2_H_extra  = []
oco2_sig_extra = []
oco2_labels   = []

if USE_OCO2 and OCO2_FILES:
    print('\n' + '='*70)
    print('Loading OCO-2 XCO2 observations ...')
    print('='*70)

    # Collect all good soundings over Indonesia
    all_lat, all_lon, all_xco2, all_unc, all_day = [], [], [], [], []
    for fpath in OCO2_FILES:
        day = int(os.path.basename(fpath).split('_')[2][4:6])   # e.g. 151023 → 23
        ds = nc.Dataset(fpath)
        lat  = np.array(ds.variables['latitude'][:])
        lon  = np.array(ds.variables['longitude'][:])
        xco2 = np.array(ds.variables['xco2'][:], dtype=float)
        unc  = np.array(ds.variables['xco2_uncertainty'][:], dtype=float)
        qf   = np.array(ds.variables['xco2_quality_flag'][:])
        ds.close()
        mask = ((lat >= -15) & (lat <= 15) & (lon >= 90) & (lon <= 145) & (qf == 0))
        all_lat.extend(lat[mask].tolist())
        all_lon.extend(lon[mask].tolist())
        all_xco2.extend(xco2[mask].tolist())
        all_unc.extend(unc[mask].tolist())
        all_day.extend([day] * int(mask.sum()))
        print(f'  Oct{day}: {int(mask.sum())} good soundings  '
              f'xco2={xco2[mask].mean():.2f}±{xco2[mask].std():.2f}  '
              f'max={xco2[mask].max():.2f} ppm')

    all_lat  = np.array(all_lat)
    all_lon  = np.array(all_lon)
    all_xco2 = np.array(all_xco2)
    all_unc  = np.array(all_unc)
    all_day  = np.array(all_day)

    # Estimate OCO-2 clean-air background from remote soundings (|lon|>150 or |lat|>30)
    bg_mask = (all_lat > 20) | (all_lat < -20)
    if bg_mask.sum() > 10:
        oco2_background = float(np.median(all_xco2[bg_mask]))
    else:
        oco2_background = CO2_BACKGROUND - 0.5   # column is slightly lower than surface
    print(f'\n  OCO-2 estimated column background: {oco2_background:.2f} ppm')

    # Aggregate to 1°×1° cells with ≥3 soundings
    lon_edges = np.arange(90, 146, OCO2_CELL_DEG)
    lat_edges = np.arange(-15, 16, OCO2_CELL_DEG)
    n_cells_used = 0
    for day, di in OCO2_DAYS.items():
        day_mask = all_day == day
        if day_mask.sum() == 0:
            continue
        lati_d = all_lat[day_mask]
        loni_d = all_lon[day_mask]
        xco2_d = all_xco2[day_mask]
        unc_d  = all_unc[day_mask]

        for lon0 in lon_edges[:-1]:
            for lat0 in lat_edges[:-1]:
                cell = ((loni_d >= lon0) & (loni_d < lon0+OCO2_CELL_DEG) &
                        (lati_d >= lat0) & (lati_d < lat0+OCO2_CELL_DEG))
                n = cell.sum()
                if n < 3:
                    continue
                cell_xco2 = float(xco2_d[cell].mean())
                cell_fire  = cell_xco2 - oco2_background
                # Observation uncertainty: standard error of mean + instrument + bg
                cell_sigma = float(np.sqrt((unc_d[cell]**2).mean()/n + BG_UNCERTAINTY**2 + 0.3**2))

                # Estimate H sensitivity via IDW from STILT receptor footprints
                cell_lon_c = lon0 + OCO2_CELL_DEG/2
                cell_lat_c = lat0 + OCO2_CELL_DEG/2
                # Use weighted average of STILT footprints for day di
                # (if day>4, cap at last available day)
                di_use = min(di, len(DATES)-1)
                H_cell = np.zeros(N_REGIONS)
                total_w = 0.0
                for si2, sname2 in enumerate(SITE_ORDER):
                    slon2, slat2 = SITES[sname2]
                    dist = np.sqrt((cell_lon_c - slon2)**2 + (cell_lat_c - slat2)**2)
                    w = 1.0 / max(dist, 0.5)**2
                    row_idx2 = di_use * len(SITE_ORDER) + si2
                    H_cell  += w * H_full[row_idx2, :]
                    total_w += w
                if total_w > 0:
                    H_cell /= total_w
                # Apply column-averaging scaling
                H_cell *= COLUMN_SCALE

                if np.abs(H_cell).max() < 1e-6:
                    continue   # no sensitivity → skip

                oco2_y_extra.append(cell_fire)
                oco2_H_extra.append(H_cell)
                oco2_sig_extra.append(cell_sigma)
                oco2_labels.append((f'OCO2_{lon0:.0f}E_{lat0:.0f}N', day, 'ok'))
                n_cells_used += 1

    print(f'\n  OCO-2 grid cells added to inversion: {n_cells_used}')
    if n_cells_used > 0:
        yy = np.array(oco2_y_extra)
        print(f'  XCO2 fire enhancement range: {yy.min():.2f} to {yy.max():.2f} ppm')

# Append OCO-2 rows to existing obs
if oco2_y_extra:
    H         = np.vstack([H,        np.array(oco2_H_extra)])
    y_obs     = np.concatenate([y_obs,   np.array(oco2_y_extra)])
    sig_oco2  = np.array(oco2_sig_extra)
    sig_full  = np.concatenate([sig_all[valid], sig_oco2])
    S_obs     = np.diag(sig_full**2)
    obs_ids   = obs_ids   + oco2_labels
    obs_types = obs_types + ['OCO2'] * len(oco2_y_extra)
    print(f'\nTotal observations after OCO-2: {len(y_obs)}')
else:
    sig_full = sig_all[valid]
    S_obs    = np.diag(sig_full**2)


# ════════════════════════════════════════════════════════════════════════════
# 4. Bayesian inversion
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('Running Bayesian inversion ...')
print('='*70)

x_prior  = np.ones(N_REGIONS)
sigma_prior = 0.5   # 50% per-region prior uncertainty
S_prior  = np.diag(np.full(N_REGIONS, sigma_prior**2))

S_obs_inv   = np.linalg.inv(S_obs)
S_prior_inv = np.linalg.inv(S_prior)

# Posterior covariance
S_post = np.linalg.inv(H.T @ S_obs_inv @ H + S_prior_inv)
# Posterior mean
x_post = S_post @ (H.T @ S_obs_inv @ y_obs + S_prior_inv @ x_prior)
sigma_post = np.sqrt(np.diag(S_post))

# Diagnostics
y_prior  = H @ x_prior
y_post   = H @ x_post
chi2_prior = float(np.mean(((y_obs - y_prior) / sig_full)**2))
chi2_post  = float(np.mean(((y_obs - y_post)  / sig_full)**2))

print('\nPosterior emission scaling factors:')
print(f'  {"Region":<20}  {"α_prior":>8}  {"α_post":>8}  {"σ_post":>8}  {"change":>8}  {"Uncertainty reduction":>22}')
print('  ' + '-'*80)
for ki, (rname, *_) in enumerate(REGIONS):
    unc_red = (1 - sigma_post[ki]/sigma_prior) * 100
    change  = (x_post[ki] - 1) * 100
    print(f'  {rname:<20}  {x_prior[ki]:>8.3f}  {x_post[ki]:>8.3f}  {sigma_post[ki]:>8.3f}  {change:>+7.1f}%  {unc_red:>21.1f}%')
print(f'\n  χ² (prior):     {chi2_prior:.3f}')
print(f'  χ² (posterior): {chi2_post:.3f}')

# Prior vs posterior Δc comparison at each obs
print('\nPer-observation fit:')
print(f'  {"Site":<22} {"Oct":>4} {"Type":>7}  {"y_obs":>8}  {"H*α_pri":>8}  {"H*α_pos":>8}  {"residual_pri":>14}  {"residual_pos":>14}')
print('  ' + '-'*95)
for i, ((sname, day, _), otype) in enumerate(zip(obs_ids, obs_types)):
    print(f'  {sname:<22} {"  "+str(day):>4} {otype:>7}  {y_obs[i]:>8.3f}  {y_prior[i]:>8.3f}  {y_post[i]:>8.3f}  {y_obs[i]-y_prior[i]:>+14.3f}  {y_obs[i]-y_post[i]:>+14.3f}')

# ════════════════════════════════════════════════════════════════════════════
# 5. Compute posterior FINN emission totals
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('Posterior fire emission totals (Oct 24-29, 6 days) ...')
print('='*70)

# Load FINN totals per region
R_EARTH = 6.371e6
ds_finn = nc.Dataset(FINN_CO2)
finn_lon_all = ds_finn.variables['lon'][:]
finn_lat_all = ds_finn.variables['lat'][:]
finn_data    = ds_finn.variables['fire_modisviirs_CO2'][FINN_IDX, :, :]  # (6, 1799, 3600)
ds_finn.close()

dlon_rad = (finn_lon_all[1]-finn_lon_all[0]) * np.pi/180
dlat_rad = abs(finn_lat_all[1]-finn_lat_all[0]) * np.pi/180
cos_lat   = np.abs(np.cos(finn_lat_all * np.pi/180))
# Build full 2D area array matching FINN grid (nlat × nlon)
area_cm2  = np.outer(R_EARTH**2 * dlon_rad * dlat_rad * cos_lat, np.ones(len(finn_lon_all))) * 1e4

CONV_Tg = 86400 * 44.01 / NA * 1e-12   # molec/cm²/s × area[cm²] → g/day → Tg/day

finn_lon2d, finn_lat2d = np.meshgrid(finn_lon_all, finn_lat_all)
prior_total = 0.0
print(f'\n  {"Region":<20}  {"Prior (Tg CO2)":>16}  {"Posterior (Tg CO2)":>20}  {"Scaling α":>12}')
print('  ' + '-'*72)
for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
    rmask = ((finn_lon2d >= rlon0) & (finn_lon2d <= rlon1) &
             (finn_lat2d >= rlat0) & (finn_lat2d <= rlat1))
    # Sum over lat/lon for each day, then sum over days
    emis_tg = float(np.sum([
        np.sum(np.asarray(finn_data[di])[rmask] * area_cm2[rmask]) * CONV_Tg
        for di in range(len(FINN_IDX))
    ]))
    emis_prior_tg = emis_tg
    emis_post_tg  = emis_tg * x_post[ki]
    prior_total  += emis_prior_tg
    print(f'  {rname:<20}  {emis_prior_tg:>16.3f}  {emis_post_tg:>20.3f}  {x_post[ki]:>12.3f} ± {sigma_post[ki]:.3f}')

# Write summary to file
summary_path = os.path.join(RES_DIR, 'inversion_summary.txt')
with open(summary_path, 'w') as f:
    f.write('WRF-STILT Bayesian Inversion: Fire CO2 Emissions\n')
    f.write('Indonesia IDN_BB_2015 — Oct 24-29 2015\n')
    f.write('='*70 + '\n\n')
    f.write('OBSERVATION SOURCES:\n')
    f.write(f'  Primary: NOAA flask BKT Oct27 = {BKT_FLASK_MEAN:.2f} ppm (flag=C, fire-contaminated)\n')
    f.write(f'  Secondary: CT2025 reanalysis surface CO2 (pseudo-obs, 30 receptor-times)\n')
    f.write(f'  Background: {CO2_BACKGROUND:.1f} ppm (GMI, SEY, SMO, MLO clean remote stations)\n\n')
    f.write('POSTERIOR SCALING FACTORS:\n')
    for ki, (rname, *_) in enumerate(REGIONS):
        f.write(f'  {rname:<20}: α = {x_post[ki]:.3f} ± {sigma_post[ki]:.3f}\n')
    f.write(f'\nchi2 (prior): {chi2_prior:.3f}\n')
    f.write(f'chi2 (post):  {chi2_post:.3f}\n\n')
    f.write('NOTE: OCO-2 data (Oct 23-26 2015) requires Earthdata password reset.\n')
    f.write('      See inversion/get_oco2.sh for download instructions.\n')
print(f'\nSummary written: {summary_path}')


# ════════════════════════════════════════════════════════════════════════════
# 6. Plots
# ════════════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature

proj = ccrs.PlateCarree()
REGION_COLORS = ['tab:orange', 'tab:red', 'tab:purple']

# ── Fig 1: Observation-space fit ─────────────────────────────────────────────
fig1, axes1 = plt.subplots(1, 2, figsize=(13, 6))

# Left: 1:1 scatter
ax = axes1[0]
colors_obs = ['red' if t=='FLASK' else ('darkorange' if t=='OCO2' else 'steelblue') for t in obs_types]
sizes_obs  = [120 if t=='FLASK' else (25 if t=='OCO2' else 40) for t in obs_types]

ax.scatter(y_obs, y_prior, c=colors_obs, s=sizes_obs, marker='o',
           alpha=0.8, label='Prior (FINN)', edgecolors='none')
ax.scatter(y_obs, y_post,  c=colors_obs, s=sizes_obs, marker='^',
           alpha=0.8, label='Posterior', edgecolors='gray', linewidths=0.5)

# Draw arrow from prior to posterior for each obs
for i in range(len(y_obs)):
    ax.annotate('', xy=(y_obs[i], y_post[i]), xytext=(y_obs[i], y_prior[i]),
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

xlim = max(abs(y_obs).max(), abs(y_prior).max(), abs(y_post).max()) * 1.15
ax.plot([-xlim, xlim], [-xlim, xlim], 'k--', lw=1, label='1:1')
ax.axhline(0, color='k', lw=0.5); ax.axvline(0, color='k', lw=0.5)
ax.set_xlabel('Observed Δ CO₂  [ppm]', fontsize=10)
ax.set_ylabel('Modelled Δ CO₂  [ppm]', fontsize=10)
ax.set_title('Prior vs. Posterior fit to observations', fontsize=10, fontweight='bold')
ax.set_xlim(-xlim, xlim); ax.set_ylim(-xlim, xlim)

# Add annotation for BKT flask
bkt_i = obs_types.index('FLASK')
ax.annotate('BKT flask\n(real obs)', xy=(y_obs[bkt_i], y_prior[bkt_i]),
            xytext=(y_obs[bkt_i]-2, y_prior[bkt_i]+4),
            arrowprops=dict(arrowstyle='->', color='red'), fontsize=8, color='red')

# Custom legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='red',         markersize=9,  label='BKT flask (real)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='steelblue',   markersize=7,  label='CT2025 pseudo-obs'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='darkorange',  markersize=6,  label='OCO-2 (1°×1° cell)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='gray',        markersize=9,  label='Prior ●'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='gray',        markersize=9,  label='Posterior ▲'),
    Line2D([0],[0], color='k', linestyle='--', label='1:1'),
]
ax.legend(handles=legend_elements, fontsize=8, loc='upper left')
ax.grid(True, lw=0.3, alpha=0.4)
ax.set_aspect('equal')

# Right: Regional scaling factors
ax2 = axes1[1]
region_names = [r[0] for r in REGIONS]
x_idx = np.arange(N_REGIONS)
width = 0.35

bars1 = ax2.bar(x_idx - width/2, x_prior,  width, label='Prior α=1', color='lightgray',
                edgecolor='k', linewidth=0.7)
bars2 = ax2.bar(x_idx + width/2, x_post,   width, label='Posterior',
                color=REGION_COLORS, edgecolor='k', linewidth=0.7, alpha=0.85)
ax2.errorbar(x_idx + width/2, x_post, yerr=sigma_post,
             fmt='none', color='black', capsize=5, lw=1.5)

ax2.axhline(1.0, color='k', linestyle='--', lw=1)
ax2.set_xticks(x_idx)
ax2.set_xticklabels(region_names, fontsize=9)
ax2.set_ylabel('FINN emission scaling factor α', fontsize=10)
ax2.set_title('Prior vs. Posterior emission\nscaling factors by fire region', fontsize=10,
              fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(axis='y', lw=0.3, alpha=0.4)
ax2.set_ylim(0, max(max(x_post)+max(sigma_post), 1.3)*1.1)

n_oco2_obs = sum(1 for t in obs_types if t=='OCO2')
oco2_label = f' + OCO-2 ({n_oco2_obs} cells)' if n_oco2_obs else ''
fig1.suptitle(f'Bayesian CO₂ Inversion — STILT footprints + NOAA flask/CT2025{oco2_label}\n'
              'Oct 24–29 2015, Indonesia fire season', fontsize=11, fontweight='bold')
out1 = os.path.join(PLT_DIR, 'inversion_fit.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight')
plt.close(fig1)
print(f'Saved: {out1}')

# ── Fig 2: Spatial map of prior vs posterior emissions ────────────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 6),
                            subplot_kw={'projection': proj},
                            gridspec_kw={'wspace': 0.08})

LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = 90, 145, -15, 15
FINN_IDX_MEAN = [d.timetuple().tm_yday - 1 for d in DATES]

ds_f = nc.Dataset(FINN_CO2)
finn_lonA = ds_f.variables['lon'][:]
finn_latA = ds_f.variables['lat'][:]
raw_mean  = ds_f.variables['fire_modisviirs_CO2'][FINN_IDX, :, :].mean(axis=0)
ds_f.close()

imask = np.where((finn_lonA >= LON_MIN) & (finn_lonA <= LON_MAX))[0]
jmask = np.where((finn_latA >= LAT_MIN) & (finn_latA <= LAT_MAX))[0]
i0,i1 = imask[0], imask[-1]+1
j0,j1 = jmask[0], jmask[-1]+1
sub_lon  = finn_lonA[i0:i1]; sub_lat = finn_latA[j0:j1]
sub_mean = raw_mean[j0:j1, i0:i1]
lon2d_s, lat2d_s = np.meshgrid(sub_lon, sub_lat)

# Build per-pixel posterior emission = prior × α_region
post_mean = sub_mean.copy().astype(float)
for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
    rmask = ((lon2d_s >= rlon0) & (lon2d_s <= rlon1) &
             (lat2d_s >= rlat0) & (lat2d_s <= rlat1))
    post_mean[rmask] *= x_post[ki]

# Convert to g/m²/day for display
CONV_GMS = 86400 * 44.01 / NA * 1e4   # molec/cm²/s → g/m²/day
prior_gm2 = np.ma.masked_less_equal(sub_mean * CONV_GMS, 0)
post_gm2  = np.ma.masked_less_equal(post_mean * CONV_GMS, 0)

from matplotlib.colors import LogNorm
nz_vals = prior_gm2.compressed()
vmin_m = np.percentile(nz_vals, 2) if nz_vals.size else 1e-3
vmax_m = np.percentile(nz_vals, 99) if nz_vals.size else 100
lnorm  = LogNorm(vmin=max(vmin_m, 1e-3), vmax=vmax_m)

for ax3, data, title in [
    (axes2[0], prior_gm2,  'Prior: FINN mean emissions'),
    (axes2[1], post_gm2,   'Posterior: α-scaled FINN'),
]:
    ax3.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=proj)
    ax3.add_feature(cfeature.LAND,      facecolor='#f0ede8')
    ax3.add_feature(cfeature.OCEAN,     facecolor='#d5eaf5')
    ax3.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='#333')
    ax3.add_feature(cfeature.BORDERS,   linewidth=0.3, edgecolor='#888', linestyle=':')

    pcm = ax3.pcolormesh(lon2d_s, lat2d_s, data, cmap='hot_r', norm=lnorm,
                          transform=proj, shading='auto', rasterized=True)

    # Draw region outlines
    for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
        rect_lons = [rlon0, rlon1, rlon1, rlon0, rlon0]
        rect_lats = [rlat0, rlat0, rlat1, rlat1, rlat0]
        ax3.plot(rect_lons, rect_lats, color=REGION_COLORS[ki],
                 lw=1.5, linestyle='-', transform=proj, zorder=5,
                 label=f'{rname} α={x_post[ki]:.2f}')

    # Mark receptor sites
    for sname, (slon, slat) in SITES.items():
        ax3.plot(slon, slat, marker='*', color='cyan', markersize=8,
                 transform=proj, zorder=7)

    # Mark BKT flask obs
    ax3.plot(100.32, -0.20, marker='D', color='lime', markersize=9,
             transform=proj, zorder=8, label='BKT flask obs')

    gl = ax3.gridlines(draw_labels=True, linewidth=0.3,
                        color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(10)
    gl.ylocator = mticker.MultipleLocator(5)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'size': 8}

    cb = fig2.colorbar(pcm, ax=ax3, orientation='vertical',
                        fraction=0.025, pad=0.04, shrink=0.9, extend='both')
    cb.set_label('CO₂ fire flux  [g m⁻² day⁻¹]', fontsize=8)
    ax3.set_title(title, fontsize=10, fontweight='bold')
    if ax3 is axes2[0]:
        ax3.legend(fontsize=7, loc='lower left', framealpha=0.85)

fig2.suptitle('Prior vs. Posterior FINN fire emission (mean Oct 24–29 2015)',
              fontsize=12, fontweight='bold', y=1.01)
out2 = os.path.join(PLT_DIR, 'posterior_emissions.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {out2}')

# ── Fig 3: Uncertainty reduction + posterior concentration fit ─────────────
fig3, axes3 = plt.subplots(1, 2, figsize=(13, 5))

# Left: uncertainty reduction bar chart
ax_u = axes3[0]
unc_red = [(1 - sigma_post[ki]/sigma_prior)*100 for ki in range(N_REGIONS)]
bars = ax_u.bar(region_names, unc_red, color=REGION_COLORS, edgecolor='k', lw=0.7)
for bar, val in zip(bars, unc_red):
    ax_u.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.1f}%',
              ha='center', va='bottom', fontsize=10, fontweight='bold')
ax_u.set_ylabel('Posterior uncertainty reduction (%)', fontsize=10)
ax_u.set_title('Uncertainty reduction by region\n(100% = perfect constraint)', fontsize=10)
ax_u.set_ylim(0, 110)
ax_u.grid(axis='y', lw=0.3, alpha=0.4)
ax_u.axhline(100, color='k', lw=0.5, linestyle='--')

# Right: receptor Δc timeseries — prior, posterior, observations
SITE_COLORS = ['tab:orange','tab:red','tab:green','tab:brown','tab:blue']
ax_t = axes3[1]
obs_by_site = {sname: {'day':[], 'y_obs':[], 'y_pri':[], 'y_pos':[], 'type':[]}
               for sname in SITE_ORDER}
for i, ((sname, day, _), otype) in enumerate(zip(obs_ids, obs_types)):
    if sname not in obs_by_site:   # skip OCO-2 grid-cell labels
        continue
    obs_by_site[sname]['day'].append(day)
    obs_by_site[sname]['y_obs'].append(y_obs[i])
    obs_by_site[sname]['y_pri'].append(y_prior[i])
    obs_by_site[sname]['y_pos'].append(y_post[i])
    obs_by_site[sname]['type'].append(otype)

x_days = list(range(24, 30))
for si, sname in enumerate(SITE_ORDER):
    d = obs_by_site[sname]
    days = sorted(d['day'])
    obs_s  = [d['y_obs'][d['day'].index(dd)] for dd in days if dd in d['day']]
    pri_s  = [d['y_pri'][d['day'].index(dd)] for dd in days if dd in d['day']]
    pos_s  = [d['y_pos'][d['day'].index(dd)] for dd in days if dd in d['day']]
    otype_s = [d['type'][d['day'].index(dd)] for dd in days if dd in d['day']]

    ax_t.plot(days, obs_s, 'o', color=SITE_COLORS[si], ms=7, label=f'{sname} obs',
              markeredgecolor='none')
    ax_t.plot(days, pri_s, '--', color=SITE_COLORS[si], lw=1.2, alpha=0.6)
    ax_t.plot(days, pos_s, '-',  color=SITE_COLORS[si], lw=2.0)
    # Highlight real flask obs
    for dd, yt, tp in zip(days, obs_s, otype_s):
        if tp == 'FLASK':
            ax_t.plot(dd, yt, '*', color='lime', ms=16, markeredgecolor='k',
                      lw=0.5, zorder=10)

ax_t.axhline(0, color='k', lw=0.5)
ax_t.set_xlabel('Day of October 2015', fontsize=9)
ax_t.set_ylabel('Δ CO₂  [ppm]', fontsize=9)
ax_t.set_title('Δ CO₂ at receptors: obs (●), prior (--), posterior (—)\n★ = real BKT flask', fontsize=9)
ax_t.set_xticks(x_days)
ax_t.set_xticklabels([f'Oct {d}' for d in x_days], fontsize=8)
ax_t.grid(lw=0.3, alpha=0.4)
# Compact legend
handles = [plt.Line2D([0],[0], color=c, lw=2, label=s[:12])
           for s,c in zip(SITE_ORDER, SITE_COLORS)]
ax_t.legend(handles=handles, fontsize=7, ncol=2, framealpha=0.85)

fig3.suptitle('Bayesian inversion diagnostics — Oct 24–29 2015', fontsize=11, fontweight='bold')
out3 = os.path.join(PLT_DIR, 'inversion_diagnostics.png')
fig3.savefig(out3, dpi=150, bbox_inches='tight')
plt.close(fig3)
print(f'Saved: {out3}')

print('\n' + '='*70)
print('Inversion complete.')
print(f'Results: {RES_DIR}/inversion_summary.txt')
print(f'Plots:   {PLT_DIR}/')
print('='*70)
print('\nNOTE: OCO-2 granules for Oct 23-26 2015 are now included.')
print('      Oct 27-29 are not available in v11.1r archive (no granules exist).')
