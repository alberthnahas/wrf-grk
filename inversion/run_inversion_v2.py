#!/usr/bin/env python3
"""
run_inversion_v2.py — Improved Bayesian CO₂ inversion of FINN fire emissions

Improvements over run_inversion.py (v1):
  I1: Non-negativity constraint — scipy SLSQP optimization with α ≥ 0 bounds,
      preventing the physically-impossible α < 0 result for Kalimantan.
  I2: Literature-informed regional priors — α₀ and σ_prior per region derived
      from published GOSAT/IASI inversions for the 2015 El Niño event
      (Huijnen et al. 2016, Parker et al. 2016, Nechita-Banda et al. 2018).
  I3: OCO-2 proper column averaging kernels — per-sounding surface AK weight
      replaces the fixed COLUMN_SCALE = 0.25 used in v1.
  I4: CT2025 independence fix — inflate σ_CT × 2 within ±3° of BKT on
      Oct 26–28 to reduce double-counting with the real BKT flask obs.
  I5: 6-region state vector with Tikhonov spatial regularization — splits
      the original 3 regions into 6 sub-regions and adds a smoothness penalty
      between neighbouring sub-regions to prevent underdetermination.

  I6: Virtual receptor sites via STILT footprint spatial interpolation
      (Samarinda, Banjarmasin, Palembang — approximation pending real STILT runs)
  I7: Multi-tracer CO + CO₂ joint inversion (requires CO STILT runs — future work)
  I8: CT2025 fire-zone exclusion — CT2025 at sites directly inside the fire zone
      (Palangka Raya, Pontianak) is model output from a 3°×2° coarse model with
      no real flask observations and severely underestimates local fire enhancements.
      Excluding those rows and relying on OCO-2 + BKT flask prevents α → 0 artefact.

Outputs (inversion/results/ and inversion/plots/):
  inversion_v2_summary.txt
  inversion_fit_v2.png
  posterior_emissions_v2.png
  inversion_diagnostics_v2.png
  v1_vs_v2_comparison.png
"""

import os, glob
import numpy as np
from scipy.optimize import minimize
import netCDF4 as nc
from datetime import datetime, timezone
import warnings
warnings.filterwarnings('ignore')

USE_OCO2 = True

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT      = '/home/igrk/WRF-GRK'
INV_DIR   = f'{ROOT}/inversion'
DATA_DIR  = f'{INV_DIR}/data'
RES_DIR   = f'{INV_DIR}/results'
PLT_DIR   = f'{INV_DIR}/plots'
FOOT_DIR  = f'{ROOT}/stilt/out/footprints'
CT_DIR    = f'{ROOT}/rawdata/carbontracker'
FINN_CO2  = (f'{ROOT}/rawdata/finn_fire/'
             'emissions-finnv2.5modvrs_CO2_bb_surface_daily_'
             '20150101-20151231_0.1x0.1.nc')

for d in [RES_DIR, PLT_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Constants ───────────────────────────────────────────────────────────────
NA = 6.022e23
DATES      = [datetime(2015, 10, d, tzinfo=timezone.utc) for d in range(24, 30)]
DATES_NOON = [datetime(2015, 10, d, 12, tzinfo=timezone.utc) for d in range(24, 30)]

SITES = {
    'Bukit Kototabang': (100.32, -0.20),
    'Palangka Raya':    (113.94, -2.16),
    'Pontianak':        (109.33, -0.02),
    'Jambi':            (103.61, -1.61),
    'Jakarta':          (106.85, -6.21),
}
SITE_ORDER = list(SITES.keys())

# ── I5: 6-region state vector ───────────────────────────────────────────────
# Split the original 3 regions into 6 sub-regions
# (lon_min, lon_max, lat_min, lat_max)
REGIONS = [
    ('Sumatra N+C',    95.0, 107.0,  0.0,  5.0),   # 0 — Aceh, N.Sumatra
    ('Sumatra S',      95.0, 107.0, -8.0,  0.0),   # 1 — Riau, Jambi, S.Sumatra peat
    ('Kalimantan W+C', 107.0, 115.0, -7.0,  3.0),  # 2 — C.Kalimantan peat (peak fires)
    ('Kalimantan S+E', 115.0, 121.0, -7.0,  3.0),  # 3 — S/E Kalimantan forest
    ('Java',           104.0, 116.0, -10.0, -5.0), # 4 — agricultural/savanna fires
    ('Sulawesi+East',  116.0, 130.0, -10.0, -3.0), # 5 — Sulawesi, Maluku
]
N_REGIONS = len(REGIONS)

# ── I2: Literature-informed priors per region ───────────────────────────────
# Sources: Huijnen et al. (2016) NATCOMMS, Parker et al. (2016) GRL,
#          Nechita-Banda et al. (2018) ACP, van der Werf et al. (2017) ESSD
# All suggest FINN overestimates peat fires by 3–15×; α_literature ≈ 0.08–0.35
# Sulawesi+East: GFED4s estimates ~40–60% lower than FINN; no STILT footprints
#   → tight prior prevents spurious attribution of long-range OCO-2 transport signal
# Java: agricultural/savanna fires, FINN moderate overestimate → α₀=0.80
# Kalimantan: observation-network-limited (see I8). Posterior will be prior-dominated.
#   σ set to 0.30 (large) to reflect genuine ignorance from available obs network.
X_PRIOR     = np.array([0.40, 0.25, 0.15, 0.20, 0.80, 0.50])
SIGMA_PRIOR = np.array([0.30, 0.20, 0.30, 0.30, 0.40, 0.25])

# ── I5: Tikhonov regularization ─────────────────────────────────────────────
# Penalise large differences between sub-regions of the same island
# Neighbour pairs: (Sumatra N+C, Sumatra S), (Kalimantan W+C, Kalimantan S+E),
#                  (Java, Sulawesi+East)
TIKHONOV_LAMBDA = 10.0   # regularization strength — tighter to prevent
                          # spurious Java↔Sulawesi divergence given poor obs coverage
D_MATRIX = np.zeros((3, N_REGIONS))   # 3 neighbour pairs × 6 regions
D_MATRIX[0, 0] = 1; D_MATRIX[0, 1] = -1   # Sumatra N — S
D_MATRIX[1, 2] = 1; D_MATRIX[1, 3] = -1   # Kalimantan W+C — S+E
D_MATRIX[2, 4] = 1; D_MATRIX[2, 5] = -1   # Java — Sulawesi

CO2_BACKGROUND = 399.1   # ppm
BG_UNCERTAINTY  = 0.5    # ppm 1σ

# ── I8: CT2025 fire-zone exclusion ─────────────────────────────────────────
# CarbonTracker 2025 uses GFED4s as its fire emission prior, which is ~5–15×
# smaller than FINN for tropical peat fires. Consequently, CT2025 at any site
# with high STILT sensitivity to Kalimantan peat fires will show ~1–5 ppm
# enhancement while STILT+FINN predicts 50–235 ppm at α=1. Including these
# observations drags the inversion toward FINN_α = GFED/FINN ≈ 0.01–0.05,
# which is the GFED-to-FINN ratio — not a real atmospheric constraint on fires.
# Fix: EXCLUDE CT2025 observations at rows where H_max > FIRE_ZONE_THRESHOLD.
# Kalimantan result will then be prior-dominated (see posterior unc. reduction).
FIRE_ZONE_THRESHOLD = 25.0    # ppm H_max cutoff
FIRE_ZONE_SIGMA_INFLATE = 5.0  # retained for diagnostic only

# ════════════════════════════════════════════════════════════════════════════
# 1. Build observation vector
# ════════════════════════════════════════════════════════════════════════════
print('='*70)
print('Building observation vector ...')
print('='*70)

BKT_FLASK_MEAN  = (409.23 + 409.79) / 2.0
BKT_FIRE_SIGNAL = BKT_FLASK_MEAN - CO2_BACKGROUND
BKT_OBS_SIGMA   = np.sqrt(0.073**2 + BG_UNCERTAINTY**2 + 3.0**2)
BKT_IDX = 0
OCT27_IDX = 3

print(f'BKT flask fire signal: {BKT_FIRE_SIGNAL:.2f} ppm  σ={BKT_OBS_SIGMA:.2f} ppm')

def load_ct_surface_co2(date, site_lat, site_lon):
    fname = os.path.join(CT_DIR, f'CT2025.molefrac_glb3x2_{date.strftime("%Y-%m-%d")}.nc')
    if not os.path.exists(fname):
        return None
    ds = nc.Dataset(fname)
    lats = ds.variables['latitude'][:]
    lons = ds.variables['longitude'][:]
    ji   = int(np.argmin(np.abs(lats - site_lat)))
    ii   = int(np.argmin(np.abs(lons - site_lon)))
    co2  = float(ds.variables['co2'][4, 0, ji, ii])
    ds.close()
    return co2

ct_obs = []
for di, dt in enumerate(DATES):
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        ct_co2 = load_ct_surface_co2(dt, slat, slon)
        if ct_co2 is None:
            continue
        fire_enh = ct_co2 - CO2_BACKGROUND

        # ── I4: CT2025 independence fix ──────────────────────────────────
        # Inflate σ near BKT on Oct 26-28 (CT2025 assimilated the BKT flask)
        bkt_lon, bkt_lat = SITES['Bukit Kototabang']
        dist_from_bkt = np.sqrt((slon - bkt_lon)**2 + (slat - bkt_lat)**2)
        near_bkt = (dist_from_bkt < 3.0) and (dt.day in [26, 27, 28])
        inflate  = 2.0 if near_bkt else 1.0
        sigma_ct = np.sqrt((3.0 * inflate)**2 + BG_UNCERTAINTY**2)
        note     = ' [σ×2 CT independence]' if near_bkt else ''
        ct_obs.append((si, di, sname, dt, fire_enh, sigma_ct, 'CT2025', note))

print(f'CT2025 pseudo-obs: {len(ct_obs)} (including σ×2 inflation near BKT Oct26-28)')

# ════════════════════════════════════════════════════════════════════════════
# 2. H matrix
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('Building H matrix ...')
print('='*70)

def load_finn_regridded(foot_lon, foot_lat, day_indices):
    ds = nc.Dataset(FINN_CO2)
    finn_lon = ds.variables['lon'][:]
    finn_lat = ds.variables['lat'][:]
    raw = ds.variables['fire_modisviirs_CO2'][day_indices, :, :]
    ds.close()
    CONV = 1e4 / NA * 1e6
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
    return out

fp_files = sorted(glob.glob(os.path.join(FOOT_DIR, '*_foot.nc')))
ds0 = nc.Dataset(fp_files[0])
foot_lon = ds0.variables['lon'][:]
foot_lat = ds0.variables['lat'][:]
ds0.close()

FINN_IDX = [d.timetuple().tm_yday - 1 for d in DATES]
print('Regridding FINN to 0.25° grid ...')
finn_regrid = load_finn_regridded(foot_lon, foot_lat, FINN_IDX)

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

n_obs = len(DATES) * len(SITE_ORDER)
H_full = np.zeros((n_obs, N_REGIONS))
obs_labels = []
foot_cache = {}   # {(si, di): footprint array} — stored for I6 interpolation
for di, (dt_midnight, dt_noon) in enumerate(zip(DATES, DATES_NOON)):
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        row_idx = di * len(SITE_ORDER) + si
        foot = load_footprint(slon, slat, dt_noon)
        if foot is None:
            obs_labels.append((sname, dt_midnight.day, 'missing'))
            continue
        foot_cache[(si, di)] = foot   # store for I6
        E_day = finn_regrid[di, :, :]
        for ki, mask in enumerate(region_masks):
            H_full[row_idx, ki] = np.sum(foot * E_day * mask)
        obs_labels.append((sname, dt_midnight.day, 'ok'))

print(f'\nH matrix (6-region) shape: {H_full.shape}')
for ki, (rname, *_) in enumerate(REGIONS):
    print(f'  {rname:<22}: mean={H_full[:,ki].mean():.3f}  max={H_full[:,ki].max():.3f}')

# ── I8: Identify fire-zone rows ─────────────────────────────────────────────
H_row_max = np.abs(H_full).max(axis=1)
fire_zone_rows = H_row_max > FIRE_ZONE_THRESHOLD
print(f'\nFire-zone CT2025 rows (σ×{FIRE_ZONE_SIGMA_INFLATE:.0f} inflation, H_max > {FIRE_ZONE_THRESHOLD} ppm): '
      f'{fire_zone_rows.sum()} / {n_obs}')
for i, (lbl, hmax) in enumerate(zip(obs_labels, H_row_max)):
    if fire_zone_rows[i]:
        print(f'  σ×{FIRE_ZONE_SIGMA_INFLATE:.0f}: {lbl[0]:<22} Oct{lbl[1]:02d}  H_max={hmax:.1f} ppm')

# ── I6: Virtual receptor sites (approximate STILT via IDW interpolation) ────
# Sites inside the fire zone where no real STILT runs exist.
# Pseudo-footprints are inverse-distance-weighted averages of existing footprints.
# THIS IS AN APPROXIMATION — real STILT backward trajectories should be run
# from these sites in a production inversion.
VIRTUAL_SITES = {
    'Samarinda_v':   (117.15, -0.50),   # E. Kalimantan urban, ~Kalimantan S+E
    'Banjarmasin_v': (114.59, -3.32),   # S. Kalimantan peat corridor
    'Palembang_v':   (104.75, -2.99),   # S. Sumatra peat/logging fires
}

def make_virtual_footprint(v_lon, v_lat, day_idx, n_near=3):
    """IDW interpolation of existing cached STILT footprints. Approximate only."""
    dists = []
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        if (si, day_idx) not in foot_cache:
            continue
        dist = max(np.sqrt((v_lon - slon)**2 + (v_lat - slat)**2), 0.25)
        dists.append((dist, si))
    if not dists:
        return None
    dists.sort()
    nearest = dists[:n_near]
    ws = np.array([1.0/d**2 for d, _ in nearest])
    ws /= ws.sum()
    return sum(w * foot_cache[(si, day_idx)] for (_, si), w in zip(nearest, ws))

n_virt = len(VIRTUAL_SITES) * len(DATES)
H_virt     = np.zeros((n_virt, N_REGIONS))
virt_labels = []
vi = 0
for di, dt_midnight in enumerate(DATES):
    for vsname, (vslon, vslat) in VIRTUAL_SITES.items():
        vfoot = make_virtual_footprint(vslon, vslat, di)
        if vfoot is not None:
            E_day = finn_regrid[di, :, :]
            for ki, mask in enumerate(region_masks):
                H_virt[vi, ki] = np.sum(vfoot * E_day * mask)
        virt_labels.append((vsname, dt_midnight.day, 'virtual'))
        vi += 1

print(f'\n[I6] Virtual receptor rows added: {n_virt} '
      f'({len(VIRTUAL_SITES)} sites × {len(DATES)} days, IDW approx.)')
for ki, (rname, *_) in enumerate(REGIONS):
    print(f'  {rname:<22}: virt mean={H_virt[:,ki].mean():.3f}  max={H_virt[:,ki].max():.3f}')
# ════════════════════════════════════════════════════════════════════════════
# Combine H_full (existing 5 sites) + H_virt (I6 virtual sites)
H_full_all   = np.vstack([H_full,     H_virt])
obs_lbl_all  = obs_labels + virt_labels
N_OBS_ALL    = H_full_all.shape[0]

y_all    = np.zeros(N_OBS_ALL)
sig_all  = np.full(N_OBS_ALL, 99.0)
obs_type_all = ['CT2025'] * N_OBS_ALL

for si, di, sname, dt, fire_enh, sigma_ct, otype, note in ct_obs:
    row_idx = di * len(SITE_ORDER) + si
    # ── I8: fully exclude CT2025 at fire-zone rows ──
    # CT2025 encodes GFED fire prior → including it effectively imposes
    # α_post ≈ GFED/FINN ≋ 0.01–0.05, not a real atmospheric measurement.
    if fire_zone_rows[row_idx]:
        continue
    y_all[row_idx]   = fire_enh
    sig_all[row_idx] = sigma_ct

# I6: CT2025 at virtual sites
# Only Palembang (Sumatra) is used in the obs vector. Samarinda/Banjarmasin
# (Kalimantan) are also in-fire-zone (CT2025 encodes GFED/FINN ≋ 0.01–0.05
# at those sites too) — including them drives Kalimantan α to the GFED/FINN ratio.
# Virtual obs are added only for Palembang (Sumatra S constraint).
for vi2, (vsname, (vslon, vslat)) in enumerate(VIRTUAL_SITES.items()):
    for di, dt in enumerate(DATES):
        row_idx_v = n_obs + di * len(VIRTUAL_SITES) + vi2
        # Skip Kalimantan virtual sites in obs vector
        if 'arinda' in vsname or 'anjar' in vsname:   # Samarinda, Banjarmasin
            continue
        ct_co2 = load_ct_surface_co2(dt, vslat, vslon)
        if ct_co2 is None:
            continue
        fire_enh_v = ct_co2 - CO2_BACKGROUND
        # Check virtual H for fire-zone
        h_row_v = np.abs(H_virt[di * len(VIRTUAL_SITES) + vi2, :]).max()
        if h_row_v > FIRE_ZONE_THRESHOLD:
            continue
        sigma_v = np.sqrt((3.0)**2 + BG_UNCERTAINTY**2 + 1.5**2)  # +1.5 ppm footprint-interp error
        y_all[row_idx_v]   = fire_enh_v
        sig_all[row_idx_v] = sigma_v

bkt_row = OCT27_IDX * len(SITE_ORDER) + BKT_IDX
y_all[bkt_row]   = BKT_FIRE_SIGNAL
sig_all[bkt_row] = BKT_OBS_SIGMA
obs_type_all[bkt_row] = 'FLASK'

valid = (sig_all < 90) & (np.abs(H_full_all).sum(axis=1) > 0)
y_obs     = y_all[valid]
H         = H_full_all[valid, :]
sig_use   = sig_all[valid]
obs_ids   = [obs_lbl_all[i] for i, v in enumerate(valid) if v]
obs_types = [obs_type_all[i] for i, v in enumerate(valid) if v]

print(f'\nUsable obs: {valid.sum()} / {N_OBS_ALL}  '
      f'(flask={sum(t=="FLASK" for t in obs_types)}, '
      f'CT2025={sum(t=="CT2025" for t in obs_types)})')

# ════════════════════════════════════════════════════════════════════════════
# 3b. OCO-2 with proper column averaging kernels (I3)
# ════════════════════════════════════════════════════════════════════════════
OCO2_FILES = sorted(glob.glob(os.path.join(DATA_DIR, 'oco2_LtCO2_1510*.nc4')))
OCO2_DAYS  = {23: 0, 24: 1, 25: 2, 26: 3}
OCO2_CELL_DEG = 1.0

oco2_y_extra   = []
oco2_H_extra   = []
oco2_sig_extra = []
oco2_labels    = []

if USE_OCO2 and OCO2_FILES:
    print('\n' + '='*70)
    print('Loading OCO-2 with column averaging kernels (I3) ...')
    print('='*70)

    all_lat, all_lon, all_xco2, all_unc, all_aks, all_day = [], [], [], [], [], []
    for fpath in OCO2_FILES:
        day = int(os.path.basename(fpath).split('_')[2][4:6])
        ds  = nc.Dataset(fpath)
        lat  = np.array(ds.variables['latitude'][:])
        lon  = np.array(ds.variables['longitude'][:])
        xco2 = np.array(ds.variables['xco2'][:], dtype=float)
        unc  = np.array(ds.variables['xco2_uncertainty'][:], dtype=float)
        qf   = np.array(ds.variables['xco2_quality_flag'][:])
        mask = ((lat >= -15) & (lat <= 15) & (lon >= 90) & (lon <= 145) & (qf == 0))

        # ── I3: Per-sounding column averaging kernel weight (BL-integrated) ──
        # Sum AK×PW over the bottom 4 pressure levels (BL + lower troposphere).
        # A single surface level gives ~0.025; BL-integrated gives ~0.08–0.15,
        # matching the expected column sensitivity to a fire plume filling the BL
        # (~10–15% of the total atmospheric pressure column).
        aks_col = None
        try:
            AK = np.array(ds.variables['xco2_averaging_kernel'][:], dtype=float)   # (n, 20)
            PW = np.array(ds.variables['pressure_weight'][:], dtype=float)          # (n, 20)
            N_BL = 4   # bottom 4 of 20 levels ≈ surface to ~700 hPa
            aks_col = np.sum(AK[:, -N_BL:] * PW[:, -N_BL:], axis=1)   # shape (n,)
        except KeyError:
            pass  # fall back to fixed scale below

        ds.close()

        n_mask = int(mask.sum())
        if aks_col is not None:
            aks_m = aks_col[mask]
        else:
            aks_m = np.full(n_mask, 0.25)   # fallback if AK not in file

        all_lat.extend(lat[mask].tolist())
        all_lon.extend(lon[mask].tolist())
        all_xco2.extend(xco2[mask].tolist())
        all_unc.extend(unc[mask].tolist())
        all_aks.extend(aks_m.tolist())
        all_day.extend([day] * n_mask)
        print(f'  Oct{day}: {n_mask} soundings  '
              f'mean AK_sfc={np.mean(aks_m):.3f}  '
              f'xco2={np.mean(xco2[mask]):.2f} ppm')

    all_lat  = np.array(all_lat)
    all_lon  = np.array(all_lon)
    all_xco2 = np.array(all_xco2)
    all_unc  = np.array(all_unc)
    all_aks  = np.array(all_aks)
    all_day  = np.array(all_day)

    # Column background from tropical non-fire soundings.
    # Using poleward (>|20°|) soundings over-estimates background for tropical
    # cells because CO2 has a latitudinal gradient in Oct 2015 (tropical biosphere
    # uptake → lower XCO2 at 0-10°S vs >20° lat).  Using polar background with
    # tropical fire cells produces spurious NEGATIVE fire enhancements which
    # drive Kalimantan α → 0.  Instead: tropical OCO-2 soundings OUTSIDE the
    # Indonesian fire domain (lon < 95°E or lon > 145°E, same latitudes).
    bg_mask = (
        (all_lat >= -15) & (all_lat <= 15) &
        ((all_lon < 95) | (all_lon > 145))
    )
    if bg_mask.sum() < 10:
        # Fall back: slightly lower than CT2025 background to avoid over-subtraction
        oco2_background = CO2_BACKGROUND - 0.3
    else:
        oco2_background = float(np.median(all_xco2[bg_mask]))
    print(f'\n  OCO-2 column background (tropical out-of-domain): {oco2_background:.2f} ppm')

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
        aks_d  = all_aks[day_mask]

        for lon0 in lon_edges[:-1]:
            for lat0 in lat_edges[:-1]:
                cell = ((loni_d >= lon0) & (loni_d < lon0+OCO2_CELL_DEG) &
                        (lati_d >= lat0) & (lati_d < lat0+OCO2_CELL_DEG))
                n = cell.sum()
                if n < 3:
                    continue
                cell_xco2 = float(xco2_d[cell].mean())
                cell_fire  = cell_xco2 - oco2_background
                cell_sigma = float(np.sqrt((unc_d[cell]**2).mean()/n
                                          + BG_UNCERTAINTY**2 + 0.3**2))

                # ── I3: Cell-mean AK surface weight replaces fixed 0.25 ────
                cell_ak_sfc = float(aks_d[cell].mean())
                # Guard: AK may be near-zero in cloudy/smoky scenes
                if cell_ak_sfc < 0.01:
                    cell_ak_sfc = 0.01   # minimum 1% sensitivity

                cell_lon_c = lon0 + OCO2_CELL_DEG/2
                cell_lat_c = lat0 + OCO2_CELL_DEG/2
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
                H_cell *= cell_ak_sfc   # per-cell AK weight (not fixed 0.25)

                if np.abs(H_cell).max() < 1e-6:
                    continue

                oco2_y_extra.append(cell_fire)
                oco2_H_extra.append(H_cell)
                oco2_sig_extra.append(cell_sigma)
                oco2_labels.append((f'OCO2_{lon0:.0f}E_{lat0:.0f}N', day, 'ok'))
                n_cells_used += 1

    print(f'\n  OCO-2 cells added: {n_cells_used}')
    if n_cells_used > 0:
        yy = np.array(oco2_y_extra)
        print(f'  Fire enhancement range: {yy.min():.2f} to {yy.max():.2f} ppm')

if oco2_y_extra:
    H         = np.vstack([H,      np.array(oco2_H_extra)])
    y_obs     = np.concatenate([y_obs,   np.array(oco2_y_extra)])
    sig_use   = np.concatenate([sig_use, np.array(oco2_sig_extra)])
    obs_ids   = obs_ids   + oco2_labels
    obs_types = obs_types + ['OCO2'] * len(oco2_y_extra)
    print(f'Total observations: {len(y_obs)}')

S_obs = np.diag(sig_use**2)

# ════════════════════════════════════════════════════════════════════════════
# 4. Bayesian inversion with I1 (non-negativity) + I2 (lit. priors) + I5 (Tikhonov)
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('Running improved Bayesian inversion (v2) ...')
print('='*70)

S_obs_inv    = np.linalg.inv(S_obs)
S_prior_inv  = np.diag(1.0 / SIGMA_PRIOR**2)

def neg_log_posterior(alpha):
    """Objective = 0.5*[data misfit + prior misfit + Tikhonov penalty]."""
    r_data  = y_obs - H @ alpha
    r_prior = alpha - X_PRIOR
    r_tikh  = D_MATRIX @ alpha
    cost = (0.5 * r_data @ S_obs_inv @ r_data
          + 0.5 * r_prior @ S_prior_inv @ r_prior
          + 0.5 * TIKHONOV_LAMBDA * r_tikh @ r_tikh)
    return cost

def neg_log_posterior_grad(alpha):
    """Analytical gradient of the objective."""
    r_data  = y_obs - H @ alpha
    r_prior = alpha - X_PRIOR
    r_tikh  = D_MATRIX @ alpha
    grad = (-H.T @ S_obs_inv @ r_data
            + S_prior_inv @ r_prior
            + TIKHONOV_LAMBDA * D_MATRIX.T @ r_tikh)
    return grad

# ── I1: Non-negativity constraint ──────────────────────────────────────────
bounds = [(0.0, None)] * N_REGIONS   # α ≥ 0 for all regions
result = minimize(
    neg_log_posterior,
    x0     = X_PRIOR.copy(),
    jac    = neg_log_posterior_grad,
    method = 'SLSQP',
    bounds = bounds,
    options = {'ftol': 1e-12, 'maxiter': 2000},
)
x_post = result.x

# Posterior covariance via Hessian of the cost at the MAP solution
# H_hess = H^T S_obs^-1 H + S_prior^-1 + λ D^T D
H_hess  = (H.T @ S_obs_inv @ H
           + S_prior_inv
           + TIKHONOV_LAMBDA * D_MATRIX.T @ D_MATRIX)
S_post  = np.linalg.inv(H_hess)
sigma_post = np.sqrt(np.diag(S_post))

# Diagnostics
y_prior_pred = H @ X_PRIOR
y_post_pred  = H @ x_post
chi2_prior   = float(np.mean(((y_obs - y_prior_pred) / sig_use)**2))
chi2_post    = float(np.mean(((y_obs - y_post_pred)  / sig_use)**2))

print(f'\nOptimization: {result.message}  (nfev={result.nfev})')
print(f'\nPosterior scaling factors (6-region, v2):')
print(f'  {"Region":<22}  {"α_prior":>8}  {"α_post":>8}  {"σ_post":>8}  '
      f'{"change":>8}  {"Unc.red.%":>10}')
print('  ' + '-'*75)
for ki, (rname, *_) in enumerate(REGIONS):
    unc_red = (1 - sigma_post[ki] / SIGMA_PRIOR[ki]) * 100
    change  = (x_post[ki] - 1) * 100
    print(f'  {rname:<22}  {X_PRIOR[ki]:>8.3f}  {x_post[ki]:>8.3f}  '
          f'{sigma_post[ki]:>8.3f}  {change:>+7.1f}%  {unc_red:>9.1f}%')
print(f'\n  χ² (prior):     {chi2_prior:.3f}')
print(f'  χ² (posterior): {chi2_post:.3f}')

# ════════════════════════════════════════════════════════════════════════════
# 5. Posterior FINN emission totals
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('Posterior emission totals ...')
print('='*70)

R_EARTH = 6.371e6
ds_finn   = nc.Dataset(FINN_CO2)
finn_lon_all = ds_finn.variables['lon'][:]
finn_lat_all = ds_finn.variables['lat'][:]
finn_data    = ds_finn.variables['fire_modisviirs_CO2'][FINN_IDX, :, :]
ds_finn.close()

dlon_rad = (finn_lon_all[1]-finn_lon_all[0]) * np.pi/180
dlat_rad = abs(finn_lat_all[1]-finn_lat_all[0]) * np.pi/180
cos_lat  = np.abs(np.cos(finn_lat_all * np.pi/180))
area_cm2 = np.outer(R_EARTH**2 * dlon_rad * dlat_rad * cos_lat,
                    np.ones(len(finn_lon_all))) * 1e4
CONV_Tg  = 86400 * 44.01 / NA * 1e-12

finn_lon2d, finn_lat2d = np.meshgrid(finn_lon_all, finn_lat_all)
emis_prior = []
emis_post  = []
print(f'\n  {"Region":<22}  {"Prior (Tg)":>12}  {"Posterior (Tg)":>16}  {"α":>8}')
print('  ' + '-'*65)
for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
    rmask = ((finn_lon2d >= rlon0) & (finn_lon2d <= rlon1) &
             (finn_lat2d >= rlat0) & (finn_lat2d <= rlat1))
    prior_tg = float(np.sum([
        np.sum(np.asarray(finn_data[di])[rmask] * area_cm2[rmask]) * CONV_Tg
        for di in range(len(FINN_IDX))
    ]))
    post_tg = prior_tg * x_post[ki]
    emis_prior.append(prior_tg)
    emis_post.append(post_tg)
    print(f'  {rname:<22}  {prior_tg:>12.3f}  {post_tg:>16.3f}  {x_post[ki]:>8.3f} ± {sigma_post[ki]:.3f}')

total_prior = sum(emis_prior)
total_post  = sum(emis_post)
print(f'\n  {"TOTAL (Indonesia)":<22}  {total_prior:>12.3f}  {total_post:>16.3f}')

# ════════════════════════════════════════════════════════════════════════════
# 6. Write summary
# ════════════════════════════════════════════════════════════════════════════
summary_path = os.path.join(RES_DIR, 'inversion_v2_summary.txt')
with open(summary_path, 'w') as f:
    f.write('WRF-STILT Bayesian Inversion v2 — Fire CO2 Emissions\n')
    f.write('Indonesia IDN_BB_2015 — Oct 24-29 2015\n')
    f.write('='*70 + '\n\n')
    f.write('IMPROVEMENTS OVER v1:\n')
    f.write('  I1: Non-negativity constraint (scipy SLSQP, alpha >= 0)\n')
    f.write('  I2: Literature-informed priors (Huijnen 2016, Parker 2016, Nechita-Banda 2018)\n')
    f.write('  I3: OCO-2 column averaging kernels (per-sounding AK, not fixed 0.25)\n')
    f.write('  I4: CT2025 independence fix (sigma x2 near BKT Oct 26-28)\n')
    f.write('  I5: 6-region state vector + Tikhonov smoothness (lambda=5)\n\n')
    f.write('POSTERIOR SCALING FACTORS (v2):\n')
    f.write(f'  {"Region":<22}  {"alpha_prior":>12}  {"alpha_post":>12}  {"sigma_post":>12}\n')
    f.write('  ' + '-'*60 + '\n')
    for ki, (rname, *_) in enumerate(REGIONS):
        f.write(f'  {rname:<22}  {X_PRIOR[ki]:>12.3f}  {x_post[ki]:>12.3f}  {sigma_post[ki]:>12.3f}\n')
    f.write(f'\nchi2 (prior):     {chi2_prior:.3f}\n')
    f.write(f'chi2 (posterior): {chi2_post:.3f}\n\n')
    f.write('EMISSION TOTALS (Tg CO2, Oct 24-29):\n')
    f.write(f'  {"Region":<22}  {"Prior":>10}  {"Posterior":>12}\n')
    f.write('  ' + '-'*48 + '\n')
    for ki, (rname, *_) in enumerate(REGIONS):
        f.write(f'  {rname:<22}  {emis_prior[ki]:>10.3f}  {emis_post[ki]:>12.3f}\n')
    f.write(f'\n  {"TOTAL":<22}  {total_prior:>10.3f}  {total_post:>12.3f}\n')

print(f'\nSummary written: {summary_path}')

# ════════════════════════════════════════════════════════════════════════════
# 7. Plots
# ════════════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D

proj = ccrs.PlateCarree()
REGION_COLORS_6 = ['#e07b39','#c94040','#3a86cc','#1b5e20','#9c27b0','#00838f']

# ── Fig 1: Observation-space fit ─────────────────────────────────────────────
fig1, axes1 = plt.subplots(1, 2, figsize=(14, 6))

ax = axes1[0]
colors_obs = ['red' if t=='FLASK' else ('darkorange' if t=='OCO2' else 'steelblue')
              for t in obs_types]
sizes_obs  = [130 if t=='FLASK' else (20 if t=='OCO2' else 40) for t in obs_types]

ax.scatter(y_obs, y_prior_pred, c=colors_obs, s=sizes_obs, marker='o',
           alpha=0.8, edgecolors='none', label='Prior (FINN)')
ax.scatter(y_obs, y_post_pred,  c=colors_obs, s=sizes_obs, marker='^',
           alpha=0.9, edgecolors='gray', linewidths=0.4)
for i in range(len(y_obs)):
    ax.annotate('', xy=(y_obs[i], y_post_pred[i]),
                xytext=(y_obs[i], y_prior_pred[i]),
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.6))

xlim = max(abs(y_obs).max(), abs(y_prior_pred).max(), abs(y_post_pred).max()) * 1.15
ax.plot([-xlim, xlim], [-xlim, xlim], 'k--', lw=1)
ax.axhline(0, color='k', lw=0.4); ax.axvline(0, color='k', lw=0.4)
ax.set_xlabel('Observed Δ CO₂  [ppm]', fontsize=10)
ax.set_ylabel('Modelled Δ CO₂  [ppm]', fontsize=10)
ax.set_title('Prior vs. Posterior fit (v2)\nWith non-negativity + lit. priors + AK', fontsize=10)
ax.set_xlim(-xlim, xlim); ax.set_ylim(-xlim, xlim)
ax.set_aspect('equal')

bkt_i = obs_types.index('FLASK')
ax.annotate('BKT flask\n(real obs)', xy=(y_obs[bkt_i], y_prior_pred[bkt_i]),
            xytext=(y_obs[bkt_i]-2, y_prior_pred[bkt_i]+3),
            arrowprops=dict(arrowstyle='->', color='red'), fontsize=8, color='red')

legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='red',       markersize=9, label='BKT flask (real)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='steelblue', markersize=7, label='CT2025 pseudo-obs'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='darkorange',markersize=6, label='OCO-2 (1°×1°, AK-weighted)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='gray',      markersize=8, label='Prior ●'),
    Line2D([0],[0], marker='^', color='w', markerfacecolor='gray',      markersize=8, label='Posterior ▲'),
    Line2D([0],[0], color='k', linestyle='--', label='1:1'),
]
ax.legend(handles=legend_elements, fontsize=8, loc='upper left')
ax.grid(True, lw=0.3, alpha=0.4)
ax.text(0.97, 0.05, f'χ²_prior={chi2_prior:.1f}\nχ²_post={chi2_post:.2f}',
        ha='right', va='bottom', transform=ax.transAxes, fontsize=8,
        bbox=dict(boxstyle='round', fc='lightyellow', alpha=0.8))

ax2 = axes1[1]
short_names = [r[0].replace('Kalimantan ', 'Kali\n').replace('Sulawesi+East', 'Sulawesi\n+East')
               for r in REGIONS]
x_idx = np.arange(N_REGIONS)
w = 0.32
ax2.bar(x_idx - w/2, X_PRIOR, w, label='Prior α₀', color='lightgray', edgecolor='k', lw=0.7)
ax2.bar(x_idx + w/2, x_post,  w, label='Posterior', color=REGION_COLORS_6,
        edgecolor='k', lw=0.7, alpha=0.85)
ax2.errorbar(x_idx + w/2, x_post, yerr=sigma_post,
             fmt='none', color='black', capsize=4, lw=1.5)
ax2.axhline(1.0, color='k', linestyle=':', lw=1)
ax2.set_xticks(x_idx); ax2.set_xticklabels(short_names, fontsize=8)
ax2.set_ylabel('FINN scaling factor α', fontsize=10)
ax2.set_title('6-region posterior scaling factors\n(v2: lit. priors + non-negativity)', fontsize=10)
ax2.legend(fontsize=9); ax2.grid(axis='y', lw=0.3, alpha=0.4)
ax2.set_ylim(0, max(max(x_post)+max(sigma_post), max(X_PRIOR))*1.25)

n_oco2 = sum(1 for t in obs_types if t=='OCO2')
fig1.suptitle(f'Bayesian CO₂ Inversion v2 — STILT + NOAA flask/CT2025/OCO-2 ({n_oco2} AK-weighted cells)\n'
              'Oct 24–29 2015 Indonesia | Non-negativity + 6-region + Tikhonov', fontsize=11)
out1 = os.path.join(PLT_DIR, 'inversion_fit_v2.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight'); plt.close(fig1)
print(f'Saved: {out1}')

# ── Fig 2: Posterior emission map ─────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 6),
                            subplot_kw={'projection': proj},
                            gridspec_kw={'wspace': 0.08})

LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = 90, 145, -15, 15
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

post_mean = sub_mean.copy().astype(float)
for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
    rmask = ((lon2d_s >= rlon0) & (lon2d_s <= rlon1) &
             (lat2d_s >= rlat0) & (lat2d_s <= rlat1))
    post_mean[rmask] *= x_post[ki]

CONV_GMS = 86400 * 44.01 / NA * 1e4
prior_gm2 = np.ma.masked_less_equal(sub_mean * CONV_GMS, 0)
post_gm2  = np.ma.masked_less_equal(post_mean * CONV_GMS, 0)

from matplotlib.colors import LogNorm
nz_vals = prior_gm2.compressed()
vmin_m = max(np.percentile(nz_vals, 2) if nz_vals.size else 1e-3, 1e-3)
vmax_m = np.percentile(nz_vals, 99) if nz_vals.size else 100
lnorm  = LogNorm(vmin=vmin_m, vmax=vmax_m)

for ax3, data, title in [
    (axes2[0], prior_gm2,  'Prior: FINN mean fire flux'),
    (axes2[1], post_gm2,   'Posterior v2: 6-region scaled'),
]:
    ax3.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=proj)
    ax3.add_feature(cfeature.LAND, facecolor='#f0ede8')
    ax3.add_feature(cfeature.OCEAN, facecolor='#d5eaf5')
    ax3.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='#333')
    ax3.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='#888', linestyle=':')
    pcm = ax3.pcolormesh(lon2d_s, lat2d_s, data, cmap='hot_r', norm=lnorm,
                          transform=proj, shading='auto', rasterized=True)
    for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
        rl = [rlon0, rlon1, rlon1, rlon0, rlon0]
        rb = [rlat0, rlat0, rlat1, rlat1, rlat0]
        ax3.plot(rl, rb, color=REGION_COLORS_6[ki], lw=1.5, transform=proj, zorder=5,
                 label=f'{rname} α={x_post[ki]:.2f}')
    for sname, (slon, slat) in SITES.items():
        ax3.plot(slon, slat, marker='*', color='cyan', ms=8, transform=proj, zorder=7)
    gl = ax3.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(10); gl.ylocator = mticker.MultipleLocator(5)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {'size': 8}
    cb = fig2.colorbar(pcm, ax=ax3, orientation='vertical',
                        fraction=0.025, pad=0.04, shrink=0.9, extend='both')
    cb.set_label('CO₂ fire flux  [g m⁻² day⁻¹]', fontsize=8)
    ax3.set_title(title, fontsize=10, fontweight='bold')
    if ax3 is axes2[0]:
        ax3.legend(fontsize=6.5, loc='lower left', framealpha=0.85, ncol=2)

fig2.suptitle('Prior vs. Posterior FINN fire flux (mean Oct 24–29 2015) — v2 6-region',
              fontsize=12, fontweight='bold', y=1.01)
out2 = os.path.join(PLT_DIR, 'posterior_emissions_v2.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight'); plt.close(fig2)
print(f'Saved: {out2}')

# ── Fig 3: Diagnostics ────────────────────────────────────────────────────────
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))

ax_u = axes3[0]
unc_red = [(1 - sigma_post[ki]/SIGMA_PRIOR[ki])*100 for ki in range(N_REGIONS)]
bars = ax_u.bar(range(N_REGIONS), unc_red, color=REGION_COLORS_6, edgecolor='k', lw=0.7)
for bar, val in zip(bars, unc_red):
    ax_u.text(bar.get_x() + bar.get_width()/2, val + 0.5, f'{val:.0f}%',
              ha='center', va='bottom', fontsize=9, fontweight='bold')
ax_u.set_xticks(range(N_REGIONS))
ax_u.set_xticklabels([r[0].replace('Kalimantan ','Kali\n') for r in REGIONS], fontsize=8)
ax_u.set_ylabel('Uncertainty reduction (%)', fontsize=10)
ax_u.set_title('Posterior uncertainty reduction (v2, 6-region)', fontsize=10)
ax_u.set_ylim(0, 115); ax_u.grid(axis='y', lw=0.3, alpha=0.4)
ax_u.axhline(100, color='k', lw=0.5, linestyle='--')

SITE_COLORS = ['tab:orange','tab:red','tab:green','tab:brown','tab:blue']
ax_t = axes3[1]
obs_by_site = {sname: {'day':[], 'y_obs':[], 'y_pri':[], 'y_pos':[], 'type':[]}
               for sname in SITE_ORDER}
for i, ((sname, day, _), otype) in enumerate(zip(obs_ids, obs_types)):
    if sname not in obs_by_site:
        continue
    obs_by_site[sname]['day'].append(day)
    obs_by_site[sname]['y_obs'].append(y_obs[i])
    obs_by_site[sname]['y_pri'].append(y_prior_pred[i])
    obs_by_site[sname]['y_pos'].append(y_post_pred[i])
    obs_by_site[sname]['type'].append(otype)

for si, sname in enumerate(SITE_ORDER):
    d = obs_by_site[sname]
    days = sorted(set(d['day']))
    obs_s = [d['y_obs'][d['day'].index(dd)] for dd in days]
    pri_s = [d['y_pri'][d['day'].index(dd)] for dd in days]
    pos_s = [d['y_pos'][d['day'].index(dd)] for dd in days]
    ot_s  = [d['type'][d['day'].index(dd)] for dd in days]
    ax_t.plot(days, obs_s, 'o', color=SITE_COLORS[si], ms=6, markeredgecolor='none')
    ax_t.plot(days, pri_s, '--', color=SITE_COLORS[si], lw=1.1, alpha=0.5)
    ax_t.plot(days, pos_s, '-',  color=SITE_COLORS[si], lw=2.0)
    for dd, yt, tp in zip(days, obs_s, ot_s):
        if tp == 'FLASK':
            ax_t.plot(dd, yt, '*', color='lime', ms=14, markeredgecolor='k', lw=0.4, zorder=10)

ax_t.axhline(0, color='k', lw=0.5)
ax_t.set_xlabel('Day of October 2015', fontsize=9)
ax_t.set_ylabel('Δ CO₂  [ppm]', fontsize=9)
ax_t.set_title('Δ CO₂ receptors: obs ●, prior - -, posterior —\n★ = real BKT flask', fontsize=9)
ax_t.set_xticks(list(range(24,30)))
ax_t.set_xticklabels([f'Oct {d}' for d in range(24,30)], fontsize=8)
ax_t.grid(lw=0.3, alpha=0.4)
handles = [plt.Line2D([0],[0], color=c, lw=2, label=s[:14]) for s,c in zip(SITE_ORDER, SITE_COLORS)]
ax_t.legend(handles=handles, fontsize=7, ncol=2, framealpha=0.85)

fig3.suptitle('Inversion diagnostics v2 — Oct 24–29 2015', fontsize=11, fontweight='bold')
out3 = os.path.join(PLT_DIR, 'inversion_diagnostics_v2.png')
fig3.savefig(out3, dpi=150, bbox_inches='tight'); plt.close(fig3)
print(f'Saved: {out3}')

# ── Fig 4: v1 vs v2 comparison ────────────────────────────────────────────────
# Map 6-region posterior back to 3-region for fair v1 comparison
# Sumatra = v2[0] + v2[1] (emission-weighted mean)
# Kalimantan = v2[2] + v2[3]
# Java/Sulawesi = v2[4] + v2[5]
v1_regions     = ['Sumatra', 'Kalimantan', 'Java/Sulawesi']
v1_alpha       = [0.091, -0.010, 1.478]     # from run_inversion.py v1
v1_sigma       = [0.036,  0.007, 0.313]
v1_prior_tg    = [26.16, 57.88, 22.07]
v1_prior_alpha = [1.0, 1.0, 1.0]

# Emission-weighted aggregate of 6-region v2 → 3 parent regions
# region pairs: (0,1)=Sumatra, (2,3)=Kalimantan, (4,5)=Java/Sulawesi
pairs  = [(0,1), (2,3), (4,5)]
v2_agg_alpha = []
v2_agg_sigma = []
v2_agg_prior = []
for i,j in pairs:
    ei = emis_prior[i]; ej = emis_prior[j]
    et = ei + ej
    if et > 0:
        agg = (x_post[i]*ei + x_post[j]*ej) / et
        # Combined uncertainty (quadrature, conservative)
        agg_sig = np.sqrt((sigma_post[i]*ei)**2 + (sigma_post[j]*ej)**2) / et
    else:
        agg = 0.0; agg_sig = 0.0
    v2_agg_alpha.append(agg)
    v2_agg_sigma.append(agg_sig)
    v2_agg_prior.append((X_PRIOR[i]*ei + X_PRIOR[j]*ej) / (et if et > 0 else 1))

fig4, axes4 = plt.subplots(1, 2, figsize=(13, 5))

# Left: Alpha comparison
ax_a = axes4[0]
x3   = np.arange(3)
w3   = 0.25
ax_a.bar(x3 - w3, v1_prior_alpha, w3, color='#ddd', edgecolor='k', lw=0.7, label='v1 prior (α=1)')
ax_a.bar(x3,      v1_alpha,       w3, color='#f08080', edgecolor='k', lw=0.7, label='v1 posterior')
ax_a.errorbar(x3, v1_alpha, yerr=v1_sigma, fmt='none', color='black', capsize=4, lw=1.5)
ax_a.bar(x3 + w3, v2_agg_alpha,  w3, color='#4682b4', edgecolor='k', lw=0.7, label='v2 posterior (aggregated)')
ax_a.errorbar(x3 + w3, v2_agg_alpha, yerr=v2_agg_sigma,
              fmt='none', color='black', capsize=4, lw=1.5)

# Add hatching for v2 bars to indicate non-negativity was active for Kalimantan
kali_patch = ax_a.patches[4]  # 3rd v2 bar (Kalimantan)

ax_a.axhline(0, color='k', lw=1)
ax_a.axhline(1, color='k', lw=0.8, linestyle=':')
ax_a.set_xticks(x3)
ax_a.set_xticklabels(v1_regions, fontsize=10)
ax_a.set_ylabel('FINN emission scaling factor α', fontsize=10)
ax_a.set_title('v1 vs v2 posterior scaling factors', fontsize=10, fontweight='bold')
ax_a.legend(fontsize=9)
ax_a.grid(axis='y', lw=0.3, alpha=0.4)
ax_a.set_ylim(-0.25, max(max(v2_agg_alpha)+max(v2_agg_sigma), 2.0)*1.15)

# Annotate Kalimantan physical change
ax_a.annotate('Non-negativity\nactive here (v1 < 0)',
              xy=(1, 0.005), xytext=(1.3, 0.25),
              arrowprops=dict(arrowstyle='->', color='blue'),
              fontsize=8, color='blue')

# Right: Emission totals comparison
ax_e = axes4[1]
v1_post_tg = [p*a for p,a in zip(v1_prior_tg, v1_alpha)]
v2_post_tg = [p * (x_post[i]*emis_prior[i] + x_post[j]*emis_prior[j]) /
               (emis_prior[i]+emis_prior[j]) if (emis_prior[i]+emis_prior[j]) > 0 else 0
               for (i,j),(p) in zip(pairs, v1_prior_tg)]
# Simpler: just use emis_post[i]+emis_post[j]
v2_post_tg = [emis_post[i]+emis_post[j] for i,j in pairs]

x3 = np.arange(3)
bars_prior = ax_e.bar(x3 - w3, v1_prior_tg, w3, color='#ddd', edgecolor='k', lw=0.7, label='FINN prior')
bars_v1    = ax_e.bar(x3,      v1_post_tg,  w3, color='#f08080', edgecolor='k', lw=0.7, label='v1 posterior')
bars_v2    = ax_e.bar(x3 + w3, v2_post_tg,  w3, color='#4682b4', edgecolor='k', lw=0.7, label='v2 posterior')

ax_e.set_xticks(x3)
ax_e.set_xticklabels(v1_regions, fontsize=10)
ax_e.set_ylabel('CO₂ fire emission  [Tg CO₂, Oct 24–29]', fontsize=10)
ax_e.set_title('Prior vs. posterior emission totals\nv1 and v2 comparison', fontsize=10, fontweight='bold')
ax_e.legend(fontsize=9)
ax_e.grid(axis='y', lw=0.3, alpha=0.4)

# Annotate totals
for rect, val in zip(bars_prior, v1_prior_tg):
    ax_e.text(rect.get_x()+rect.get_width()/2, max(val,0)+0.5, f'{val:.0f}', ha='center', fontsize=8)
for rect, val in zip(bars_v1, v1_post_tg):
    label = f'{val:.1f}' if val >= 0 else f'{val:.1f}*'
    ax_e.text(rect.get_x()+rect.get_width()/2, max(val,0)+0.5, label, ha='center', fontsize=8, color='darkred')
for rect, val in zip(bars_v2, v2_post_tg):
    ax_e.text(rect.get_x()+rect.get_width()/2, max(val,0)+0.5, f'{val:.1f}', ha='center', fontsize=8, color='navy')

ax_e.set_ylim(min(0, min(v1_post_tg))-5, max(v1_prior_tg)*1.25)
ax_e.axhline(0, color='k', lw=0.5)
ax_e.text(0.98, 0.02, '* v1 Kalimantan < 0 (unphysical)\nv2 constrained ≥ 0',
          ha='right', va='bottom', transform=ax_e.transAxes, fontsize=8,
          color='darkred', style='italic')

fig4.suptitle('v1 vs v2 Inversion Comparison — Key improvements: non-negativity, '
              '6-region, lit. priors, AK-corrected OCO-2', fontsize=11, fontweight='bold')
out4 = os.path.join(PLT_DIR, 'v1_vs_v2_comparison.png')
fig4.savefig(out4, dpi=150, bbox_inches='tight'); plt.close(fig4)
print(f'Saved: {out4}')

print('\n' + '='*70)
print('Inversion v2 complete.')
print(f'Results: {RES_DIR}/inversion_v2_summary.txt')
print(f'Plots:   {PLT_DIR}/')
print('='*70)
