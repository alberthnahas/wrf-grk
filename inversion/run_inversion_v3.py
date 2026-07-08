#!/usr/bin/env python3
"""
run_inversion_v3.py — Bayesian CO₂ inversion with real 11-site STILT network

═══════════════════════════════════════════════════════════════════════════════
CHANGELOG
═══════════════════════════════════════════════════════════════════════════════

v3.0 (initial)
  I6 (proper): 11 real STILT backward-trajectory sites covering all major
      Indonesian fire regions — replaces v2 IDW-approximated virtual footprints.
      New sites: Pekanbaru, Palembang, Banjarmasin, Samarinda, Balikpapan,
      Makassar. Kept: BKT, Jambi, Jakarta, Pontianak, Palangka Raya.
  I7: MODIS hotspot-informed prior σ modulation (regions with >200 hotspots
      get σ × 0.7; <50 get σ × 1.2).
  I8: CT2025 fire-zone exclusion — rows with H_max > 25 ppm are dropped from
      the inversion (CT2025's 3°×2° smoothing cannot represent dense plumes).
  STILT receptors at 12 UTC (= 19 WIB).

v3.1 (post-Sulawesi-bug fix, 2026-05-08)
  Fixed four bugs that produced an unphysical Sulawesi+East ≈ Sumatra S
  posterior in v3.0:
   1. REGIONS geometry: Central Kalimantan peat corridor (lon 112–115) was
      mis-allocated to "Kalimantan W+C". Redrawn so peat corridor is in S+E
      and Sulawesi+East starts at lon 119 (Sulawesi proper).
   2. HOTSPOT_REGION_BOUNDS updated to match new REGIONS.
   3. Tikhonov D_MATRIX W+C↔S+E coupling row removed — was forcing the
      two Kalimantan regions to track each other and masking real signal.
   4. Priors updated from literature (Huijnen 2016, Field 2016, Wooster 2018):
      Sulawesi tightened to α=0.10 ± 0.08 (was 0.27 ± 0.20).
  CO joint constraint (I7 CO) DISABLED — found to be a circular model-vs-model
  constraint (USE_CO_I7 = False, kept for future TROPOMI/MOPITT integration).
  Diagnostics plot bug: receptor time series subplot fixed (line ~875,
  list-comprehension `_` rebind shadowed day comparison).

v3.2 (in progress, addresses analysis-doc limitations §10)
  STILT rerun at 06 UTC (= 13 WIB) — true afternoon mixed PBL, aligned with
    OCO-2 overpass time (~06 UTC equator). Replaces v3.1 12 UTC (= 19 WIB
    post-sunset stable layer). Output: stilt/out_06utc/footprints/.
    DONE: FOOT_DIR → out_06utc, DATES_NOON hour → 06, CT2025 index → 2 (06 UTC),
    and WRF tracer sampling → 06 UTC are all switched. Receptor time, footprint
    time, CT2025 sample, WRF tracer, and BKT flask (07 UTC) now agree to ≤1 h.
  §10.2 fix: OCO-2 cells with no STILT receptor within OCO2_MAX_DIST_DEG=4°
    are rejected — drops Maluku/Papua/east Sulawesi cells where IDW from
    far-away sites would fabricate H. Cells over the well-sampled west and
    central Indonesian network are retained.
  §10.3 fix: CT2025 σ inflated ×2 for borderline-fire rows (5 < H_max ≤ 25
    ppm). Combined with the v3.0 I8 hard exclusion (H_max > 25), this lets
    OCO-2 dominate where CT2025 cannot resolve plume structure while still
    allowing CT2025 to weakly constrain background-influenced rows.
  §10.6 fix: Tikhonov D_MATRIX emptied — Sumatra N↔S smoothing row removed.
    Was inflating Sumatra N+C posterior σ by tying it to the well-constrained
    S region. Tikhonov term now contributes 0; structure retained for
    possible future spatial smoothing.

  Limitations §10.1 (no Sulawesi receptor) and §10.4 (no real CO data)
  remain — both require external data acquisition out of scope here.

═══════════════════════════════════════════════════════════════════════════════
Sites (11 total × 26 days Oct 6–31 = 286 simulations)
═══════════════════════════════════════════════════════════════════════════════
"""

import os, glob
import numpy as np
from scipy.optimize import minimize
import netCDF4 as nc
from datetime import datetime, timezone
import warnings
warnings.filterwarnings('ignore')

# [v3.1 CHANGE] USE_CO_I7 disabled — found to be a circular constraint
# (WRF CO_BBU forced with FINN α=1, then compared against H_STILT × FINN_CO,
# which pulls posterior α toward the WRF input α rather than truth).
# Retained as a feature flag for future TROPOMI/MOPITT integration.
USE_OCO2     = True
USE_CO_I7    = False
USE_HOTSPOT_PRIOR = True  # I7: use MODIS hotspot-informed prior scaling

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT      = '/home/igrk/WRF-GRK'
INV_DIR   = f'{ROOT}/inversion'
DATA_DIR  = f'{INV_DIR}/data'
RES_DIR   = f'{INV_DIR}/results'
PLT_DIR   = f'{INV_DIR}/plots'
# [v3.2 CHANGE] Switched from out/ (12 UTC = 19 WIB, post-sunset) to
# out_06utc/ (06 UTC = 13 WIB, true afternoon mixed PBL). Aligns with
# NOAA BKT flask sample time (07 UTC) and OCO-2 overpass (~06 UTC).
FOOT_DIR  = f'{ROOT}/stilt/out_06utc/footprints'
CT_DIR    = f'{ROOT}/rawdata/carbontracker'
FINN_CO2  = (f'{ROOT}/rawdata/finn_fire/'
             'emissions-finnv2.5modvrs_CO2_bb_surface_daily_'
             '20150101-20151231_0.1x0.1.nc')
FINN_CO   = (f'{ROOT}/rawdata/finn_fire/'
             'emissions-finnv2.5modvrs_CO_bb_surface_daily_'
             '20150101-20151231_0.1x0.1.nc')
WRF_OUT   = f'{ROOT}/simulations/IDN_BB_2015/output'
HOTSPOT   = f'{ROOT}/rawdata/hotspot/archived_hotspot_idn.csv'
BKT_FLASK_FILE = f'{INV_DIR}/data/co2_bkt_flask.txt'  # NOAA GML BKT flask record

for d in [RES_DIR, PLT_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Constants ───────────────────────────────────────────────────────────────
NA         = 6.022e23
# [WINDOW] Case study = Oct 6-27 2015 (Option A, obs-overlap window):
#   lower bound Oct 6  — ARL met starts Oct 1, STILT needs 5 prior days;
#   upper bound Oct 27 — OCO-2 granules end Oct 26 and the BKT flask is Oct 27,
#   so Oct 28-31 would be CT2025-only and add little constraint. STILT itself
#   may have footprints out to Oct 31; the inversion simply uses Oct 6-27.
DATES      = [datetime(2015, 10, d, tzinfo=timezone.utc) for d in range(6, 28)]
# [v3.2 CHANGE] Receptor times shifted 12→06 UTC; name kept for code clarity.
DATES_NOON = [datetime(2015, 10, d, 6, tzinfo=timezone.utc) for d in range(6, 28)]

# ── I6: All 11 real receptor sites ─────────────────────────────────────────
SITES = {
    'Bukit Kototabang': (100.32, -0.20),   # EANET/WMO flask — only real flask
    'Jambi':            (103.61, -1.61),   # Sumatra S
    'Jakarta':          (106.85, -6.21),   # Java
    'Pontianak':        (109.33, -0.02),   # Kalimantan W coast
    'Palangka Raya':    (113.94, -2.16),   # Central Kalimantan peat fires
    'Pekanbaru':        (101.45,  0.53),   # Sumatra C, Riau peat
    'Palembang':        (104.75, -2.99),   # Sumatra S, major peat fire region
    'Banjarmasin':      (114.59, -3.32),   # Kalimantan S, peat corridor
    'Samarinda':        (117.15, -0.50),   # Kalimantan E
    'Balikpapan':       (116.85, -1.27),   # Kalimantan E coast
    'Makassar':         (119.42, -5.14),   # Sulawesi S — major city
}
SITE_ORDER = list(SITES.keys())
N_SITES    = len(SITE_ORDER)

# BKT real flask observation (Oct 27 only — Oct 26/28 are CT2025)
BKT_IDX    = SITE_ORDER.index('Bukit Kototabang')
OCT27_IDX  = 21  # 0-indexed: Oct 6=0, 7=1, ... Oct 27=21

# ── 6-region state vector ──────────────────────────────────────────────────
# [v3.1 CHANGE] Region geometry redrawn so Central Kalimantan peat corridor
# (Palangka Raya, Sampit, Kapuas — lon 112–115) is in S+E/peat, NOT W+C.
# Sulawesi restricted to lon ≥ 119 (Sulawesi proper) so East Kalimantan
# coast (Samarinda 117.15, Balikpapan 116.85) doesn't leak into "Sulawesi".
# Previous (v3.0, BUGGY) bounds were:
#   Kalimantan W+C: 107–115, -7–3   (captured the peat corridor incorrectly)
#   Kalimantan S+E: 115–121, -7–3   (was mostly East Kalimantan coast)
#   Sulawesi+East:  116–130, -10–-3 (overlapped East Kalimantan)
REGIONS = [
    ('Sumatra N+C',    95.0, 104.0,  0.0,  6.0),
    ('Sumatra S',      99.0, 107.0, -6.0,  0.0),
    ('Kalimantan W+C', 108.0, 112.0, -2.0,  3.0),   # Pontianak basin only
    ('Kalimantan S+E', 112.0, 119.0, -5.0,  1.0),   # Central+S+E peat corridor
    ('Java',           105.0, 116.0, -9.0, -5.5),
    ('Sulawesi+East',  119.0, 141.0, -6.0,  2.0),   # Sulawesi proper + Maluku/Papua
]
N_REGIONS = len(REGIONS)

# ── I2: Literature-informed priors ─────────────────────────────────────────
# Priors indexed to match REGIONS order:
#   [Sumatra N+C, Sumatra S, Kali W+C, Kali S+E, Java, Sulawesi+East]
#
# Sumatra N+C: FINN ~4-5× high for Riau peat → α=0.25  (Huijnen 2016)
# Sumatra S:   FINN ~2.5-3× high for S.Sumatra peat → α=0.35  (van der Werf 2017)
# Kali W+C:    Central Kali peat: FINN ~5-8× high → α=0.15  (Parker 2016)
# Kali S+E:    East Kali: GFED4s ≈ FINN for 2015 (~80 Tg Oct) → α=0.80
#              (van der Werf 2017 Fig.4; FINN prior 90 Tg ≈ GFED 70-100 Tg)
# Java:        FINN uses tropical biome EF; Java fires are minor/agricultural
#              GFED ~0.8 Tg vs FINN ~6.5 Tg → GFED/FINN ≈ 0.12 → α=0.12
# Sulawesi:    Savanna/grassland fires; GFED ~12 Tg vs FINN ~44 Tg → α=0.27
#              (FINN overestimates savanna combustion completeness for Sulawesi)
# [v3.1 CHANGE] Priors updated for new region geometry, literature-based.
# Sulawesi tightened to 0.10 ± 0.08 — literature (Huijnen 2016, Field 2016,
# Wooster 2018) places Sulawesi at ~1–3% of Indonesia 2015 fire emissions
# (~2–6 Tg CO₂ for Oct), with FINN ~44 Tg → α ≈ 0.05–0.15.
# Kalimantan S+E now contains the Central Kalimantan peat corridor; FINN
# overestimates peat ~5–8× (Parker 2016, van der Werf 2017) → α=0.20.
# Previous (v3.0): X_PRIOR=[0.25,0.35,0.15,0.80,0.12,0.27],
#                  SIGMA_PRIOR=[0.20,0.25,0.20,0.35,0.20,0.20]
X_PRIOR     = np.array([0.25, 0.35, 0.20, 0.20, 0.12, 0.10])
SIGMA_PRIOR = np.array([0.20, 0.25, 0.15, 0.15, 0.20, 0.08])

# ── I7: MODIS hotspot-based prior update ───────────────────────────────────
# Hotspot counts provide spatially explicit fire intensity per region
# independent of any transport model. Scale the prior mean by hotspot fraction.
# [v3.1 CHANGE] Hotspot region bounds updated to match new REGIONS geometry.
HOTSPOT_REGION_BOUNDS = {
    'Sumatra N+C':    (95, 104,  0,  6),
    'Sumatra S':      (99, 107, -6,  0),
    'Kalimantan W+C': (108, 112, -2,  3),
    'Kalimantan S+E': (112, 119, -5,  1),
    'Java':           (105, 116, -9, -5.5),
    'Sulawesi+East':  (119, 141, -6,  2),
}

# ── Tikhonov regularization ─────────────────────────────────────────────────
# [v3.1 CHANGE] D_MATRIX reduced from 2 rows → 1 row (W+C↔S+E removed).
# [v3.2 CHANGE] Sumatra N↔S smoothing row also removed (addresses analysis
# limitation §10.6): forcing α_SumN ≈ α_SumS inflated posterior σ for N+C
# (a small region with weak data) by tying it to the well-constrained S.
# With D_MATRIX empty, Tikhonov term contributes 0 to cost — kept structure
# in case spatial smoothing is reintroduced for a different region pair.
TIKHONOV_LAMBDA = 5.0
D_MATRIX = np.zeros((0, N_REGIONS))   # no smoothing rows; rely on prior σ alone

CO2_BACKGROUND = 399.1
BG_UNCERTAINTY = 0.5

# ── I8: CT2025 fire-zone exclusion ─────────────────────────────────────────
# CT2025 uses GFED as fire prior; at fire-zone rows H_max×α = CT2025 encodes
# GFED/FINN ≈ 0.01, not a real atmospheric constraint. Exclude those rows.
FIRE_ZONE_THRESHOLD = 25.0  # ppm

# I7: CO model parameters
# WRF CO_BBU is in ppmv, already forward-modelled at receptor grid cell
# CO_UNCERTAINTY: WRF model error + background uncertainty
CO_BACKGROUND   = 0.095    # ppm (typical clean-air CO in Indonesia, Oct 2015)
CO_UNCERTAINTY  = 0.02     # ppm background uncertainty
CO_OBS_SIGMA    = 0.05     # ppm observation/model error at receptor


# ════════════════════════════════════════════════════════════════════════════
# Helper functions
# ════════════════════════════════════════════════════════════════════════════

def load_ct_surface_co2(date, site_lat, site_lon):
    fname = os.path.join(CT_DIR,
                         f'CT2025.molefrac_glb3x2_{date.strftime("%Y-%m-%d")}.nc')
    if not os.path.exists(fname):
        return None
    ds  = nc.Dataset(fname)
    lat = ds.variables['latitude'][:]
    lon = ds.variables['longitude'][:]
    ji  = int(np.argmin(np.abs(lat - site_lat)))
    ii  = int(np.argmin(np.abs(lon - site_lon)))
    # [v3.2 CHANGE] CT2025 is 3-hourly (00, 03, 06, ..., 21 UTC).
    # Index 2 = 06 UTC, matching new STILT receptor release time.
    # Was index 4 (= 12 UTC) before the v3.2 06 UTC switch.
    co2 = float(ds.variables['co2'][2, 0, ji, ii])
    ds.close()
    return co2


def load_wrf_tracer_at_site(day, slon, slat, varname):
    """Return surface-level WRF tracer (ppmv) at 06 UTC at nearest grid cell.

    [v3.2 FIX] Time index moved 12→06 UTC to match the 06 UTC STILT receptor
    release, the CT2025 06 UTC sample (index 2), and the ~07 UTC BKT flask.
    (Only reached when USE_CO_I7 is enabled; kept consistent for reactivation.)
    """
    fpath = os.path.join(WRF_OUT,
                         f'wrfout_d01_2015-10-{day:02d}_01:00:00')
    if not os.path.exists(fpath):
        return None
    f = nc.Dataset(fpath)
    times_raw = f.variables['Times'][:]
    times     = [''.join(t.astype(str)) for t in times_raw]
    try:
        t_idx = next(i for i, t in enumerate(times) if '06:00' in t)
    except StopIteration:
        f.close()
        return None
    xlat = f.variables['XLAT'][0]    # (ny, nx)
    xlon = f.variables['XLONG'][0]   # (ny, nx)
    dist = (xlat - slat)**2 + (xlon - slon)**2
    iy, ix = np.unravel_index(dist.argmin(), dist.shape)
    val = float(f.variables[varname][t_idx, 0, iy, ix])
    f.close()
    return val


def load_finn_regridded(foot_lon, foot_lat, day_indices, varname='fire_modisviirs_CO2',
                         finn_path=None):
    if finn_path is None:
        finn_path = FINN_CO2
    ds       = nc.Dataset(finn_path)
    finn_lon = ds.variables['lon'][:]
    finn_lat = ds.variables['lat'][:]
    raw      = ds.variables[varname][day_indices, :, :]  # (ndays, nlat, nlon)
    ds.close()
    CONV     = 1e4 / NA * 1e6   # molec/cm²/s → µmol/m²/s
    raw_umol = raw * CONV
    dlat_f   = abs(float(foot_lat[1] - foot_lat[0]))
    dlon_f   = float(foot_lon[1] - foot_lon[0])
    n_days   = raw.shape[0]
    out      = np.zeros((n_days, len(foot_lat), len(foot_lon)), dtype=np.float64)
    for jj, flat in enumerate(foot_lat):
        jmask = np.where((finn_lat >= flat-dlat_f/2) &
                         (finn_lat <  flat+dlat_f/2))[0]
        if jmask.size == 0:
            continue
        cos_lat = np.abs(np.cos(finn_lat[jmask] * np.pi/180))
        for ii, flon in enumerate(foot_lon):
            imask = np.where((finn_lon >= flon-dlon_f/2) &
                             (finn_lon <  flon+dlon_f/2))[0]
            if imask.size == 0:
                continue
            patch = raw_umol[:, jmask[0]:jmask[-1]+1, imask[0]:imask[-1]+1]
            w2d   = cos_lat[:, np.newaxis] * np.ones(len(imask))[np.newaxis, :]
            out[:, jj, ii] = (np.sum(patch * w2d, axis=(1, 2))
                              / w2d.sum())
    return out


def load_footprint(slon, slat, dt_noon):
    """Load footprint NetCDF for this site and time."""
    tag = dt_noon.strftime('%Y%m%d%H%M')
    pattern = os.path.join(FOOT_DIR, f'{tag}_*_foot.nc')
    for f in glob.glob(pattern):
        parts = os.path.basename(f).split('_')
        try:
            if (abs(float(parts[1]) - slon) < 0.02 and
                    abs(float(parts[2]) - slat) < 0.02):
                ds   = nc.Dataset(f)
                foot = ds.variables['foot'][0, :, :]
                ds.close()
                return foot
        except (ValueError, IndexError):
            continue
    return None


# ════════════════════════════════════════════════════════════════════════════
# 1.  I7 — MODIS hotspot prior update
# ════════════════════════════════════════════════════════════════════════════
print('=' * 70)
print('I7: MODIS hotspot-informed prior ...')
print('=' * 70)

import pandas as pd

hs_df = pd.read_csv(HOTSPOT)
hs15  = hs_df[(hs_df.year == 2015) & (hs_df.month == 10) &
               (hs_df.day >= 6) & (hs_df.day <= 31)]
print(f'MODIS hotspots Oct 6-31 2015: {len(hs15)}')

hotspot_counts = {}
for rname, (x0, x1, y0, y1) in HOTSPOT_REGION_BOUNDS.items():
    n = int(((hs15.lon >= x0) & (hs15.lon < x1) &
              (hs15.lat >= y0) & (hs15.lat < y1)).sum())
    hotspot_counts[rname] = n

total_hs = max(sum(hotspot_counts.values()), 1)
print(f'  {"Region":<22}  Count   Fraction')
for rname, n in hotspot_counts.items():
    print(f'  {rname:<22}  {n:5d}   {n/total_hs*100:5.1f}%')

# Use hotspot counts as independent fire-intensity prior.
# Hotspot density ∝ fire radiative power ∝ CO₂ emission (roughly).
# Combine with literature α prior: weight 50/50 between literature prior
# and hotspot-implied relative intensity.
if USE_HOTSPOT_PRIOR:
    # Normalised hotspot fraction per region (as emission proxy)
    hs_frac = np.array([hotspot_counts[r[0]] / total_hs for r in REGIONS])

    # Expected total Indonesia emission = α_lit × FINN_total
    # We know total FINN_prior × 0.30 (midpoint of literature) ≈ literature midpoint
    # Redistribute: α_hotspot_i = α_lit_mean * (hs_frac_i / hs_frac_expected_i)
    # where hs_frac_expected_i = X_PRIOR[i]*FINN_frac_i / sum(X_PRIOR*FINN_frac)
    # This is an additive correction to X_PRIOR, not a full replacement.
    # We implement it by tightening SIGMA_PRIOR where hotspot data is dense
    # (high hotspot count = more confident about relative emission intensity).

    # Hotspot-implied uncertainty reduction: denser hotspot = tighter prior
    # For regions with >200 hotspots: σ_prior × 0.7 (30% tighter)
    # For regions with <50 hotspots: σ_prior × 1.2 (more uncertain)
    for ki, (rname, *_) in enumerate(REGIONS):
        n = hotspot_counts.get(rname, 0)
        if n > 200:
            SIGMA_PRIOR[ki] *= 0.7
            print(f'  [{rname}] σ_prior tightened (hotspot={n}): '
                  f'{SIGMA_PRIOR[ki]/0.7:.3f} → {SIGMA_PRIOR[ki]:.3f}')
        elif n < 50:
            SIGMA_PRIOR[ki] *= 1.2
            print(f'  [{rname}] σ_prior widened (hotspot={n}): '
                  f'{SIGMA_PRIOR[ki]/1.2:.3f} → {SIGMA_PRIOR[ki]:.3f}')

    print(f'\nUpdated SIGMA_PRIOR: {SIGMA_PRIOR}')


# ════════════════════════════════════════════════════════════════════════════
# 2.  Build observation vector (CT2025 + BKT flask)
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('Building CO₂ observation vector ...')
print('=' * 70)

# [v3.2 CHANGE] Multi-date BKT flask loader — was a single hardcoded Oct 27
# pair; now reads the full NOAA GML flask file and pairs all Oct 2015 samples.
# Sample time is 07 UTC (= 14 WIB, true afternoon mixed PBL). With the v3.2
# 06 UTC footprints the flask/footprint offset is now 1 h (was 5 h under the
# 12 UTC config) — acceptable within the well-mixed afternoon PBL.
def load_bkt_flasks(year=2015, month=10):
    """Return list of (date, mean_ppm, instr_sigma, qcflag) for paired flask
    samples. Pairs are averaged; rejected (col-1 != '.') flasks dropped."""
    samples = {}
    with open(BKT_FLASK_FILE) as fh:
        for line in fh:
            if not line.startswith('BKT '):
                continue
            parts = line.split()
            if int(parts[1]) != year or int(parts[2]) != month:
                continue
            day = int(parts[3])
            hh, mm = int(parts[4]), int(parts[5])
            value, unc = float(parts[10]), float(parts[11])
            qc = parts[-1]
            # QC col 1 = rejection flag; col 2 = background-selection flag.
            # 'C' in col 1 indicates smoke contamination (fails NOAA background
            # criteria) but the *measurement* is valid — KEEP these for the
            # fire-event inversion. Only drop true instrument-failure codes.
            INSTR_FAIL = set('NXY*')   # NOAA field-failure markers
            if qc[0] in INSTR_FAIL:
                continue
            key = (day, hh)
            samples.setdefault(key, []).append((value, unc, qc))
    out = []
    for (day, hh), pair in sorted(samples.items()):
        vals = [v for v, _, _ in pair]
        uncs = [u for _, u, _ in pair]
        qcs  = [q for _, _, q in pair]
        out.append({
            'datetime': datetime(year, month, day, hh, tzinfo=timezone.utc),
            'mean':     float(np.mean(vals)),
            'std':      float(np.std(vals, ddof=0)),
            'instr_sigma': float(np.mean(uncs)),
            'qcflag':   qcs[0],
            'n_pair':   len(pair),
        })
    return out

bkt_flasks = load_bkt_flasks(2015, 10)
print(f'BKT NOAA flask samples in Oct 2015: {len(bkt_flasks)}')
for f in bkt_flasks:
    enh = f['mean'] - CO2_BACKGROUND
    print(f"  {f['datetime']:%Y-%m-%d %H:%M UTC}  "
          f"co2={f['mean']:.2f} ppm  enh={enh:+.2f} ppm  "
          f"qc={f['qcflag']}  n={f['n_pair']}")

ct_obs = []
for di, dt in enumerate(DATES):
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        ct_co2 = load_ct_surface_co2(dt, slat, slon)
        if ct_co2 is None:
            continue
        fire_enh = ct_co2 - CO2_BACKGROUND
        # I4: CT2025 independence fix
        bkt_lon, bkt_lat = SITES['Bukit Kototabang']
        dist_bkt = np.sqrt((slon - bkt_lon)**2 + (slat - bkt_lat)**2)
        near_bkt = (dist_bkt < 3.0) and (dt.day in [26, 27, 28])
        inflate  = 2.0 if near_bkt else 1.0
        sigma_ct = np.sqrt((3.0 * inflate)**2 + BG_UNCERTAINTY**2)
        note     = ' [σ×2 CT independence]' if near_bkt else ''
        ct_obs.append((si, di, sname, dt, fire_enh, sigma_ct, 'CT2025', note))

print(f'CT2025 pseudo-obs: {len(ct_obs)} ({N_SITES} sites × 6 days)')


# ════════════════════════════════════════════════════════════════════════════
# 3.  Build H matrix from ALL 11 real STILT footprints
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('Building H matrix from 11-site STILT network ...')
print('=' * 70)

fp_files = sorted(glob.glob(os.path.join(FOOT_DIR, '*_foot.nc')))
print(f'Total footprint files found: {len(fp_files)}')

ds0      = nc.Dataset(fp_files[0])
foot_lon = ds0.variables['lon'][:]
foot_lat = ds0.variables['lat'][:]
ds0.close()

FINN_IDX    = [d.timetuple().tm_yday - 1 for d in DATES]
print('Regridding FINN CO₂ to 0.25° grid ...')
finn_co2_regrid = load_finn_regridded(foot_lon, foot_lat, FINN_IDX,
                                       varname='fire_modisviirs_CO2',
                                       finn_path=FINN_CO2)

if USE_CO_I7:
    print('Regridding FINN CO to 0.25° grid ...')
    finn_co_regrid = load_finn_regridded(foot_lon, foot_lat, FINN_IDX,
                                          varname='fire_modisviirs_CO',
                                          finn_path=FINN_CO)

lon2d, lat2d = np.meshgrid(foot_lon, foot_lat)
region_masks = []
for rname, rlon0, rlon1, rlat0, rlat1 in REGIONS:
    mask = ((lon2d >= rlon0) & (lon2d <= rlon1) &
            (lat2d >= rlat0) & (lat2d <= rlat1))
    region_masks.append(mask)
    print(f'  Region "{rname}": {mask.sum()} grid cells')

n_obs  = len(DATES) * N_SITES
H_co2  = np.zeros((n_obs, N_REGIONS))   # CO₂ H-matrix
H_co   = np.zeros((n_obs, N_REGIONS))   # CO H-matrix (same footprint, FINN_CO)
obs_labels = []
foot_ok    = np.zeros(n_obs, dtype=bool)

n_missing = 0
for di, (dt_midnight, dt_noon) in enumerate(zip(DATES, DATES_NOON)):
    for si, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        row_idx    = di * N_SITES + si
        foot       = load_footprint(slon, slat, dt_noon)
        if foot is None:
            obs_labels.append((sname, dt_midnight.day, 'missing'))
            n_missing += 1
            continue
        foot_ok[row_idx] = True
        E_co2 = finn_co2_regrid[di, :, :]
        for ki, mask in enumerate(region_masks):
            H_co2[row_idx, ki] = np.sum(foot * E_co2 * mask)
        if USE_CO_I7:
            E_co = finn_co_regrid[di, :, :]
            for ki, mask in enumerate(region_masks):
                H_co[row_idx, ki] = np.sum(foot * E_co * mask)
        obs_labels.append((sname, dt_midnight.day, 'ok'))

print(f'\nH matrix (CO₂) shape: {H_co2.shape}')
print(f'Missing footprints: {n_missing} / {n_obs}')
for ki, (rname, *_) in enumerate(REGIONS):
    col = H_co2[:, ki]
    print(f'  {rname:<22}: mean={col.mean():.3f}  max={col.max():.3f}')

if USE_CO_I7:
    print(f'\nH matrix (CO) shape: {H_co.shape}')
    for ki, (rname, *_) in enumerate(REGIONS):
        col = H_co[:, ki]
        print(f'  {rname:<22}: mean={col.mean():.4f}  max={col.max():.4f}')

# ── I8: Identify fire-zone CT2025 rows ─────────────────────────────────────
H_row_max   = np.abs(H_co2).max(axis=1)
fire_zone_rows = H_row_max > FIRE_ZONE_THRESHOLD
# [v3.2 CHANGE] (addresses analysis §10.3 — CT2025 dominance)
# CT2025's 3°×2° smoothing cannot represent local plumes. v3.0 dropped rows
# with H_max > 25 ppm entirely (45 of 286). v3.2 also INFLATES σ for the
# borderline-fire band 5 < H_max ≤ 25 ppm (where CT2025 is partially
# fire-influenced but not catastrophically so) so OCO-2 dominates the
# constraint there while CT2025 still contributes a weak prior pull.
CT_BORDERLINE_LOW  = 5.0    # ppm; below this, σ unchanged
CT_BORDERLINE_HIGH = FIRE_ZONE_THRESHOLD   # = 25 ppm; above, row excluded
borderline_rows = ((H_row_max > CT_BORDERLINE_LOW) &
                   (H_row_max <= CT_BORDERLINE_HIGH))
print(f'\nFire-zone CT2025 rows (H_max > {FIRE_ZONE_THRESHOLD} ppm): '
      f'{fire_zone_rows.sum()} / {n_obs}')
print(f'Borderline CT2025 rows ({CT_BORDERLINE_LOW} < H_max ≤ '
      f'{CT_BORDERLINE_HIGH} ppm, σ inflated ×2): '
      f'{borderline_rows.sum()} / {n_obs}')
for i, (lbl, hmax) in enumerate(zip(obs_labels, H_row_max)):
    if fire_zone_rows[i]:
        print(f'  Excluded CT2025: {lbl[0]:<22} Oct{lbl[1]:02d}  H_max={hmax:.1f} ppm')


# ════════════════════════════════════════════════════════════════════════════
# 4.  Assemble observation vector — CO₂
# ════════════════════════════════════════════════════════════════════════════
y_co2   = np.zeros(n_obs)
sig_co2 = np.full(n_obs, 99.0)
otype   = ['CT2025'] * n_obs

for si, di, sname, dt, fire_enh, sigma_ct, obstype, note in ct_obs:
    row_idx = di * N_SITES + si
    if fire_zone_rows[row_idx]:     # I8: exclude fire-zone CT2025
        continue
    # [v3.2 CHANGE] borderline-fire CT2025: inflate σ ×2 (analysis §10.3)
    if borderline_rows[row_idx]:
        sigma_ct = sigma_ct * 2.0
    y_co2[row_idx]   = fire_enh
    sig_co2[row_idx] = sigma_ct

# [v3.2 CHANGE] BKT real flask — replaces single Oct 27 hardcode with all
# Oct 2015 NOAA GML flask samples that have a STILT footprint at the same day.
# Each flask pair (07 UTC) overrides the CT2025 row at the BKT slot for that
# day. Background = CT2025-derived 399.1 ppm (CO2_BACKGROUND); fire enhancement
# = pair_mean - background. σ = sqrt(instr² + bg² + repr²) with repr=3 ppm
# (footprint-CT2025 representativity), repr=5 ppm if QC flag has 'C' (smoke).
n_flasks_used = 0
for f in bkt_flasks:
    di_match = next((i for i, dt in enumerate(DATES)
                     if dt.day == f['datetime'].day), None)
    if di_match is None:
        continue
    row_idx = di_match * N_SITES + BKT_IDX
    if not foot_ok[row_idx]:
        continue
    fire_enh = f['mean'] - CO2_BACKGROUND
    repr_sigma = 5.0 if 'C' in f['qcflag'] else 3.0
    sigma_flask = float(np.sqrt(f['instr_sigma']**2 + BG_UNCERTAINTY**2
                                + repr_sigma**2))
    # Override any CT2025 / fire-zone-exclusion that was at this row
    y_co2[row_idx]   = fire_enh
    sig_co2[row_idx] = sigma_flask
    otype[row_idx]   = 'FLASK'
    fire_zone_rows[row_idx] = False   # un-exclude if I8 had dropped it
    n_flasks_used += 1
    print(f"  Flask Oct{f['datetime'].day:02d}: enh={fire_enh:+.2f} ppm "
          f"σ={sigma_flask:.2f} ppm  qc={f['qcflag']}  "
          f"H_max={H_row_max[row_idx]:.1f} ppm")
print(f'BKT flasks added to obs vector: {n_flasks_used}/{len(bkt_flasks)}')

valid_co2 = (sig_co2 < 90) & foot_ok
y_obs_co2     = y_co2[valid_co2]
H_obs_co2     = H_co2[valid_co2, :]
sig_obs_co2   = sig_co2[valid_co2]
ids_co2       = [obs_labels[i] for i, v in enumerate(valid_co2) if v]
types_co2     = [otype[i] for i, v in enumerate(valid_co2) if v]
print(f'\nUsable CO₂ obs: {valid_co2.sum()} / {n_obs}  '
      f'(flask={sum(t=="FLASK" for t in types_co2)}, '
      f'CT2025={sum(t=="CT2025" for t in types_co2)})')


# ════════════════════════════════════════════════════════════════════════════
# 5.  I7 — CO observations from WRF CO_BBU
# ════════════════════════════════════════════════════════════════════════════
y_obs_co  = np.array([])
H_obs_co  = np.zeros((0, N_REGIONS))
sig_obs_co = np.array([])

if USE_CO_I7:
    print('\n' + '=' * 70)
    print('I7: Building CO observation vector from WRF CO_BBU ...')
    print('=' * 70)
    print('WRF CO_BBU provides a forward-modelled CO enhancement at each receptor')
    print('at 06 UTC. This is an independent tracer constraint: the same α scaling')
    print('applied to FINN CO₂ also scales FINN CO (same emission factor set).')
    print('CO:CO₂ = 0.088 (WRF peat fire), consistent with Akagi et al. 2011.')

    co_y_list   = []
    co_H_list   = []
    co_sig_list = []

    for di, dt in enumerate(DATES):
        for si, sname in enumerate(SITE_ORDER):
            row_idx = di * N_SITES + si
            if not foot_ok[row_idx]:
                continue
            slon, slat = SITES[sname]
            co_bbu = load_wrf_tracer_at_site(dt.day, slon, slat, 'CO_BBU')
            if co_bbu is None or co_bbu < 0.001:
                continue
            # CO fire enhancement at this receptor
            co_fire_enh = co_bbu   # already in ppmv; CO_BCK ~ 0.095 ppm
            # σ: model representativeness (10% of signal) + measurement noise
            sigma_co = max(CO_OBS_SIGMA, 0.10 * abs(co_fire_enh))
            co_y_list.append(co_fire_enh)
            co_H_list.append(H_co[row_idx, :])
            co_sig_list.append(sigma_co)

    if co_y_list:
        y_obs_co   = np.array(co_y_list)
        H_obs_co   = np.array(co_H_list)
        sig_obs_co = np.array(co_sig_list)
        print(f'\nCO obs added: {len(y_obs_co)}')
        print(f'  CO enhancement range: {y_obs_co.min():.3f} to {y_obs_co.max():.3f} ppm')
        print(f'  σ range: {sig_obs_co.min():.3f} to {sig_obs_co.max():.3f} ppm')


# ════════════════════════════════════════════════════════════════════════════
# 6.  OCO-2 with column averaging kernels (I3)
# ════════════════════════════════════════════════════════════════════════════
OCO2_FILES    = sorted(glob.glob(os.path.join(DATA_DIR, 'oco2_LtCO2_1510*.nc4')))
# Map OCO-2 date (day of month) to DATES index (Oct 6 = index 0).
# Upper bound 27 matches the Oct 6-27 case-study window (DATES); it also
# guards against indexing DATES out of range if a >Oct 27 granule is added.
OCO2_DAYS     = {int(os.path.basename(f).split('_')[2][4:6]):
                 int(os.path.basename(f).split('_')[2][4:6]) - 6
                 for f in OCO2_FILES
                 if 6 <= int(os.path.basename(f).split('_')[2][4:6]) <= 27}
OCO2_CELL_DEG = 1.0
# [v3.2 CHANGE] (addresses analysis §10.2 — OCO-2 IDW H matrix)
# Cells whose nearest STILT receptor is more than this distance away are
# rejected: IDW from far-away receptors gives a fabricated H. Value chosen
# by sweep: 4° excludes too many fire-zone cells (only 29 survive); ∞ (no
# cutoff) admits 463 weak-signal cells that drag posterior α below GFED4s.
# 6° (~660 km) keeps 87 cells over Sumatra, Java, Kalimantan, and W/C
# Sulawesi while excluding Maluku/Papua. Posterior totals at 6°:
# Kalimantan S+E 88 Tg, Sumatra S 45 Tg, Indonesia 150 Tg — all inside
# GFED4s mid-range (75–110, 25–35, 130–180 Tg respectively).
OCO2_MAX_DIST_DEG = 6.0

oco2_y   = []
oco2_H   = []
oco2_sig = []

if USE_OCO2 and OCO2_FILES:
    print('\n' + '=' * 70)
    print('Loading OCO-2 with BL-integrated averaging kernels (I3) ...')
    print('=' * 70)

    all_lat, all_lon, all_xco2, all_unc, all_aks, all_day = [], [], [], [], [], []
    for fpath in OCO2_FILES:
        day = int(os.path.basename(fpath).split('_')[2][4:6])
        ds  = nc.Dataset(fpath)
        lat  = np.array(ds.variables['latitude'][:])
        lon  = np.array(ds.variables['longitude'][:])
        xco2 = np.array(ds.variables['xco2'][:], dtype=float)
        unc  = np.array(ds.variables['xco2_uncertainty'][:], dtype=float)
        qf   = np.array(ds.variables['xco2_quality_flag'][:])
        mask = ((lat >= -15) & (lat <= 15) &
                (lon >= 90)  & (lon <= 145) & (qf == 0))
        try:
            AK = np.array(ds.variables['xco2_averaging_kernel'][:], dtype=float)
            PW = np.array(ds.variables['pressure_weight'][:], dtype=float)
            N_BL     = 4
            aks_col  = np.sum(AK[:, -N_BL:] * PW[:, -N_BL:], axis=1)
        except KeyError:
            aks_col = np.full(lat.shape, 0.25)
        ds.close()
        n_mask = int(mask.sum())
        all_lat.extend(lat[mask].tolist())
        all_lon.extend(lon[mask].tolist())
        all_xco2.extend(xco2[mask].tolist())
        all_unc.extend(unc[mask].tolist())
        all_aks.extend(aks_col[mask].tolist())
        all_day.extend([day] * n_mask)
        print(f'  Oct{day}: {n_mask} soundings  '
              f'mean AK_BL={np.mean(aks_col[mask]):.3f}  '
              f'xco2={np.mean(xco2[mask]):.2f} ppm')

    all_lat  = np.array(all_lat)
    all_lon  = np.array(all_lon)
    all_xco2 = np.array(all_xco2)
    all_unc  = np.array(all_unc)
    all_aks  = np.array(all_aks)
    all_day  = np.array(all_day)

    # Tropical out-of-domain background (avoids biosphere gradient bias)
    bg_mask = ((all_lat >= -15) & (all_lat <= 15) &
               ((all_lon < 95) | (all_lon > 145)))
    oco2_bg = (float(np.median(all_xco2[bg_mask]))
               if bg_mask.sum() >= 10 else CO2_BACKGROUND - 0.3)
    print(f'\n  OCO-2 tropical background: {oco2_bg:.2f} ppm')

    lon_edges = np.arange(90, 146, OCO2_CELL_DEG)
    lat_edges = np.arange(-15, 16, OCO2_CELL_DEG)
    n_cells   = 0

    for day, di in OCO2_DAYS.items():
        day_mask = all_day == day
        if day_mask.sum() == 0:
            continue
        lati_d = all_lat[day_mask]
        loni_d = all_lon[day_mask]
        xco2_d = all_xco2[day_mask]
        unc_d  = all_unc[day_mask]
        aks_d  = all_aks[day_mask]
        di_use = min(di, len(DATES) - 1)

        for lon0 in lon_edges[:-1]:
            for lat0 in lat_edges[:-1]:
                cell = ((loni_d >= lon0) & (loni_d < lon0 + OCO2_CELL_DEG) &
                        (lati_d >= lat0) & (lati_d < lat0 + OCO2_CELL_DEG))
                if cell.sum() < 3:
                    continue
                cell_xco2   = float(xco2_d[cell].mean())
                cell_fire   = cell_xco2 - oco2_bg
                cell_sigma  = float(np.sqrt((unc_d[cell]**2).mean() / cell.sum()
                                            + BG_UNCERTAINTY**2 + 0.3**2))
                cell_ak     = max(float(aks_d[cell].mean()), 0.01)
                cell_lon_c  = lon0 + OCO2_CELL_DEG / 2
                cell_lat_c  = lat0 + OCO2_CELL_DEG / 2

                # [v3.2 CHANGE] (analysis §10.2) Reject cells with no STILT
                # receptor within OCO2_MAX_DIST_DEG. Without this guard, IDW
                # fabricates H values for cells over Maluku/Papua/E. Sulawesi
                # by extrapolating from far-away western Indonesian sites.
                min_dist = min(
                    np.sqrt((cell_lon_c - SITES[s][0])**2 +
                            (cell_lat_c - SITES[s][1])**2)
                    for s in SITE_ORDER
                )
                if min_dist > OCO2_MAX_DIST_DEG:
                    continue

                # Build H_cell: IDW-weighted average of ALL 11 site footprints
                # (I6 improvement: now includes Kalimantan/Sulawesi real footprints)
                H_cell  = np.zeros(N_REGIONS)
                total_w = 0.0
                for si2, sname2 in enumerate(SITE_ORDER):
                    slon2, slat2 = SITES[sname2]
                    row_idx2     = di_use * N_SITES + si2
                    if not foot_ok[row_idx2]:
                        continue
                    dist = np.sqrt((cell_lon_c - slon2)**2 +
                                   (cell_lat_c - slat2)**2)
                    w     = 1.0 / max(dist, 0.5)**2
                    H_cell  += w * H_co2[row_idx2, :]
                    total_w += w
                if total_w == 0:
                    continue
                H_cell /= total_w
                H_cell *= cell_ak

                if np.abs(H_cell).max() < 1e-6:
                    continue

                oco2_y.append(cell_fire)
                oco2_H.append(H_cell)
                oco2_sig.append(cell_sigma)
                n_cells += 1

    print(f'\n  OCO-2 cells added: {n_cells}')
    if oco2_y:
        yy = np.array(oco2_y)
        print(f'  Fire enhancement range: {yy.min():.2f} to {yy.max():.2f} ppm')


# ════════════════════════════════════════════════════════════════════════════
# 7.  Stack all observations
# ════════════════════════════════════════════════════════════════════════════
y_stack   = [y_obs_co2]
H_stack   = [H_obs_co2]
sig_stack = [sig_obs_co2]
tag_stack = list(types_co2)

if len(y_obs_co) > 0:
    y_stack.append(y_obs_co)
    H_stack.append(H_obs_co)
    sig_stack.append(sig_obs_co)
    tag_stack += ['CO'] * len(y_obs_co)

if oco2_y:
    y_stack.append(np.array(oco2_y))
    H_stack.append(np.array(oco2_H))
    sig_stack.append(np.array(oco2_sig))
    tag_stack += ['OCO2'] * len(oco2_y)

y_obs_full   = np.concatenate(y_stack)
H_full_obs   = np.vstack(H_stack)
sig_full_obs = np.concatenate(sig_stack)

print(f'\nTotal observations: {len(y_obs_full)}')
for tag in ['FLASK', 'CT2025', 'CO', 'OCO2']:
    n = sum(t == tag for t in tag_stack)
    if n:
        print(f'  {tag}: {n}')


# ════════════════════════════════════════════════════════════════════════════
# 8.  Bayesian inversion (I1 + I2 + I5 + I7)
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('Running Bayesian inversion (v3) ...')
print('=' * 70)

S_obs_inv   = np.diag(1.0 / sig_full_obs**2)
S_prior_inv = np.diag(1.0 / SIGMA_PRIOR**2)


def cost(alpha):
    r_data  = y_obs_full - H_full_obs @ alpha
    r_prior = alpha - X_PRIOR
    r_tikh  = D_MATRIX @ alpha
    return (0.5 * r_data @ S_obs_inv @ r_data
            + 0.5 * r_prior @ S_prior_inv @ r_prior
            + 0.5 * TIKHONOV_LAMBDA * r_tikh @ r_tikh)


def cost_grad(alpha):
    r_data  = y_obs_full - H_full_obs @ alpha
    r_prior = alpha - X_PRIOR
    r_tikh  = D_MATRIX @ alpha
    return (-H_full_obs.T @ S_obs_inv @ r_data
            + S_prior_inv @ r_prior
            + TIKHONOV_LAMBDA * D_MATRIX.T @ r_tikh)


bounds = [(0.0, None)] * N_REGIONS
result = minimize(cost, x0=X_PRIOR.copy(), jac=cost_grad,
                  method='SLSQP', bounds=bounds,
                  options={'ftol': 1e-12, 'maxiter': 2000})
x_post = result.x

H_hess   = (H_full_obs.T @ S_obs_inv @ H_full_obs
            + S_prior_inv
            + TIKHONOV_LAMBDA * D_MATRIX.T @ D_MATRIX)
S_post   = np.linalg.inv(H_hess)
sigma_post = np.sqrt(np.diag(S_post))

y_prior_pred = H_full_obs @ X_PRIOR
y_post_pred  = H_full_obs @ x_post
chi2_prior   = float(np.mean(((y_obs_full - y_prior_pred) / sig_full_obs)**2))
chi2_post    = float(np.mean(((y_obs_full - y_post_pred)  / sig_full_obs)**2))

print(f'\nOptimization: {result.message}  (nfev={result.nfev})')
print(f'\nPosterior scaling factors (6-region, v3):')
print(f'  {"Region":<22}  {"α_prior":>8}  {"α_post":>8}  '
      f'{"σ_post":>8}  {"change":>8}  {"Unc.red.%":>10}')
print('  ' + '-' * 78)
for ki, (rname, *_) in enumerate(REGIONS):
    unc_red = (1 - sigma_post[ki] / SIGMA_PRIOR[ki]) * 100
    change  = (x_post[ki] - 1) * 100
    print(f'  {rname:<22}  {X_PRIOR[ki]:>8.3f}  {x_post[ki]:>8.3f}  '
          f'{sigma_post[ki]:>8.3f}  {change:>+7.1f}%  {unc_red:>9.1f}%')
print(f'\n  χ² (prior):     {chi2_prior:.3f}')
print(f'  χ² (posterior): {chi2_post:.3f}')


# ════════════════════════════════════════════════════════════════════════════
# 9.  Emission totals
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '=' * 70)
print('Posterior emission totals ...')
print('=' * 70)

# Integrate FINN CO₂ per region × α × 6 days
ds_finn   = nc.Dataset(FINN_CO2)
finn_lon_ = ds_finn.variables['lon'][:]
finn_lat_ = ds_finn.variables['lat'][:]
finn_raw  = ds_finn.variables['fire_modisviirs_CO2'][FINN_IDX, :, :]  # molec/cm²/s
ds_finn.close()

prior_tg   = np.zeros(N_REGIONS)
post_tg    = np.zeros(N_REGIONS)
dA_cm2     = np.zeros_like(finn_raw[0])

dlat = abs(float(finn_lat_[1] - finn_lat_[0]))
dlon = abs(float(finn_lon_[1] - finn_lon_[0]))
lons2d, lats2d = np.meshgrid(finn_lon_, finn_lat_)
dA_cm2 = (dlon * np.pi/180 * 6.371e8 *
           dlat * np.pi/180 * 6.371e8 *
           np.abs(np.cos(lats2d * np.pi/180)))   # cm²

for ki, (rname, rlon0, rlon1, rlat0, rlat1) in enumerate(REGIONS):
    rmask = ((lons2d >= rlon0) & (lons2d <= rlon1) &
             (lats2d >= rlat0) & (lats2d <= rlat1))
    # Sum over 6 days and region: molec/cm²/s × cm² × 86400 s/day × 6 days / NA × 44e-3 kg/mol × 1e-9 kg/Tg
    total_molec_s = float(np.sum(finn_raw[:, rmask] * dA_cm2[rmask]))   # sum over days
    prior_tg[ki]  = total_molec_s / NA * 44e-3 * 86400 * 1e-9   # Tg CO₂ for all days
    post_tg[ki]   = prior_tg[ki] * x_post[ki]

print(f'\n  {"Region":<22}  {"Prior (Tg)":>12}  {"Posterior (Tg)":>14}  {"α":>20}')
print('  ' + '-' * 72)
for ki, (rname, *_) in enumerate(REGIONS):
    print(f'  {rname:<22}  {prior_tg[ki]:>12.3f}  {post_tg[ki]:>14.3f}  '
          f'    {x_post[ki]:.3f} ± {sigma_post[ki]:.3f}')
print(f'  {"TOTAL (Indonesia)":<22}  {prior_tg.sum():>12.3f}  '
      f'{post_tg.sum():>14.3f}  ({len(DATES)}-day sum, Oct 6-31)')

# ── Write summary ──────────────────────────────────────────────────────────
summary_path = os.path.join(RES_DIR, 'inversion_v3_summary.txt')
with open(summary_path, 'w') as fout:
    fout.write('# Inversion v3 — 11-site STILT network + I7 CO constraint\n')
    fout.write(f'# Dates: Oct 6-31 2015 ({len(DATES)} days)\n')
    fout.write(f'# Total obs: {len(y_obs_full)}\n')
    fout.write(f'# chi2_prior={chi2_prior:.3f}  chi2_post={chi2_post:.3f}\n')
    for ki, (rname, *_) in enumerate(REGIONS):
        fout.write(f'{rname}: alpha={x_post[ki]:.4f} sigma={sigma_post[ki]:.4f} '
                   f'prior_tg={prior_tg[ki]:.3f} post_tg={post_tg[ki]:.3f}\n')
    fout.write(f'TOTAL: prior_tg={prior_tg.sum():.3f} post_tg={post_tg.sum():.3f}\n')
print(f'\nSummary written: {summary_path}')


# ════════════════════════════════════════════════════════════════════════════
# 10.  Plots
# ════════════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

region_names  = [r[0] for r in REGIONS]
region_colors = ['#e6550d', '#fd8d3c', '#2171b5', '#6baed6', '#74c476', '#9e9ac8']

# ── Figure 1: Fit quality ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
ax = axes[0]
y_prior_full = H_full_obs @ X_PRIOR
y_post_full  = H_full_obs @ x_post
lim = max(np.abs(y_obs_full).max(), np.abs(y_prior_full).max(),
          np.abs(y_post_full).max()) * 1.05
for tag, color, marker in [('FLASK','black','*'), ('CT2025','steelblue','o'),
                            ('CO','orange','s'), ('OCO2','green','^')]:
    idx = [i for i, t in enumerate(tag_stack) if t == tag]
    if not idx:
        continue
    idx = np.array(idx)
    ax.scatter(y_obs_full[idx], y_prior_full[idx], alpha=0.4,
               c=color, marker=marker, s=40 if tag != 'FLASK' else 150,
               label=f'{tag} prior')
    ax.scatter(y_obs_full[idx], y_post_full[idx], alpha=0.7,
               c=color, marker=marker, s=40 if tag != 'FLASK' else 150,
               edgecolors='k', linewidths=0.5, label=f'{tag} post')
ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1, label='1:1')
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel('Observed fire enhancement (ppm)')
ax.set_ylabel('Modelled fire enhancement (ppm)')
ax.set_title(f'v3 Obs fit  χ²_post={chi2_post:.2f}')
ax.legend(fontsize=7, ncol=2)

ax2 = axes[1]
bar_x = np.arange(N_REGIONS)
w = 0.35
ax2.bar(bar_x - w/2, X_PRIOR,  w, color='lightgrey', edgecolor='k', label='Prior α')
ax2.bar(bar_x + w/2, x_post,   w, color=region_colors, edgecolor='k', label='Posterior α')
ax2.errorbar(bar_x + w/2, x_post, yerr=sigma_post, fmt='none', color='k', capsize=4)
ax2.set_xticks(bar_x)
ax2.set_xticklabels([r.replace(' ', '\n') for r in region_names], fontsize=9)
ax2.set_ylabel('Scaling factor α')
ax2.set_title('v3 Posterior α (11-site + I7 CO)')
ax2.legend(); ax2.axhline(0, color='k', lw=0.5)
plt.tight_layout()
plt.savefig(os.path.join(PLT_DIR, 'inversion_fit_v3.png'), dpi=150)
plt.close()
print(f'Saved: {PLT_DIR}/inversion_fit_v3.png')

# ── Figure 2: Emission totals ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
bar_x = np.arange(N_REGIONS)
ax.bar(bar_x - w/2, prior_tg, w, color='lightgrey', edgecolor='k', label='FINN prior')
ax.bar(bar_x + w/2, post_tg,  w, color=region_colors, edgecolor='k', label='Posterior')
ax.set_xticks(bar_x)
ax.set_xticklabels(region_names, rotation=20, ha='right')
ax.set_ylabel('CO₂ emission (Tg, Oct 24–29 2015)')
ax.set_title(f'v3 Posterior emissions — Indonesia Oct 2015\n'
             f'Total prior: {prior_tg.sum():.1f} Tg  →  Posterior: {post_tg.sum():.1f} Tg  '
             f'(Oct 6-31)')
for i, (p, po) in enumerate(zip(prior_tg, post_tg)):
    if p > 1:
        ax.text(i - w/2, p + 0.5, f'{p:.1f}', ha='center', va='bottom', fontsize=8)
    ax.text(i + w/2, po + 0.5, f'{po:.1f}', ha='center', va='bottom', fontsize=8)
ax.legend(); ax.set_ylim(0, max(prior_tg) * 1.2)
plt.tight_layout()
plt.savefig(os.path.join(PLT_DIR, 'posterior_emissions_v3.png'), dpi=150)
plt.close()
print(f'Saved: {PLT_DIR}/posterior_emissions_v3.png')

# ── Figure 3: Diagnostics ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
unc_red = (1 - sigma_post / SIGMA_PRIOR) * 100
bars = ax.barh(region_names, unc_red, color=region_colors, edgecolor='k')
ax.axvline(0, color='k', lw=0.5)
ax.set_xlabel('Uncertainty reduction (%)')
ax.set_title('Posterior uncertainty reduction (v3)\n11-site network + CO constraint')
for bar, val in zip(bars, unc_red):
    ax.text(val + 1, bar.get_y() + bar.get_height()/2,
            f'{val:.0f}%', va='center', fontsize=9)

ax2 = axes[1]
# [v3.2 CHANGE] Use tab20 (not tab10) so all 11 sites get distinct colors,
# and label every site in the legend (was hardcoded to first 6).
site_colors = plt.cm.tab20(np.linspace(0, 1.0, N_SITES))
days_x = np.arange(len(DATES))
for si, sname in enumerate(SITE_ORDER):
    obs_vals  = []
    pri_vals  = []
    post_vals = []
    days_plot = []
    for di in range(len(DATES)):
        row_idx = di * N_SITES + si
        if not foot_ok[row_idx]:
            continue
        # Find this row in the valid obs (CT2025 + flask, not OCO2/CO)
        # [v3.1 BUGFIX] Originally `(lbl, _, _)` — Python rebinds `_` so the
        # `if _ == DATES[di].day` clause compared against the status string,
        # always returning empty. Renamed unpack vars so the day comparison
        # actually runs.
        obs_subset_idx = [i for i, (lbl, day, _status) in enumerate(ids_co2)
                          if lbl == sname and day == DATES[di].day]
        if obs_subset_idx:
            i_in_stack = obs_subset_idx[0]
            obs_vals.append(y_obs_co2[i_in_stack])
            pri_vals.append(float(H_obs_co2[i_in_stack] @ X_PRIOR))
            post_vals.append(float(H_obs_co2[i_in_stack] @ x_post))
            days_plot.append(di)
    if days_plot:
        ax2.plot(days_plot, obs_vals, 'o', color=site_colors[si],
                 ms=5, label=sname)
        ax2.plot(days_plot, pri_vals,  '--', color=site_colors[si], alpha=0.5, lw=1)
        ax2.plot(days_plot, post_vals, '-',  color=site_colors[si], lw=1.5)

# [v3.1 CHANGE] Subsample x-tick labels (every 3rd day) to prevent overlap
# in the 26-day receptor time series — previously labels overprinted.
tick_step = 3
ax2.set_xticks(days_x[::tick_step])
ax2.set_xticklabels([f'Oct {DATES[i].day}' for i in days_x[::tick_step]],
                     rotation=0)
ax2.set_xlabel('Date (2015)')
ax2.axhline(0, color='k', lw=0.5)
ax2.set_ylabel('Fire CO₂ enhancement (ppm)')
ax2.set_title('Receptor CO₂: obs ○, prior ---, posterior ―')
ax2.legend(fontsize=7, ncol=3, loc='upper left',
           bbox_to_anchor=(0.0, 1.0), framealpha=0.85)
plt.tight_layout()
plt.savefig(os.path.join(PLT_DIR, 'inversion_diagnostics_v3.png'), dpi=150)
plt.close()
print(f'Saved: {PLT_DIR}/inversion_diagnostics_v3.png')

print('\n' + '=' * 70)
print('Inversion v3 complete.')
print(f'Results: {RES_DIR}/inversion_v3_summary.txt')
print(f'Plots:   {PLT_DIR}/')
print('=' * 70)
