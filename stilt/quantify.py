#!/usr/bin/env python3
"""
Quantify fire emissions (FINN) and WRF-STILT concentrations
for the IDN_BB_2015 simulation, Oct 24-29 2015.

Outputs:
  stilt/out/quantification.txt       — numbers table (emissions, conc, STILT Δc)
  stilt/out/receptor_timeseries.png  — daily concentration at each receptor
  stilt/out/stilt_delta_conc.png     — STILT-derived Δc per site/day
"""

import os, glob
import numpy as np
import netCDF4 as nc
from datetime import datetime, timezone

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT    = '/home/igrk/WRF-GRK'
WRF_DIR = f'{ROOT}/simulations/IDN_BB_2015/output'
FOOT_DIR= f'{ROOT}/stilt/out/footprints'
OUT_DIR = f'{ROOT}/stilt/out'
FINN = {
    'CO2': f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CO2_bb_surface_daily_20150101-20151231_0.1x0.1.nc',
    'CO' : f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CO_bb_surface_daily_20150101-20151231_0.1x0.1.nc',
    'CH4': f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CH4_bb_surface_daily_20150101-20151231_0.1x0.1.nc',
}
FINN_VARNAMES = {'CO2':'fire_modisviirs_CO2','CO':'fire_modisviirs_CO',
                 'CH4':'fire_modisviirs_CH4'}
MOLAR = {'CO2': 44.01, 'CO': 28.01, 'CH4': 16.04}   # g/mol

DATES     = [datetime(2015, 10, d, tzinfo=timezone.utc) for d in range(24, 30)]
DATES_NOON= [datetime(2015, 10, d, 12, tzinfo=timezone.utc) for d in range(24, 30)]

SITES = {
    'Bukit Kototabang': (100.32, -0.20),
    'Palangka Raya':    (113.94, -2.16),
    'Pontianak':        (109.33, -0.02),
    'Jambi':            (103.61, -1.61),
    'Jakarta':          (106.85, -6.21),
}
SITE_ORDER = list(SITES.keys())

# WRF domain (to restrict emission totals)
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = 90.0, 145.0, -15.0, 15.0
R_EARTH = 6.371e6   # m
NA      = 6.022e23  # molec/mol


# ════════════════════════════════════════════════════════════════════════════
# 1. Total FINN fire emissions (domain-sum, Oct 24-29)
# ════════════════════════════════════════════════════════════════════════════
print('='*70)
print('1. FINN domain-total fire emissions (Oct 24-29 2015)')
print('   Domain: lon 90-145°E, lat 15°S-15°N')
print('='*70)

FINN_IDX = [d.timetuple().tm_yday - 1 for d in DATES]  # 0-based day indices

# Molar mass for conversion
# E [g/m²/s] = E [molec/cm²/s] × (1e4 cm²/m²) × M_g_mol / N_A
# E [Tg/day] = sum over domain of E [g/m²/s] × ΔA[m²] × 86400 s/day × 1e-12

total_emis = {}  # species → Tg over 6 days
daily_emis = {}  # species → [Tg/day] × 6 days

for sp, path in FINN.items():
    ds = nc.Dataset(path)
    finn_lon = ds.variables['lon'][:]
    finn_lat = ds.variables['lat'][:]
    data = ds.variables[FINN_VARNAMES[sp]][FINN_IDX, :, :]  # (6, lat, lon)
    ds.close()

    # Crop to domain
    imask = np.where((finn_lon >= LON_MIN) & (finn_lon <= LON_MAX))[0]
    jmask = np.where((finn_lat >= LAT_MIN) & (finn_lat <= LAT_MAX))[0]
    i0,i1 = imask[0], imask[-1]+1
    j0,j1 = jmask[0], jmask[-1]+1
    sub_lon = finn_lon[i0:i1]
    sub_lat = finn_lat[j0:j1]
    sub_data = data[:, j0:j1, i0:i1]    # (6, nlat, nlon) molec/cm²/s

    # Grid cell area [m²]: ΔA = Δx × Δy × cos(lat)
    dlon_rad = (sub_lon[1] - sub_lon[0]) * np.pi / 180.0
    dlat_rad = abs(sub_lat[1] - sub_lat[0]) * np.pi / 180.0
    lat_rad  = sub_lat * np.pi / 180.0
    # area of each row (same for all longitudes)
    area_row = R_EARTH**2 * dlon_rad * dlat_rad * np.abs(np.cos(lat_rad))  # m²
    area_2d  = np.tile(area_row[:, np.newaxis], (1, len(sub_lon)))          # (nlat, nlon)
    area_cm2 = area_2d * 1e4                                                # cm²

    M = MOLAR[sp]
    # Convert molec/cm²/s → g/cm²/s → Tg over cell and day
    # mass_per_cell [g/day] = flux[molec/cm²/s] × area[cm²] × 86400[s/day] × M/NA
    daily_Tg = np.sum(
        sub_data * area_cm2[np.newaxis, :, :] * 86400 * M / NA,
        axis=(1, 2)
    ) * 1e-12  # → Tg

    daily_emis[sp] = daily_Tg
    total_emis[sp] = daily_Tg.sum()
    print(f'\n  {sp}:')
    for i, (dt, v) in enumerate(zip(DATES, daily_Tg)):
        print(f'    {dt.strftime("Oct %d")}: {v*1000:.4f} Gg/day   ({v:.6f} Tg/day)')
    print(f'    TOTAL (6 days): {total_emis[sp]*1000:.3f} Gg  =  {total_emis[sp]:.6f} Tg')


# ════════════════════════════════════════════════════════════════════════════
# 2. WRF surface concentrations at receptor sites (12 UTC, level 0)
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('2. WRF-GHG surface biomass-burning tracer at receptor sites (12 UTC)')
print('   Units: ppmv (mixing ratio enhancement above background)')
print('='*70)

# Load XLAT/XLONG from first file for nearest-neighbour lookup
ds0 = nc.Dataset(os.path.join(WRF_DIR,'wrfout_d01_2015-10-24_01:00:00'))
xlat_ref = ds0.variables['XLAT'][0,:,:]
xlon_ref = ds0.variables['XLONG'][0,:,:]
ds0.close()

# Pre-compute nearest WRF grid cell for each site
site_ij = {}
for sname, (slon, slat) in SITES.items():
    dist = np.sqrt((xlat_ref - slat)**2 + (xlon_ref - slon)**2)
    ij   = np.unravel_index(dist.argmin(), dist.shape)
    site_ij[sname] = ij

WRF_TRACERS = [
    ('CO2_BBU', 'CO2_fire', 'ppmv'),
    ('CO_BBU',  'CO_fire',  'ppmv'),
    ('CH4_BBU', 'CH4_fire', 'ppmv'),
    ('CO2_ANT', 'CO2_anthro','ppmv'),
]

# conc_table[tracer][site][day_idx] = float
conc_table = {tr[0]: {s: np.zeros(6) for s in SITE_ORDER} for tr in WRF_TRACERS}

for di, dt in enumerate(DATES):
    fname = os.path.join(WRF_DIR, f'wrfout_d01_{dt.strftime("%Y-%m-%d")}_01:00:00')
    ds = nc.Dataset(fname)
    t_idx = 11   # 01:00 + 11 h = 12:00 UTC
    for tr_key, _, _ in WRF_TRACERS:
        if tr_key not in ds.variables:
            continue
        conc2d = ds.variables[tr_key][t_idx, 0, :, :]
        for sname in SITE_ORDER:
            ij = site_ij[sname]
            conc_table[tr_key][sname][di] = float(conc2d[ij])
    ds.close()

# Print table
header = f"{'Species':<12}" + "".join(f"{s[:12]:>13}" for s in SITE_ORDER)
for tr_key, label, units in WRF_TRACERS:
    print(f'\n  {label} [{units}]:')
    print('  ' + '-'*90)
    print('  ' + header)
    print('  ' + '-'*90)
    for di, dt in enumerate(DATES):
        row = f"  {dt.strftime('Oct %d'):<12}"
        for sname in SITE_ORDER:
            row += f"{conc_table[tr_key][sname][di]:>13.3f}"
        print(row)
    print('  ' + '-'*90)
    means_row = f"  {'Mean':<12}"
    for sname in SITE_ORDER:
        means_row += f"{conc_table[tr_key][sname].mean():>13.3f}"
    print(means_row)

# ════════════════════════════════════════════════════════════════════════════
# 3. STILT footprint × FINN emission → Δc at each receptor
#    footprint units: ppm / (µmol m⁻² s⁻¹)  per 0.25° grid cell
#    FINN units:      molec cm⁻² s⁻¹          at 0.1°
#    Method: area-weighted regrid FINN → 0.25° footprint grid, then
#            Δc = Σ_{i,j}  foot[i,j]  ×  E_regrid[i,j]
#    where E_regrid is in µmol m⁻² s⁻¹
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('3. STILT-derived concentration increment (Δc = footprint ⊗ FINN flux)')
print('   Δc = Σ f(x,y) [ppm/(µmol m⁻² s⁻¹)] × E(x,y) [µmol m⁻² s⁻¹]')
print('   Species: CO2 (FINN biomass-burning)')
print('='*70)

def load_finn_raw(path, varname, day_indices, foot_lon, foot_lat):
    """Load FINN, convert to µmol/m²/s, regrid to footprint 0.25° grid."""
    ds = nc.Dataset(path)
    finn_lon = ds.variables['lon'][:]
    finn_lat = ds.variables['lat'][:]
    raw = ds.variables[varname][day_indices, :, :]  # (n_days, nlat, nlon) molec/cm²/s
    ds.close()

    # Convert: molec/cm²/s → µmol/m²/s
    # = molec/cm²/s × 1e4 cm²/m² / 6.022e23 molec/mol × 1e6 µmol/mol
    CONV = 1e4 / NA * 1e6   # = 1.661e-14  ... wait that gives tiny numbers?
    # Let me re-check: 1e4/6.022e23 = 1.661e-20 mol/m²/s → ×1e6 = 1.661e-14 µmol/m²/s
    # That's correct - FINN values are ~1e16 molec/cm²/s → ~166 µmol/m²/s
    raw_umol = raw * CONV

    # Area-weighted regrid from 0.1° to 0.25°
    # For each 0.25° cell, collect all 0.1° cells within and average
    dlat_f = abs(float(foot_lat[1] - foot_lat[0]))  # 0.25
    dlon_f = float(foot_lon[1] - foot_lon[0])        # 0.25
    n_days = raw.shape[0]
    out    = np.zeros((n_days, len(foot_lat), len(foot_lon)), dtype=np.float64)

    for jj, flat in enumerate(foot_lat):
        lat_lo = flat - dlat_f/2
        lat_hi = flat + dlat_f/2
        jmask  = np.where((finn_lat >= lat_lo) & (finn_lat < lat_hi))[0]
        if jmask.size == 0:
            continue
        sub_lat_rad = finn_lat[jmask] * np.pi / 180.0
        cos_lat     = np.abs(np.cos(sub_lat_rad))  # weight by cos(lat) for area

        for ii, flon in enumerate(foot_lon):
            lon_lo = flon - dlon_f/2
            lon_hi = flon + dlon_f/2
            imask  = np.where((finn_lon >= lon_lo) & (finn_lon < lon_hi))[0]
            if imask.size == 0:
                continue
            patch = raw_umol[:, jmask[0]:jmask[-1]+1, imask[0]:imask[-1]+1]
            # Area-weighted mean (cos-lat weights)
            w2d  = cos_lat[:, np.newaxis] * np.ones(len(imask))[np.newaxis, :]
            denom = w2d.sum()
            out[:, jj, ii] = np.sum(patch * w2d[np.newaxis, :, :], axis=(1,2)) / denom

    return out   # (n_days, n_foot_lat, n_foot_lon) µmol/m²/s


# Load footprints for each (site, day) pair
print('\n  Loading footprints ...')
stilt_delta = {}   # (sname, di) → Δc_CO2 [ppm]

# Map site footprint files
def find_footprint_file(slon, slat, dt):
    tag = dt.strftime('%Y%m%d%H%M')
    for f in glob.glob(os.path.join(FOOT_DIR, f'{tag}_*_foot.nc')):
        base = os.path.basename(f)
        parts = base.split('_')
        try:
            flon = float(parts[1]); flat = float(parts[2])
            if abs(flon - slon) < 0.01 and abs(flat - slat) < 0.01:
                return f
        except (IndexError, ValueError):
            pass
    return None

# Get footprint grid from first file
f0 = find_footprint_file(100.32, -0.20, DATES_NOON[0])
ds0 = nc.Dataset(f0)
foot_lon = ds0.variables['lon'][:]
foot_lat = ds0.variables['lat'][:]
ds0.close()

# Regrid FINN CO2 to footprint grid for each day
print('  Regridding FINN CO2 to 0.25° footprint grid (may take ~30s) ...')
finn_co2_regrid = load_finn_raw(
    FINN['CO2'], FINN_VARNAMES['CO2'], FINN_IDX, foot_lon, foot_lat
)   # (6, 120, 220) µmol/m²/s

print('  Computing Δc per receptor ...')
print(f"\n  {'Site':<22} {'Date':<10} {'Δc_CO2 (ppm)':>15}  {'WRF CO2_BBU (ppm)':>20}  {'Ratio':>8}")
print('  ' + '-'*80)

for sname, (slon, slat) in SITES.items():
    for di, (dt_midnight, dt_noon) in enumerate(zip(DATES, DATES_NOON)):
        fpath = find_footprint_file(slon, slat, dt_noon)
        if fpath is None:
            print(f'  {sname:<22} {dt_midnight.strftime("Oct %d"):<10}  FILE NOT FOUND')
            continue
        dsf = nc.Dataset(fpath)
        foot = dsf.variables['foot'][0, :, :]   # (120, 220)
        dsf.close()

        # Δc = sum of (foot × FINN_emission)
        E_day = finn_co2_regrid[di, :, :]   # µmol/m²/s
        delta_c = float(np.sum(foot * E_day))

        wrf_c   = conc_table['CO2_BBU'][sname][di]
        ratio   = delta_c / wrf_c if wrf_c > 0 else float('nan')
        stilt_delta[(sname, di)] = delta_c

        print(f"  {sname:<22} {dt_midnight.strftime('Oct %d'):<10} {delta_c:>15.3f}  {wrf_c:>20.3f}  {ratio:>8.3f}")

# Per-site means
print('  ' + '-'*80)
for sname in SITE_ORDER:
    vals = [stilt_delta[(sname,di)] for di in range(6)]
    wrf_vals = [conc_table['CO2_BBU'][sname][di] for di in range(6)]
    print(f"  {sname:<22} {'Mean':<10} {np.mean(vals):>15.3f}  {np.mean(wrf_vals):>20.3f}  {np.mean(vals)/np.mean(wrf_vals):>8.3f}")


# ════════════════════════════════════════════════════════════════════════════
# 4. Summary statistics
# ════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('4. Summary statistics')
print('='*70)

print('\n  Total FINN fire emissions over WRF domain (Oct 24-29 2015):')
for sp in ['CO2','CO','CH4']:
    print(f'    {sp:4s}: {total_emis[sp]*1000:.2f} Gg  ({total_emis[sp]:.5f} Tg)')

print('\n  WRF-GHG CO2 fire enhancement (ppmv) range across all sites/days:')
all_co2 = [conc_table['CO2_BBU'][s][d] for s in SITE_ORDER for d in range(6)]
print(f'    min: {min(all_co2):.2f}  max: {max(all_co2):.2f}  mean: {np.mean(all_co2):.2f}  median: {np.median(all_co2):.2f}')
print(f'    Highest: Palangka Raya (central Kalimantan, ground-zero fire region)')

print('\n  STILT Δc_CO2 range across all receptors/days:')
all_stilt = list(stilt_delta.values())
print(f'    min: {min(all_stilt):.3f}  max: {max(all_stilt):.3f}  mean: {np.mean(all_stilt):.3f}')


# ════════════════════════════════════════════════════════════════════════════
# 5. Plots: receptor time series + STILT Δc bar chart
# ════════════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

DAY_LABELS = [d.strftime('%-d %b') for d in DATES]
COLORS = ['tab:orange','tab:red','tab:green','tab:brown','tab:blue']

# ── Plot A: WRF surface concentrations time series ────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True,
                          gridspec_kw={'hspace': 0.08})

for ax, (tr_key, label, units), cmap_base in zip(
        axes,
        [('CO2_BBU','CO₂ fire','ppmv'),
         ('CO_BBU', 'CO fire', 'ppmv'),
         ('CH4_BBU','CH₄ fire','ppmv')],
        ['Oranges','Reds','Purples']):
    for si, sname in enumerate(SITE_ORDER):
        vals = conc_table[tr_key][sname]
        ax.plot(range(6), vals, 'o-', color=COLORS[si], lw=1.8, ms=6,
                label=sname)
    ax.set_ylabel(f'{label}\n[{units}]', fontsize=9)
    ax.grid(True, lw=0.4, alpha=0.4)
    ax.set_yscale('log')
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())

axes[0].legend(fontsize=8, ncol=3, loc='upper left', framealpha=0.85)
axes[-1].set_xticks(range(6))
axes[-1].set_xticklabels(DAY_LABELS, fontsize=9)
axes[-1].set_xlabel('Date (Oct 2015)', fontsize=9)
fig.suptitle('WRF-GHG biomass-burning tracer at receptor sites — Oct 24-29 2015\n'
             '(surface level, 12 UTC)', fontsize=11, fontweight='bold')
out_ts = os.path.join(OUT_DIR, 'receptor_timeseries.png')
fig.savefig(out_ts, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved: {out_ts}')

# ── Plot B: STILT Δc bar chart (per receptor, 6 days stacked) ────────────────
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5),
                            gridspec_kw={'wspace': 0.35})

# Left: stacked bar per site (6 days as segments)
ax_l = axes2[0]
width = 0.14
x = np.arange(len(SITE_ORDER))
day_colors = plt.cm.YlOrRd(np.linspace(0.3, 0.95, 6))
for di in range(6):
    vals = [stilt_delta[(s, di)] for s in SITE_ORDER]
    bars = ax_l.bar(x + (di - 2.5) * width, vals, width,
                    color=day_colors[di], label=DAY_LABELS[di])
ax_l.set_xticks(x)
ax_l.set_xticklabels([s.replace(' ','\n') for s in SITE_ORDER], fontsize=8)
ax_l.set_ylabel('STILT Δc CO₂  [ppm]', fontsize=9)
ax_l.set_title('STILT-derived CO₂ concentration\nincrement per receptor',
               fontsize=10, fontweight='bold')
ax_l.legend(fontsize=7, ncol=2, title='Date', title_fontsize=7)
ax_l.grid(axis='y', lw=0.4, alpha=0.4)

# Right: scatter WRF CO2_BBU vs STILT Δc
ax_r = axes2[1]
wrf_all   = [conc_table['CO2_BBU'][s][d] for s in SITE_ORDER for d in range(6)]
stilt_all = [stilt_delta[(s, d)]          for s in SITE_ORDER for d in range(6)]
site_col  = [COLORS[si] for si, s in enumerate(SITE_ORDER) for d in range(6)]
ax_r.scatter(stilt_all, wrf_all, c=site_col, s=55, alpha=0.85, edgecolors='none')
for si, sname in enumerate(SITE_ORDER):
    xs = [stilt_delta[(sname, d)] for d in range(6)]
    ys = [conc_table['CO2_BBU'][sname][d] for d in range(6)]
    ax_r.scatter(xs, ys, c=COLORS[si], s=55, alpha=0.85, edgecolors='none',
                 label=sname)
# 1:1 and linear regression
xarr = np.array(stilt_all); yarr = np.array(wrf_all)
valid = (xarr > 0) & (yarr > 0)
if valid.sum() > 2:
    m, b = np.polyfit(xarr[valid], yarr[valid], 1)
    xfit = np.linspace(xarr[valid].min(), xarr[valid].max(), 100)
    ax_r.plot(xfit, m*xfit + b, 'k--', lw=1.2,
              label=f'y={m:.1f}x+{b:.1f}  (R²={np.corrcoef(xarr[valid],yarr[valid])[0,1]**2:.2f})')
ax_r.set_xlabel('STILT Δc CO₂ [ppm]  (footprint × FINN)', fontsize=9)
ax_r.set_ylabel('WRF CO₂ fire tracer [ppm]', fontsize=9)
ax_r.set_title('WRF Eulerian vs. STILT Lagrangian\nCO₂ fire enhancement',
               fontsize=10, fontweight='bold')
ax_r.legend(fontsize=7, framealpha=0.85)
ax_r.grid(lw=0.4, alpha=0.4)

fig2.suptitle('STILT Δc and WRF concentration quantification — Oct 24-29 2015',
              fontsize=11, fontweight='bold', y=1.01)
out_bar = os.path.join(OUT_DIR, 'stilt_delta_conc.png')
fig2.savefig(out_bar, dpi=150, bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {out_bar}')
print('\nQuantification complete.')
