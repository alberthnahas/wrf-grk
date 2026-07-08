#!/usr/bin/env python3
"""
Plot fire emissions (FINN) and WRF-GHG surface concentrations for the
IDN_BB_2015 biomass-burning simulation, Oct 24–29 2015.

Outputs:
  stilt/out/emission_finn.png        — FINN CO2 fire flux, 6-day panel
  stilt/out/concentration_wrf.png    — WRF surface CO2/CO/CH4 BBU, 6-day panels
  stilt/out/emission_vs_footprint.png — FINN mean CO2 emission + STILT footprint
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
from datetime import datetime, timezone

# ── Paths ─────────────────────────────────────────────────────────────────────
ROOT       = '/home/igrk/WRF-GRK'
FINN_CO2   = f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CO2_bb_surface_daily_20150101-20151231_0.1x0.1.nc'
FINN_CO    = f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CO_bb_surface_daily_20150101-20151231_0.1x0.1.nc'
FINN_CH4   = f'{ROOT}/rawdata/finn_fire/emissions-finnv2.5modvrs_CH4_bb_surface_daily_20150101-20151231_0.1x0.1.nc'
WRF_DIR    = f'{ROOT}/simulations/IDN_BB_2015/output'
FOOT_DIR   = f'{ROOT}/stilt/out/footprints'
OUT_DIR    = f'{ROOT}/stilt/out'
os.makedirs(OUT_DIR, exist_ok=True)

# ── Config ────────────────────────────────────────────────────────────────────
DATES      = [datetime(2015, 10, d, tzinfo=timezone.utc) for d in range(24, 30)]
DATES_NOON = [datetime(2015, 10, d, 12, tzinfo=timezone.utc) for d in range(24, 30)]
DATE_LABELS = [d.strftime('Oct %-d') for d in DATES]

# WRF domain extent
LON_MIN, LON_MAX = 90.0, 145.0
LAT_MIN, LAT_MAX = -15.0, 15.0

# Receptor sites
SITES = {
    'Bukit Kototabang': (100.32, -0.20),
    'Palangka Raya':    (113.94, -2.16),
    'Pontianak':        (109.33, -0.02),
    'Jambi':            (103.61, -1.61),
    'Jakarta':          (106.85, -6.21),
}

proj = ccrs.PlateCarree()


# ── Helper: WRF land/coast features ──────────────────────────────────────────
def add_map_features(ax, linewidth=0.4):
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=proj)
    ax.add_feature(cfeature.LAND,      facecolor='#f0ede8', linewidth=0)
    ax.add_feature(cfeature.OCEAN,     facecolor='#d5eaf5', linewidth=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=linewidth, edgecolor='#333')
    ax.add_feature(cfeature.BORDERS,   linewidth=linewidth*0.7,
                   edgecolor='#777', linestyle=':')


def add_sites(ax, ms=6, color='blue'):
    for sname, (slon, slat) in SITES.items():
        ax.plot(slon, slat, marker='*', color=color, markersize=ms,
                transform=proj, zorder=6)


def sparse_gridlines(ax, dx=15, dy=10):
    gl = ax.gridlines(draw_labels=False, linewidth=0.3,
                      color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(dx)
    gl.ylocator = mticker.MultipleLocator(dy)
    return gl


# ── Load FINN on WRF-domain sub-grid ─────────────────────────────────────────
def load_finn(path, varname, date_indices):
    """Returns (lon1d, lat1d, data3d[ndays, nlat, nlon]) cropped to WRF domain."""
    ds = nc.Dataset(path)
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    imask = np.where((lons >= LON_MIN - 0.5) & (lons <= LON_MAX + 0.5))[0]
    jmask = np.where((lats >= LAT_MIN - 0.5) & (lats <= LAT_MAX + 0.5))[0]
    i0, i1 = imask[0], imask[-1] + 1
    j0, j1 = jmask[0], jmask[-1] + 1
    data = ds.variables[varname][date_indices, j0:j1, i0:i1]
    data = np.ma.masked_values(data, 0.0)
    ds.close()
    return lons[i0:i1], lats[j0:j1], data


# Get Oct 24-29 indices in FINN (day-of-year 297-302, 0-indexed = 296-301)
FINN_IDX = [d.timetuple().tm_yday - 1 for d in DATES]   # 0-based

# Convert molecules/cm² /s → g CO2 m⁻² day⁻¹
# 1 molecule CO2 = 44 / 6.022e23 g
# 1 cm² = 1e-4 m²  →  1 mol/cm²/s = 1e4 mol/m²/s
# flux [g/m²/day] = flux [molec/cm²/s] × 86400 s/day × 44 g/mol / 6.022e23 × 1e4
MOLEC_CM2_S_TO_G_M2_DAY = 86400 * 44 / 6.022e23 * 1e4   # CO2
MOLEC_CM2_S_TO_G_M2_DAY_CO  = 86400 * 28 / 6.022e23 * 1e4
MOLEC_CM2_S_TO_G_M2_DAY_CH4 = 86400 * 16 / 6.022e23 * 1e4

print('Loading FINN emissions ...')
finn_lon, finn_lat, finn_co2 = load_finn(
    FINN_CO2, 'fire_modisviirs_CO2', FINN_IDX)
finn_co2_g = finn_co2 * MOLEC_CM2_S_TO_G_M2_DAY    # g CO2 m⁻² day⁻¹

_, _, finn_co  = load_finn(FINN_CO,  'fire_modisviirs_CO',  FINN_IDX)
finn_co_g  = finn_co  * MOLEC_CM2_S_TO_G_M2_DAY_CO
_, _, finn_ch4 = load_finn(FINN_CH4, 'fire_modisviirs_CH4', FINN_IDX)
finn_ch4_g = finn_ch4 * MOLEC_CM2_S_TO_G_M2_DAY_CH4


# ── Load WRF concentrations at surface (level 0), noon UTC ───────────────────
def load_wrf_noon(date):
    """Return (xlat, xlon, CO2_BBU, CO_BBU, CH4_BBU) at surface noon."""
    fname = os.path.join(
        WRF_DIR,
        f'wrfout_d01_{date.strftime("%Y-%m-%d")}_01:00:00'
    )
    if not os.path.exists(fname):
        return (None,) * 5
    ds = nc.Dataset(fname)
    # Times: 24 frames from 01:00 to 00:00 next day → index 11 = 12:00
    t_idx = 11   # 01 + 11h = 12:00 UTC
    xlat  = ds.variables['XLAT'][t_idx, :, :]
    xlon  = ds.variables['XLONG'][t_idx, :, :]
    co2   = ds.variables['CO2_BBU'][t_idx, 0, :, :]   # lev 0 = surface
    co    = ds.variables['CO_BBU' ][t_idx, 0, :, :]
    ch4   = ds.variables['CH4_BBU'][t_idx, 0, :, :]
    ds.close()
    return xlat, xlon, co2, co, ch4


print('Loading WRF surface concentrations ...')
wrf_data = {}
for dt in DATES:
    wrf_data[dt.day] = load_wrf_noon(dt)
# Get common XLAT/XLONG from first valid entry
xlat_ref = xlon_ref = None
for v in wrf_data.values():
    if v[0] is not None:
        xlat_ref, xlon_ref = v[0], v[1]
        break


# ── Load mean STILT footprint ─────────────────────────────────────────────────
def load_mean_footprint():
    files = glob.glob(os.path.join(FOOT_DIR, '*_foot.nc'))
    foot_sum = None
    count = 0
    for f in files:
        ds = nc.Dataset(f)
        foot = ds.variables['foot'][0, :, :]
        lon  = ds.variables['lon'][:]
        lat  = ds.variables['lat'][:]
        if foot_sum is None:
            foot_sum = np.zeros_like(foot, dtype=np.float64)
            flon, flat = lon, lat
        foot_sum += foot
        count += 1
        ds.close()
    return flon, flat, foot_sum / count if count else foot_sum

print('Loading STILT mean footprint ...')
foot_lon, foot_lat, foot_mean = load_mean_footprint()


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1 — FINN fire emission flux (CO2), 6-day panel
# ═══════════════════════════════════════════════════════════════════════════════
print('Plotting Fig 1: FINN emissions ...')

fig1, axes1 = plt.subplots(
    1, 6, figsize=(22, 4.5),
    subplot_kw={'projection': proj},
    gridspec_kw={'wspace': 0.04},
)

# Determine colour scale from nonzero values across all 6 days
nz = finn_co2_g.compressed()
nz = nz[nz > 0]
vmin_f = np.percentile(nz, 5)  if nz.size else 1e-3
vmax_f = np.percentile(nz, 99) if nz.size else 10.0
norm_f = LogNorm(vmin=max(vmin_f, 1e-4), vmax=vmax_f)

for ci, (dt, label) in enumerate(zip(DATES, DATE_LABELS)):
    ax = axes1[ci]
    add_map_features(ax)
    sparse_gridlines(ax)

    data_day = finn_co2_g[ci, :, :]
    lon2d, lat2d = np.meshgrid(finn_lon, finn_lat)
    pcm = ax.pcolormesh(
        lon2d, lat2d, data_day,
        cmap='hot_r', norm=norm_f,
        transform=proj, shading='auto', rasterized=True,
    )
    add_sites(ax, ms=5, color='dodgerblue')
    ax.set_title(label, fontsize=9, fontweight='bold')

    if ci == 0:
        ax.set_ylabel('Lat', fontsize=8)

cbar_ax1 = fig1.add_axes([0.15, 0.04, 0.70, 0.025])
sm1 = plt.cm.ScalarMappable(cmap='hot_r', norm=norm_f)
sm1.set_array([])
cb1 = fig1.colorbar(sm1, cax=cbar_ax1, orientation='horizontal', extend='both')
cb1.set_label('FINN CO₂ fire emission  [g CO₂ m⁻² day⁻¹]', fontsize=9)

fig1.suptitle(
    'FINN v2.5 daily biomass-burning CO₂ emission flux — Oct 24–29 2015',
    fontsize=11, y=1.01, fontweight='bold',
)
out1 = os.path.join(OUT_DIR, 'emission_finn.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight')
plt.close(fig1)
print(f'Saved: {out1}')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2 — WRF surface concentrations (3 species × 6 days)
# ═══════════════════════════════════════════════════════════════════════════════
print('Plotting Fig 2: WRF surface concentrations ...')

SPECIES = [
    ('CO2_BBU', 'CO₂ fire (ppmv)', 'OrRd'),
    ('CO_BBU',  'CO fire (ppmv)',   'YlOrBr'),
    ('CH4_BBU', 'CH₄ fire (ppmv)',  'PuRd'),
]
NROWS_C = 3
NCOLS_C = 6

fig2, axes2 = plt.subplots(
    NROWS_C, NCOLS_C,
    figsize=(NCOLS_C * 3.5, NROWS_C * 3.0),
    subplot_kw={'projection': proj},
    gridspec_kw={'hspace': 0.06, 'wspace': 0.04},
)

# Determine per-species colour limits from all days combined
sp_norms = {}
sp_keys  = ['co2', 'co', 'ch4']
for si, (sp_key, sp_label, sp_cmap) in enumerate(SPECIES):
    vals = []
    for dt in DATES:
        v = wrf_data[dt.day]
        if v[0] is None:
            continue
        arr = v[si + 2]   # co2=2, co=3, ch4=4
        nonz = arr[arr > 0]
        if nonz.size:
            vals.append(nonz)
    if vals:
        combined = np.concatenate(vals)
        vmin_c = 0.0
        vmax_c = np.percentile(combined, 99)
    else:
        vmin_c, vmax_c = 0, 1
    sp_norms[si] = mcolors.Normalize(vmin=vmin_c, vmax=vmax_c)

for ri, (sp_key, sp_label, sp_cmap) in enumerate(SPECIES):
    norm_c = sp_norms[ri]
    for ci, dt in enumerate(DATES):
        ax = axes2[ri, ci]
        add_map_features(ax)
        sparse_gridlines(ax)

        v = wrf_data[dt.day]
        if v[0] is not None:
            conc = v[ri + 2]
            pcm = ax.pcolormesh(
                xlon_ref, xlat_ref, conc,
                cmap=sp_cmap, norm=norm_c,
                transform=proj, shading='auto', rasterized=True,
            )
        add_sites(ax, ms=5, color='blue')

        if ri == 0:
            ax.set_title(DATE_LABELS[ci], fontsize=8, fontweight='bold')
        if ci == 0:
            ax.set_ylabel(sp_label, fontsize=8, labelpad=4)

    # Per-row colorbar on right
    cax_r = fig2.add_axes([0.92, 0.69 - ri * 0.32, 0.008, 0.26])
    sm_r  = plt.cm.ScalarMappable(cmap=sp_cmap, norm=norm_c)
    sm_r.set_array([])
    cb_r  = fig2.colorbar(sm_r, cax=cax_r, orientation='vertical', extend='max')
    cb_r.set_label(sp_label, fontsize=8, rotation=270, labelpad=14)

fig2.suptitle(
    'WRF-GHG surface biomass-burning tracer concentrations (level 0, 12 UTC)'
    ' — Oct 24–29 2015',
    fontsize=11, y=1.002, fontweight='bold',
)
out2 = os.path.join(OUT_DIR, 'concentration_wrf.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {out2}')


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Mean FINN CO2 emission overlaid with STILT mean footprint
# ═══════════════════════════════════════════════════════════════════════════════
print('Plotting Fig 3: Emission vs footprint overlay ...')

# Mean emission over the 6 days
finn_co2_mean = finn_co2_g.mean(axis=0)
finn_co2_mean_masked = np.ma.masked_less_equal(finn_co2_mean, 0)
lon2d, lat2d = np.meshgrid(finn_lon, finn_lat)

# Mean WRF CO2_BBU surface concentration
wrf_co2_list = [wrf_data[dt.day][2] for dt in DATES if wrf_data[dt.day][0] is not None]
wrf_co2_mean = np.mean(wrf_co2_list, axis=0) if wrf_co2_list else None

fig3, axes3 = plt.subplots(
    1, 2, figsize=(16, 6),
    subplot_kw={'projection': proj},
    gridspec_kw={'wspace': 0.1},
)

# ── Left panel: FINN emission + footprint contour ────────────────────────────
ax3a = axes3[0]
add_map_features(ax3a, linewidth=0.5)

nz_m = finn_co2_mean_masked.compressed()
vmax_m = np.percentile(nz_m, 99) if nz_m.size else 1
pcm3a = ax3a.pcolormesh(
    lon2d, lat2d, finn_co2_mean_masked,
    cmap='hot_r',
    norm=LogNorm(vmin=max(nz_m.min(), 1e-4) if nz_m.size else 1e-4, vmax=vmax_m),
    transform=proj, shading='auto', rasterized=True,
)

# STILT footprint as contours
foot_log = np.where(foot_mean > 0, np.log10(foot_mean), np.nan)
foot_lon2d, foot_lat2d = np.meshgrid(foot_lon, foot_lat)
levels = [-3.5, -3.0, -2.5, -2.0, -1.5]
cs = ax3a.contour(
    foot_lon2d, foot_lat2d, foot_log,
    levels=levels, colors='cyan', linewidths=[0.6, 0.8, 1.0, 1.2, 1.5],
    transform=proj, zorder=5,
)
ax3a.clabel(cs, fmt='%.1f', fontsize=7, inline=True)

add_sites(ax3a, ms=8, color='lime')

gl3a = ax3a.gridlines(draw_labels=True, linewidth=0.4,
                      color='gray', alpha=0.5, linestyle='--')
gl3a.xlocator = mticker.MultipleLocator(10)
gl3a.ylocator = mticker.MultipleLocator(5)
gl3a.top_labels = gl3a.right_labels = False
gl3a.xlabel_style = gl3a.ylabel_style = {'size': 8}

cb3a = fig3.colorbar(pcm3a, ax=ax3a, orientation='vertical',
                     fraction=0.025, pad=0.04, shrink=0.9, extend='both')
cb3a.set_label('FINN CO₂ mean emission  [g m⁻² day⁻¹]', fontsize=9)
ax3a.set_title('Fire emission flux + STILT footprint (cyan contours, log₁₀)',
               fontsize=10, fontweight='bold')

# ── Right panel: WRF CO2_BBU mean + footprint contour ────────────────────────
ax3b = axes3[1]
add_map_features(ax3b, linewidth=0.5)

if wrf_co2_mean is not None:
    nz_w = wrf_co2_mean[wrf_co2_mean > 0]
    vmax_w = np.percentile(nz_w, 99) if nz_w.size else 1
    pcm3b = ax3b.pcolormesh(
        xlon_ref, xlat_ref, wrf_co2_mean,
        cmap='OrRd',
        norm=mcolors.Normalize(vmin=0, vmax=vmax_w),
        transform=proj, shading='auto', rasterized=True,
    )
    cb3b = fig3.colorbar(pcm3b, ax=ax3b, orientation='vertical',
                         fraction=0.025, pad=0.04, shrink=0.9, extend='max')
    cb3b.set_label('WRF CO₂ biomass-burning mean conc. (ppmv)', fontsize=9)

# Overlay footprint contours
ax3b.contour(
    foot_lon2d, foot_lat2d, foot_log,
    levels=levels, colors='blue', linewidths=[0.6, 0.8, 1.0, 1.2, 1.5],
    transform=proj, zorder=5,
)
add_sites(ax3b, ms=8, color='blue')

gl3b = ax3b.gridlines(draw_labels=True, linewidth=0.4,
                      color='gray', alpha=0.5, linestyle='--')
gl3b.xlocator = mticker.MultipleLocator(10)
gl3b.ylocator = mticker.MultipleLocator(5)
gl3b.top_labels = gl3b.right_labels = False
gl3b.xlabel_style = gl3b.ylabel_style = {'size': 8}

# Label receptor sites
for sname, (slon, slat) in SITES.items():
    ax3b.text(slon + 0.6, slat + 0.6, sname, fontsize=7, transform=proj,
              zorder=6, bbox=dict(boxstyle='round,pad=0.2', fc='white',
                                  alpha=0.75, ec='none'))

ax3b.set_title('WRF CO₂ fire tracer + STILT footprint (blue contours, log₁₀)',
               fontsize=10, fontweight='bold')

fig3.suptitle(
    'Fire source vs. atmospheric sensitivity — Oct 24–29 2015 mean',
    fontsize=12, y=1.01, fontweight='bold',
)
out3 = os.path.join(OUT_DIR, 'emission_vs_footprint.png')
fig3.savefig(out3, dpi=150, bbox_inches='tight')
plt.close(fig3)
print(f'Saved: {out3}')
print('All done.')
