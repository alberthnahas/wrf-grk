#!/usr/bin/env python3
"""
Plot STILT backward-trajectory footprints for the IDN_BB_2015 simulation.

Produces two figures:
  1. stilt/out/footprint_grid.png  — 5-column × 6-row panel (site × day)
  2. stilt/out/footprint_mean.png  — single map of mean footprint across all
     30 receptor-times

Footprint units: ppm (µmol⁻¹ m² s), i.e. atmospheric influence per unit
surface flux.  Values are log10-scaled for display (zeros masked).
"""

import os
import glob
import re
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc

# ── Configuration ────────────────────────────────────────────────────────────
FOOT_DIR  = '/home/igrk/WRF-GRK/stilt/out/footprints'
OUT_DIR   = '/home/igrk/WRF-GRK/stilt/out'
os.makedirs(OUT_DIR, exist_ok=True)

# Receptor sites: label, lon, lat
SITES = {
    'Bukit Kototabang': (100.32, -0.20),
    'Palangka Raya':    (113.94, -2.16),
    'Pontianak':        (109.33, -0.02),
    'Jambi':            (103.61, -1.61),
    'Jakarta':          (106.85, -6.21),
}
SITE_ORDER = list(SITES.keys())

# Receptor dates (Oct 24–29 2015 noon UTC)
DATES = [datetime(2015, 10, d, 12, tzinfo=timezone.utc) for d in range(24, 30)]

# Map extent (WRF domain)
LON_MIN, LON_MAX = 90.0, 145.0
LAT_MIN, LAT_MAX = -15.0, 15.0

# Colormap: white-yellow-orange-red for influence
CMAP = plt.cm.YlOrRd
NORM_VMIN, NORM_VMAX = -4, 0   # log10 range (10⁻⁴ … 10⁰ ppm µmol⁻¹ m² s)


# ── Helper: load footprint ───────────────────────────────────────────────────
def find_footprint(site_lon, site_lat, dt):
    """Return (lon1d, lat1d, foot2d) for a given site and datetime."""
    tag = dt.strftime('%Y%m%d%H%M')
    pattern = os.path.join(FOOT_DIR, f'{tag}_{site_lon}_*{site_lat}*_foot.nc')
    files = glob.glob(pattern)
    if not files:
        # Try both sign variants and floating-point repr differences
        pattern2 = os.path.join(FOOT_DIR, f'{tag}_*_foot.nc')
        for f in glob.glob(pattern2):
            base = os.path.basename(f)
            # base: 201510241200_100.32_-0.2_10_foot.nc
            parts = base.split('_')
            try:
                flon = float(parts[1])
                flat = float(parts[2])
                if abs(flon - site_lon) < 0.01 and abs(flat - site_lat) < 0.01:
                    files = [f]
                    break
            except (IndexError, ValueError):
                continue
    if not files:
        return None, None, None
    ds = nc.Dataset(files[0])
    lon  = ds['lon'][:]
    lat  = ds['lat'][:]
    foot = ds['foot'][0, :, :]   # (lat, lon)
    ds.close()
    return lon, lat, foot


def log_foot(foot):
    """Mask zeros and return log10 array."""
    f = np.where(foot > 0, foot, np.nan)
    return np.log10(f)


# ── Figure 1: 6-row × 5-col panel (date × site) ─────────────────────────────
print('Building panel figure ...')

NROWS = len(DATES)     # 6
NCOLS = len(SITE_ORDER)  # 5

proj = ccrs.PlateCarree()
fig, axes = plt.subplots(
    NROWS, NCOLS,
    figsize=(NCOLS * 4, NROWS * 3.2),
    subplot_kw={'projection': proj},
    gridspec_kw={'hspace': 0.08, 'wspace': 0.05},
)

# Collect all log-foot arrays to set a common colour scale
all_lf = []
data_cache = {}
for ri, dt in enumerate(DATES):
    for ci, sname in enumerate(SITE_ORDER):
        slon, slat = SITES[sname]
        lon, lat, foot = find_footprint(slon, slat, dt)
        data_cache[(ri, ci)] = (lon, lat, foot)
        if foot is not None:
            lf = log_foot(foot)
            all_lf.append(lf[np.isfinite(lf)])

if all_lf:
    combined = np.concatenate(all_lf)
    vmin = max(NORM_VMIN, np.percentile(combined, 2))
    vmax = min(NORM_VMAX, np.percentile(combined, 99))
else:
    vmin, vmax = NORM_VMIN, NORM_VMAX

norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

for ri, dt in enumerate(DATES):
    for ci, sname in enumerate(SITE_ORDER):
        ax = axes[ri, ci]
        slon, slat = SITES[sname]
        lon, lat, foot = data_cache[(ri, ci)]

        ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=proj)
        ax.add_feature(cfeature.LAND,        facecolor='#f5f5f0', linewidth=0)
        ax.add_feature(cfeature.OCEAN,       facecolor='#d0e8f5', linewidth=0)
        ax.add_feature(cfeature.COASTLINE,   linewidth=0.4, edgecolor='#444')
        ax.add_feature(cfeature.BORDERS,     linewidth=0.3, edgecolor='#888',
                       linestyle=':')

        if foot is not None:
            lf = log_foot(foot)
            pcm = ax.pcolormesh(
                lon, lat, lf,
                cmap=CMAP, norm=norm,
                transform=proj, shading='auto',
                rasterized=True,
            )

        # Receptor marker
        ax.plot(slon, slat, marker='*', color='blue', markersize=6,
                transform=proj, zorder=5)

        # Gridlines — only outer edges labelled
        gl = ax.gridlines(draw_labels=False, linewidth=0.3,
                          color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = mticker.MultipleLocator(15)
        gl.ylocator = mticker.MultipleLocator(10)

        # Column header (site name) on top row
        if ri == 0:
            ax.set_title(sname, fontsize=8, pad=4, fontweight='bold')

        # Row label (date) on left column
        if ci == 0:
            ax.set_ylabel(dt.strftime('%b %-d'), fontsize=8, labelpad=4)

# Shared colorbar
cbar_ax = fig.add_axes([0.15, 0.02, 0.70, 0.015])
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cb.set_label('log₁₀ footprint  [ppm µmol⁻¹ m² s]', fontsize=9)

fig.suptitle(
    'WRF-STILT 5-day backward footprints — Indonesia fire season Oct 2015',
    fontsize=11, y=0.995, fontweight='bold',
)

out_panel = os.path.join(OUT_DIR, 'footprint_grid.png')
fig.savefig(out_panel, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out_panel}')


# ── Figure 2: mean footprint across all 30 simulations ───────────────────────
print('Building mean footprint map ...')

# Accumulate on common grid (from first file)
lon_ref = lat_ref = foot_sum = count = None
for (ri, ci), (lon, lat, foot) in data_cache.items():
    if foot is None:
        continue
    if foot_sum is None:
        lon_ref   = lon
        lat_ref   = lat
        foot_sum  = np.zeros_like(foot, dtype=np.float64)
        count     = 0
    foot_sum += foot
    count    += 1

foot_mean = foot_sum / count if count else foot_sum

fig2, ax2 = plt.subplots(
    1, 1, figsize=(10, 6),
    subplot_kw={'projection': proj},
)
ax2.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=proj)
ax2.add_feature(cfeature.LAND,      facecolor='#f5f5f0', linewidth=0)
ax2.add_feature(cfeature.OCEAN,     facecolor='#d0e8f5', linewidth=0)
ax2.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='#444')
ax2.add_feature(cfeature.BORDERS,   linewidth=0.4, edgecolor='#888',
                linestyle=':')

lf_mean = log_foot(foot_mean)
pcm2 = ax2.pcolormesh(
    lon_ref, lat_ref, lf_mean,
    cmap=CMAP, norm=norm,
    transform=proj, shading='auto', rasterized=True,
)

# Mark all receptor sites
for sname, (slon, slat) in SITES.items():
    ax2.plot(slon, slat, marker='*', color='blue', markersize=9,
             transform=proj, zorder=5)
    ax2.text(slon + 0.5, slat + 0.5, sname, fontsize=7,
             transform=proj, zorder=5,
             bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7, ec='none'))

gl2 = ax2.gridlines(draw_labels=True, linewidth=0.4,
                    color='gray', alpha=0.5, linestyle='--')
gl2.xlocator = mticker.MultipleLocator(10)
gl2.ylocator = mticker.MultipleLocator(5)
gl2.top_labels = False
gl2.right_labels = False
gl2.xlabel_style = {'size': 8}
gl2.ylabel_style = {'size': 8}

cb2 = fig2.colorbar(pcm2, ax=ax2, orientation='vertical',
                    fraction=0.025, pad=0.04, shrink=0.85)
cb2.set_label('log₁₀ footprint  [ppm µmol⁻¹ m² s]', fontsize=9)

ax2.set_title(
    f'Mean WRF-STILT footprint — 5 sites × Oct 24–29 2015 (n={count})',
    fontsize=11, fontweight='bold',
)

out_mean = os.path.join(OUT_DIR, 'footprint_mean.png')
fig2.savefig(out_mean, dpi=150, bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {out_mean}')
print('Done.')
