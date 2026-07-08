#!/usr/bin/env python3
"""Plot prior vs. posterior CO₂ emissions as region boxes on an Indonesia map."""
import os
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT      = '/home/igrk/WRF-GRK'
FINN_CO2  = (f'{ROOT}/rawdata/finn_fire/'
             'emissions-finnv2.5modvrs_CO2_bb_surface_daily_'
             '20150101-20151231_0.1x0.1.nc')
PLT_DIR   = f'{ROOT}/inversion/plots'

REGIONS = [
    ('Sumatra N+C',    95.0, 104.0,  0.0,  6.0),
    ('Sumatra S',      99.0, 107.0, -6.0,  0.0),
    ('Kalimantan W+C', 108.0, 112.0, -2.0,  3.0),
    ('Kalimantan S+E', 112.0, 119.0, -5.0,  1.0),
    ('Java',           105.0, 116.0, -9.0, -5.5),
    ('Sulawesi+East',  119.0, 141.0, -6.0,  2.0),
]

# v3.2 final: 06 UTC STILT footprints + 4 NOAA BKT flasks
PRIOR_TG = np.array([3.089, 126.416, 7.982, 215.506, 5.912, 112.668])
POST_TG  = np.array([0.866, 58.363, 1.952, 136.533, 1.095, 15.726])
ALPHA    = np.array([0.280, 0.462, 0.245, 0.634, 0.185, 0.140])

SITES = {
    'BKT':  (100.32, -0.20),  'Jambi': (103.61, -1.61),
    'Jakarta': (106.85, -6.21), 'Pontianak': (109.33, -0.02),
    'Palangka Raya': (113.94, -2.16), 'Pekanbaru': (101.45, 0.53),
    'Palembang': (104.75, -2.99), 'Banjarmasin': (114.59, -3.32),
    'Samarinda': (117.15, -0.50), 'Balikpapan': (116.85, -1.27),
    'Makassar': (119.42, -5.14),
}

# ── Load and aggregate FINN to a common 0.1° grid for background shading ───
print('Loading FINN CO₂ for background ...')
ds = nc.Dataset(FINN_CO2)
lon = ds.variables['lon'][:]
lat = ds.variables['lat'][:]
# Oct 6-31 = day-of-year 278..303 (0-indexed: 278-1=277..303-1=302)
day_idx = list(range(278, 304))
flux = ds.variables['fire_modisviirs_CO2'][day_idx, :, :].sum(axis=0)  # molec/cm²/s sum
ds.close()

# Restrict to domain
lon_mask = (lon >= 93) & (lon <= 145)
lat_mask = (lat >= -12) & (lat <= 8)
lon_d = lon[lon_mask]; lat_d = lat[lat_mask]
flux_d = flux[np.ix_(lat_mask, lon_mask)]

# Convert molec/cm²/s × 86400 → daily molec/cm² (then sum already done over days)
# For visualization, just use log of raw integrated flux (relative intensity)
flux_d = np.where(flux_d <= 0, np.nan, flux_d)

# Distribute prior/posterior emissions spatially using FINN as the spatial pattern
# For each region: spatial pattern × (region_total / sum(spatial_pattern in region))
prior_grid = np.zeros_like(flux_d)
post_grid  = np.zeros_like(flux_d)
LON2D, LAT2D = np.meshgrid(lon_d, lat_d)

for ki, (rname, x0, x1, y0, y1) in enumerate(REGIONS):
    rmask = (LON2D >= x0) & (LON2D <= x1) & (LAT2D >= y0) & (LAT2D <= y1)
    raw = np.where(np.isnan(flux_d), 0, flux_d) * rmask
    rsum = raw.sum()
    if rsum > 0:
        # Each cell gets a fraction of the regional Tg
        prior_grid[rmask] += (raw[rmask] / rsum) * PRIOR_TG[ki]
        post_grid[rmask]  += (raw[rmask] / rsum) * POST_TG[ki]

prior_grid = np.where(prior_grid <= 0, np.nan, prior_grid)
post_grid  = np.where(post_grid  <= 0, np.nan, post_grid)

# ── Plot ────────────────────────────────────────────────────────────────────
region_colors = ['#e6550d', '#fd8d3c', '#2171b5', '#6baed6', '#74c476', '#9e9ac8']

proj = ccrs.PlateCarree()
fig, axes = plt.subplots(1, 2, figsize=(18, 7),
                         subplot_kw={'projection': proj})

vmax = np.nanmax(prior_grid)
vmin = vmax / 1e4

for ax, grid, title, total in [
    (axes[0], prior_grid, f'FINN PRIOR — Total: {PRIOR_TG.sum():.1f} Tg CO₂', PRIOR_TG),
    (axes[1], post_grid,  f'POSTERIOR (v3.1) — Total: {POST_TG.sum():.1f} Tg CO₂', POST_TG),
]:
    # Land/ocean background
    ax.add_feature(cfeature.OCEAN.with_scale('50m'),
                   facecolor='#dceaf5', zorder=0)
    ax.add_feature(cfeature.LAND.with_scale('50m'),
                   facecolor='#f4ede0', zorder=0)

    im = ax.pcolormesh(lon_d, lat_d, grid, norm=LogNorm(vmin=vmin, vmax=vmax),
                       cmap='YlOrRd', shading='auto', transform=proj, zorder=2)

    # Coastlines and country borders on top of the data
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'),
                   linewidth=0.6, edgecolor='black', zorder=3)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'),
                   linewidth=0.4, edgecolor='gray', linestyle=':', zorder=3)

    # Region boxes with labels
    for ki, (rname, x0, x1, y0, y1) in enumerate(REGIONS):
        rect = mpatches.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                   linewidth=2.2, edgecolor=region_colors[ki],
                                   facecolor='none', linestyle='-',
                                   transform=proj, zorder=4)
        ax.add_patch(rect)
        tg = total[ki]
        label = f'{rname}\n{tg:.1f} Tg'
        ax.text(x1 - 0.3, y1 - 0.3, label,
                ha='right', va='top', fontsize=8, fontweight='bold',
                color=region_colors[ki], transform=proj, zorder=5,
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                          edgecolor=region_colors[ki], alpha=0.9))

    # Receptor sites
    for sname, (slon, slat) in SITES.items():
        ax.plot(slon, slat, marker='^', color='black', markersize=7,
                markerfacecolor='yellow', markeredgewidth=1,
                transform=proj, zorder=6)
        ax.annotate(sname, (slon, slat), xytext=(3, 3),
                    textcoords='offset points', fontsize=6.5, zorder=6)

    ax.set_extent([93, 145, -12, 8], crs=proj)
    ax.set_title(title, fontsize=12, fontweight='bold')
    gl = ax.gridlines(draw_labels=True, alpha=0.3, linestyle=':',
                      linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

# Single colorbar
cbar = fig.colorbar(im, ax=axes, orientation='horizontal',
                    pad=0.08, shrink=0.55, aspect=40)
cbar.set_label('CO₂ emission per 0.1° cell  (Tg, Oct 6–31 2015, log scale)')

fig.suptitle('Indonesia 2015 Fire CO₂ Emissions — Spatial Distribution by Region\n'
             'Prior (FINN) vs. Posterior (v3.1 Bayesian inversion, 11-site STILT + OCO-2)',
             fontsize=13, fontweight='bold', y=1.02)

out = os.path.join(PLT_DIR, 'spatial_prior_vs_posterior_v3.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f'Saved: {out}')
