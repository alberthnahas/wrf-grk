#!/usr/bin/env python3
"""STILT footprint analysis plots for v3.2 (06 UTC release time).

Creates two diagnostic figures:
  1. footprint_examples_v3.png — Per-site mean footprint over Oct 6-31, log-scale
     map with coastlines.
  2. footprint_diagnostics_v3.png — (a) Time series of footprint integral per
     site, (b) per-site footprint hotspot count and max value table.
"""
import os, glob
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT     = '/home/igrk/WRF-GRK'
FOOT_DIR = f'{ROOT}/stilt/out_06utc/footprints'
PLT_DIR  = f'{ROOT}/inversion/plots'

SITES = {
    'Bukit Kototabang': (100.32, -0.20),
    'Jambi':            (103.61, -1.61),
    'Jakarta':          (106.85, -6.21),
    'Pontianak':        (109.33, -0.02),
    'Palangka Raya':    (113.94, -2.16),
    'Pekanbaru':        (101.45,  0.53),
    'Palembang':        (104.75, -2.99),
    'Banjarmasin':      (114.59, -3.32),
    'Samarinda':        (117.15, -0.50),
    'Balikpapan':       (116.85, -1.27),
    'Makassar':         (119.42, -5.14),
}
SITE_ORDER = list(SITES.keys())

def load_foot(slon, slat, day):
    tag = f'201510{day:02d}0600'
    pattern = os.path.join(FOOT_DIR, f'{tag}_*_foot.nc')
    for f in glob.glob(pattern):
        parts = os.path.basename(f).split('_')
        try:
            if (abs(float(parts[1]) - slon) < 0.02 and
                    abs(float(parts[2]) - slat) < 0.02):
                with nc.Dataset(f) as ds:
                    return (ds['foot'][0].data, ds['lon'][:].data,
                            ds['lat'][:].data)
        except (ValueError, IndexError):
            continue
    return None

# ── Aggregate: mean footprint per site over all 26 days ────────────────────
print('Loading footprints ...')
foot_means = {}
foot_lon = foot_lat = None
foot_integrals = {s: [] for s in SITE_ORDER}
for sname in SITE_ORDER:
    slon, slat = SITES[sname]
    stack = []
    for day in range(6, 32):
        r = load_foot(slon, slat, day)
        if r is None:
            continue
        f, lon, lat = r
        if foot_lon is None:
            foot_lon, foot_lat = lon, lat
        stack.append(f)
        foot_integrals[sname].append((day, float(f.sum())))
    foot_means[sname] = np.mean(stack, axis=0) if stack else None
    print(f'  {sname:<22}: {len(stack)} footprints, '
          f'mean integral = {np.mean([x[1] for x in foot_integrals[sname]]):.2f} '
          f'ppm·s·m²/µmol')

# ─────────────────────────────────────────────────────────────────────────
# Figure 1 — Per-site mean footprint maps (3×4 grid)
# ─────────────────────────────────────────────────────────────────────────
proj = ccrs.PlateCarree()
fig, axes = plt.subplots(3, 4, figsize=(20, 11),
                         subplot_kw={'projection': proj})
axes = axes.ravel()

vmin, vmax = 1e-4, 1e1
for i, sname in enumerate(SITE_ORDER):
    ax = axes[i]
    ax.add_feature(cfeature.OCEAN.with_scale('50m'),
                   facecolor='#dceaf5', zorder=0)
    ax.add_feature(cfeature.LAND.with_scale('50m'),
                   facecolor='#f4ede0', zorder=0)
    fmean = foot_means.get(sname)
    if fmean is not None and fmean.max() > 0:
        plot_data = np.where(fmean > 0, fmean, np.nan)
        im = ax.pcolormesh(foot_lon, foot_lat, plot_data,
                           norm=LogNorm(vmin=vmin, vmax=vmax),
                           cmap='magma_r', shading='auto',
                           transform=proj, zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'),
                   linewidth=0.4, edgecolor='black', zorder=3)

    slon, slat = SITES[sname]
    ax.plot(slon, slat, marker='*', color='cyan',
            markersize=14, markeredgecolor='black',
            markeredgewidth=1, transform=proj, zorder=5)
    ax.set_extent([90, 145, -12, 8], crs=proj)
    ax.set_title(f'{sname}\n({slon:.2f}°E, {slat:.2f}°N)', fontsize=10)
    gl = ax.gridlines(draw_labels=False, alpha=0.2, linewidth=0.3)

# Hide last unused subplot
for j in range(len(SITE_ORDER), len(axes)):
    axes[j].set_axis_off()

cbar = fig.colorbar(im, ax=axes.tolist(), orientation='horizontal',
                    pad=0.04, shrink=0.5, aspect=40)
cbar.set_label('Mean footprint sensitivity  (ppm / (µmol m⁻² s⁻¹), log scale)\n'
               'Oct 6–31 2015 mean, release at 06 UTC (= 13 WIB) — afternoon mixed PBL',
               fontsize=10)
fig.suptitle('STILT v3.2 — Mean backward footprints, 11-site Indonesian network',
             fontsize=14, fontweight='bold', y=0.98)
out1 = os.path.join(PLT_DIR, 'stilt_footprint_means_v3.png')
plt.savefig(out1, dpi=130, bbox_inches='tight')
plt.close()
print(f'Saved: {out1}')

# ─────────────────────────────────────────────────────────────────────────
# Figure 2 — Footprint integral time series per site
# ─────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

ax = axes[0]
colors = plt.cm.tab20(np.linspace(0, 1, len(SITE_ORDER)))
for i, sname in enumerate(SITE_ORDER):
    days = [d for d, _ in foot_integrals[sname]]
    vals = [v for _, v in foot_integrals[sname]]
    ax.plot(days, vals, '-o', color=colors[i], label=sname,
            ms=4, lw=1.2)
ax.set_xlabel('Date in October 2015')
ax.set_ylabel('∫ footprint dA  (ppm·s·m²/µmol  ≈  total upwind influence)')
ax.set_title('Footprint integral by receptor and day\n'
             '(higher = more land-surface coverage in the past 5 days)')
ax.legend(fontsize=7, ncol=2, loc='upper right')
ax.grid(alpha=0.3)

# Right panel: bar chart of mean integral + fire-zone H_max per site
ax2 = axes[1]
mean_int = [np.mean([v for _, v in foot_integrals[s]]) for s in SITE_ORDER]
ypos = np.arange(len(SITE_ORDER))
bars = ax2.barh(ypos, mean_int, color=colors)
ax2.set_yticks(ypos)
ax2.set_yticklabels(SITE_ORDER)
ax2.invert_yaxis()
ax2.set_xlabel('Mean footprint integral over Oct 6–31  (ppm·s·m²/µmol)')
ax2.set_title('Mean per-site total upwind influence\n'
              '(sites listed in receptor order)')
ax2.grid(alpha=0.3, axis='x')
for b, v in zip(bars, mean_int):
    ax2.text(v + 0.05, b.get_y() + b.get_height()/2,
             f'{v:.2f}', va='center', fontsize=8)

plt.tight_layout()
out2 = os.path.join(PLT_DIR, 'stilt_footprint_integrals_v3.png')
plt.savefig(out2, dpi=140, bbox_inches='tight')
plt.close()
print(f'Saved: {out2}')

# ─────────────────────────────────────────────────────────────────────────
# Figure 3 — BKT-day-by-day footprint comparison
# Highlights what each NOAA flask sample is sensitive to.
# ─────────────────────────────────────────────────────────────────────────
flask_days = [6, 13, 20, 27]
fig, axes = plt.subplots(1, 4, figsize=(22, 6),
                         subplot_kw={'projection': proj})
slon, slat = SITES['Bukit Kototabang']
for i, day in enumerate(flask_days):
    ax = axes[i]
    ax.add_feature(cfeature.OCEAN.with_scale('50m'),
                   facecolor='#dceaf5', zorder=0)
    ax.add_feature(cfeature.LAND.with_scale('50m'),
                   facecolor='#f4ede0', zorder=0)
    r = load_foot(slon, slat, day)
    if r is not None:
        f, lon, lat = r
        plot_data = np.where(f > 0, f, np.nan)
        im = ax.pcolormesh(lon, lat, plot_data,
                           norm=LogNorm(vmin=vmin, vmax=vmax),
                           cmap='magma_r', shading='auto',
                           transform=proj, zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'),
                   linewidth=0.5, edgecolor='black', zorder=3)
    ax.plot(slon, slat, marker='*', color='cyan',
            markersize=18, markeredgecolor='black',
            markeredgewidth=1.2, transform=proj, zorder=5)
    ax.set_extent([90, 130, -10, 8], crs=proj)
    enh_lookup = {6: '+5.24', 13: '+7.17', 20: '+5.60', 27: '+10.41'}
    qc = 'C..' if day == 27 else '...'
    ax.set_title(f'BKT — Oct {day}  06 UTC\n'
                 f'flask enh = {enh_lookup[day]} ppm  (qc={qc})',
                 fontsize=11)

cbar = fig.colorbar(im, ax=axes.tolist(), orientation='horizontal',
                    pad=0.06, shrink=0.6, aspect=45)
cbar.set_label('Footprint sensitivity (ppm / (µmol m⁻² s⁻¹), log scale)',
               fontsize=10)
fig.suptitle('BKT footprints on each NOAA flask sampling day\n'
             'Where each flask measurement "looks" upwind',
             fontsize=13, fontweight='bold', y=1.02)
out3 = os.path.join(PLT_DIR, 'stilt_bkt_flask_footprints_v3.png')
plt.savefig(out3, dpi=140, bbox_inches='tight')
plt.close()
print(f'Saved: {out3}')
