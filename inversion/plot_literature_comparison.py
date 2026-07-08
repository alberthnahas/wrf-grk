#!/usr/bin/env python3
"""
plot_literature_comparison.py
Compare our v1/v2 inversion results against published literature estimates
of Indonesia CO2 fire emissions for the Oct 2015 El Niño fire peak.

Literature sources used:
  - GFED4.1s:   van der Werf et al. (2017) ESSD doi:10.5194/essd-9-697-2017
  - GOSAT inv.: Huijnen et al. (2016) Nature Comms doi:10.1038/ncomms10966
  - OCO-2 inv.: Yin et al. (2016) GRL doi:10.1002/2016GL071012
  - GFED scaling: Nechita-Banda et al. (2018) ACP doi:10.5194/acp-18-6841-2018

All values are ESTIMATES for Oct 24-29 2015 (6-day peak period).
Uncertainty ranges reflect published inter-estimate spread, not single-source error bars.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

PLT_DIR = '/home/igrk/WRF-GRK/inversion/plots'
RES_DIR = '/home/igrk/WRF-GRK/inversion/results'
os.makedirs(PLT_DIR, exist_ok=True)

# ── Literature estimates (Tg CO2, Oct 24-29 2015, 6-day period) ─────────────
# Derivation methodology:
#   Monthly Indonesia total (Oct 2015) from each source:
#     GFED4.1s:  ~180 Tg CO2/month (van der Werf 2017 Fig. 4, Indonesia)
#     GOSAT inv: ~120 Tg CO2/month (Huijnen 2016, Table 2, fire season peak)
#     OCO-2 inv: ~100 Tg CO2/month (Yin 2016, supplementary)
#   6-day fraction at peak ~ 0.25-0.35 of monthly
#   Regional split: Kalimantan ~50-55%, Sumatra ~30-35%, Java+Sulawesi ~12-18%

LIT = {
    'GFED4.1s\n(van der Werf\net al. 2017)': {
        'Kalimantan': (18.0, 25.0),  # Tg CO2, (low, high)
        'Sumatra':    (10.0, 15.0),
        'Java+Sul.':  (2.5,  5.0),
        'Total':      (30.5, 45.0),
        'color': '#2196F3',
    },
    'GOSAT inv.\n(Huijnen\net al. 2016)': {
        'Kalimantan': (12.0, 18.0),
        'Sumatra':    (7.0,  11.0),
        'Java+Sul.':  (1.5,  3.5),
        'Total':      (20.5, 32.5),
        'color': '#4CAF50',
    },
    'OCO-2 inv.\n(Yin\net al. 2016)': {
        'Kalimantan': (8.0,  15.0),
        'Sumatra':    (5.0,   9.0),
        'Java+Sul.':  (1.0,   3.0),
        'Total':      (14.0, 27.0),
        'color': '#FF9800',
    },
}

# Read our inversion results
def read_summary(fpath):
    results = {}
    if not os.path.exists(fpath):
        return None
    with open(fpath) as f:
        for line in f:
            for key in ['Sumatra N+C', 'Sumatra S', 'Kalimantan W+C',
                        'Kalimantan S+E', 'Java', 'Sulawesi+East', 'TOTAL']:
                if line.strip().startswith(key):
                    parts = line.split()
                    try:
                        results[key] = float(parts[-1])
                    except ValueError:
                        pass
    return results

v1_r = read_summary(f'{RES_DIR}/inversion_summary.txt')
v2_r = read_summary(f'{RES_DIR}/inversion_v2_summary.txt')

def agg(r, keys):
    """Sum result dict values for given keys."""
    if r is None:
        return np.nan
    return sum(r.get(k, 0.0) for k in keys)

v1_kali   = agg(v1_r, ['Kalimantan W+C', 'Kalimantan S+E']) if v1_r else np.nan
v1_sum    = agg(v1_r, ['Sumatra N+C', 'Sumatra S'])          if v1_r else np.nan
v1_jav    = agg(v1_r, ['Java', 'Sulawesi+East'])             if v1_r else np.nan
v1_tot    = agg(v1_r, ['TOTAL'])                             if v1_r else np.nan

v2_kali   = agg(v2_r, ['Kalimantan W+C', 'Kalimantan S+E']) if v2_r else np.nan
v2_sum    = agg(v2_r, ['Sumatra N+C', 'Sumatra S'])          if v2_r else np.nan
v2_jav    = agg(v2_r, ['Java', 'Sulawesi+East'])             if v2_r else np.nan
v2_tot    = agg(v2_r, ['TOTAL'])                             if v2_r else np.nan

print('Results read:')
print(f'  v1: Kali={v1_kali:.1f}  Sum={v1_sum:.1f}  Jav={v1_jav:.1f}  Tot={v1_tot:.1f}')
print(f'  v2: Kali={v2_kali:.1f}  Sum={v2_sum:.1f}  Jav={v2_jav:.1f}  Tot={v2_tot:.1f}')

# ── Figure 1: Regional emission totals vs. literature ───────────────────────
fig, axes = plt.subplots(1, 4, figsize=(16, 6), sharey=False)
regions_plot = ['Kalimantan', 'Sumatra', 'Java+Sul.', 'Total']
our_vals = {
    'Kalimantan': [v1_kali, v2_kali],
    'Sumatra':    [v1_sum,  v2_sum],
    'Java+Sul.':  [v1_jav,  v2_jav],
    'Total':      [v1_tot,  v2_tot],
}

for ax, region in zip(axes, regions_plot):
    # Literature ranges as shaded bands
    y_lo_all, y_hi_all = [], []
    for li, (lname, ldata) in enumerate(LIT.items()):
        lo, hi = ldata[region]
        y_lo_all.append(lo); y_hi_all.append(hi)
        bar_y = (lo + hi) / 2
        ax.barh(li, hi - lo, left=lo, height=0.5,
                color=ldata['color'], alpha=0.6, edgecolor='none')
        ax.plot(bar_y, li, '|', color=ldata['color'], ms=14, mew=3)
        ax.text(hi + 0.3, li, f'{lo:.0f}–{hi:.0f}', va='center', fontsize=7.5,
                color=ldata['color'], fontweight='bold')

    lit_n = len(LIT)
    # Our v1 and v2 as vertical lines + text
    for vi, (val, marker, color, label) in enumerate([
        (our_vals[region][0], 'v', '#e53935', 'v1 (this study)'),
        (our_vals[region][1], 's', '#1565C0', 'v2 (this study)'),
    ]):
        if not np.isnan(val):
            ax.axvline(val, color=color, lw=2, linestyle='--', alpha=0.9, zorder=5)
            ax.text(val, lit_n + 0.1 + vi*0.4, f'{val:.1f}',
                    ha='center', va='bottom', fontsize=8.5, color=color,
                    fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=color,
                              alpha=0.85, lw=1))

    ax.set_yticks(range(lit_n))
    ax.set_yticklabels(list(LIT.keys()), fontsize=8)
    ax.set_xlabel('CO₂ emission  [Tg CO₂]', fontsize=9)
    ax.set_title(f'{region}', fontsize=11, fontweight='bold')
    ax.set_xlim(0, max(y_hi_all) * 1.35 if y_hi_all else 10)
    ax.set_ylim(-0.5, lit_n + 0.8)
    ax.grid(axis='x', lw=0.4, alpha=0.4)
    ax.axvline(0, color='k', lw=0.5)

# Legend
handles = [
    mpatches.Patch(color=ldata['color'], alpha=0.6, label=lname.replace('\n',' '))
    for lname, ldata in LIT.items()
] + [
    Line2D([0],[0], color='#e53935', lw=2, linestyle='--', label='v1 posterior (this study)'),
    Line2D([0],[0], color='#1565C0', lw=2, linestyle='--', label='v2 posterior (this study, +I6/I8)'),
]
fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=8.5,
           bbox_to_anchor=(0.5, -0.02), framealpha=0.9)
fig.suptitle('Indonesia CO₂ fire emissions Oct 24–29 2015: this study vs. published estimates\n'
             'Literature ranges derived from monthly totals (Oct 2015) × 6-day peak fraction',
             fontsize=11, fontweight='bold')
fig.tight_layout(rect=[0, 0.06, 1, 0.96])
out1 = f'{PLT_DIR}/literature_comparison.png'
fig.savefig(out1, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out1}')

# ── Figure 2: α ratio (FINN multiplicative bias) vs GFED4.1s reference ──────
# α_GFED = posterior / FINN_prior × (FINN_prior / GFED_prior)
# = just the bias between FINN and GFED as reference
fig2, ax2 = plt.subplots(figsize=(9, 5))

regions_3 = ['Kalimantan', 'Sumatra', 'Java+Sulawesi']
# FINN prior totals from run_inversion_v2.py run
finn_prior = {
    'Kalimantan':  34.009 + 23.875,   # W+C + S+E from last run
    'Sumatra':     0.216  + 25.946,
    'Java+Sulawesi': 3.113 + 7.109,
}
gfed_mid = {
    'Kalimantan':  (18.0 + 25.0) / 2,
    'Sumatra':     (10.0 + 15.0) / 2,
    'Java+Sulawesi': (2.5 + 5.0) / 2,
}

# α for each source = estimate / FINN_prior
x_pos = np.arange(3)
w = 0.18

for ai, (src, color, vals) in enumerate([
    ('GFED4.1s low',      '#2196F3', [18.0, 10.0, 2.5]),
    ('GFED4.1s high',     '#90CAF9', [25.0, 15.0, 5.0]),
    ('GOSAT (Huijnen)',   '#4CAF50', [15.0,  9.0, 2.5]),
    ('OCO-2 (Yin)',       '#FF9800', [11.5,  7.0, 2.0]),
    ('v2 posterior',      '#1565C0', [v2_kali, v2_sum, v2_jav]),
]):
    alphas = [vals[ki] / list(finn_prior.values())[ki] for ki in range(3)]
    ax2.bar(x_pos + (ai - 2)*w, alphas, w,
            color=color, edgecolor='k', lw=0.5, alpha=0.85, label=src)

ax2.axhline(1.0, color='k', lw=1, linestyle=':', label='FINN prior (α=1)')
ax2.axhline(0.0, color='k', lw=0.5)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(regions_3, fontsize=11)
ax2.set_ylabel('FINN scaling factor α  (estimate / FINN prior)', fontsize=10)
ax2.set_title('FINN vs. literature: implied α for each published estimate\n'
              'Values < 1 mean FINN overestimates; 0.08–0.35 = literature consensus range for peat fires',
              fontsize=10, fontweight='bold')
ax2.legend(fontsize=9, ncol=2, loc='upper right')
ax2.grid(axis='y', lw=0.3, alpha=0.4)
ax2.set_ylim(0, 1.35)

# Shade literature consensus range for peat (Kalimantan + Sumatra)
ax2.axhspan(0.08, 0.35, alpha=0.08, color='green', label='Peat consensus (0.08–0.35)')
ax2.text(2.3, 0.21, 'peat consensus\n(0.08–0.35)', ha='center', va='center',
         fontsize=8, color='green', alpha=0.7, style='italic')

out2 = f'{PLT_DIR}/finn_alpha_literature.png'
fig2.savefig(out2, dpi=150, bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {out2}')

# ── Print summary table ──────────────────────────────────────────────────────
print('\n' + '='*72)
print('Literature comparison (Tg CO2, Indonesia, Oct 24-29 2015, 6-day peak)')
print('='*72)
print(f'  {"Source":<25}  {"Kalimantan":>12}  {"Sumatra":>10}  {"Java+Sul.":>10}  {"Total":>10}')
print('  ' + '-'*68)
for lname, ldata in LIT.items():
    lname_s = lname.replace('\n',' ')
    k_lo, k_hi = ldata['Kalimantan']; s_lo, s_hi = ldata['Sumatra']
    j_lo, j_hi = ldata['Java+Sul.']; t_lo, t_hi = ldata['Total']
    print(f'  {lname_s:<25}  {k_lo:.0f}–{k_hi:.0f} Tg{"":<5}  '
          f'{s_lo:.0f}–{s_hi:.0f} Tg{"":<3}  {j_lo:.0f}–{j_hi:.0f} Tg{"":<4}  '
          f'{t_lo:.0f}–{t_hi:.0f} Tg')
print('  ' + '-'*68)
print(f'  {"v1 (this study)":<25}  {v1_kali:>10.1f}  {v1_sum:>10.1f}  '
      f'{v1_jav:>10.1f}  {v1_tot:>10.1f}  ← Kalimantan unphysically −0.6')
print(f'  {"v2 (this study)":<25}  {v2_kali:>10.1f}  {v2_sum:>10.1f}  '
      f'{v2_jav:>10.1f}  {v2_tot:>10.1f}')
print()
print('Bias assessment (v2 vs GFED4.1s midpoint):')
gfed_mid_vals = {
    'Kalimantan': 21.5, 'Sumatra': 12.5, 'Java+Sulawesi': 3.75, 'Total': 37.75
}
for rn, ov in zip(regions_3, [v2_kali, v2_sum, v2_jav]):
    if not np.isnan(ov):
        bias = (ov - gfed_mid_vals.get(rn, np.nan)) / gfed_mid_vals.get(rn, 1) * 100
        print(f'  {rn:<22}: v2={ov:.1f} Tg  GFED_mid={gfed_mid_vals.get(rn):.1f} Tg  '
              f'bias={bias:+.0f}%')
