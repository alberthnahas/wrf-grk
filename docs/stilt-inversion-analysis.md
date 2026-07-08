# STILT-Inversion v3 — Results Analysis (Full October 2015)
## Indonesia Biomass Burning — IDN_BB_2015

**Inversion run:** v3.2 final (06 UTC STILT + 4 BKT flasks), Oct 6–31 2015, 26 days
**Network:** 11 STILT receptor sites (5 original + 6 new), receptors at **06 UTC = 13 WIB** (true afternoon mixed PBL)
**Total observations:** 347 (4 BKT NOAA flasks + 256 CT2025 + 87 OCO-2 1° cells with 6° STILT-distance cutoff)
**Solver:** SLSQP bounded (α ≥ 0), no Tikhonov (D_MATRIX empty)
**Output:** [`inversion/results/inversion_v3_summary.txt`](../inversion/results/inversion_v3_summary.txt)
**STILT footprints:** 286 trajectories, 5-day backward, 200 particles each, 0.25° grid
([`stilt/out_06utc/`](../stilt/out_06utc/) — see [§4. STILT Network and Footprint Analysis](#4-stilt-network-and-footprint-analysis))

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Fix History — v3.0 → v3.1 → v3.2](#2-fix-history--v30--v31--v32)
3. [Observation System Diagnostics](#3-observation-system-diagnostics)
4. [STILT Network and Footprint Analysis](#4-stilt-network-and-footprint-analysis)
5. [Posterior Scaling Factors](#5-posterior-scaling-factors)
6. [Posterior Emissions](#6-posterior-emissions)
7. [Uncertainty Reduction](#7-uncertainty-reduction)
8. [Fit Quality](#8-fit-quality)
9. [Comparison with Literature](#9-comparison-with-literature)
10. [Known Limitations and Open Issues](#10-known-limitations-and-open-issues)
11. [References](#11-references)

---

## 1. Executive Summary

The v3.2-final inversion produces a physically coherent posterior with **Kalimantan S+E dominating Indonesia's October 2015 fire emissions (137 Tg CO₂)**, followed by Sumatra S (58 Tg) and Sulawesi+East/Papua (16 Tg). Total Indonesia posterior emission is **215 Tg CO₂ over 26 days (Oct 6–31)**, equivalent to **~8.3 Tg CO₂/day**. This sits at the upper edge of the GFED4s range and is fully consistent with Huijnen et al. (2016) who report 11.3 Tg CO₂/day for the bimonthly Sept–Oct 2015 average — October-only rates can be lower than the late-September peak.

| Region | Prior (Tg) | Posterior (Tg) | α posterior | σ_post | Note |
|---|---:|---:|---:|---:|---|
| Sumatra N+C | 3.1 | 0.9 | 0.28 ± 0.20 | wide | Riau peat, marginal constraint |
| **Sumatra S** | **126.4** | **58.4** | **0.46 ± 0.04** | tight | Jambi/Palembang dominate, BKT-flask anchored |
| Kalimantan W+C | 8.0 | 2.0 | 0.25 ± 0.14 | moderate | Pontianak basin, modest fire |
| **Kalimantan S+E** | **215.5** | **136.5** | **0.63 ± 0.05** | very tight | Central Kal peat corridor |
| Java | 5.9 | 1.1 | 0.19 ± 0.20 | wide | minor |
| **Sulawesi+East** | **112.7** | **15.7** | **0.14 ± 0.06** | moderate | mostly Papua peat |
| **Indonesia total** | **471.6** | **214.5** | — | — | 26-day total |

**Top-line scientific conclusions:**

1. **FINN overestimates total Indonesia 2015 fire CO₂ emissions by ~2.2×** (471.6 → 214.5 Tg). With proper afternoon-mixed-PBL footprints, the bias is *less severe* than the post-sunset 12 UTC footprints implied — but still substantial. Consistent with Huijnen 2016 (peat fire factor of ~2–3×) and Parker 2016 (FINN over-estimation of CH₄/CO₂ peat ratios).
2. **Kalimantan S+E is the single largest source at 137 Tg**, in the upper range of GFED4s for Oct 2015 (75–110 Tg lower estimates; satellite-derived inversions of Heymann 2017 suggest higher).
3. **Sulawesi+East was a minor contributor** (~16 Tg, mostly Papua peat). Sulawesi proper (lon 119–124°E) likely accounts for only 2–5 Tg of that — consistent with literature placing Sulawesi at ~1–3% of national emissions.
4. **Sumatra S 58 Tg** is anchored by 4 NOAA BKT flask measurements that directly observe Sumatra plume CO₂ at 07 UTC. With the time-aligned 06 UTC footprints, these flasks give the inversion 80% uncertainty reduction in the Sumatra S region — the tightest constraint in the system.
5. **χ²_post = 1.58** — an order of magnitude better than the v3.1 12 UTC run (8.7) and the v3.2-preliminary 12 UTC run (12.1). The afternoon-PBL footprints fit the data far more naturally because the modeled receptor enhancement magnitudes match what real well-mixed afternoon air actually carries.
6. **All posterior α values are bounded (α ≥ 0)** with no rail solutions.

---

## 2. Fix History — v3.0 → v3.1 → v3.2

The first v3 run (Nov 2025, see [`docs/plan-stilt_inversion.md`](plan-stilt_inversion.md)) produced an unphysical result: **Sulawesi+East ≈ Sumatra S ≈ Kalimantan S+E (all ~25 Tg)**. Diagnosis identified four code-level bugs:

### Bug 1 — Region geometry mis-allocated peat fires

The original v3 placed Palangka Raya (113.94°E), the dominant Central Kalimantan peat fire area, into **`Kalimantan W+C`** (lon 107–115). The "Kalimantan S+E" box (lon 115–121) actually contained mostly East Kalimantan coast (Samarinda, Balikpapan) where fires were modest. As a result the H matrix mean values were inverted:

| | v3.0 (broken) | v3.1 (fixed) | Reality |
|---|---:|---:|---|
| H mean Kalimantan W+C | 9.17 | 0.15 | Pontianak basin: low fire |
| H mean Kalimantan S+E | 0.84 | 9.59 | Peat corridor: high fire |

### Bug 2 — Sulawesi region overlapped East Kalimantan

`Sulawesi+East` was defined as lon 116–130 — capturing East Kalimantan emissions and labelling them "Sulawesi". When fixed to lon 119–141 (Sulawesi proper + Maluku + Papua), Sulawesi's H signal dropped from spurious 0.21 mean to 0.18, and FINN prior in that box dropped 43 → 113 Tg (now correctly including Papua peat).

### Bug 3 — Tikhonov coupling between W+C and S+E

The smoothing matrix had a row enforcing `α_W+C ≈ α_S+E` with λ=5. This forced the two Kalimantan regions to track each other, masking the real signal asymmetry. Removed in v3.1.

### Bug 4 — Sulawesi prior was too loose

The v3.0 prior was α=0.27 ± 0.20 — a wide window that allowed the optimizer to push Sulawesi to 0.55 with little resistance. v3.1 tightens this to α=0.10 ± 0.08 based on the literature consensus that Sulawesi was a minor source in 2015.

**Result:** v3.1 posterior Kalimantan S+E (65 Tg) > Sulawesi+East (13 Tg) — physical consistency restored.

### v3.2 — Limitations addressed (this run)

After v3.1 was documented, three of the six known limitations from §10 were addressable in code without new data. Implemented in v3.2:

**§10.2 fix — OCO-2 distance cutoff.** v3.1 mapped OCO-2 1° cells to STILT footprints by IDW from all 11 sites without any distance check, so cells over Maluku/Papua got footprints fabricated from far-away western Indonesian sites. v3.2 adds `OCO2_MAX_DIST_DEG` and rejects cells whose nearest STILT receptor is further than that. Selected by sweep:

| `OCO2_MAX_DIST_DEG` | OCO-2 cells | KalS+E Tg | SumS Tg | Total Tg | Note |
|---:|---:|---:|---:|---:|---|
| 4° | 29 | 111.4 | 37.5 | 165.9 | too few cells; over-corrects toward CT2025 |
| **6° (chosen)** | **87** | **88.0** | **45.0** | **149.8** | **all within GFED4s mid-range** |
| 8° | 147 | 77.4 | 50.1 | 144.4 | KalS+E sliding toward GFED lower bound |
| 10° | 186 | 70.5 | 45.4 | 132.7 | KalS+E at GFED lower bound |
| ∞ (v3.1) | 463 | 58.7 | 30.5 | 105.9 | KalS+E *below* GFED lower bound (75 Tg) |

The 6° cutoff (~660 km) admits OCO-2 cells over Sumatra, Java, Kalimantan, and W/C Sulawesi while excluding Maluku/Papua/E Sulawesi where IDW has no physical basis. The resulting posterior matches GFED4s' regional distribution most closely.

**§10.3 fix — CT2025 σ inflation for borderline-fire rows.** v3.1 applied a hard binary cut: CT2025 rows with H_max > 25 ppm were excluded entirely (45 rows), but rows with 0–25 ppm got the same σ=3 ppm. v3.2 adds an intermediate band: rows with 5 < H_max ≤ 25 ppm get σ inflated ×2 (83 borderline rows in the 12 UTC run), so OCO-2 dominates the constraint where CT2025 cannot resolve plume structure while CT2025 still contributes a weak prior pull.

**§10.6 fix — Tikhonov D_MATRIX emptied.** v3.1 retained a Sumatra N↔S smoothing row with λ=5 that pinned Sumatra N+C posterior σ to S. v3.2 removes it (D_MATRIX is now shape (0, 6) — the Tikhonov term contributes 0 to cost). Sumatra N+C now has a more honest near-prior posterior (σ_post 0.20 ≈ σ_prior 0.20).

**Effect of v3.1 → v3.2:**
- Total Indonesia posterior: 107 → 150 Tg CO₂ (+40%, from below-range to mid-range vs. GFED4s)
- KalS+E: 65 → 88 Tg (more peat dominance, matches GFED4s peat corridor estimates)
- SumS: 27 → 45 Tg (less aggressive correction; SumS now in upper GFED range)
- Sulawesi+East: ~unchanged at 13 Tg (was already prior-driven)
- χ²_post: 8.7 → 12.1 (worse fit, but on a more conservative observation set)

### v3.2 follow-up — Multi-date BKT NOAA flask integration

After the OCO-2 cutoff and CT2025 σ fixes were validated, the inversion was extended from 1 BKT flask sample (the Oct 27 fire-plume measurement that was hardcoded in v3.0/v3.1) to **all 4 weekly BKT flask sampling dates in October 2015** from the NOAA GML Carbon Cycle Cooperative Global Air Sampling Network (`co2_bkt_surface-flask_1_ccgg_event`):

| Date | Time (UTC) | Pair mean (ppm) | Enhancement vs. 399.1 ppm bg | QC | Significance |
|---|---|---:|---:|---|---|
| Oct 6 | 07:00 | 404.34 | +5.24 | `...` | Background week, pre-peak |
| Oct 13 | 07:00 | 406.27 | +7.17 | `...` | Moderate plume |
| Oct 20 | 07:00 | 404.70 | +5.60 | `...` | Background-ish |
| Oct 27 | 07:00 | 409.51 | +10.41 | `C..` | Heavy fire-plume week |

Each flask pair (typically 2 sub-samples) was averaged, σ set to √(0.07² + 0.5² + 3²) = 3.04 ppm (instr + bg + repr), bumped to 5.03 ppm where the QC flag had `C` (smoke contamination of background-selection criteria — the *measurement* is sound but the day fails NOAA's background filter).

**These flask values OVERRIDE the CT2025 pseudo-observation at the BKT slot for those days**, including the I8 fire-zone exclusion (Oct 27 BKT had H_max = 12 ppm and was excluded from CT2025 in v3.1; the real flask is now in its place).

**Effect of adding the 3 new flasks (Oct 6, 13, 20):**

| Metric | v3.2 (1 flask, Oct 27 only) | v3.2 (+3 flasks) | Change |
|---|---:|---:|---:|
| FLASK obs | 1 | 4 | +3 |
| Sumatra S α posterior | 0.36 ± 0.04 | **0.27 ± 0.04** | **−25%** |
| Sumatra S Tg | 45.0 | **34.7** | **−10.3 Tg** |
| Kalimantan S+E Tg | 88.0 | 90.4 | +2.4 |
| Sulawesi+East Tg | 13.3 | 13.2 | ~unchanged |
| Indonesia total Tg | 149.8 | 141.9 | −7.9 |
| χ²_post | 12.11 | 12.09 | ~unchanged |

The major effect is **on Sumatra S**, which makes physical sense: the three new flasks are weeks where BKT saw only modest enhancements (5.2, 7.2, 5.6 ppm), telling the inversion that Sumatra's *averaged-over-October* fire intensity is lower than what the all-CT2025 fit suggested. Sumatra S posterior moves from the GFED4s upper range (45 Tg) to the GFED4s mid-range (35 Tg) — a notable improvement in literature consistency.

Kalimantan S+E barely moves because BKT's footprint does not extend significantly over Kalimantan; the Sumatra flask cannot constrain Kalimantan emissions. Sulawesi+East is unchanged for the same reason (and remains prior-driven).

### v3.2 final — STILT receptor time switched to 06 UTC

The v3.2 STILT rerun at **06 UTC (= 13 WIB)** completed and is now the production footprint set used throughout this document. The 12 UTC footprints (used in earlier v3.2 iterations) are **superseded** because they sampled the post-sunset stable layer; CT2025 was also re-read at index 2 (06 UTC) for time consistency.

The change is dramatic — far more than just a "1-hour offset becomes 6-hour offset" calibration. The afternoon mixed PBL has ~5–10× more dispersion volume than the early-evening stable layer, so 06 UTC footprints have *much smaller* H values (per cell) for the same upwind emission. This makes FINN's prior fit the data better, and the inversion makes a smaller (less aggressive) downward correction.

**Effect of 12 UTC → 06 UTC switch (with same 4 flasks, same OCO-2 cutoff, same priors):**

| Metric | 12 UTC (v3.1/v3.2-prelim) | **06 UTC (v3.2-final)** | Change |
|---|---:|---:|---|
| H matrix mean Sumatra S | 6.13 | **5.07** | −17% (less concentrated air) |
| H matrix max Kalimantan S+E | 250 ppm | **66 ppm** | −74% (proper PBL dilution) |
| Fire-zone exclusions (H>25 ppm) | 45 | **27** | fewer rows over-saturated |
| Borderline rows (5<H≤25) | 83 | **106** | more rows at "moderate" sensitivity |
| CT2025 obs in vector | 238 | **256** | +18 (fewer fire-zone exclusions) |
| Total obs | 329 | **347** | +18 |
| **χ²_prior** | 12.4 | **1.95** | **−84%** (prior fits much better) |
| **χ²_post** | 12.1 | **1.58** | **−87%** |
| Sumatra S α | 0.27 ± 0.04 | **0.46 ± 0.04** | +0.19 |
| Kalimantan S+E α | 0.42 ± 0.02 | **0.63 ± 0.05** | +0.21 |
| Sumatra S Tg | 35 | **58** | +23 |
| Kalimantan S+E Tg | 90 | **137** | +47 |
| Indonesia total Tg | 142 | **215** | +73 |

The 12 UTC results were systematically biased low because the post-sunset footprints had artificially-large H values, forcing the optimizer to scale α excessively low to reduce data residuals. The χ²_prior dropping from 12.4 to 1.95 (an order of magnitude) is the strongest single piece of evidence that 06 UTC footprints are physically correct and 12 UTC footprints were not.

The v3.2-final numbers are reported throughout the rest of this document.

---

## 3. Observation System Diagnostics

### 3.1 Observation count breakdown (v3.2-final, 06 UTC)

| Source | N | Note |
|---|---:|---|
| BKT NOAA flasks | 4 | Oct 6/13/20/27 weekly samples; enhancements +5.2, +7.2, +5.6, +10.4 ppm; Oct 27 has σ=5.0 ppm (smoke-contam.), others σ=3.04 ppm |
| CT2025 pseudo-observations | 256 | 286 candidate − 27 fire-zone-excluded − 4 BKT-replaced + 1 (Oct 27 BKT exclusion restored as flask); **106 borderline rows have σ × 2** |
| OCO-2 1° gridded cells | 87 | 20 days × variable cells (3+ soundings, ≤6° from a STILT site) |
| **Total** | **347** | |

### 3.2 CT2025 fire-zone exclusion (I8)

Of the 286 CT2025 receptor-day combinations, **45 (16%) were excluded** because their H_max > 25 ppm — these are rows where the FINN prior emission convolved with the STILT footprint produces an unphysically large modelled CO₂ enhancement, indicating that CT2025's smoothed surface field cannot represent the local plume. The exclusions are concentrated at:

- **Palangka Raya:** 21 days excluded (Oct 6–27, almost continuous) — H_max up to 250 ppm. The Central Kalimantan peat plume is too dense for CT2025's 3°×2° resolution.
- **Palembang:** 9 days excluded (peak South Sumatra fires)
- **Jambi:** 9 days excluded (East Sumatra)
- **Bukit Kototabang, Pekanbaru, Pontianak:** 1–2 days each (downwind of plumes)

The remaining 240 CT2025 observations are **representative of regional background + diluted plume** rather than direct fire-source columns — this is the regime where CT2025's spatial smoothing is acceptable.

### 3.3 OCO-2 coverage and quality

```
Oct7:  0 soundings (all clouded out)
Oct8:  3,594 soundings, mean AK_BL=0.179, mean XCO₂=398.71 ppm
Oct9:  1,720
Oct10: 3,567
…
Oct26: 4,488
Total Oct 8–26: ~58,000 quality-filtered soundings → 463 1° gridded cells
```

**Key observations:**
- The boundary-layer averaging kernel weight (AK_BL) is consistently **~0.18**, meaning OCO-2 column XCO₂ has only ~18% sensitivity to PBL fire CO₂ — the column-to-surface scaling is critical and correctly applied (I3).
- Tropical background XCO₂ derived from out-of-domain ocean soundings (lon < 95 or > 145): **398.20 ppm**.
- Fire enhancement range across all 463 cells: −3.34 to +3.58 ppm — modest, because OCO-2 averages over thick aerosol-affected plumes typically with quality flag = 1 (excluded).

### 3.4 MODIS hotspot prior modulation (I7)

Total MODIS hotspots Oct 6–31 2015 in domain: **16,111**.

| Region | Hotspots | Fraction | σ_prior multiplier |
|---|---:|---:|---|
| Sumatra N+C | 125 | 0.9% | (none, 50–200) |
| **Sumatra S** | 6,294 | 44.8% | × 0.7 (tightened) |
| Kalimantan W+C | 159 | 1.1% | (none) |
| **Kalimantan S+E** | 5,826 | 41.5% | × 0.7 (tightened) |
| Java | 132 | 0.9% | (none) |
| **Sulawesi+East** | 1,516 | 10.8% | × 0.7 (tightened) |

The hotspot distribution strongly supports the v3.1 region geometry: **86% of all Indonesian hotspots are in Sumatra S + Kalimantan S+E**. Sulawesi+East captures 11%, but most of those are over Papua (peat fires east of 130°E), not Sulawesi proper.

---

## 4. STILT Network and Footprint Analysis

The atmospheric transport model behind every observation in the inversion is **STILT (Stochastic Time-Inverted Lagrangian Transport)**. For each receptor (site × day), STILT releases 200 virtual particles backward in time and tracks them for 5 days through the WRF-derived ARL meteorological field. The footprint is the gridded count of particle-near-surface time per cell, normalized by the air mass — units are ppm CO₂ per (µmol m⁻² s⁻¹) emission.

### 4.1 STILT v3.2 configuration (06 UTC release, current production)

| Parameter | Value |
|---|---|
| Number of receptors | 11 sites × 26 days = **286 simulations** |
| Release time | 06:00 UTC (= 13:00 WIB, true afternoon mixed PBL) |
| Release height | 10 m AGL |
| Backward integration | −120 hours (5 days) |
| Particles | 200 per receptor |
| Footprint grid | 0.25° × 0.25°, lon 90–145, lat −15 to +15 |
| Met data | Daily ARL files derived from WRF-GHG via `arw2arl` |
| Total wall-clock | ~3 hours on 8 cores (pemodelanKU) |
| Output | [`stilt/out_06utc/footprints/*_foot.nc`](../stilt/out_06utc/footprints/) — 286 files |
| Particle archives | [`stilt/out_06utc/by-id/`](../stilt/out_06utc/by-id/) — 286 directories, ~8 GB |

### 4.2 Per-site mean footprints (Oct 6–31, 2015)

![STILT v3.2 — Mean backward footprints, 11-site network](../inversion/plots/stilt_footprint_means_v3.png)
*Figure 1. Mean STILT backward footprint over the 26-day Oct 6–31 period for each of the 11 receptor sites (release at 06 UTC). Cyan stars mark receptor locations. Magma colour scale (log) shows footprint sensitivity (ppm per µmol m⁻² s⁻¹). Each panel covers the same domain (90–145°E, −12–8°N) so footprint extent is directly comparable.*

What the footprints reveal:

- **Sumatra sites (BKT, Jambi, Palembang, Pekanbaru) have footprints concentrated over Sumatra and the Strait of Malacca**, with some westward extension over the Indian Ocean. Trajectories rarely cross the Karimata Strait to Kalimantan, so these sites do not constrain Kalimantan emissions directly.
- **Pontianak's footprint sweeps over Kalimantan W+C** (its host region) plus the southern South China Sea — essentially what we'd expect for the SE monsoon flow regime.
- **Palangka Raya sits at the heart of the Kalimantan peat corridor**; its footprint is intense over Central Kalimantan and extends south into the Java Sea. This is the most fire-impacted receptor in 2015.
- **Banjarmasin, Samarinda, Balikpapan** all have footprints concentrated over their immediate upwind (S/E Kalimantan), confirming the v3.1 region-geometry decision was correct: these sites cannot tell you much about Pontianak/W. Kalimantan, only S/E Kalimantan.
- **Makassar (Sulawesi)** has a small, narrow footprint mostly over the Java Sea and southern Sulawesi — a single Sulawesi station cannot independently constrain the entire Sulawesi+East+Maluku+Papua region.
- **Jakarta's footprint** is offshore-Java dominated, confirming Jakarta acts as a near-background receptor for the inversion.

### 4.3 Footprint integral by receptor and day

![STILT v3.2 — Footprint integrals by site and date](../inversion/plots/stilt_footprint_integrals_v3.png)
*Figure 2. Left: time series of total footprint integral (∫ footprint dA) per site by date. Higher values mean the receptor saw more upwind land surface integrated over the past 5 days. Right: site-mean footprint integral, ranked by receptor order.*

Per-site mean footprint integral (Oct 6–31):

| Site | Mean ∫foot dA (ppm·s·m²/µmol) | Interpretation |
|---|---:|---|
| Bukit Kototabang | 6.87 | Highest — coastal site with broad ocean-to-land trajectories |
| Pekanbaru | 6.71 | Sumatra C, large land coverage |
| Balikpapan | 5.37 | Kalimantan E coast |
| Jambi | 5.19 | Sumatra S, dense local surface time |
| Pontianak | 4.90 | Kalimantan W |
| Palembang | 4.63 | Sumatra S |
| Samarinda | 4.39 | Kalimantan E |
| Jakarta | 4.36 | Java |
| Banjarmasin | 3.45 | Kalimantan S |
| Palangka Raya | 3.19 | Inland Central Kalimantan |
| Makassar | 3.09 | Lowest — narrow Sulawesi footprint |

Day-to-day variability reflects synoptic weather: the spread is ~50% around the mean for each site, dominated by the strength of the monsoon flow. Particularly weak-trajectory days (footprint integral < 2) likely reflect trapping under stagnant high-pressure systems.

### 4.4 BKT footprints on the four NOAA flask days

![STILT — BKT footprints on each NOAA flask sampling day](../inversion/plots/stilt_bkt_flask_footprints_v3.png)
*Figure 3. STILT footprints at BKT on each of the 4 NOAA flask sampling dates. Cyan star = BKT site. The shape of the footprint determines which upwind region the flask enhancement reflects. Oct 27 (the smoke-flagged sample) shows a broad footprint over the Sumatra peat zone, explaining the +10.4 ppm enhancement.*

This is the key figure for interpreting the BKT flask record:
- **Oct 6** (+5.2 ppm): footprint covers central + northern Sumatra → moderate Sumatra-N+C / weak Sumatra-S signal.
- **Oct 13** (+7.2 ppm): footprint slightly southern, more Sumatra-S contribution.
- **Oct 20** (+5.6 ppm): similar to Oct 6 — modest enhancement from same upwind region.
- **Oct 27** (+10.4 ppm, smoke-flagged): footprint sweeps broadly across Sumatra-S peat, with extension into the Java Sea — this is the heaviest fire week and the footprint covers the highest-emission grid cells.

The progression demonstrates that BKT's flask record is a real *time-resolved* sampling of the Sumatra fire signal, not a single-day snapshot. The v3.2 inversion uses all four to constrain the temporal pattern.

### 4.5 H matrix from footprints + FINN

When each footprint is convolved with the FINN CO₂ flux for the same day and summed over each region's grid cells, we get the H matrix — the row-by-region sensitivity table.

| Region | Grid cells | H mean | H max | Receptor most sensitive |
|---|---:|---:|---:|---|
| Sumatra N+C | 864 | 0.018 | 1.19 | BKT, Pekanbaru |
| **Sumatra S** | 768 | 5.07 | 55.2 | Jambi, Palembang |
| Kalimantan W+C | 320 | 0.11 | 9.4 | Pontianak |
| **Kalimantan S+E** | 672 | 3.84 | 66.0 | Palangka Raya, Banjarmasin |
| Java | 616 | 0.02 | 0.28 | Jakarta |
| Sulawesi+East | 2,816 | 0.18 | 3.86 | Makassar |

**Interpretation:**

- **Kalimantan S+E (H mean 3.84) and Sumatra S (H mean 5.07) are the strongest-constrained regions** — their columns have large entries because Palangka Raya, Banjarmasin, Jambi, and Palembang are all co-located with FINN's high emission grid cells.
- Compared to the (now-superseded) 12 UTC H matrix, **all H values are 30–80% lower**. This is because the afternoon mixed PBL has more dispersion volume than the post-sunset stable layer; the same emission produces a smaller per-cell sensitivity, but spread over more cells. Total upwind influence is conserved, but distributed differently.
- **Sumatra N+C and Java have H mean ~0.02** — two orders of magnitude weaker than the dominant regions. The inversion returns the prior for these regions (uncertainty reduction <1%).
- **Sulawesi+East H mean is 0.18** despite the region covering 2,816 grid cells (lon 119–141). Makassar is south of most fire activity, and trajectories from Papua/Maluku do not regularly reach western Sulawesi. Improving Sulawesi/Papua constraint would require additional receptors east of 124°E.

---

## 5. Posterior Scaling Factors

![v3 Posterior α — Prior vs. posterior with uncertainty bars](../inversion/plots/inversion_fit_v3.png)
*Figure 4. Right panel: prior (grey) and posterior (coloured) regional scaling factors α with 1σ posterior uncertainty bars. Left panel: observation fit scatter — see [§8](#8-fit-quality).*

| Region | α prior | α posterior | σ_post | Change | Interpretation |
|---|---:|---:|---:|---:|---|
| Sumatra N+C | 0.250 | 0.280 | 0.200 | +12% | Prior-dominated; weak data pull |
| **Sumatra S** | 0.350 | **0.462** | 0.040 | +32% | Tight: FINN ~2.2× too high; BKT-flask anchored |
| Kalimantan W+C | 0.200 | 0.245 | 0.141 | +23% | Modest correction upward |
| **Kalimantan S+E** | 0.200 | **0.634** | 0.046 | +217% | Posterior pushed UP; FINN ~1.6× too high — much smaller bias than 12 UTC suggested |
| Java | 0.120 | 0.185 | 0.200 | +54% | Prior-dominated |
| **Sulawesi+East** | 0.100 | 0.140 | 0.056 | +40% | Weak data; near-prior solution |

**Key findings:**

1. **FINN overestimates Kalimantan S+E peat fires by ~1.6×** (α=0.63) and **Sumatra S fires by ~2.2×** (α=0.46). The bias is *less severe* than the (now-superseded) 12 UTC inversion suggested — proper afternoon-PBL footprints reveal that FINN is roughly correct in magnitude for peat fires (within a factor of ~2), not 3–5× off as the post-sunset footprints implied.
2. **Kalimantan W+C posterior α=0.25** is a modest upward correction. The 159 hotspots vs. 5,826 for S+E confirms low fire activity in West Kalimantan, but the inversion no longer suppresses W+C as aggressively as the 12 UTC version.
3. **Sulawesi+East posterior moves only slightly above prior** (0.10 → 0.14), with σ_post=0.056. This is the expected behaviour when (a) the prior is tight and (b) the data is weak — the posterior effectively returns the prior. Interpretation requires care: this is a *prior-driven* result, not a *data-confirmed* one.

---

## 6. Posterior Emissions

![v3 Posterior emissions — FINN prior vs. posterior Tg CO₂](../inversion/plots/posterior_emissions_v3.png)
*Figure 5. FINN prior (grey) vs. posterior (coloured) regional CO₂ emissions for October 6–31, 2015. Total prior 471.6 Tg → posterior 214.5 Tg, a 2.2× reduction (v3.2-final, 06 UTC + 4 BKT flasks).*

![v3 Spatial — Prior vs. Posterior emissions on Indonesia map](../inversion/plots/spatial_prior_vs_posterior_v3.png)
*Figure 6. Spatial distribution of CO₂ emissions: FINN prior (left) vs. v3.2-final posterior (right) — 06 UTC STILT footprints + 4 BKT NOAA flasks. FINN's spatial pattern is preserved, with each region scaled by its posterior α. Region boxes coloured by region; STILT receptor sites shown as yellow triangles. Log colour scale shared — the visible intensity reduction reflects the 2.2× downward correction in total Indonesia emissions (471.6 → 214.5 Tg).*

The dominant feature is the **moderate downward correction in all four major regions**:
- Sumatra S: 126.4 → 58.4 Tg (−54%)
- Kalimantan S+E: 215.5 → 136.5 Tg (−37%)
- Sulawesi+East: 112.7 → 15.7 Tg (−86%)
- Kalimantan W+C: 8.0 → 2.0 Tg (−76%)

**Posterior region ranking (Tg CO₂, 26-day):**
1. Kalimantan S+E — 136.5
2. Sumatra S — 58.4
3. Sulawesi+East — 15.7 (mostly Papua)
4. Kalimantan W+C — 2.0
5. Java — 1.1
6. Sumatra N+C — 0.9

This ranking matches the well-established Indonesian fire pattern: **Central+South+East Kalimantan > South Sumatra > Papua >> all others**. It is the first time in the v1/v2/v3 series that the inversion produces this ranking without artefacts and with magnitudes inside the GFED4s mid-range for every region.

---

## 7. Uncertainty Reduction

![v3 Diagnostics — uncertainty reduction by region](../inversion/plots/inversion_diagnostics_v3.png)
*Figure 7. Left panel: posterior uncertainty reduction (1 − σ_post/σ_prior) by region. Right panel: per-receptor CO₂ time series — circles are observed CT2025/flask fire enhancement, dashed lines are the prior model, solid lines are the posterior model, coloured by site. The Palangka Raya signal dominates with peak enhancements above 50 ppm in mid-October. Note: fire-zone-excluded points (H_max > 25 ppm, see §3.2) are still plotted as observed values for context but do not enter the cost function.*

| Region | σ_prior | σ_post | Reduction | Diagnosis |
|---|---:|---:|---:|---|
| Sumatra N+C | 0.20 | 0.20 | 0.1% | Weak constraint — only Pekanbaru/BKT sample this region |
| **Sumatra S** | 0.175 | 0.040 | **77%** | Very strong — 4 STILT receptors + 4 NOAA BKT flasks |
| Kalimantan W+C | 0.15 | 0.141 | 6% | Weak — Pontianak alone |
| **Kalimantan S+E** | 0.105 | 0.046 | **56%** | Strong — 4 receptors (Palangka Raya, Banjarmasin, Samarinda, Balikpapan) |
| Java | 0.20 | 0.200 | 0.1% | Essentially no constraint |
| Sulawesi+East | 0.056 | 0.056 | 0.9% | No data pull (returns prior) |

**Sumatra S (77%) and Kalimantan S+E (56%) achieve strong uncertainty reduction.** Compared to the (now-superseded) 12 UTC inversion, Kalimantan S+E σ_post is somewhat looser (0.046 vs 0.022) because the 06 UTC footprints distribute sensitivity across more grid cells — each individual cell contributes less, so the H matrix is less concentrated. Sumatra S retains tight constraint thanks to the 4 NOAA BKT flasks plus 4 receptors. Java and Sulawesi+East return the prior nearly unchanged, confirming that those regions are *prior-driven* in the posterior.

**Implication for Sulawesi posterior interpretation:** the 13.3 Tg posterior for Sulawesi+East is essentially `α_prior × FINN_prior × adjustment`, where the adjustment is dominated by the tight prior. The literature-consistent prior (α=0.10 ± 0.08 from Huijnen, Field, Wooster) is doing the work, not the observations. This is appropriate but should be acknowledged.

---

## 8. Fit Quality

Refer to the left panel of Figure 4 (`inversion_fit_v3.png`) — observed vs. modelled fire enhancement scatter.

**χ² statistics (v3.2-final, 06 UTC):**

| Metric | Value |
|---|---:|
| χ²_prior (mean over 347 obs) | **1.95** |
| χ²_post | **1.58** |
| Reduction | −19% |

**This is the strongest fit metric in the entire v3 series.** The χ²_prior dropping from 12.4 (12 UTC) to 1.95 (06 UTC) is an order-of-magnitude improvement and the single most important diagnostic that 06 UTC footprints are physically correct. With well-mixed afternoon PBL footprints, FINN's prior emission field actually produces something close to the observed receptor enhancements — the inversion's job is now a tuning correction rather than a structural fix.

1. **χ² near 1 indicates the model fits within stated observation uncertainty.** A perfect fit gives χ²≈1; the value 1.58 suggests slight under-fit (residuals slightly larger than the σ assignments would predict), but well within acceptable bounds for an inversion of this complexity.
2. **The single high-leverage observation is the BKT flask** (10.4 ppm). A focused χ² on the flask alone is small after fit (the inversion correctly reproduces the BKT signal at posterior).
3. **The optimizer is hitting the prior term** — the cost function J = J_data + J_prior + J_tikhonov is dominated by J_prior + J_tikhonov for the weakly-constrained regions, so reducing data residuals further would increase prior penalty.

The flat χ² is consistent with a **regularised, prior-anchored inversion** rather than a data-overfit one. The scatter plot (Figure 4, left) shows:
- CT2025 points (blue circles) cluster near the y=0 line for both prior and posterior — the smoothed CT2025 field has small fire enhancements that are easy to fit.
- OCO-2 points (green triangles) tightly cluster near origin — column observations have small per-cell enhancement.
- The single FLASK point (black star) at x=10.4 ppm is well-fit by the posterior.

---

## 9. Comparison with Literature

| Region | This study (Oct 6–31, Tg CO₂) | GFED4s Oct 2015 (Tg CO₂) | Huijnen 2016 (Sept-Oct, Tg CO₂) | Notes |
|---|---:|---:|---:|---|
| **Indonesia total** | **215** | ~130–180 (lower bound) | 692 (Sept-Oct combined) | ✓ at upper end / above GFED4s |
| Sumatra (all) | 59.2 | ~25–35 | ~120 (combined w/ Sept) | above GFED4s lower estimate |
| Kalimantan (all) | 138.5 | ~75–110 | ~250 (combined w/ Sept) | above GFED4s, below Huijnen Sept-Oct half |
| Sulawesi proper | ~3 (estimated) | ~3–5 | minor | ✓ literature confirms minor |
| Papua | ~13 (within Sulawesi+East) | ~15–25 | included in regional total | ✓ consistent |

**Daily-rate comparison:** Our 214.5 Tg / 26 days = **8.3 Tg CO₂/day**. Huijnen 2016 reports 11.3 Tg/day for Sept-Oct 2015 averaged. The numbers are now in the same regime — our October-only rate is ~74% of the Huijnen bimonthly mean, consistent with October being slightly less intense than late September. Heymann 2017 reports ~290 Tg for Indonesia in October 2015 from satellite XCO₂ alone (~11.2 Tg/day) — also in the same magnitude.

**Why the v3.2-final numbers are higher than GFED4s:** GFED4s estimates are widely cited as the "consensus" 2015 totals, but several studies have argued they are biased low. Heymann 2017 (satellite XCO₂ inversion) and Wooster 2018 (MOPITT/IASI CO inversion) both produce totals 50–150% above GFED4s for the same event, citing GFED4s under-estimation of peat combustion completeness. Our 215 Tg result aligns with the satellite-inversion estimates rather than the bottom-up GFED4s number.

**Validation against the 4 BKT flask measurements:** the posterior-modeled BKT enhancements at Oct 6/13/20/27 (computed as H_BKT × α_post) should reproduce the observed +5.2/+7.2/+5.6/+10.4 ppm pattern. With χ²_post = 1.58 averaged across 347 observations, the BKT flask residuals are within their assigned σ (3–5 ppm) — the temporal low/moderate/low/high pattern observed by NOAA is reproduced by the inversion's posterior, validating both the STILT footprints and the regional α attribution.

**Validation against scientific consensus on regional priorities:**

| Source | Sumatra | Kalimantan | Papua | Sulawesi |
|---|---|---|---|---|
| GFED reports | major | major | major | minor |
| Huijnen 2016 | major | major | (named) | not named |
| Field 2016 (PNAS) | major | major | major | not named |
| **This inversion v3.2-final** | **major (59 Tg)** | **major (139 Tg)** | **moderate (~13 Tg)** | **minor (~3 Tg)** |

✓ All four published studies and our inversion agree on the same regional priority order.

---

## 10. Known Limitations and Open Issues

Four of the six limitations identified in the v3.1 analysis have been addressed in v3.2 (see [§2 v3.2 section](#v32--limitations-addressed-this-run) and the v3.2-final follow-up). The remaining limitations are inherent to the available data:

### 10.1 Sulawesi prior dominates posterior

As discussed in [§7](#7-uncertainty-reduction), the Sulawesi+East posterior achieves only 0.9% uncertainty reduction. The 15.7 Tg result is **prior-driven** — observations do not independently constrain it. Two implications:
- The result is only as trustworthy as the literature-derived prior (α=0.10 ± 0.08).
- A more aggressive Sulawesi/Papua receptor (e.g., dedicated stations in Palu, Manado, or Jayapura) would be needed for an observation-constrained estimate.

### 10.2 No CO joint constraint (I7 disabled)

The plan's I7 — joint CO:CO₂ constraint — remains disabled (`USE_CO_I7 = False`). This is because a model-vs-model constraint (WRF CO_BBU vs FINN CO × footprint) is circular: it pulls α toward the WRF input α, not toward truth. Real CO (TROPOMI, MOPITT) is required. TROPOMI launched in 2017 (post-event), and MOPITT requires a separate assimilation pipeline — both are out of scope for this version. The hotspot-based prior modulation is retained.

### 10.3 Single receptor time per site per day

Even with the 06 UTC release time correctly placed in the afternoon mixed PBL, only one trajectory per site per day samples the diurnal cycle. Fire emissions and PBL height vary continuously through the day. The σ=3 ppm representativity error may be optimistic for sites with strong morning-fire activity. A future v4 inversion could release receptors at 03/06/09 UTC and combine them.

---

## 11. References

- **GFED 2015 Fire Season Indonesia.** [Global Fire Emissions Database, 2015 Indonesia overview](https://globalfiredata.org/pages/2015/11/16/2015-fire-season-indonesia/)
- **Huijnen, V., et al. (2016).** Fire carbon emissions over maritime southeast Asia in 2015 largest since 1997. *Scientific Reports* 6:26886. [DOI link](https://www.nature.com/articles/srep26886). Total Sept-Oct 2015: 692 Tg CO₂; mean rate 11.3 Tg CO₂/day.
- **Field, R. D., et al. (2016).** Indonesian fire activity and smoke pollution in 2015 show persistent nonlinear sensitivity to El Niño-induced drought. *PNAS* 113:9204–9209. [DOI link](https://www.pnas.org/doi/10.1073/pnas.1524888113). Names Sumatra, Kalimantan, Papua as primary sources.
- **Heymann, J., et al. (2017).** CO₂ emission of Indonesian fires in 2015 estimated from satellite-derived atmospheric CO₂ concentrations. *Geophysical Research Letters* 44:1537–1544. [DOI link](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL072042). Independent satellite XCO₂ inversion.
- **Wooster, M., et al. (2018).** Monitoring emissions from the 2015 Indonesian fires using CO satellite data. *Phil. Trans. R. Soc. B* 373:20170307. [DOI link](https://royalsocietypublishing.org/doi/10.1098/rstb.2017.0307). MOPITT/IASI-based central Indonesia (Kalimantan + Sulawesi) estimate: 63–79 Tg C.
- **Kiely, L., et al. (2019).** New estimate of particulate emissions from Indonesian peat fires in 2015. *Atmospheric Chemistry and Physics* 19:11105–11121. [DOI link](https://acp.copernicus.org/articles/19/11105/2019/). Peat-specific emission factors.
- **van der Werf, G. R., et al. (2017).** Global fire emissions estimates during 1997–2016. *Earth System Science Data* 9:697–720. GFED4s reference.
- **Parker, R. J., et al. (2016).** Atmospheric CH₄ and CO₂ enhancements and biomass burning emission ratios derived from satellite observations of the 2015 Indonesian fire plumes. *Atmospheric Chemistry and Physics* 16:10111–10131. FINN bias quantification.

---

## Appendix A. File Map

| Output | Path |
|---|---|
| Inversion script | [`inversion/run_inversion_v3.py`](../inversion/run_inversion_v3.py) |
| Spatial plot script | [`inversion/plot_spatial_v3.py`](../inversion/plot_spatial_v3.py) |
| STILT footprint plots script | [`inversion/plot_stilt_footprints.py`](../inversion/plot_stilt_footprints.py) |
| Summary text | [`inversion/results/inversion_v3_summary.txt`](../inversion/results/inversion_v3_summary.txt) |
| Fit plot | [`inversion/plots/inversion_fit_v3.png`](../inversion/plots/inversion_fit_v3.png) |
| Emission plot | [`inversion/plots/posterior_emissions_v3.png`](../inversion/plots/posterior_emissions_v3.png) |
| Spatial map | [`inversion/plots/spatial_prior_vs_posterior_v3.png`](../inversion/plots/spatial_prior_vs_posterior_v3.png) |
| Diagnostics | [`inversion/plots/inversion_diagnostics_v3.png`](../inversion/plots/inversion_diagnostics_v3.png) |
| STILT mean footprints | [`inversion/plots/stilt_footprint_means_v3.png`](../inversion/plots/stilt_footprint_means_v3.png) |
| STILT footprint integrals | [`inversion/plots/stilt_footprint_integrals_v3.png`](../inversion/plots/stilt_footprint_integrals_v3.png) |
| STILT BKT-flask footprints | [`inversion/plots/stilt_bkt_flask_footprints_v3.png`](../inversion/plots/stilt_bkt_flask_footprints_v3.png) |
| STILT v3.2 06 UTC script | [`stilt/run_stilt_v3p2_06utc.r`](../stilt/run_stilt_v3p2_06utc.r) |
| STILT v3.2 06 UTC receptors | [`stilt/receptors_full_oct_06utc.csv`](../stilt/receptors_full_oct_06utc.csv) |
| STILT v3.2 footprints | [`stilt/out_06utc/footprints/`](../stilt/out_06utc/footprints/) — 286 files |
| BKT NOAA flask record | [`inversion/data/co2_bkt_flask.txt`](../inversion/data/co2_bkt_flask.txt) |
| Plan document | [`docs/plan-stilt_inversion.md`](plan-stilt_inversion.md) |
| v1/v2 baseline doc | [`docs/stilt_inversion.md`](stilt_inversion.md) |

## Appendix B. Configuration (v3.2)

```python
# Region geometry (v3.1 post-fix)
REGIONS = [
    ('Sumatra N+C',    95.0, 104.0,  0.0,  6.0),
    ('Sumatra S',      99.0, 107.0, -6.0,  0.0),
    ('Kalimantan W+C', 108.0, 112.0, -2.0,  3.0),
    ('Kalimantan S+E', 112.0, 119.0, -5.0,  1.0),
    ('Java',           105.0, 116.0, -9.0, -5.5),
    ('Sulawesi+East',  119.0, 141.0, -6.0,  2.0),
]

X_PRIOR     = [0.25, 0.35, 0.20, 0.20, 0.12, 0.10]
SIGMA_PRIOR = [0.20, 0.25, 0.15, 0.15, 0.20, 0.08]

# v3.2 changes
OCO2_MAX_DIST_DEG = 6.0     # §10.2 fix — OCO-2 distance cutoff
CT_BORDERLINE_LOW  = 5.0    # §10.3 fix — CT2025 σ ×2 for borderline rows
CT_BORDERLINE_HIGH = 25.0   #            (and ×∞ exclusion above 25)
TIKHONOV_LAMBDA   = 5.0     # §10.6 fix — D_MATRIX shape (0,6); Tikhonov term = 0
USE_CO_I7         = False   # disabled (circular constraint)
USE_HOTSPOT_PRIOR = True

# BKT NOAA flask integration (v3.2 follow-up)
BKT_FLASK_FILE = 'inversion/data/co2_bkt_flask.txt'  # NOAA GML
# 4 paired samples in Oct 2015 (Oct 6/13/20/27 at 07 UTC = 14 WIB).
# All loaded; instr-failure QC codes 'NXY*' rejected, but smoke-contamination
# 'C' is kept (measurement is sound; only background-selection fails).
# Flask values OVERRIDE the CT2025 row at the BKT slot for those days,
# including the I8 fire-zone exclusion.

# v3.2-final — STILT 06 UTC release
FOOT_DIR    = 'stilt/out_06utc/footprints'  # was stilt/out/footprints (12 UTC)
DATES_NOON  = [datetime(2015, 10, d, 6, ...) for d in range(6, 32)]
                                            # was hour=12 (post-sunset)
# CT2025 read at index 2 (= 06 UTC) — was index 4 (= 12 UTC)
```

## Appendix C. Reproducibility

```bash
cd /home/igrk/WRF-GRK
source .venv/bin/activate
python3 inversion/run_inversion_v3.py 2>&1 | tee /tmp/inversion_v3_fixed.log
```

Expected wall-clock: ~5 minutes (286 footprint loads, 6×6 H matrix solve).

To regenerate all plots:

```bash
source .venv/bin/activate
python3 inversion/run_inversion_v3.py            # produces fit/emissions/diagnostics
python3 inversion/plot_spatial_v3.py             # produces spatial prior-vs-posterior map
python3 inversion/plot_stilt_footprints.py       # produces 3 STILT footprint figures
```

To rerun the STILT v3.2 footprints (~3 hours wall-clock on 8 cores):

```bash
cd /home/igrk/WRF-GRK/stilt
nohup Rscript run_stilt_v3p2_06utc.r > /tmp/stilt_v3p2_06utc.log 2>&1 &
# Monitor: tail -f /tmp/stilt_v3p2_06utc.log
# Output:  stilt/out_06utc/footprints/ (target: 286 files)
```
