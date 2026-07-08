# Plan: WRF-STILT-Inversion v3 — Full October 2015
## Indonesia Biomass Burning — IDN_BB_2015

**Document scope:** Detailed scientific and computational plan for inversion version 3,
covering the full October 2015 El Niño fire event with a real 11-site STILT network,
CO tracer joint constraint, and OCO-2 satellite assimilation.

**Status at writing:** STILT full-October run in progress (PID 3896721, Oct 8 of 26 days).

---

## Table of Contents

1. [Motivation and Deficiencies of v1/v2](#1-motivation-and-deficiencies-of-v1v2)
2. [Domain, Period, and Data Inventory](#2-domain-period-and-data-inventory)
3. [STILT Network Design (I6)](#3-stilt-network-design-i6)
4. [STILT Run Configuration](#4-stilt-run-configuration)
5. [Observation System](#5-observation-system)
6. [FINN Prior Emissions](#6-finn-prior-emissions)
7. [Bayesian Inversion Framework](#7-bayesian-inversion-framework)
8. [Improvements over v2 (I1–I8)](#8-improvements-over-v2-i1i8)
9. [Expected Results and Validation Targets](#9-expected-results-and-validation-targets)
10. [File Map and Execution Sequence](#10-file-map-and-execution-sequence)
11. [Known Risks and Mitigations](#11-known-risks-and-mitigations)

---

## 1. Motivation and Deficiencies of v1/v2

### 1.1 v1 Results

v1 (`run_inversion.py`) used 5 receptor sites × 6 days (Oct 24–29), 3 aggregated regions
(Sumatra, Kalimantan, Java/Sulawesi), and a simple analytical Bayesian solve. Key results:

| Region | α posterior | Tg CO₂ posterior | Physical plausibility |
|---|---:|---:|---|
| Sumatra | 0.091 | 2.4 Tg | Plausible but probably too low |
| Kalimantan | −0.010 | ≈ 0 Tg | **Unphysical** — fires did burn |
| Java/Sulawesi | 1.478 | 32.6 Tg | Sulawesi >> Kalimantan is wrong |

The Kalimantan near-zero result is a numerical artifact: FINN's prior emissions at Palangka Raya
are 100–600× too high relative to what CT2025 and OCO-2 observe. The optimizer collapses
α_Kalimantan to zero to eliminate the massive outliers at those receptor sites. The 5-site network
had no real footprint sampling Sulawesi or eastern Kalimantan, so those regions were poorly
constrained by IDW extrapolation.

### 1.2 v2 Results

v2 (`run_inversion_v2.py`) added: bounded optimization (α ≥ 0), per-sounding OCO-2 averaging
kernels, literature-informed priors, and virtual receptor sites interpolated by IDW. Key results:

| Region | α posterior | Tg CO₂ (6 days) |
|---|---:|---:|
| Sumatra S | 0.31 | ~5 Tg |
| Sumatra N+C | 0.12 | ~0.8 Tg |
| Java/Bali | 0.15 | ~0.4 Tg |
| Kalimantan W+C | 0.08 | ~1.2 Tg |
| **Kalimantan S+E** | **0.07** | **~0.9 Tg** |
| **Sulawesi+East** | **0.82** | **~6.6 Tg** |

**Sulawesi (6.6 Tg) >> Kalimantan S+E (0.9 Tg)** — physically wrong. The 2015 El Niño drove
intense peat fires in South and Central Kalimantan that were among the largest on record. GFED4.1s
attributes ~18–25 Tg CO₂ to Kalimantan and only ~2–4 Tg to Sulawesi for October 2015. The IDW
virtual sites placed over Sulawesi and east Kalimantan had footprints that were entirely fabricated
from the distance-weighted blend of the 5 western Indonesian sites — no real trajectory ever
reached Sulawesi in those footprints. The inversion assigned a high Sulawesi scaling factor
because the CT2025 "observations" at those virtual sites were high but the (fabricated) H matrix
was low, forcing a large α to compensate.

### 1.3 Root Cause

**The IDW virtual footprint approach has no scientific validity.** STILT backward trajectories
are physically driven by transport — they cannot be interpolated between sites because trajectories
from Sulawesi go through a completely different airflow regime than trajectories from Sumatra.
The fix is simple: **run real STILT trajectories at Sulawesi and east Kalimantan sites.** That
is the entire purpose of v3.

---

## 2. Domain, Period, and Data Inventory

### 2.1 Spatial Domain

| Parameter | Value |
|---|---|
| WRF domain | 27 km, 231 × 127 cells, 90–145 °E, 15 °S–15 °N |
| STILT footprint grid | 0.25°, 220 × 120 cells, 90.125–144.875 °E, 14.875 °S–14.875 °N |
| FINN emission grid | 0.1° × 0.1°, same geographic extent |
| CT2025 grid | 3° × 2° global, interpolated to site locations |
| OCO-2 assimilation grid | 1° × 1° cells, Indonesia only |

### 2.2 Temporal Period

| Layer | Period | Notes |
|---|---|---|
| WRF-GHG simulation | Oct 1–31, 2015 | 31 daily wrfout files at 01:00 UTC |
| ARL met files | Oct 1–31, 2015 | 31 × 240 MB files in `stilt/met/` — **all complete** |
| STILT receptors | Oct 6–31, 2015 | Oct 1–5 excluded: STILT needs 5 ARL days prior to receptor |
| CT2025 reanalysis | Oct 1–31, 2015 | **All 31 files already downloaded** |
| FINN emissions | Full year 2015 | Daily 0.1° CO₂, CO, CH₄ files available |
| OCO-2 granules | Oct 7–26, 2015 | 20 granules **downloaded** (server has none for Oct 27–31) |
| MODIS hotspots | Full year 2015 | `rawdata/hotspot/archived_hotspot_idn.csv` |

**Why not Oct 1–5?** The STILT backward trajectory runs 120 hours (5 days) before the receptor
time. A receptor on Oct 1 would need ARL met data back to Sep 26. We only have WRF-derived ARL
from Oct 1 onward. Oct 6 is the earliest receptor date for which the full 5-day backward
trajectory stays within the available ARL window.

**Why not Oct 27–31 OCO-2?** v11.1r archive on NASA GES DISC contains no granules for these
dates over Indonesia — likely cloud/geometry coverage gaps. The 20 granules for Oct 7–26 provide
good temporal coverage of the peak fire period.

### 2.3 Data Status

| Dataset | Path | Count | Status |
|---|---|---|---|
| WRF wrfout | `simulations/IDN_BB_2015/output/` | 31 files | ✅ |
| ARL met | `stilt/met/*.arl` | 31 files | ✅ |
| CT2025 | `rawdata/carbontracker/CT2025.molefrac_glb3x2_2015-10-*.nc` | 31 files | ✅ |
| FINN CO₂ | `rawdata/finn_fire/emissions-finnv2.5modvrs_CO2_bb_*.nc` | 1 file (full year) | ✅ |
| FINN CO | `rawdata/finn_fire/emissions-finnv2.5modvrs_CO_bb_*.nc` | 1 file (full year) | ✅ |
| OCO-2 | `inversion/data/oco2_LtCO2_1510*.nc4` | 20 files | ✅ |
| MODIS hotspots | `rawdata/hotspot/archived_hotspot_idn.csv` | ~80,000 rows | ✅ |

---

## 3. STILT Network Design (I6)

### 3.1 Original 5-Site Network (v1/v2)

| Site | Code | Lon | Lat | Region sampled |
|---|---|---|---|---|
| Bukit Kototabang | BKT | 100.32 °E | 0.20 °S | West Sumatra |
| Jambi | JAM | 103.61 °E | 1.61 °S | East Sumatra |
| Jakarta | JKT | 106.85 °E | 6.21 °S | Java (background) |
| Pontianak | PON | 109.33 °E | 0.02 °S | West Kalimantan |
| Palangka Raya | PLK | 113.94 °E | 2.16 °S | Central Kalimantan |

**Network gap:** All 5 sites are in the western half of the domain (100–114 °E). No site covers:
- East Kalimantan (117–118 °E) — largest peat fire area in 2015
- South Kalimantan (114–116 °E, south) — second largest
- Sulawesi (119–125 °E) — moderate fire activity (~2–4 Tg)
- Maluku/Papua (130–141 °E) — minor fires

### 3.2 Six New Sites Added (v3)

Selected to fill the geographic gap and provide direct footprint coverage of the
East Kalimantan and Sulawesi fire zones:

| Site | Code | Lon | Lat | Scientific rationale |
|---|---|---|---|---|
| Pekanbaru | PKU | 101.45 °E | 0.53 °N | Riau province — central to Sumatra fire region; BMKG monitoring station |
| Palembang | PLM | 104.75 °E | 2.99 °S | South Sumatra peat fires; second most fire-affected Sumatra province |
| Banjarmasin | BDJ | 114.59 °E | 3.32 °S | South Kalimantan; IATA airport with long-term records; upwind of major 2015 peat fires |
| Samarinda | SRI | 117.15 °E | 0.50 °S | East Kalimantan capital; BMKG station directly in the 2015 East Kalimantan fire zone |
| Balikpapan | BPN | 116.85 °E | 1.27 °S | East Kalimantan coast; ~80 km south of Samarinda; additional east Kalimantan constraint |
| Makassar | UPG | 119.42 °E | 5.14 °S | South Sulawesi capital; allows the inversion to independently constrain Sulawesi |

### 3.3 Combined 11-Site Network

The full network spans 100–119 °E with balanced east-west coverage. The 11 sites
provide real backward trajectory footprints for all six inversion regions, eliminating
the need for IDW virtual sites entirely.

**Simulation count:** 11 sites × 26 days (Oct 6–31) = **286 simulations**  
**Plus existing:** 30 (original, Oct 24–29) + 66 (new sites + overlap, Oct 24–29) = 96 already done  
**Net new simulations:** 286 − 66 (Oct 24–29 already done) = 220 new runs  

---

## 4. STILT Run Configuration

### 4.1 Key Parameters

| Parameter | Value | Rationale |
|---|---|---|
| `n_hours` | −120 (5 days backward) | Captures regional-scale fire influence; longer would leave WRF domain |
| `numpar` | 200 particles | Balance between precision and computation time |
| `n_cores` | 8 | All available on pemodelanKU |
| Receptor height (`zagl`) | 10 m AGL | Surface in-situ equivalent |
| Receptor time | 12:00 UTC (20:00 WIB) | Afternoon well-mixed boundary layer |
| Footprint grid | 0.25°, 220 × 120 | Matches FINN 0.1° downsampled, OCO-2 1° oversample |
| Met format | `%Y-%m-%d.arl` | Daily WRF-derived ARL files |
| `met_path` | `stilt/met/` | |

### 4.2 STILT Scripts

| Script | Receptor file | Simulations | Status |
|---|---|---|---|
| `stilt/run_stilt.r` | `receptors.csv` (30 rows) | 5 sites × Oct 24–29 | ✅ Done |
| `stilt/run_stilt_new.r` | `receptors_new.csv` (36 rows) | 6 new sites × Oct 24–29 | ✅ Done |
| `stilt/run_stilt_full_oct.r` | `receptors_full_oct.csv` (286 rows) | 11 sites × Oct 6–31 | 🔄 Running |

`run_stilt_full_oct.r` will skip simulations for which `by-id/` output already exists, so the
66 Oct 24–29 footprints already computed will not be rerun.

### 4.3 Footprint Output

Each simulation writes:
- `stilt/out/footprints/{timestamp}_{lon}_{lat}_{zagl}_foot.nc` — NetCDF, shape (1, 120, 220),
  units `ppm / (µmol m⁻² s⁻¹)`, variable `foot`
- `stilt/out/by-id/{simulation_id}/` — particle trajectory files

Expected final count: **286 footprints** (220 new + 66 existing).

---

## 5. Observation System

### 5.1 CarbonTracker CT2025 Pseudo-Observations

**Source:** NOAA CarbonTracker 2025 global reanalysis  
**Files:** `rawdata/carbontracker/CT2025.molefrac_glb3x2_2015-10-{01..31}.nc`  
**Variable:** `co2` [mol mol⁻¹], 34 pressure levels, 8 time steps day⁻¹ (3-hourly)  
**Level used:** Level 0 (~981 hPa surface)  
**Time step selected:** Index 3 = 10:30 UTC ≈ afternoon well-mixed local time

CT2025 assimilated real NOAA flask, aircraft (HIPPO, ORCAS), and tower observations into an
ensemble Kalman filter. For Indonesia, the nearest real observations are from BKT (Sumatra) and
sparse Pacific aircraft profiles. The CT2025 surface field in Indonesia is therefore primarily
model-inferred with ~3° × 2° resolution, providing a physically consistent but spatially smoothed
representation of fire CO₂ enhancements.

**In v3, CT2025 provides:**
- 286 "observed" daily noon CO₂ values (one per receptor site per day)
- Background reference for fire enhancement calculation
- Independent constraint for days when OCO-2 has no overpass

**Uncertainty:** σ = 3.04 ppm (representativity), increased to 5 ppm for days/sites where
CT2025 itself shows fire influence (enhancement > 3 ppm) to avoid double-counting with OCO-2.

**BKT independence fix (I4):** CT2025 assimilated the BKT flask. On Oct 27 (the flask day),
the BKT CT2025 observation is excluded from the inversion and replaced by the real flask value.

### 5.2 BKT Flask (Real In-Situ)

**Source:** NOAA GML Carbon Cycle Cooperative Global Air Sampling Network  
**Site:** Bukit Kototabang, West Sumatra (100.32 °E, 0.20 °S, 864 m ASL)  
**Observation:** Oct 27, 2015, 07:00 UTC — two flask samples, averaged: **409.51 ppm**  
**Flag:** `C` (fire/smoke contamination, per NOAA QC)  
**Background subtracted:** 409.51 − 399.1 = **10.41 ppm fire enhancement**  
**Uncertainty:** σ = √(0.07² + 0.5² + 3.0²) = 3.04 ppm (instrument + background + representativity)

This is the single most valuable observation in the system — a real measurement of the
fire plume passing over Sumatra on the day of peak fire activity. Its weight in the inversion
is commensurate with its uncertainty (3 ppm), comparable to a single CT2025 pseudo-observation.

### 5.3 OCO-2 Column CO₂ (I3)

**Source:** OCO-2 Lite FP v11.1r (B11100Ar), NASA GES DISC  
**Files:** `inversion/data/oco2_LtCO2_1510{07..26}_B11100Ar_*.nc4` (20 granules)  
**Variable:** `xco2` [ppm], quality flag `xco2_quality_flag = 0`  
**Averaging kernel:** `xco2_averaging_kernel` (20 pressure levels, per sounding)  
**Footprint:** `warn_level < 15`, land or ocean glint  

**Processing pipeline:**
1. Filter to Indonesia bounding box (90–145 °E, 15 °S–15 °N) with quality flag = 0
2. Compute tropical background: median XCO₂ of soundings over ocean (15 °S–15 °N, west of 60 °E)
   with no fire influence → `oco2_bg` ≈ 398.6–399.1 ppm depending on date
3. Compute fire enhancement: `xco2_fire = xco2 - oco2_bg`
4. Grid to 1° × 1° cells, minimum 3 soundings per cell
5. For each 1° cell, construct H_cell via IDW from all 11 real STILT footprints (not fabricated)

**Column-to-surface scaling (I3 improvement from v1):**
Rather than a fixed factor of 0.25, v3 uses per-sounding pressure-weighted averaging kernels:
$$A_{eff} = \frac{\sum_l AK_l \cdot \Delta p_l \cdot \mathbf{1}[\text{PBL}]}{\sum_l \Delta p_l \cdot \mathbf{1}[\text{PBL}]}$$
where the sum is over pressure levels within the planetary boundary layer (below 850 hPa).
This correctly accounts for the fact that fire plumes confined to the PBL have reduced XCO₂
sensitivity compared to a surface flask.

**Fire zone exclusion (I8):** OCO-2 grid cells with `xco2_fire > 25 ppm` are excluded from
the inversion. Such anomalously high XCO₂ values arise from dense smoke aerosol that biases
the OCO-2 retrieval toward lower altitudes and inflates XCO₂. This threshold (25 ppm) is
conservative — OCO-2 XCO₂ biases > 3 ppm are well documented in thick smoke, and 25 ppm
is well above the typical 1–7 ppm fire signal seen in clean retrievals.

**Expected observation count:** ~200–400 1° grid cells over 20 days (variable by cloud cover)

### 5.4 CO Tracer Joint Constraint (I7)

**Rationale:** CO is co-emitted with CO₂ in biomass burning. The CO:CO₂ molar emission ratio
for Indonesian peat fires is approximately 0.08–0.12 mol CO / mol CO₂ (smoldering combustion).
If the inversion finds α_Kalimantan ≈ 0 (implying no fire CO₂), but WRF simulates very high
CO_BBU at the Kalimantan sites, this is physically inconsistent and suggests model transport
error rather than zero emissions. The CO constraint adds a model-data consistency check.

**Implementation:**
- WRF `CO_BBU` (fire CO tracer) is sampled at all 11 receptor sites for all 26 days
- FINN CO flux is used to build H_CO with the same footprints as H_CO₂
- Predicted CO enhancement: `Δco_model = H_CO · α` (assuming same α scales both CO and CO₂)
- Observed CO: WRF `CO_BBU` at receptor (used as a proxy for "observed" CO in absence of
  real CO network data in Indonesia)
- Additional constraint equation: `Δco_wrf ≈ H_CO · α` with σ_CO = 50 ppb

**Note:** This is a model-to-model constraint, not an independent satellite CO observation
(TROPOMI CO launched 2017; MOPITT CO not used because the vertical sensitivity profile
and retrieval operator would require additional implementation). The CO constraint primarily
prevents the inversion from finding α = 0 in regions where WRF shows genuine CO transport,
acting as a physical consistency regularizer.

### 5.5 MODIS Hotspot Prior Modulation (I7 — hotspot component)

MODIS active fire detections (Terra + Aqua, Collection 6.1) from `rawdata/hotspot/archived_hotspot_idn.csv`
are aggregated by region and day. The hotspot counts modulate the prior uncertainty `σ_prior`:

| Hotspot count (region, Oct 6–31) | σ_prior multiplier | Effect |
|---|---|---|
| > 200 | × 0.70 | Tighten prior — more confident fires are real |
| 50–200 | × 1.00 | No change |
| < 50 | × 1.20 | Loosen prior — allow larger correction if few fires detected |

**Reasoning:** FINN derives fire locations from MODIS FRP. Regions with many MODIS detections
are genuine high-fire areas and the inversion should not push α far from the prior without
strong observational evidence. Regions with few detections either have genuinely low fire activity
(prior is probably low anyway) or are consistently cloud-obscured (prior is uncertain).
The hotspot modulation is applied per-region per-day before assembling the full S_prior matrix.

---

## 6. FINN Prior Emissions

### 6.1 Files

| Tracer | File |
|---|---|
| CO₂ | `rawdata/finn_fire/emissions-finnv2.5modvrs_CO2_bb_surface_daily_20150101-20151231_0.1x0.1.nc` |
| CO | `rawdata/finn_fire/emissions-finnv2.5modvrs_CO_bb_surface_daily_20150101-20151231_0.1x0.1.nc` |
| CH₄ | `rawdata/finn_fire/emissions-finnv2.5modvrs_CH4_bb_surface_daily_20150101-20151231_0.1x0.1.nc` |

**Variable:** `fire_modisviirs_CO2` [molecules cm⁻² s⁻¹], dimensions `(time, lat, lon)` with
`time` = day-of-year (0-indexed, so Oct 1 = index 273 in 2015).

### 6.2 Inversion Regions (6 regions)

| Index | Region | Lon bounds | Lat bounds | Prior α | Prior σ |
|---|---|---|---|---|---|
| 0 | Sumatra South | 104–107 °E | 5 °S–2 °N | 0.40 | 0.30 |
| 1 | Sumatra North+Central | 95–104 °E | 2 °S–5 °N | 0.25 | 0.20 |
| 2 | Java/Bali | 105–116 °E | 9 °S–6 °S | 0.15 | 0.30 |
| 3 | Kalimantan West+Central | 108–116 °E | 4 °S–2 °N | 0.20 | 0.30 |
| 4 | Kalimantan South+East | 114–120 °E | 5 °S–2 °N | 0.80 | 0.40 |
| 5 | Sulawesi+East | 118–135 °E | 8 °S–3 °N | 0.50 | 0.25 |

**Prior α rationale:** Literature-based from GFED4.1s vs FINN ratios for Indonesian peat fires
during 2015 (Huijnen et al. 2016; van der Werf et al. 2017). FINN consistently overestimates
peat fire CO₂ by 3–10× relative to GFED for this event. Kalimantan S+E has the highest prior
(0.80) because the site network (Samarinda, Balikpapan) now provides direct constraint, reducing
dependence on the prior. Sulawesi prior (0.50) is intermediate — real fires exist but FINN may
be overestimating as it often does for savanna-type fires.

### 6.3 Prior Total Emissions (Oct 6–31, 26 days)

Computed as `sum(FINN_flux[day, region_mask] * dA_cm2) * 86400 / N_A * 44e-3 * 1e-9 [Tg CO₂]`
for each of the 26 days. Approximate totals:

| Region | Approx FINN prior (26-day, Tg CO₂) | Approx posterior target |
|---|---|---|
| Sumatra South | ~60 Tg | ~24 Tg |
| Sumatra N+C | ~25 Tg | ~6 Tg |
| Java/Bali | ~5 Tg | ~0.75 Tg |
| Kalimantan W+C | ~80 Tg | ~16 Tg |
| Kalimantan S+E | ~150 Tg | ~120 Tg |
| Sulawesi+East | ~35 Tg | ~17 Tg |
| **Indonesia total** | **~355 Tg** | **~184 Tg** |

*Posterior targets are approximate, based on applying prior α values and GFED4.1s literature.*

---

## 7. Bayesian Inversion Framework

### 7.1 State Vector

$$\boldsymbol{\alpha} = [\alpha_0, \alpha_1, \alpha_2, \alpha_3, \alpha_4, \alpha_5]^\top \in \mathbb{R}^6$$

Each $\alpha_k$ is a dimensionless regional scaling factor applied to the FINN daily emissions
in region $k$. Posterior $\alpha_k$ is constrained to $[0, 10]$ via SLSQP bounded optimization.

### 7.2 Observation Vector

$$\mathbf{y} \in \mathbb{R}^{n_{obs}}$$

assembled in order:
1. CT2025 pseudo-observations: 286 values (11 sites × 26 days)
2. BKT flask: 1 value (Oct 27)
3. OCO-2 1° grid cells: ~200–400 values (variable by day/coverage)
4. CO consistency: 286 values (one per receptor-day, as WRF CO_BBU)

**Approximate total:** ~800–1000 observations

### 7.3 H Matrix (Jacobian)

Shape: $(n_{obs}, 6)$. Each row $i$ gives the sensitivity of observation $i$ to each regional
scaling factor.

**For CT2025 / flask (surface receptor):**
$$H_{ik} = \sum_{(x,y) \in \text{region}_k} F_i(x,y) \cdot E_k(x,y) \cdot \Delta A$$
where $F_i$ is the STILT footprint [ppm / (µmol m⁻² s⁻¹)] for receptor $i$, $E_k$ is FINN flux
[µmol m⁻² s⁻¹] in region $k$, and $\Delta A$ is grid cell area [m²].

**For OCO-2 grid cells (column):**
$$H_{ik}^{OCO2} = A_{eff} \cdot \frac{\sum_{s \in \text{cell}} w_s \cdot H_{s,k}^{STILT}}{\sum_{s} w_s}$$
where $w_s = 1/d_s^2$ (IDW from all 11 real STILT sites), $A_{eff}$ is the effective PBL
averaging kernel weight (computed per-granule), and $H_{s,k}^{STILT}$ is the surface H value
for the nearest STILT receptor.

**For CO consistency:**
$$H_{ik}^{CO} = r_{CO:CO2} \cdot H_{ik}$$
where $r_{CO:CO2}$ is the FINN molar CO:CO₂ emission ratio for region $k$ (computed daily
from the two FINN files).

### 7.4 Error Covariance

$$\mathbf{S}_{obs} = \text{diag}(\sigma_1^2, \ldots, \sigma_{n_{obs}}^2)$$

Off-diagonal terms are neglected (observation errors treated as independent).

| Observation type | σ |
|---|---|
| CT2025 (low-fire cell, enhancement < 1 ppm) | 3.0 ppm |
| CT2025 (moderate fire, 1–5 ppm) | 5.0 ppm |
| BKT flask | 3.04 ppm |
| OCO-2 1° cell | $\sqrt{\bar{\sigma}_{XCO2}^2 / N_{cell} + \sigma_{bg}^2 + 0.3^2}$ ≈ 0.3–0.9 ppm |
| CO consistency | 50 ppb |

$$\mathbf{S}_{prior} = \text{diag}(\sigma_{prior,0}^2, \ldots, \sigma_{prior,5}^2)$$
with $\sigma_{prior,k}$ from the table in §6.2, modulated by MODIS hotspot counts (§5.5).

### 7.5 Posterior Solution (SLSQP)

Minimise the cost function subject to $\alpha_k \geq 0$:

$$J(\boldsymbol{\alpha}) = \frac{1}{2}(\mathbf{y} - \mathbf{H}\boldsymbol{\alpha})^\top \mathbf{S}_{obs}^{-1} (\mathbf{y} - \mathbf{H}\boldsymbol{\alpha}) + \frac{1}{2}(\boldsymbol{\alpha} - \boldsymbol{\alpha}_0)^\top \mathbf{S}_{prior}^{-1} (\boldsymbol{\alpha} - \boldsymbol{\alpha}_0) + \lambda \|\boldsymbol{\alpha}\|^2$$

where $\lambda = 5.0$ is the Tikhonov regularization parameter. The additional Tikhonov term
prevents runaway scaling factors when a region has very low H values (few footprints) but
non-zero observations.

**Analytical posterior covariance** (for uncertainty reporting):
$$\mathbf{S}_{post} = \left(\mathbf{H}^\top \mathbf{S}_{obs}^{-1} \mathbf{H} + \mathbf{S}_{prior}^{-1} + \lambda \mathbf{I}\right)^{-1}$$

---

## 8. Improvements over v2 (I1–I8)

| ID | Improvement | v2 status | v3 status |
|---|---|---|---|
| **I1** | Bounded optimization (α ≥ 0, SLSQP) | ✅ Added in v2 | ✅ Retained |
| **I2** | Literature-informed 6-region prior α | ✅ Added in v2 | ✅ Retained + updated |
| **I3** | OCO-2 per-sounding averaging kernel | ✅ Added in v2 | ✅ Retained |
| **I4** | CT2025/BKT independence fix | ✅ Added in v2 | ✅ Retained |
| **I5** | Tikhonov regularization (λ=5) | ✅ Added in v2 | ✅ Retained |
| **I6** | **Real STILT at 6 Kalimantan/Sulawesi sites** | ❌ IDW virtual | ✅ **New in v3** |
| **I7** | **CO tracer + MODIS hotspot constraints** | ❌ Not present | ✅ **New in v3** |
| **I8** | Fire-zone OCO-2 exclusion (XCO₂ > 25 ppm) | ❌ Not present | ✅ **New in v3** |

### I6 Detail

The six new sites (Pekanbaru, Palembang, Banjarmasin, Samarinda, Balikpapan, Makassar) were
chosen to provide direct footprint coverage of each inversion region. Specifically:
- Samarinda + Balikpapan directly sample Kalimantan S+E (region 4) — the region with the
  largest FINN prior and the most physically important correction
- Makassar directly samples Sulawesi+East (region 5) — previously only constrained by IDW
- Banjarmasin samples the boundary between Kalimantan W+C and S+E

### I7 Detail

The CO joint constraint exploits the fact that CO and CO₂ are co-emitted in fixed ratios by
fire. A solution that sets α_Kalimantan = 0 while WRF CO_BBU shows 100+ ppb CO enhancement
at Samarinda is physically inconsistent — the inversion should not be free to choose it.
The CO constraint rows add equations of the form `observed_CO ≈ H_CO · α`, which strongly
penalise solutions where fire CO₂ is near zero but WRF CO is high.

The MODIS hotspot modulation addresses a subtler issue: regions with dense hotspot detection
(e.g., South Sumatra in 2015 with 1,200+ fires per day) are clearly real fire zones and
the prior should not be widened unnecessarily. Conversely, regions where MODIS barely detected
fires may have FINN priors that are more uncertain — the prior uncertainty is widened to allow
the observations to pull the scaling factor lower.

---

## 9. Expected Results and Validation Targets

### 9.1 Scientific Expectations

Based on GFED4.1s (van der Werf et al. 2017) and previous studies of the 2015 Indonesian fires
(Huijnen et al. 2016, Field et al. 2016, Yin et al. 2016), the v3 inversion should produce:

| Region | Expected posterior α | Expected Tg CO₂ (26-day) | Literature reference |
|---|---|---|---|
| Sumatra South | 0.30–0.45 | 18–27 Tg | GFED4.1s ~22 Tg |
| Sumatra N+C | 0.15–0.30 | 4–7 Tg | GFED4.1s ~6 Tg |
| Java/Bali | 0.10–0.20 | 0.5–1.0 Tg | GFED4.1s ~0.8 Tg |
| Kalimantan W+C | 0.15–0.35 | 12–28 Tg | GFED4.1s ~20 Tg |
| Kalimantan S+E | 0.40–0.80 | 60–120 Tg | GFED4.1s ~70–100 Tg |
| Sulawesi+East | 0.20–0.50 | 7–17 Tg | GFED4.1s ~12 Tg |
| **Total** | | **~100–180 Tg** | **GFED4.1s ~130 Tg** |

**Critical physical check:** Kalimantan S+E must be >> Sulawesi+East in the posterior.
The 2015 event was driven by extreme peat fires in South and Central Kalimantan (the
"year without rain" drought; SPI < −2.5 for the whole season). Any posterior where
Sulawesi > Kalimantan signals a structural problem in the inversion setup.

### 9.2 Diagnostic Plots Produced

The inversion script generates these outputs in `inversion/plots/`:

1. **`inversion_fit_v3.png`** — 4-panel: (a) prior vs observed Δc scatter, (b) posterior vs
   observed scatter, (c) regional bar chart (prior/posterior Tg), (d) 2D map of posterior
   emission spatial distribution
2. **`posterior_emissions_v3.png`** — Indonesia map with 6 regions colored by posterior Tg CO₂,
   overlaid with MODIS active fire count for Oct 2015
3. **`inversion_diagnostics_v3.png`** — uncertainty reduction by region, χ² decomposition
   by observation type
4. **`footprint_grid_v3.png`** — 11×26 grid of all footprints (too large to show individually;
   binned to regional means)

### 9.3 Validation Checks

Before trusting results, these sanity checks must pass:

1. **χ² reduction:** Posterior χ² < prior χ² / 10 (inversion is actually improving fit)
2. **Non-negativity:** All posterior α ≥ 0 (guaranteed by SLSQP bounds)
3. **Kalimantan >> Sulawesi:** α₄ × FINN_KalE > α₅ × FINN_Sulawesi
4. **BKT posterior fit:** Posterior Δc at BKT Oct 27 within 2σ of 10.4 ppm observation
5. **OCO-2 residuals:** Mean OCO-2 residual < 1 ppm (no systematic bias)
6. **Prior uncertainty usage:** Posterior σ_k < prior σ_k for all regions (inversion is
   informative, not just returning prior)

---

## 10. File Map and Execution Sequence

### 10.1 Directory Structure

```
WRF-GRK/
├── rawdata/
│   ├── carbontracker/CT2025.molefrac_glb3x2_2015-10-*.nc  (31 files)
│   ├── finn_fire/emissions-finnv2.5modvrs_CO2_bb_*.nc
│   ├── finn_fire/emissions-finnv2.5modvrs_CO_bb_*.nc
│   └── hotspot/archived_hotspot_idn.csv
├── stilt/
│   ├── met/
│   │   ├── 2015-10-{01..31}.arl   (31 ARL files, all complete)
│   │   └── wrfout_links/          (hard-linked colon-free wrfout .nc)
│   ├── out/
│   │   ├── footprints/            (target: 286 *_foot.nc)
│   │   └── by-id/                 (particle outputs, ~286 dirs)
│   ├── receptors_full_oct.csv     (286 rows, Oct 6–31 × 11 sites)
│   ├── run_stilt_full_oct.r       (RUNNING — PID 3896721)
│   └── exe/arw2arl                (WRF→ARL converter)
├── inversion/
│   ├── data/oco2_LtCO2_1510*.nc4  (20 files, Oct 7–26)
│   ├── run_inversion_v3.py        (complete, syntax validated)
│   └── plots/                     (output directory)
├── simulations/IDN_BB_2015/output/wrfout_d01_2015-10-*.nc  (31 files)
└── docs/
    ├── stilt_inversion.md          (v1/v2 results documentation)
    └── plan-stilt_inversion.md     (this file)
```

### 10.2 Execution Sequence

```bash
# Step 1: Wait for STILT to complete (~3 hours remaining)
# Monitor: tail -f /tmp/stilt_full_oct.log
# Target: ls stilt/out/footprints/ | wc -l → 286

# Step 2: Verify all 286 footprints are present and valid
source .venv/bin/activate
python3 - <<'EOF'
import glob, netCDF4 as nc
fps = sorted(glob.glob('stilt/out/footprints/*_foot.nc'))
print(f'{len(fps)} footprints found')
bad = []
for f in fps:
    try:
        with nc.Dataset(f) as ds:
            v = ds['foot'][:]
            if v.max() == 0: bad.append(f)
    except Exception as e:
        bad.append(f'{f}: {e}')
print(f'{len(bad)} bad footprints: {bad[:5]}')
EOF

# Step 3: Run the v3 inversion
source .venv/bin/activate
python3 inversion/run_inversion_v3.py 2>&1 | tee /tmp/inversion_v3.log

# Step 4: Review results
# Check: inversion/plots/inversion_fit_v3.png
# Check: /tmp/inversion_v3.log for posterior α values
# Verify: Kalimantan S+E posterior Tg > Sulawesi posterior Tg
```

### 10.3 Estimated Run Times

| Step | Estimated time |
|---|---|
| STILT full October (220 remaining simulations) | ~3–4 hours (8 cores) |
| Footprint validation | ~2 minutes |
| v3 inversion (including H matrix build) | ~5–10 minutes |
| Plot generation | ~2 minutes |

---

## 11. Known Risks and Mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| STILT fails for Oct 6–18 (early days, new ARL files) | Medium | ARL files tested and confirmed 240 MB each; arw2arl ran clean |
| OCO-2 cells < 3 soundings (thin coverage) | Low | Only 20 of 26 days have OCO-2; the 6 days without are covered by CT2025 |
| Kalimantan S+E still comes out low due to WRF transport error | Medium | CO constraint (I7) prevents α→0 where WRF shows CO; if still low, increase CO weight σ_CO→20 ppb |
| BKT posterior underprediction (same as v2) | High | Known issue — single flask outweighed by 286 CT2025 values; will report residual and note limitation |
| FINN CO and CO₂ unit mismatch causing wrong CO:CO₂ ratio | Low | Both files are `molecules cm⁻² s⁻¹`; ratio computed directly from files, not hardcoded |
| STILT R process killed (SIGPIPE from log monitoring) | Low | Running with `nohup`, no pipe to monitoring process |
| Memory error in H matrix (286×6 × large grid) | Very low | H matrix is 286×6 floats (~14 KB), not a memory concern |
| OCO-2 fire zone threshold too tight (25 ppm) | Low | Only affects direct-overpass cells; background cells retained |
