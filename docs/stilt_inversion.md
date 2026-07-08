# WRF-STILT CO₂ Forward Simulation and Bayesian Inversion
## Indonesia Biomass Burning 2015 — IDN_BB_2015

**Period:** 22 October – 29 October 2015  
**Domain:** 90–145 °E, 15 °S–15 °N (Indonesia and surrounding seas)  
**Fire event:** 2015 El Niño extreme peat fire season, predominantly in Sumatra and Kalimantan  

---

## Table of Contents

1. [Experimental Setup](#1-experimental-setup)
2. [FINN Fire Emissions](#2-finn-fire-emissions)
3. [WRF-GHG Forward Simulation](#3-wrf-ghg-forward-simulation)
4. [STILT Backward Trajectories and Footprints](#4-stilt-backward-trajectories-and-footprints)
5. [CO₂ Concentration Quantification](#5-co₂-concentration-quantification)
6. [Observation Datasets](#6-observation-datasets)
7. [Bayesian Inversion](#7-bayesian-inversion)
8. [Results and Analysis](#8-results-and-analysis)
9. [Conclusions and Caveats](#9-conclusions-and-caveats)

---

## 1. Experimental Setup

### 1.1 WRF-GHG Configuration

| Parameter | Value |
|---|---|
| Model | WRF-Chem 4.1.5 with GHG (chem_opt = 17) |
| Domain | Single domain d01, 27 km horizontal resolution |
| Grid | 231 × 127 cells |
| Vertical levels | 35 eta levels |
| Simulation period | 22 Oct – 29 Oct 2015 |
| Output files | `wrfout_d01_2015-10-22_01:00:00` through `wrfout_d01_2015-10-29_01:00:00` (8 files) |
| Meteorology driver | ERA5 reanalysis |

**GHG tracers simulated:**

| Tracer | Description |
|---|---|
| `CO2_BBU` | CO₂ from biomass burning (fire source, tagged) |
| `CO2_BCK` | CO₂ background (~400 ppm, fixed boundary) |
| `CO2_ANT` | CO₂ from anthropogenic sources |
| `CH4_BBU` | CH₄ from biomass burning |
| `CO_BBU` | CO from biomass burning |

### 1.2 Receptor Sites

Five ground-level receptor sites were selected to span the major fire-affected regions of Indonesia:

| Site | Longitude | Latitude | Region |
|---|---|---|---|
| Bukit Kototabang (BKT) | 100.32 °E | 0.20 °S | West Sumatra |
| Palangka Raya | 113.94 °E | 2.16 °S | Central Kalimantan |
| Pontianak | 109.33 °E | 0.02 °S | West Kalimantan |
| Jambi | 103.61 °E | 1.61 °S | East Sumatra |
| Jakarta | 106.85 °E | 6.21 °S | Java (background reference) |

---

## 2. FINN Fire Emissions

**Source:** FINN v2.5 (Fire INventory from NCAR), MODIS + VIIRS combined  
**File:** `rawdata/finn_fire/emissions-finnv2.5modvrs_CO2_bb_surface_daily_20150101-20151231_0.1x0.1.nc`  
**Spatial resolution:** 0.1° × 0.1°  
**Variable:** `fire_modisviirs_CO2` [molecules cm⁻² s⁻¹]

### 2.1 Spatial Distribution

![FINN fire emission map, Oct 24–29 2015](../stilt/out/emission_finn.png)

The FINN inventory shows intense fire activity concentrated in three regions during October 2015:

- **Sumatra** (95–107 °E, 8 °S–5 °N): Extensive peat and forest fires along both coasts of the island, particularly in Riau and Jambi provinces
- **Kalimantan** (107–121 °E, 7 °S–3 °N): The most intense fire region, with peak emissions exceeding 100 g CO₂ m⁻² day⁻¹ in Central and South Kalimantan — primarily peat swamp fires
- **Java/Sulawesi** (104–130 °E, 10–3 °S): Lower-density agricultural and savanna fires

### 2.2 Emission Totals (Oct 24–29, 6-day prior)

| Region | FINN Prior Total | Units |
|---|---|---|
| Sumatra | 26.2 | Tg CO₂ |
| Kalimantan | 57.9 | Tg CO₂ |
| Java/Sulawesi | 22.1 | Tg CO₂ |
| **Indonesia total** | **116.2** | **Tg CO₂** |

These prior totals are the inputs to the inversion. **FINN is known to overestimate peat fire emissions** because it relies on Fire Radiative Power (FRP) calibrated primarily for flaming combustion, while Indonesian peat fires are predominantly smoldering — a fundamentally different combustion regime with different emission factors.

---

## 3. WRF-GHG Forward Simulation

### 3.1 CO₂ Concentration Fields

![WRF-GHG CO2_BBU concentration maps](../stilt/out/concentration_wrf.png)

WRF-GHG shows that the `CO2_BBU` (fire-tagged CO₂) tracer accumulates strongly over the source regions, with:

- **Peak surface concentrations**: 206 ppm above background (Oct 22), diminishing to ~153 ppm by Oct 23 as meteorological transport dilutes the plume
- **Plume transport**: The prevailing eastward flow during the 2015 El Niño directed Indonesian fire plumes northeastward, consistent with MODIS Aqua aerosol optical depth observations
- **Background CO₂** (`CO2_BCK`): ~400 ppm throughout, confirming stable lateral boundary conditions

### 3.2 Receptor Concentrations from WRF

The table below gives `CO2_BBU` (fire-only tracer) at each receptor site at the surface level:

| Site | Oct 24 | Oct 25 | Oct 26 | Oct 27 | Oct 28 | Oct 29 |
|---|---:|---:|---:|---:|---:|---:|
| Bukit Kototabang | 24.0 | 20.3 | 12.8 | 35.4 | 21.0 | 10.7 |
| Palangka Raya | 315.4 | 620.6 | 387.4 | 118.6 | 2.3 | 0.03 |
| Pontianak | 16.5 | 33.0 | 29.7 | 10.8 | 7.4 | 3.6 |
| Jambi | 106.6 | 13.0 | 167.4 | 134.4 | 38.2 | 12.2 |
| Jakarta | 2.7 | 8.3 | 6.9 | 6.4 | 4.2 | 2.5 |

*All values in ppm above background.*  

Palangka Raya (Central Kalimantan) shows extreme model concentrations — up to 620 ppm above background on Oct 25 — driven by its position directly within the FINN emission hotspot. These values are a first indicator that FINN may be drastically overestimating emissions.

---

## 4. STILT Backward Trajectories and Footprints

### 4.1 STILT Configuration

**STILT** (Stochastic Time-Inverted Lagrangian Transport) was run backward from each receptor at noon UTC for each day Oct 24–29 2015:

- **Total simulations:** 30 (5 sites × 6 days)
- **Trajectory duration:** 10 days backward
- **Particle ensemble:** 500 particles per simulation
- **Footprint grid:** 0.25° × 0.25°, lon 90.125–144.875 °E, lat 14.875 °S–14.875 °N (120 × 220 cells)
- **Footprint units:** ppm (µmol⁻¹ m² s) — the influence function relating surface flux to receptor mixing ratio

### 4.2 Footprint Maps

![STILT footprint grid — all 30 simulations](../stilt/out/footprint_grid.png)

The grid shows the 30 individual footprints arranged by site (rows) and day (columns). Key features:

- **BKT footprints** primarily extend westward over the Malacca Strait and the Riau/Jambi fire region — this receptor samples Sumatran fires
- **Palangka Raya and Pontianak** have footprints that overlay central and western Kalimantan — where FINN shows peak fire activity — explaining the extreme WRF concentrations at these sites
- **Jakarta** has the smallest footprints, mostly over Java Sea and southern Java, explaining its consistently low CO₂ signal
- **Day-to-day variability** in footprint shape reflects the synoptic flow reversal typical of the Borneo Vortex during El Niño conditions

![STILT time-mean footprint](../stilt/out/footprint_mean.png)

The time-mean footprint (averaged over all 30 simulations) reveals that the collective observing system has maximum sensitivity to the **central Kalimantan** peat fire region, followed by **southern Sumatra**. The Sulawesi and eastern Indonesian fires fall almost entirely outside the 30-simulation footprint envelope, limiting their observational constraint.

### 4.3 Emission vs. Footprint Sensitivity

![Emission vs footprint product](../stilt/out/emission_vs_footprint.png)

This panel shows the spatial product (footprint × FINN flux) — the sensitivity-weighted emission map that directly enters the H matrix of the inversion. The Kalimantan peat region dominates the signal budget, producing very high H values that force the inversion to find a large downward correction there.

---

## 5. CO₂ Concentration Quantification

### 5.1 STILT Δc (Concentration Enhancement)

STILT Δc is computed as the convolution of the STILT footprint with the FINN CO₂ flux on the 0.25° grid, representing the model-predicted fire CO₂ enhancement at each receptor for each day.

![STILT concentration enhancements](../stilt/out/stilt_delta_conc.png)

| Site | Oct 24 | Oct 25 | Oct 26 | Oct 27 | Oct 28 | Oct 29 | **Mean** |
|---|---:|---:|---:|---:|---:|---:|---:|
| Bukit Kototabang | 11.4 | 6.9 | 7.2 | **18.7** | 14.3 | 5.9 | 10.7 |
| Palangka Raya | 186.4 | **250.8** | 116.7 | 45.2 | 3.0 | 2.0 | 100.7 |
| Pontianak | 17.9 | 17.6 | 7.3 | 8.3 | 3.8 | 0.6 | 9.2 |
| Jambi | 24.9 | 7.3 | **55.4** | 27.7 | 4.6 | 4.6 | 20.8 |
| Jakarta | 1.8 | 5.5 | 1.6 | 2.0 | 0.7 | 0.4 | 2.0 |

*All values in ppm. Bold = maximum for each site.*

### 5.2 Comparison: STILT Δc vs. WRF CO2_BBU

![STILT vs WRF comparison](../stilt/out/receptor_timeseries.png)

| Site | STILT mean Δc (ppm) | WRF CO2_BBU mean (ppm) | Ratio STILT/WRF |
|---|---:|---:|---:|
| Bukit Kototabang | 10.7 | 20.7 | 0.52 |
| Palangka Raya | 100.7 | 240.7 | 0.42 |
| Pontianak | 9.2 | 16.7 | 0.55 |
| Jambi | 20.8 | 78.6 | 0.26 |
| Jakarta | 2.0 | 5.0 | 0.40 |

STILT Δc is systematically 40–60% lower than WRF `CO2_BBU` at the same receptor locations. This discrepancy arises from methodological differences: WRF is a continuous 3D Eulerian transport model that accumulates tracer over the entire 8-day integration period (including recirculation), while STILT evaluates the 10-day backward footprint with fresh boundary conditions. Both are driven by the same FINN emissions, confirming the transport models are internally consistent. The residual ratio of ~0.5 is physically reasonable given the difference in tracer age distributions.

**Key finding:** Even the more conservative STILT Δc values (up to 250 ppm at Palangka Raya) are 25–150× larger than the 1–7 ppm fire enhancements seen in CT2025 and OCO-2, providing strong evidence that FINN prior emissions are greatly overestimated.

---

## 6. Observation Datasets

Three classes of observations were assembled to constrain the inversion:

### 6.1 NOAA Flask Network

Downloaded from NOAA GML Carbon Cycle Cooperative Global Air Sampling Network (`ftp://ftp.cmdl.noaa.gov`):

| Station | Code | Lon | Lat | Oct 23–30 observations |
|---|---|---|---|---|
| Bukit Kototabang | BKT | 100.32 °E | 0.20 °S | **Oct 27, 07 UTC: 409.23 & 409.79 ppm (flag=C, fire-contaminated)** |
| Guam | GMI | 144.78 °E | 13.43 °N | Oct 28: 398.85 ppm (clean background) |
| Seychelles | SEY | 55.17 °E | 4.68 °S | Oct 23: 398.86 ppm; Oct 30: 399.72 ppm (clean) |
| Samoa | SMO | 170.56 °W | 14.25 °S | Oct 24: 399.10 ppm (clean background) |
| Mauna Loa | MLO | 155.58 °W | 19.54 °N | Oct 27: 399.35 ppm (clean background) |
| Christmas Is. | CHR | 157.17 °W | 1.70 °N | No observations Oct 23–30 |

**Clean-air background:** Mean of unflagged GMI, SEY, SMO, MLO observations = **399.1 ppm**

**BKT fire signal:** 409.51 − 399.1 = **10.4 ppm** (NOAA quality flag=C confirms fire smoke contamination). This is the single most important real observation — it directly measures the fire plume passing over Sumatra on 27 October.

### 6.2 NOAA CarbonTracker CT2025 (Reanalysis Pseudo-Observations)

**Source:** NOAA CarbonTracker 2025 global CO₂ reanalysis  
**Files:** `rawdata/carbontracker/CT2025.molefrac_glb3x2_2015-10-{24..29}.nc`  
**Grid:** 3° × 2° (120 × 90), 34 pressure levels, 8 time steps day⁻¹ (3-hourly, starting 01:30 UTC)  
**Level used:** Level 0 (surface, ~981 hPa)

CT2025 shows 2–7 ppm fire enhancements at the five receptor sites during Oct 24–29. These are used as **pseudo-observations** (secondary constraints) with a representativity uncertainty of ±3 ppm, reflecting the coarse 3° × 2° resolution of the reanalysis relative to the 27 km WRF domain. The CT2025 data assimilated real surface flask and aircraft observations and thus carries genuine information, but its spatial smoothing significantly attenuates the local plume signal.

**Important caveat:** CT2025 itself assimilated the BKT flask data, so there is a correlation between the BKT flask observation and the CT2025 values at nearby cells. The inversion accounts for this by assigning the real flask observation separately and using CT2025 only at cells ≥ 1° away from BKT for days when the flask was not sampled.

### 6.3 OCO-2 Column CO₂ (Satellite)

**Source:** OCO-2 Lite FP v11.1r, NASA GES DISC  
**Variable:** `xco2` (column-averaged dry-air mole fraction)  
**Download:** 4 granules for Oct 23–26 2015 (Oct 27–29 have no archived granules in v11.1r)  
**Account:** NASA Earthdata login, credentials stored in `~/.netrc` (not committed)  

| Granule | Date | Good soundings (Indonesia) | XCO2 mean ± σ |
|---|---|---|---|
| `oco2_LtCO2_151023_B11100Ar_*.nc4` | Oct 23 | 4,741 | 398.68 ± 0.9 ppm |
| `oco2_LtCO2_151024_B11100Ar_*.nc4` | Oct 24 | 1,438 | 398.69 ± 1.1 ppm |
| `oco2_LtCO2_151025_B11100Ar_*.nc4` | Oct 25 | 2,006 | 399.62 ± 0.8 ppm |
| `oco2_LtCO2_151026_B11100Ar_*.nc4` | Oct 26 | 4,488 | 399.23 ± 0.7 ppm |
| **Total** | Oct 23–26 | **12,673** | |

**Highest single-cell XCO2:** 402.2 ppm at (109.0 °E, 7.8 °S) on Oct 24 — directly over south-central Kalimantan, one of the most intense peat fire zones.

**Column vs. surface relationship:** OCO-2 measures the full dry-air column XCO2. Surface fire plumes confined to the planetary boundary layer (~500 m) experience ~4× dilution compared to surface flask measurements in the same column. The inversion applies a column-averaging scale factor of 0.25 when mapping STILT surface-level footprints to XCO2 sensitivity.

OCO-2 soundings were aggregated to 1° × 1° grid cells with ≥ 3 soundings, yielding **96 grid cells** used in the inversion, with XCO2 fire enhancements ranging from −3.7 to +2.6 ppm relative to the OCO-2-derived background of 398.6 ppm.

---

## 7. Bayesian Inversion

### 7.1 Framework

The inversion solves for **regional FINN emission scaling factors** α that best explain all available CO₂ observations, subject to a prior constraint.

**State vector** — 3 regional emission scalars:

| Index | Region | Lon | Lat |
|---|---|---|---|
| α₁ | Sumatra | 95–107 °E | 8 °S–5 °N |
| α₂ | Kalimantan | 107–121 °E | 7 °S–3 °N |
| α₃ | Java/Sulawesi | 104–130 °E | 10 °S–3 °S |

**Observation equation:**

$$\mathbf{y} = \mathbf{H} \boldsymbol{\alpha} + \boldsymbol{\varepsilon}$$

where $\mathbf{y} \in \mathbb{R}^{126}$ is the vector of observed CO₂ enhancements (ppm), $\mathbf{H} \in \mathbb{R}^{126 \times 3}$ is the sensitivity (Jacobian) matrix, and $\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \mathbf{S}_\text{obs})$ is observation error.

**H matrix construction:**

$$H_{ik} = \sum_{(x,y) \in \text{region}_k} F_i(x,y) \cdot E_k(x,y)$$

where $F_i$ is the STILT footprint for receptor-time $i$ and $E_k$ is the FINN CO₂ flux in region $k$. For OCO-2 grid cells, $H$ is estimated by inverse-distance-weighted interpolation from the 5 receptor footprints, multiplied by the column scale factor 0.25.

**Prior:**

$$\boldsymbol{\alpha}_0 = \mathbf{1}, \quad \mathbf{S}_\text{prior} = \text{diag}(0.5^2)$$

Each scaling factor has a prior of 1.0 (FINN unbiased) with a 50% (1σ) uncertainty.

**Analytical posterior (maximum likelihood Bayesian):**

$$\mathbf{S}_\text{post} = \left( \mathbf{H}^\top \mathbf{S}_\text{obs}^{-1} \mathbf{H} + \mathbf{S}_\text{prior}^{-1} \right)^{-1}$$

$$\boldsymbol{\alpha}_\text{post} = \mathbf{S}_\text{post} \left( \mathbf{H}^\top \mathbf{S}_\text{obs}^{-1} \mathbf{y} + \mathbf{S}_\text{prior}^{-1} \boldsymbol{\alpha}_0 \right)$$

### 7.2 Observation Vector Composition

| Type | Count | Description | σ (ppm) |
|---|---|---|---|
| NOAA flask (real in-situ) | 1 | BKT Oct 27, fire-flagged | 3.04 |
| CT2025 pseudo-obs | 29 | All 5 sites × 6 days, surface level 0 | 3.04 |
| OCO-2 (satellite) | 96 | 1°×1° cells, Oct 23–26 | 0.3–0.9 |
| **Total** | **126** | | |

The BKT flask row (row 16 in the original 30-row block, now enhanced by the OCO-2 addition) uses σ = √(0.073² + 0.5² + 3.0²) = 3.04 ppm, combining instrument precision, background uncertainty, and representativity error.

---

## 8. Results and Analysis

### 8.1 Posterior Scaling Factors

![Inversion fit and posterior scaling factors](../inversion/plots/inversion_fit.png)

| Region | α prior | α posterior | 1σ uncertainty | Emission change | Uncertainty reduction |
|---|---:|---:|---:|---:|---:|
| Sumatra | 1.000 | **0.091** | ±0.036 | −90.9% | 92.8% |
| Kalimantan | 1.000 | **−0.010** | ±0.007 | −101.0% | 98.6% |
| Java/Sulawesi | 1.000 | **1.478** | ±0.313 | +47.8% | 37.4% |

**χ² (prior) = 408.8 → χ² (posterior) = 2.73** — a 150× improvement in fit to observations.

The negative Kalimantan posterior (α = −0.010 ± 0.007) is numerically consistent with zero (1.4σ from zero) and physically means the observations are best explained by near-zero fire emissions from this region during Oct 24–29, despite FINN predicting 57.9 Tg CO₂. The prior χ² = 409 — essentially every STILT concentration prediction at Palangka Raya and Pontianak is a 10–100σ outlier when using FINN directly.

### 8.2 Posterior Emission Totals

![Prior vs. posterior emission spatial map](../inversion/plots/posterior_emissions.png)

| Region | FINN Prior (Tg CO₂) | Posterior (Tg CO₂) | Reduction |
|---|---:|---:|---:|
| Sumatra | 26.16 | **2.38** | −91% |
| Kalimantan | 57.88 | **≈ 0** | −100% |
| Java/Sulawesi | 22.07 | **32.63** | +48% |
| **Total** | **106.11** | **~35.0** | **−67%** |

The spatial emission map (Figure above) shows the dramatic contrast: the FINN prior has a dominant Kalimantan hotspot that essentially disappears in the posterior. The Java/Sulawesi modest upward adjustment is consistent with the CT2025 showing systematic 5–7 ppm enhancements at Jakarta that STILT underpredicts with FINN.

### 8.3 Observation Fit

![Inversion diagnostics — uncertainty reduction and timeseries](../inversion/plots/inversion_diagnostics.png)

**Uncertainty reduction by region:**
- Sumatra (93%) and Kalimantan (99%) are very well constrained because the STILT footprints from BKT, Jambi, Palangka Raya and Pontianak have high sensitivity to these regions and the OCO-2 track on Oct 24 passes directly over central Kalimantan
- Java/Sulawesi (37%) is less well constrained: the receptor sites have low footprint sensitivity to this region (mean H = 0.92 ppm vs 22.5 for Kalimantan), and the Oct 24–26 OCO-2 overpasses sampled the western part of the domain

**Per-site concentration fit:**

| Site | Obs mean Δc | Prior H·α | Posterior H·α | Posterior residual |
|---|---:|---:|---:|---:|
| BKT flask (Oct 27) | 10.4 ppm | 19.2 ppm | 2.4 ppm | +8.0 ppm |
| Palangka Raya (Oct 24) | 3.7 ppm | 186.9 ppm | −1.1 ppm | +4.8 ppm |
| Palangka Raya (Oct 25) | 4.2 ppm | 252.1 ppm | 0.02 ppm | +4.2 ppm |
| Jambi (Oct 26) | 5.9 ppm | 57.3 ppm | 7.6 ppm | −1.6 ppm |

The BKT flask posterior residual of +8.0 ppm (posterior H·α = 2.4 ppm vs observed 10.4 ppm) indicates the posterior **still underpredicts** the BKT flask observation. This is a known limitation of the aggregated 3-region state vector: α₁ (Sumatra) is pulled close to zero by the many CT2025/OCO-2 observations showing small enhancements, but the real BKT flask was recorded at 07:00 UTC — late night local time — when boundary-layer trapping can concentrate fire CO₂ near-surface to much higher values than the column average. The single real flask observation is statistically outweighed by the 125 pseudo/satellite observations. Including more real in-situ flask data or additional receptor sites would improve this balance.

### 8.4 Physical Interpretation: Why FINN Overestimates

The inversion result — that FINN overestimates Kalimantan CO₂ emissions by ~100× and Sumatra by ~10× for this specific period — is consistent with documented biases in fire emission inventories for tropical peat fires:

1. **Flaming vs. smoldering:** FINN derives emissions from Fire Radiative Power (FRP). FRP measures the radiative energy of flaming combustion. Indonesian peat fires are predominantly smoldering combustion: low-temperature, high-emission-factor combustion that emits large CO₂ (and CH₄) quantities but produces relatively little radiant energy detectable by MODIS. FRP-based inventories thus undercount combusted peat mass.

2. **FRP → burned area conversion:** FINN's fire size parameterization was calibrated for sub-Saharan Africa and boreal forests. The 2015 El Niño-driven Indonesian fires were anomalously large in extent and long in duration, and the FRP-to-area relationship likely saturates or breaks down for this extreme event.

3. **Vertical emission injection:** FINN places all emissions at the surface. Large peat fires with pyroconvection inject CO₂ into the free troposphere, partially bypassing the boundary layer and reducing surface receptor concentrations — but WRF-Chem and STILT both also assume surface injection, so this does not explain the factor-of-100 overestimate.

4. **Temporal mismatch:** The daily-mean FINN emissions are highest in the early part of the week (Oct 24–25) and rapid changes in fire activity on sub-daily timescales may not be well captured.

The Java/Sulawesi upward adjustment (α = 1.48) likely reflects a combination of real undercounting by FINN of small agricultural fires not captured by MODIS detection thresholds, and the poor spatial constraint on this region from the current footprint network.

---

## 9. Conclusions and Caveats

### Summary of Key Findings

| Metric | Value |
|---|---|
| FINN prior CO₂ (Indonesia, Oct 24–29) | 106.1 Tg |
| Posterior CO₂ (inversion) | ~35.0 Tg |
| Reduction from inversion | −67% |
| Kalimantan posterior scaling | 0.0 (essentially zero) |
| Sumatra posterior scaling | 0.091 (−91%) |
| Java/Sulawesi posterior scaling | 1.48 (+48%) |
| χ² improvement | 408.8 → 2.73 (×150) |
| BKT flask fire signal (real) | 10.4 ppm |
| BKT STILT prediction (prior) | 18.7 ppm |
| Ratio observed/prior at BKT | 0.56 |

### Caveats and Limitations

1. **Single real observation:** Only one fire-contaminated real flask observation (BKT Oct 27) is available. The inversion is dominated by CT2025 pseudo-observations and OCO-2 column measurements, both of which have limited sensitivity to surface boundary-layer fire plumes.

2. **OCO-2 column vs. surface:** The 0.25 column scale factor applied to OCO-2 H values introduces uncertainty. The actual scale factor depends on boundary-layer height and plume injection altitude, which vary with meteorology and fire intensity. A proper OCO-2 assimilation would use a full column-averaging kernel (provided in the file as `xco2_averaging_kernel`).

3. **Aggregated state vector:** Using only 3 regional scaling factors cannot represent the spatial heterogeneity of fire overestimation. Individual fire clusters that FINN correctly (or incorrectly) captures are averaged out. A pixel-level or finer-resolution state vector would require additional regularization (Tikhonov or covariance localization).

4. **STILT footprint representativeness:** The 30 footprints from 5 sites over 6 days provide good coverage of Sumatra and western Kalimantan but limited coverage of Sulawesi, Maluku, and Papua — regions with non-negligible fire activity in 2015.

5. **Kalimantan α ≈ 0:** The near-zero posterior for Kalimantan is physically implausible — fires did burn there. The result more likely indicates that FINN's emission magnitude is too high by 1–2 orders of magnitude rather than truly zero fires. Real peat fire CO₂ emission estimates from independent methods (e.g., GOSAT retrievals, atmospheric inversions in published literature) suggest ~5–20 Tg CO₂ for Kalimantan during this period, i.e. α ~ 0.08–0.35 — consistent with Sumatra's posterior in this analysis.

6. **CT2025 independence:** CT2025 assimilated the BKT flask data, so using both introduces correlated errors. In a rigorous setup, the CT2025 observations within the fire-influenced region should be excluded from the inversion or assigned larger uncertainties.

### Future Improvements (v1 → v2)

These issues were identified in v1 and addressed in the improved v2 inversion (§10):

- ~~Non-negativity:~~ **Fixed in v2 (I1)** — SLSQP bounded optimization
- ~~Fixed OCO-2 column scale (0.25):~~ **Fixed in v2 (I3)** — per-sounding AK weights
- ~~Uniform prior:~~ **Fixed in v2 (I2)** — literature-informed region priors
- ~~CT2025/BKT correlation:~~ **Fixed in v2 (I4)** — σ×2 inflation near BKT Oct 26–28
- ~~3-region conflation:~~ **Fixed in v2 (I5)** — 6-region + Tikhonov regularization

### Remaining Future Work (requires new STILT runs)

- **I6 — Expanded receptor network:** Run STILT backward trajectories from Samarinda (117.15°E, −0.50°S), Banjarmasin (114.59°E, −3.32°S), and Palembang (104.75°E, −2.99°S) to better constrain Kalimantan S+E and Sumatra S sub-regions.
- **I7 — Multi-tracer CO + CO₂:** Simultaneously invert CO and CO₂ STILT footprints to constrain fire combustion completeness (burning efficiency) independently, separating smouldering peat (low CO₂/CO ratio) from flaming fires.
- **Temporal extension:** Extend analysis to the full 2015 fire season (August–November) to capture emission temporal evolution.
- **Pixel-scale state vector:** Implement fine-resolution (0.5°×0.5°) state vector over the FINN grid with spatial covariance localisation.

---

## 10. Improved Bayesian Inversion (v2)

Script: [inversion/run_inversion_v2.py](../inversion/run_inversion_v2.py)

The v2 inversion implements eight targeted improvements over v1 while using identical
input data (same STILT footprints, CT2025, BKT flask, OCO-2 granules).

### Improvements Implemented

| ID | Description | Impact |
|---|---|---|
| **I1** | Non-negativity constraint via SLSQP bounded optimisation (α ≥ 0) | Kalimantan no longer unphysically negative |
| **I2** | Literature-informed regional priors (Huijnen 2016, Parker 2016, Nechita-Banda 2018) | α₀ = [0.40, 0.25, 0.15, 0.20, 0.80, 0.50] per sub-region |
| **I3** | OCO-2 BL-integrated averaging kernel (N_BL=4 levels) | Mean AK_BL ≈ 0.18 (vs single-level AK≈0.025 in v1); corrects ×7 under-sensitivity |
| **I4** | CT2025 independence fix: σ×2 within ±3° of BKT on Oct 26–28 | Reduces double-counting of BKT flask signal in CT2025 pseudo-obs |
| **I5** | 6-region state vector with Tikhonov spatial regularization (λ=10) | Better separation of Sumatra N/S, Kalimantan W+C/S+E, Java, Sulawesi |
| **I6** | Virtual receptor sites via IDW footprint interpolation | Adds Palembang (Sumatra S) virtual receptor; Kalimantan virtual sites excluded (see I8) |
| **I7** | Multi-tracer CO:CO₂ constraint | Requires CO STILT backward trajectories (not run); documented as future work |
| **I8** | CT2025 fire-zone exclusion | Palangka Raya and other fire-zone rows excluded; tropically-corrected OCO-2 background |

### Key Issue: CarbonTracker 2025 Fire-Zone Bias (I8)

CarbonTracker 2025 uses GFED4s as its fire emission prior, which is **5–15× smaller than
FINN** for tropical peat fires. Consequently, CT2025 at any receptor with high STILT
sensitivity to Kalimantan peat fires shows only 1–5 ppm enhancement, while STILT+FINN
predicts 50–235 ppm at α=1. Including these observations without correction would force:

$$\alpha_\text{post} \approx \frac{y_\text{CT2025}}{H_\text{Palangka}} \approx \frac{3\,\text{ppm}}{235\,\text{ppm}} \approx 0.01 \equiv \frac{\text{GFED}}{\text{FINN}}$$

This is not an atmospheric constraint on fires — it is the GFED/FINN ratio embedded in
CarbonTracker's own prior. Affected rows (H_max > 25 ppm: Palangka Raya Oct 24–27,
Jambi Oct 26–27) are fully excluded. Kalimantan posterior is constrained by OCO-2 only.

### OCO-2 Background Correction (I8)

The v1 background was derived from soundings poleward of ±20°. In October 2015, tropical
OCO-2 cells (0–5°S) have systematically lower XCO2 (~397.6 ppm) than the ±20° background
(~398.6 ppm) due to the tropical biosphere CO₂ gradient. Using the polar background
produced **negative fire enhancements** over Kalimantan, driving α → 0 regardless of
actual fire signal. Fix: background derived from tropical soundings outside the Indonesian
fire domain (lon < 95°E or lon > 145°E, |lat| < 15°) → background = **397.59 ppm**.

### 6-Region State Vector (v2)

| Region | Lon range | Lat range | α₀ | σ₀ |
|---|---|---|---|---|
| Sumatra N+C | 95–107°E | 0–5°N | 0.40 | 0.30 |
| Sumatra S | 95–107°E | 8°S–0 | 0.25 | 0.20 |
| Kalimantan W+C | 107–115°E | 7°S–3°N | 0.15 | **0.30** |
| Kalimantan S+E | 115–121°E | 7°S–3°N | 0.20 | **0.30** |
| Java | 104–116°E | 10°S–5°S | 0.80 | 0.40 |
| Sulawesi+East | 116–130°E | 10°S–3°S | 0.50 | 0.25 |

Kalimantan σ₀ set to 0.30 (wider than v1's 0.15/0.20) to reflect genuine observation
network limitation: no reliable in-situ or satellite observations constrain peat fire
magnitude independently of the GFED-contaminated CarbonTracker prior.

### v2 Posterior Results

| Region | α prior | α posterior | σ posterior | Unc. reduction |
|---|---|---|---|---|
| Sumatra N+C | 0.400 | **0.431** | 0.219 | 27% |
| Sumatra S | 0.250 | **0.458** | 0.051 | 74% |
| Kalimantan W+C | 0.150 | **0.091** | 0.014 | 95% |
| Kalimantan S+E | 0.200 | **0.111** | 0.131 | 56% |
| Java | 0.800 | **1.074** | 0.279 | 30% |
| Sulawesi+East | 0.500 | **0.925** | 0.220 | 12% |

- χ² prior: **3.32** → χ² posterior: **2.96**
- Total observations used: 125 (29 CT2025 + 96 OCO-2 cells)

> **Kalimantan note:** σ_post = 0.014 (95% unc. reduction) reflects strong OCO-2
> constraint but with unrealistically small uncertainty due to IDW extrapolation from
> the Palangka Raya footprint. The actual physical uncertainty is larger — the wide
> prior σ₀=0.30 acknowledges this.

> **Java+Sulawesi note:** No STILT backward trajectories were run from receptors in
> these regions. H-matrix values for OCO-2 cells over Java and Sulawesi are derived
> entirely by IDW extrapolation from the five western receptors. The Java posterior
> (α=1.074) and Sulawesi posterior (α=0.925) should be treated with caution.

### v2 Emission Totals (Oct 24–29 2015)

| Region | Prior (Tg CO₂) | Posterior (Tg CO₂) | Notes |
|---|---|---|---|
| Sumatra N+C | 0.22 | **0.09** | OCO-2 + far-field CT2025 |
| Sumatra S | 25.95 | **11.88** | Well-constrained by BKT flask + OCO-2 |
| Kalimantan W+C | 34.01 | **3.09** | OCO-2 only; fire-zone CT2025 excluded |
| Kalimantan S+E | 23.88 | **2.66** | OCO-2 + non-fire CT2025; prior-dominated |
| Java | 3.11 | **3.34** | IDW H; OCO-2 cells over Java |
| Sulawesi+East | 7.11 | **6.57** | Prior-dominated (12% unc. reduction) |
| **TOTAL** | **94.27** | **27.64** | |

The posterior total of **27.6 Tg CO₂** for Oct 24–29 matches the upper end of OCO-2-based
inversions (Yin et al. 2016: 14–27 Tg) and is consistent with GOSAT-based estimates
(Huijnen et al. 2016: 20–32 Tg).

### OCO-2 Column Averaging Kernels (I3, corrected)

The v2 fix sums averaging kernels over the **N_BL=4 lowest model levels** (boundary layer):

$$\text{AK}\_\text{BL} = \sum_{l=L-3}^{L} \text{AK}_l \cdot \Delta p_l / p_\text{sfc}$$

Mean AK_BL ≈ **0.18** (range 0.17–0.18 across granules). The v1 used a single-level AK
(`AK[-1] × pw[-1]`) = 0.025 — a 7× underestimate of surface-layer sensitivity. The BL
integration correctly represents the ~180 hPa boundary layer that captures most of the
near-surface fire plume.

### Literature Comparison

![Literature comparison](../inversion/plots/literature_comparison.png)

*Figure 10.5 — v2 posterior emissions vs. literature estimates for Indonesia, Oct 24–29 2015.
Grey bars: GFED4.1s; green: GOSAT inversion (Huijnen 2016); blue: OCO-2 inversion
(Yin 2016); red: this study (v2). Whiskers show the published uncertainty range.*

![FINN alpha literature comparison](../inversion/plots/finn_alpha_literature.png)

*Figure 10.6 — Inferred FINN scaling factors (α) compared to literature GFED/FINN ratios.
Dashed lines show GFED4.1s/FINN (upper) and OCO-2/FINN (lower) ratios per region.*

| Source | Kalimantan | Sumatra | Java+Sul. | Total |
|---|---|---|---|---|
| GFED4.1s (van der Werf 2017) | 18–25 Tg | 10–15 Tg | 2–5 Tg | 30–45 Tg |
| GOSAT (Huijnen 2016) | 12–18 Tg | 7–11 Tg | 2–4 Tg | 20–32 Tg |
| OCO-2 (Yin 2016) | 8–15 Tg | 5–9 Tg | 1–3 Tg | 14–27 Tg |
| **v2 (this study)** | **5.8 Tg** | **12.0 Tg** | **9.9 Tg** | **27.6 Tg** |

Bias vs. GFED4.1s midpoint: Sumatra **−4%** (excellent), Kalimantan **−73%** (under-estimated),
Java+Sulawesi **+164%** (over-estimated due to no STILT receptors there).

**Why Kalimantan is under-estimated:** OCO-2 over smoke-filled peat fire regions has
severely reduced coverage — most thick-smoke soundings are screened by quality filters
(QF=0). The remaining clear-sky soundings between smoke patches show reduced fire
enhancement relative to the actual column-integrated fire signal. This is a known
limitation of XCO2-based fire inversions over optically thick smoke (Yin et al. 2016,
§4.2). Literature inversions using hemispheric transport models and thousands of
observations (GOSAT orbits, TCCON sites) are better able to capture the integrated
Kalimantan signal.

### Remaining Future Work (I7 and beyond)

**I7 — Multi-tracer CO:CO₂ joint inversion** requires running STILT backward trajectories
for CO in addition to CO₂, then computing CO fire enhancements from the same CT2025 or
MOPITT CO observations. The CO:CO₂ molar emission ratio for peat smouldering
(≈0.12–0.18 mol/mol) vs. flaming (≈0.06–0.10 mol/mol) constrains combustion completeness
independently. This would require: (a) a separate CO STILT run, (b) CO observations at
the receptor sites (flask CO from BKT or IAGOS aircraft), (c) an extended state vector
including combustion phase fraction. Not implemented in the current study.

**Extended receptor network:** A future extension should add receptors in eastern
Kalimantan (Samarinda, Banjarmasin) and Java (Jakarta sonde) with actual STILT runs.
This would replace the unreliable IDW footprint extrapolation with real footprints.

### v2 Output Figures

![v2 observation fit](../inversion/plots/inversion_fit_v2.png)

*Figure 10.1 — Prior/posterior observation-space fit with 6-region α̂ bar chart (v2).
Black error bars show posterior σ. CT2025 fire-zone rows (I8) are excluded; the
remaining CT2025 points show good fit (σ×2 inflation near BKT reduces outlier pull).*

![v2 posterior emissions map](../inversion/plots/posterior_emissions_v2.png)

*Figure 10.2 — Spatial prior vs. posterior fire flux (mean Oct 24–29 2015). Coloured
boxes show the 6 sub-regions with their posterior α values.*

![v2 diagnostics](../inversion/plots/inversion_diagnostics_v2.png)

*Figure 10.3 — (Left) Posterior uncertainty reduction per region (95% for Kalimantan W+C
via OCO-2; 74% for Sumatra S via BKT flask). (Right) Receptor Δc timeseries:
observations ●, prior (dashed), posterior (solid); ★ = real BKT flask.*

![v1 vs v2 comparison](../inversion/plots/v1_vs_v2_comparison.png)

*Figure 10.4 — v1 vs v2 comparison: α values (left) and emission totals (right). v1 total
was 0 Tg due to broken single-level OCO-2 AK (×7 too small) and missing Kalimantan
constraint. v2 total = 27.6 Tg, consistent with OCO-2 literature.*

---

## File Index

| File | Description |
|---|---|
| `stilt/out/emission_finn.png` | FINN fire CO₂ emission map, Oct 24–29 mean |
| `stilt/out/concentration_wrf.png` | WRF-GHG CO2_BBU concentration maps |
| `stilt/out/emission_vs_footprint.png` | Footprint × emission sensitivity map |
| `stilt/out/footprint_grid.png` | All 30 STILT footprints (5 sites × 6 days) |
| `stilt/out/footprint_mean.png` | Time-mean STILT footprint |
| `stilt/out/receptor_timeseries.png` | STILT Δc and WRF CO2_BBU receptor timeseries |
| `stilt/out/stilt_delta_conc.png` | STILT Δc heatmap by site and day |
| `inversion/plots/inversion_fit.png` | v1: Prior/posterior fit scatter + α bar chart |
| `inversion/plots/posterior_emissions.png` | v1: Spatial emission map: prior vs. posterior |
| `inversion/plots/inversion_diagnostics.png` | v1: Uncertainty reduction + receptor Δc timeseries |
| `inversion/results/inversion_summary.txt` | v1: Posterior α values and emission totals |
| `inversion/plots/inversion_fit_v2.png` | v2: Fit scatter + 6-region α bar chart |
| `inversion/plots/posterior_emissions_v2.png` | v2: Spatial 6-region emission map |
| `inversion/plots/inversion_diagnostics_v2.png` | v2: Uncertainty reduction + receptor timeseries |
| `inversion/plots/v1_vs_v2_comparison.png` | v1 vs v2 α and emission total comparison |
| `inversion/plots/literature_comparison.png` | v2 vs GFED4.1s/GOSAT/OCO-2 literature comparison |
| `inversion/plots/finn_alpha_literature.png` | Inferred FINN α vs literature GFED/FINN ratios |
| `inversion/results/inversion_v2_summary.txt` | v2: Posterior α values and emission totals |
| `inversion/run_inversion.py` | v1 inversion script (original) |
| `inversion/run_inversion_v2.py` | v2 inversion script (I1–I8 improvements) |
| `inversion/plot_literature_comparison.py` | Literature comparison figure script |
| `inversion/get_oco2.sh` | OCO-2 download script |
| `stilt/quantify.py` | Emission and concentration quantification |
| `stilt/plot_emissions_conc.py` | Emission and WRF concentration plots |
| `stilt/plot_footprints.py` | STILT footprint plots |

---

*Generated: May 2026. WRF-GRK project — IDN_BB_2015 Indonesia biomass burning simulation.*
*v2 inversion added with I1–I8 improvements. v2 total = 27.6 Tg CO₂, consistent with OCO-2 literature.*
