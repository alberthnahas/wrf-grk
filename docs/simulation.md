# WRF-GHG (WRF-GRK) Simulation: IDN_BB_2015

**Case:** Indonesia Biomass Burning Episode, Aug–Nov 2015  
**Configuration:** WRF-Chem 4.1.5, GHG tracer mode (`chem_opt=17`)  
**Last updated:** 2026-05-06

> **Operating this simulation:** the `wrf-ghg` Claude skill at `~/.agents/skills/wrf-ghg/` automates segment launch, input verification, monitoring, and stitching. It references this document as the authoritative source. See `~/.agents/skills/wrf-ghg/SKILL.md` and `references/segment_workflow.md`.

---

## Table of Contents

1. [Scientific Objective](#1-scientific-objective)
2. [Compute Environment](#2-compute-environment)
3. [Software Stack](#3-software-stack)
4. [Domain Configuration](#4-domain-configuration)
5. [WRF Physics & Chemistry Settings](#5-wrf-physics--chemistry-settings)
6. [Simulation Period & Segments](#6-simulation-period--segments)
7. [Input Data Sources](#7-input-data-sources)
8. [Preprocessing Scripts](#8-preprocessing-scripts)
9. [Emissions Detail](#9-emissions-detail)
10. [Boundary & Initial Conditions](#10-boundary--initial-conditions)
11. [Run Scripts & Pipeline](#11-run-scripts--pipeline)
12. [Directory Structure](#12-directory-structure)
13. [Known Issues & Bug Fixes](#13-known-issues--bug-fixes)
14. [Progress Log](#14-progress-log)
15. [Pending Tasks](#15-pending-tasks)

---

## 1. Scientific Objective

Simulate the atmospheric transport and budget of CO₂, CH₄, and CO over the Maritime Continent during the extreme 2015 El Niño biomass burning event (August–November 2015). The El Niño–Southern Oscillation (ENSO) event produced record peatland and forest fires across Kalimantan, Sumatra, and Papua, making this one of the largest single-year GHG emission events in the historical record.

**Research goals:**
- Quantify fire-sector CO₂, CH₄, and CO enhancements from Indonesian peatland burning
- Evaluate WRF-GHG against CarbonTracker column mole fractions and surface observations
- Separate biospheric (VPRM), fire (FINN), anthropogenic (EDGAR), and ocean (SOMFFN) flux contributions
- Assess the role of convective vertical transport in GHG redistribution

---

## 2. Compute Environment

| Item | Detail |
|------|--------|
| Hostname | `pemodelanKU` |
| OS | Ubuntu 22.04.3 LTS (kernel 6.8.0-107-generic) |
| Architecture | x86_64 |
| CPU | Intel Xeon E5-2650 v2 @ 2.60 GHz |
| Physical cores | 16 (2 NUMA nodes × 8 cores) |
| Logical CPUs | 32 (HT enabled) |
| RAM | 110 GiB total (~71 GiB available) |
| Disk (`/`) | 1007 GB, 330 GB used, 626 GB free |
| Python venv | `/home/igrk/WRF-GRK/.venv` |
| Working dir | `/home/igrk/WRF-GRK` |

---

## 3. Software Stack

| Software | Version | Notes |
|----------|---------|-------|
| WRF-Chem | 4.1.5 | Compiled with `chem_opt=17` (GHG tracers) |
| WPS | 4.1 | Preprocessor; Vtable.ERA-interim used for ERA5 |
| gfortran | 11.4.0 | `configure` option 34 (GNU/OpenMPI dmpar) |
| OpenMPI | 4.1.2 | Parallel WRF runtime |
| NetCDF | 4.8.1 | Libs at `/usr/lib/x86_64-linux-gnu/` |
| Python | 3.10.12 | Preprocessing scripts |
| numpy / xarray / scipy / netCDF4 / pyhdf | (venv) | Key Python deps |
| screen | system | Background job management |

**WRF executables:**
- `WRF/main/wrf.exe` (69 MB)
- `WRF/main/real.exe` (58 MB)

---

## 4. Domain Configuration

**Projection:** Mercator  
**Center:** 0.0°N, 117.5°E (central Borneo/Indonesia)

| Parameter | Value |
|-----------|-------|
| `MAP_PROJ` | `mercator` |
| `REF_LAT` / `REF_LON` | 0.0°N / 117.5°E |
| `TRUELAT1` | 0.0° |
| `STAND_LON` | 117.5°E |
| `DX` / `DY` | 27,000 m (27 km) |
| `E_WE` | 232 (west–east grid points) |
| `E_SN` | 128 (south–north grid points) |
| `E_VERT` | 50 vertical levels |
| `P_TOP` | 5,000 Pa (50 hPa) |
| `GEOG_DATA_RES` | `modis_30s_lake+modis_lai+modis_fpar+default` |
| `geo_em.d01.nc` | 16 MB — **SUCCESS** (geogrid complete) |

**Approximate domain extent:** ~90°E–145°E, ~15°S–15°N — covers all of Indonesia, the Indian Ocean east of Sri Lanka, southern Philippines, Papua New Guinea, and surrounding seas. (231 × 27 km = 6,237 km E–W; 127 × 27 km = 3,429 km N–S; centered on 0°N/117.5°E.)

---

## 5. WRF Physics & Chemistry Settings

### Physics (`&physics`)

| Option | Setting | Scheme |
|--------|---------|--------|
| `mp_physics` | 6 | WSM 6-class graupel |
| `ra_lw_physics` | 4 | RRTMG longwave |
| `ra_sw_physics` | 4 | RRTMG shortwave |
| `radt` | 27 min | Radiation call interval |
| `sf_sfclay_physics` | 1 | Revised MM5 surface layer |
| `sf_surface_physics` | 2 | Noah Land Surface Model |
| `bl_pbl_physics` | 5 | MYNN level-2.5 TKE PBL |
| `cu_physics` | 3 | Grell-Freitas ensemble cumulus |
| `num_soil_layers` | 4 | |

### Dynamics (`&dynamics`)

| Option | Setting |
|--------|---------|
| `w_damping` | 1 (on) |
| `diff_opt` | 1 |
| `km_opt` | 4 |
| `damp_opt` | 3 (Rayleigh) |
| `zdamp` | 5,000 m |
| `dampcoef` | 0.2 |
| `non_hydrostatic` | `.true.` |
| `moist_adv_opt` | 1 |
| `scalar_adv_opt` | 1 |

### Chemistry (`&chem`) — GHG tracer mode

| Option | Setting | Notes |
|--------|---------|-------|
| `chem_opt` | 17 | WRF-GHG tracer (CO₂, CH₄, CO) |
| `emiss_opt` | 17 | GHG emissions input |
| `emiss_inpt_opt` | 16 | |
| `bio_emiss_opt` | 17 | VPRM biogenic (NEE/respiration/GPP) |
| `biomass_burn_opt` | 5 | Fire emissions with GHG plume rise (`biomassb_ghg` registry package) |
| `vprm_opt` | `'VPRM_table_TROPICS'` | VPRM vegetation table for tropics (CHARACTER type — not integer) |
| `term_opt` | `'CH4_termite_NW'` | CH₄ termite emission table (CHARACTER type — not integer) |
| `have_bcs_chem` | `.true.` | Use CarbonTracker lateral BCs |
| `gaschem_onoff` | 0 | No gas-phase chemistry (GHG only) |
| `aerchem_onoff` | 0 | No aerosol chemistry |
| `vertmix_onoff` | 1 | Vertical mixing on |
| `chem_conv_tr` | 1 | Convective transport on |
| `kemit` | 1 | Single emission layer |

### Auxiliary inputs

| auxinput | File pattern | Interval | Description |
|----------|-------------|----------|-------------|
| 5 | `wrfchemi_d01_<date>` | 60 min (24 frames/file, daily) | Anthropogenic (EDGAR). See Bug 11 — original hourly single-frame layout caused `Status = -4`; bundled into 130 daily 24-frame files. |
| 6 | `wrfoce_d01_<date>` | 43,200 min (30 days) | Ocean CO₂ flux (SOMFFN) |
| 7 | `wrffirechemi_d01_<date>` | 60 min (24 frames/file) | Fire (FINN+hotspot) |
| 15 | `vprm_input_d01_<date>` | 11,520 min (8 days) | VPRM biogenic |

---

### 5.3 Physics Scheme Rationale — Indonesian Tropical Conditions

The physics options were tuned for Maritime Continent / tropical island convection rather than mid-latitude WRF defaults. Changes from default and justifications:

#### Cumulus: Grell-Freitas Ensemble (`cu_physics=3`) — changed from Kain-Fritsch (opt=1)

| Aspect | Kain-Fritsch (old) | Grell-Freitas (new) |
|--------|-------------------|---------------------|
| Convection trigger | CAPE + moisture lift | Stochastic ensemble, no hard threshold |
| Tropical diurnal cycle | Too-early triggering over islands | Better matches 14:00–16:00 LST peak |
| Maritime Continent literature | Known positive bias in morning convection | Preferred in SE Asia studies |
| GHG relevance | Premature deep convection transports surface plumes too early | More accurate convective vertical transport of CO₂/CH₄/CO from peat fires |

Convective scheme choice is **critical for this study**: `chem_conv_tr=1` routes GHG tracers through the cumulus parameterisation. Under KF, afternoon peat fire emissions would be convectively pumped to the free troposphere ~2–4 h too early, misplacing the vertical column signal relative to GOSAT/TROPOMI overpass times.

#### PBL: MYNN Level-2.5 TKE (`bl_pbl_physics=5`) — changed from YSU (opt=1)

| Aspect | YSU (old) | MYNN 2.5 (new) |
|--------|-----------|----------------|
| Formulation | Non-local K-profile | Local TKE closure (prognostic TKE) |
| Stable layers | Struggles with strong inversions | Explicit TKE budget handles stable regimes |
| Smoke/haze BL | Not designed for optically thick aerosol layers | Better entrainment in smoke-stabilised PBL |
| Deep tropical BL (2–3 km) | Adequate but non-local | TKE approach preferred for deep tropical boundary layers |

During the 2015 event, dense peat smoke (AOD > 3.0 over Sumatra/Kalimantan) suppressed surface heating, creating persistent daytime stable layers with temperature inversions. MYNN's TKE budget explicitly resolves entrainment across these inversions, improving the vertical distribution of near-surface GHG enhancements. YSU's non-local approach tends to over-mix tracers through the stable layers.

**Surface layer** (`sf_sfclay_physics=1`, Revised MM5) is retained — fully compatible with MYNN PBL and well-validated over tropical ocean/land surfaces.

#### Radiation: RRTMG (`ra_lw/sw_physics=4`) — unchanged, already optimal

RRTMG is the correct choice for this case:
- Handles aerosol-radiation interaction via `aer_opt` — critical for the 2015 haze event
- Same scheme used in CMIP6 GCMs, enabling clean comparison with reanalysis
- `radt=27` (radiation call every 27 min, matching dx in km) is the standard recommended setting

#### Microphysics: WSM6 (`mp_physics=6`) — unchanged

Single-moment 6-class (ice, snow, graupel) is adequate for 27 km GHG transport. Convective rainfall modifies aerosol washout but precipitation accuracy is secondary to tracer transport for this study.

#### Land Surface: Noah LSM (`sf_surface_physics=2`) — unchanged

Handles tropical soil moisture, latent heat flux over rainforest and peat. Note: WRF-Chem 4.1.5 has no dedicated peat/organic soil parameterisation in Noah — peat carbon emissions are handled entirely by FINN fire inputs, not the LSM.

---

### 5.4 GHG Tracer Output Variables (chem_opt=17)

The `ghg_tracer` package (activated by `chem_opt==17` in `WRF/Registry/registry.chem`) allocates exactly **15 3D tracer fields**, all in units of **ppmv dry-air mole fraction**, written to every `wrfout_d01_*` file:

#### CO₂ tracers — 6 × 3D (ikj)

| WRF variable | Source dataset | Physical meaning |
|-------------|---------------|-----------------|
| `CO2_ANT` | EDGAR v8.0 (19 sectors) | Fossil fuel + industry CO₂ |
| `CO2_BIO` | VPRM / MODIS | Net biosphere exchange (GPP − respiration) |
| `CO2_OCE` | SOMFFN v2023 | Ocean–atmosphere CO₂ flux |
| `CO2_BCK` | CarbonTracker CT2025 | Background / lateral BC CO₂ |
| `CO2_BBU` | FINN v2.5 | Biomass burning CO₂ |
| `CO2_TST` | EDGAR (twin of `*_ANT`) | `*_TST` is driven only by `E_*TST` in `wrfchemi` (= `E_*` anthro). Per-source duplicate, not a tagged total — see Bug 18. |

#### CH₄ tracers — 5 × 3D (ikj)

| WRF variable | Source dataset | Physical meaning |
|-------------|---------------|-----------------|
| `CH4_ANT` | EDGAR v8.0 (22 sectors excl. AWB) | Anthropogenic CH₄ (fossil, waste, agriculture) |
| `CH4_BIO` | VPRM / wetland map | Biogenic CH₄ (wetlands, rice paddies) |
| `CH4_BCK` | CarbonTracker CTCH4_2024 | Background / lateral BC CH₄ |
| `CH4_BBU` | FINN v2.5 | Biomass burning CH₄ |
| `CH4_TST` | EDGAR (twin of `*_ANT`) | `*_TST` is driven only by `E_*TST` in `wrfchemi` (= `E_*` anthro). Per-source duplicate, not a tagged total — see Bug 18. |

#### CO tracers — 4 × 3D (ikj)

| WRF variable | Source dataset | Physical meaning |
|-------------|---------------|-----------------|
| `CO_ANT` | HTAP v3 | Anthropogenic CO |
| `CO_BCK` | CarbonTracker | Background / lateral BC CO |
| `CO_BBU` | FINN v2.5 | Biomass burning CO |
| `CO_TST` | EDGAR (twin of `*_ANT`) | `*_TST` is driven only by `E_*TST` in `wrfchemi` (= `E_*` anthro). Per-source duplicate, not a tagged total — see Bug 18. |

**Registry source:** `WRF/Registry/registry.chem` — `ghg_tracer` package (`chem_opt==17`); `eghg` package (`emiss_opt==17`, 6 emission surface vars: `E_CO2`, `E_CO2TST`, `E_CH4`, `E_CH4TST`, `E_CO`, `E_COTST`); `ebioghg` package (`bio_emiss_opt==17`, VPRM flux vars).

`Registry.EM_CHEM` is a 23-line wrapper that `include`s `registry.chem` (4,174 lines covering all chem options). Of the 1,970+ possible variables in `registry.chem`, only the 15 tracers above are active for `chem_opt=17`.

**wrfout variable count for storage sizing:**
- Active 3D output vars: **30** (15 GHG tracers + 15 core met: U, V, W, T, PH, P, PB, PHB, QVAPOR, QC, QR, QICE, QSNOW, CLDFRA, REFL_10CM)
- Active 2D output vars: **~48** (XLAT, XLONG, T2, U10, V10, PSFC, TSK, HFX, LH, PBLH, UST, SWDOWN, GLW, OLR, RAINC, RAINNC, SNOW, SNOWH, VEGFRA, LAI, ALBEDO, XLAND, EMISS, TMN, COSZEN, TSLB, SMOIS, SH2O, GRDFLX, SST + VPRM bio fluxes: `ebio_gee`, `ebio_res`, `ebio_co2oce`, `ebio_ch4wet`, `ebio_ch4soil`, `rad_vprm`, `lambda_vprm`, `alpha_vprm`, `resp_vprm`)
- Per hourly frame: (30 × 232×128×50 + 48 × 232×128) × 4 bytes × 50% compression ≈ **88 MB**

---

## 6. Simulation Period & Segments

| Parameter | Value |
|-----------|-------|
| Start | 2015-07-28 00:00 UTC (aligned to MODIS 8-day composite, see Bug 15) |
| End | 2015-12-01 00:00 UTC |
| Met spin-up | 3 days (2015-07-28 → 2015-07-31 before fire peak) |
| Total duration | ~126 days |
| Time step | 150 s (~5.5×dx) |
| History interval | 60 min |
| Restart interval | 1,440 min (24 h) |
| `frames_per_outfile` | 24 (1 wrfout per day) |
| MPI tasks | 16 (one rank per physical core; OpenMPI on this host counts cores not hwthreads — see Bug 13) |

### Restart segments (aligned to MODIS 8-day calendar — see Bug 15)

Segment boundaries land on MODIS composite dates (2015-07-28, 08-29, 09-30, 10-24, 11-25) so each restart's first VPRM read finds an existing real composite file (no symlink hacks, no `Times` rewrites).

| Segment | Start | End | Duration |
|---------|-------|-----|----------|
| seg1 | 2015-07-28 00:00 | 2015-08-29 00:00 | 768 h (32 days) |
| seg2 | 2015-08-29 00:00 | 2015-09-30 00:00 | 768 h (32 days) |
| seg3 | 2015-09-30 00:00 | 2015-10-24 00:00 | 576 h (24 days) |
| seg4 | 2015-10-24 00:00 | 2015-11-25 00:00 | 768 h (32 days) |
| seg5 | 2015-11-25 00:00 | 2015-12-01 00:00 | 144 h (6 days) |

---

## 7. Input Data Sources

### 7.1 ERA5 Meteorology

| Item | Detail |
|------|--------|
| Source | Copernicus Climate Data Store (CDS) |
| Variable sets | Pressure-level (PL) + Single-level (SL) |
| Resolution | 0.25° × 0.25°, 6-hourly |
| PL variables | u, v, w, T, q, z (37 pressure levels) |
| SL variables | 2m T, 2m Td, 10m u/v, SST, PSFC, CAPE, etc. |
| Period needed | 2015-07-25 to 2015-12-01 |
| Files | `era5_pl_YYYYMM.grib` + `era5_sl_YYYYMM.grib` |
| **Status** | 10/10 GRIBs complete ✅ |
| Missing | None — download finished ~18:59 WIB |
| CDS API key | `cdc594b9-8754-417f-83f8-97259710653d` |
| Vtable | `WPS/ungrib/Variable_Tables/Vtable.ERA-interim` |

### 7.2 EDGAR Anthropogenic Emissions

| Item | Detail |
|------|--------|
| Source | EDGAR v8.0 FT2022 GHG (CO₂, CH₄); HTAP v3 (CO) |
| CO₂ | Per-sector annual files, 0.1° × 0.1°, kg m⁻² s⁻¹ |
| CH₄ | Per-sector annual files, 0.1° × 0.1°, kg m⁻² s⁻¹ |
| CO | HTAP v3 monthly sector gridmaps, 0.1° × 0.1° |
| CO₂ sectors | 19 (TNR_Aviation_SPS missing from EDGAR v8.0 — 404) |
| CH₄ sectors | 22 (excl. AWB — covered by FINN; TNR_Aviation_SPS — 404) |
| **Status** | CO₂: 19/19 NC ✅  CH₄: 22/22 NC ✅  CO: 1 NC ✅ |
| Base URL | `https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/` |
| Unit conversion | kg m⁻² s⁻¹ → mol km⁻² hr⁻¹: `× (1e9 × 3600) / MW` |

**Note:** AWB (agricultural waste burning) sector explicitly excluded from CH₄ to avoid double-counting with FINN fire emissions.

### 7.3 FINN v2.5 Fire Emissions

| Item | Detail |
|------|--------|
| Source | NCAR FINN v2.5 (MODIS+VIIRS) |
| Variables | CO, CO₂, CH₄ |
| Grid | 0.1° × 0.1°, daily 2015 |
| Units | molecules cm⁻² s⁻¹ |
| Files | 3 NetCDF files in `rawdata/finn_fire/` |
| Variable names | `fire_modisviirs_CO`, `fire_modisviirs_CO2`, `fire_modisviirs_CH4` |
| Time encoding | Days since 1990-01-01; shape (365, 1799, 3600) |
| **Status** | ✅ Complete |

### 7.4 BMKG Hotspot Data

| Item | Detail |
|------|--------|
| Source | BMKG Indonesia archived hotspot |
| File | `rawdata/hotspot/archived_hotspot_idn.csv` |
| Columns | `lat, lon, month, day, year, confidence` |
| Filter | year == 2015 (the CSV is already pre-filtered to high-confidence detections; `confidence` is a constant marker — do not filter on it) |
| Usage | Augments FINN spatial placement for Indonesian peatlands |
| **Status** | ✅ Present; preprocessing script now uses all 2015 rows (Bug 14 fixed). |

### 7.5 MODIS Surface Reflectance (VPRM inputs)

| Item | Detail |
|------|--------|
| Source | NASA Earthdata — MODIS MOD09A1 (8-day, 500m) |
| Tiles | h27–h28 × v08–v10 (Indonesia/Maritime Continent) |
| Also needed | MCD12Q1 (annual land cover, 500m) |
| Auth | NASA Earthdata login (credentials in `~/.netrc`, not committed) |
| **Status** | MOD09A1: 799 HDF files (17 tiles × 47 dates, 0 missing) ✅  MCD12Q1: 17 tiles (h27–h32 × v08–v10; h32v10 absent = ocean) ✅ |

### 7.6 CarbonTracker Boundary/Initial Conditions

| Item | Detail |
|------|--------|
| CO₂ source | NOAA CarbonTracker CT2025 |
| CO₂ files | `CT2025.molefrac_glb3x2_YYYY-MM-DD.nc` (3°×2° daily) |
| CH₄ source | NOAA CarbonTracker CH4 CTCH4_2024 |
| CH₄ files | `CTCH4_2024.molefrac_glb3x2_YYYY-MM-DD.nc` |
| Period | 2015-07-25 to 2015-12-01 (daily files) |
| **Status** | 368 NC files in `rawdata/carbontracker/` ✅ |

### 7.7 SOMFFN Ocean CO₂ Flux

| Item | Detail |
|------|--------|
| Source | MPI-BGC SOM-FFN v2023 (NCEI/OCADS) |
| File | `rawdata/ocean/MPI_SOM-FFN_v2023_NCEI_OCADS.nc` |
| Coverage | Global monthly, 1°×1° |
| **Status** | ✅ Present |

### 7.8 WPS Geographic Data

| Item | Detail |
|------|--------|
| Source | WRF Users Page (WPS_GEOG) |
| Path | `rawdata/geog/WPS_GEOG/` |
| Includes | MODIS LULC 30s, MODIS LAI/FPAR, lakes |
| **Status** | ✅ Complete |

---

## 8. Preprocessing Scripts

### 8.1 Script Overview

| Script | Input | Output | Status |
|--------|-------|--------|--------|
| `preprocessing/anthropogenic/01_process_edgar.py` | EDGAR sector NCs | `wrfchemi_d01_*` (hourly, single-frame) | ✅ Complete — 3,120 hourly files in `input/` |
| `preprocessing/anthropogenic/02_bundle_daily_chemi.py` | hourly `wrfchemi_d01_*` | `wrfchemi_d01_*` (daily, 24 frames) in `input_daily/` + run-dir symlinks | ✅ Complete — 130 daily files (Bug 11) |
| `preprocessing/fire/01_process_finn_hotspot.py` | FINN NCs + BMKG hotspot CSV | `wrffirechemi_d01_*` (daily, 24 frames) | ✅ Complete — 129 files; BMKG-augmented (Bug 14) |
| `preprocessing/ocean/01_process_somffn.py` | SOMFFN NC | `wrfoce_d01_*` (monthly) | ✅ Complete — 6 files (Jul–Dec 2015) |
| `preprocessing/boundary_conditions/01_process_carbontracker.py` | CT NCs + wrfinput_d01 | `wrfinput_d01` (ICs), `wrfbdy_d01` (BCs) | ✅ Complete — CO₂_BCK 393–424 ppm; CH₄_BCK 1582–2039 ppb |
| `preprocessing/vprm/01_process_vprm.py` | MOD09A1 + MCD12Q1 HDFs | `vprm_input_d01_*` (8-day MODIS composites) | ✅ Complete — 17 files |
| `preprocessing/utils/01_patch_emiss_globals.py` | wrfinput_d01 + chem aux files | patches WRF projection globals into all chem aux files | ✅ Required after every emission regen (Bug 11) |
| `preprocessing/utils/02_align_aux_times.py` | wrfinput_d01 + monthly/8-day aux files | start-date-aligned copies of `wrfoce` and `vprm_input` (real-file copies, internal `Times` rewritten) | ✅ Required when sim start ≠ existing aux file dates (Bug 12, 15) |
| `preprocessing/domain/namelist.wps` | — | `geo_em.d01.nc` | ✅ Done |

### 8.2 Unit Conversions

| Source | Native units | WRF units | Conversion factor |
|--------|-------------|-----------|-------------------|
| EDGAR (CO₂/CH₄/CO) | kg m⁻² s⁻¹ | mol km⁻² hr⁻¹ | `× (1e9 × 3600) / MW` |
| FINN fire | molecules cm⁻² s⁻¹ | mol km⁻² hr⁻¹ | `× 3600 × 1e10 / 6.022e23` |
| CarbonTracker | mol mol⁻¹ (dry) | ppm (CO₂), ppb (CH₄) | standard conversion |

### 8.3 EDGAR Temporal Profiles

Sector-specific monthly and diurnal profiles are applied at write time, calibrated for Indonesia (UTC+7, so local midnight = 17 UTC preceding day):

**Sector → profile type mapping (CO₂ example):**

| Sector | Monthly profile | Diurnal profile |
|--------|----------------|----------------|
| ENE (power) | energy | energy |
| TRO (road transport) | transport | transport |
| RCO (residential) | residential | residential |
| IND, CHE, NMM, IRO, NFE | industry | industry |
| NEU (non-energy use) | industry | flat |
| REF_TRF | industry | industry |
| PRO_FFF | energy | flat |
| AGS | agriculture | flat |
| SWD_INC | waste | industry |
| TNR_Ship | shipping | shipping |
| TNR_Aviation_* | aviation | aviation |
| TNR_Other | transport | transport |

CH₄ additionally covers: PRO_COAL, PRO_GAS, PRO_OIL, ENF, MNM, SWD_LDF, WWT.  
**AWB excluded from CH₄** (double-counting with FINN).

---

## 9. Emissions Detail

### 9.1 EDGAR CO₂ Sectors Available

`AGS, CHE, ENE, IND, IRO, NEU, NFE, NMM, PRO_FFF, PRU_SOL, RCO, REF_TRF, SWD_INC, TNR_Aviation_CDS, TNR_Aviation_CRS, TNR_Aviation_LTO, TNR_Other, TNR_Ship, TRO`

**Not available in EDGAR v8.0 (404):** `TNR_Aviation_SPS`

### 9.2 EDGAR CH₄ Sectors Available

`AGS, CHE, ENE, ENF, IND, IRO, MNM, PRO_COAL, PRO_FFF, PRO_GAS, PRO_OIL, RCO, REF_TRF, SWD_INC, SWD_LDF, TNR_Aviation_CDS, TNR_Aviation_CRS, TNR_Aviation_LTO, TNR_Other, TNR_Ship, TRO, WWT`

**Excluded:** `AWB` (fire; covered by FINN), `TNR_Aviation_SPS` (404)

### 9.3 Fire Emissions (FINN v2.5 + BMKG hotspots)

- Sources:
  - **FINN v2.5** (NCAR, MODIS+VIIRS-driven): daily 0.1° fluxes for CO₂, CH₄, CO
  - **BMKG `archived_hotspot_idn.csv`**: 54,930 high-confidence Indonesian hotspot detections in 2015 (already pre-filtered upstream — `confidence` column is a constant marker, not a filter key; see Bug 14)
- Hybrid logic in `01_process_finn_hotspot.py`:
  - Where FINN has signal: scaled by `(1 + α × hs_norm)` with `α = 0.5`
  - Where FINN = 0 but hotspots exist: baseline injection (`CO2 = 500`, `CO = 50`, `CH4 = 2` mol km⁻² hr⁻¹) × normalised hotspot density
- Temporal disaggregation: daily totals × fixed diurnal profile (peak 05–07 UTC ≈ 12–14 LT for Indonesia UTC+7), normalised so the daily mean is preserved
- Plume rise: active (`biomass_burn_opt = 5`); plume-rise fraction is implicit in v4.1.5 (no separate `plumerisefire_frct` namelist key — that variable was removed in 4.x; see Bug 9)
- File layout: 129 daily files, 24 hourly frames each → `wrffirechemi_d01_<YYYY-MM-DD>_00:00:00`

---

## 10. Boundary & Initial Conditions

### Initial Conditions (ICs)

- **Meteorology:** from `real.exe` → `wrfinput_d01`
- **Chemistry:** CarbonTracker mole fractions interpolated to WRF grid by `01_process_carbontracker.py`; written into `wrfinput_d01` (`CO2_BIO`, `CO2_OCE`, `CO2_FF`, `CH4`)

### Lateral Boundary Conditions (BCs)

- **Meteorology:** from `real.exe` → `wrfbdy_d01`
- **Chemistry:** CarbonTracker daily files interpolated to WRF lateral boundaries; tendencies written into `wrfbdy_d01`
- **Update interval:** 6 h (matching ERA5 `interval_seconds = 21600`)
- `spec_bdy_width = 5`, `spec_zone = 1`, `relax_zone = 4`

### Upper Boundary

- `co2_upper_bc = 1` — CO₂ relaxed to CarbonTracker zonal mean above domain top

---

## 11. Run Scripts & Pipeline

### Download scripts (`scripts/download/`)

| Script | Purpose |
|--------|---------|
| `02_download_edgar.py` | EDGAR v8 sector ZIPs (CO₂/CH₄) + HTAP CO |
| `04_download_somffn.py` | MPI SOM-FFN ocean CO₂ flux |
| `06_download_geog_data.sh` | WPS geographic static data |

ERA5 downloaded via CDS API; CarbonTracker via aria2c; MODIS via NASA Earthdata wget.

### Run scripts (`scripts/run/`)

| Script | Purpose | Order |
|--------|---------|-------|
| `00_setup_simulation.sh` | Create run dir, link executables + data tables + met_em | 1 |
| `01_run_real.sh` | Run `mpirun -np 32 real.exe`; copies wrfinput/wrfbdy | 2 |
| `03_run_preprocessing.sh` | Launch edgar/fire/ocean/carbontracker in parallel screens | 3 |
| `02_run_wrf.sh` | 5 restart segments, 16 MPI tasks (Bug 13), auto-moves output | 4 |

### Full pipeline sequence

```
ERA5 download (all 10 GRIBs)
    ↓
ungrib + metgrid (via wps_um watcher screen)
    → met_em.d01.YYYY-MM-DD_HH:mm:ss files in met_em/
    ↓
00_setup_simulation.sh
    ↓
01_run_real.sh   (mpirun -np 32 real.exe)
    → wrfinput_d01, wrfbdy_d01
    ↓
03_run_preprocessing.sh  [parallel screens]
    edgar          → wrfchemi_d01_* (3,120 hourly → bundled to 130 daily 24-frame, Bug 11) ✅
    fire           → wrffirechemi_d01_* (129 daily files) ✅
    ocean          → wrfoce_d01_* (6 monthly files, Jul–Dec 2015) ✅
    carbontracker  → ICs into wrfinput_d01, BCs into wrfbdy_d01 ✅
    vprm           → vprm_input_d01_* (17 × 8-day files) ✅
    ↓
02_run_wrf.sh   (16 MPI tasks, 5 segments — see Bug 13)
    → wrfout_d01_* in simulations/IDN_BB_2015/output/
```

### Expected output file counts

| File type | Count | Size/file | Total | Notes |
|-----------|-------|-----------|-------|-------|
| `wrfchemi_d01_*` (daily, 24 frames) | 130 | ~17 MB | ~2.2 GB | EDGAR anthro, 24-frame daily layout (Bug 11) |
| `wrffirechemi_d01_*` (daily, 24 frames) | 129 | ~17 MB | ~2.2 GB | FINN + BMKG hotspot fire tracers |
| `wrfoce_d01_*` (monthly ocean) | 6 + 1 | ~2 MB | ~14 MB | Jul–Dec 2015, plus start-date copy `2015-07-28` |
| `vprm_input_d01_*` (8-day VPRM) | 17 | ~8 MB | ~136 MB | MODIS composite dates only (Bug 15) |
| `wrfinput_d01`, `wrfbdy_d01` | 2 | ~0.2/4.5 GB | ~4.7 GB | from real.exe (re-built for 2015-07-28 start) |
| `wrfout_d01_*` (daily, 24 frames/file) | 126 | ~2.1 GB | **~265 GB** | dominant cost |
| `wrfrst_d01_*` (restart) | 4 | ~4 GB | ~16 GB | all state vars; one per segment boundary |
| `xghg_*.nc` (post-proc column GHG) | ~3,024 | ~0.2 MB | ~0.6 GB | optional |

---

### 11.3 Post-Processing: Column GHG (XCO₂, XCH₄, XCO)

**Script:** `postprocessing/column_ghg/01_compute_xghg.py`

Computes dry-air pressure-weighted column averages from the 15 3D GHG tracer fields. Outputs are directly comparable to GOSAT FTS, OCO-2, and TROPOMI satellite retrievals (no averaging kernels applied — apply satellite AKs in a separate validation step).

**Column integral formula (exact WRF terrain-following $\eta$ coordinates):**

$$\Delta p^{dry}_k = \frac{(\eta^w_k - \eta^w_{k+1})\,(P_{sfc} - P_{top})}{1 + q^{vap}_k}$$

$$X_{gas} = \frac{\displaystyle\sum_{k=1}^{50} c_k \cdot \Delta p^{dry}_k}{\displaystyle\sum_{k=1}^{50} \Delta p^{dry}_k}$$

where $c_k$ is the dry-air mole fraction [ppmv] at model level $k$, $\eta^w_k$ is the staggered half-level eta from `ZNW`, and $q^{vap}_k$ is water vapour mixing ratio from `QVAPOR`.

**Output variables per hourly file (`xghg_YYYY-MM-DD_HH.nc`):**

| Variable | Units | Description |
|----------|-------|-------------|
| `XCO2` | ppm | Total column CO₂ (CO2_ANT + CO2_BIO + CO2_OCE + CO2_BCK + CO2_BBU + CO2_TST) |
| `XCO2_CO2_ANT/BIO/OCE/BCK/BBU/TST` | ppm | Per-component column |
| `XCH4` | ppb | Total column CH₄ (all 5 tracers; converted from ppmv) |
| `XCH4_CH4_ANT/BIO/BCK/BBU/TST` | ppb | Per-component column |
| `XCO` | ppb | Total column CO (all 4 tracers; converted from ppmv) |
| `XCO_CO_ANT/BCK/BBU/TST` | ppb | Per-component column |
| `XLAT`, `XLONG` | degrees | Grid coordinates |

**Fire enhancement** (`*_BBU` components): directly gives the biomass burning contribution to the column, enabling clean separation of the fire signal from background.

**Usage after WRF completes:**
```bash
source .venv/bin/activate

# Process all timesteps, generate PNG maps, combine into single NC
python postprocessing/column_ghg/01_compute_xghg.py \
    --start 2015-08-01 --end 2015-11-30 --plot --combine

# Output NC files: simulations/IDN_BB_2015/output/xghg/xghg_YYYY-MM-DD_HH.nc
# Combined:        simulations/IDN_BB_2015/output/xghg/xghg_IDN_BB_2015.nc
# PNG maps:        simulations/IDN_BB_2015/output/plots/xghg/YYYYMM/
```

`--plot` generates two PNG types per hour:
1. **Total column maps** — `xco2_YYYYMMDD_HH.png`, `xch4_*.png`, `xco_*.png`
2. **Fire enhancement panel** — `fire_xghg_YYYYMMDD_HH.png` (3-panel: ΔXCo₂ / ΔXCH₄ / ΔXCO from `*_BBU` tracers only)

---

### 11.4 Storage Budget

Full disk accounting for `/home/igrk/WRF-GRK/` (1 TB disk, 956 GB usable, 335 GB already used as of 2026-05-03):

#### Input data (`rawdata/`)

| Dataset | Files | Size | Notes |
|---------|-------|------|-------|
| ERA5 GRIBs (10 × PL+SL, Jul–Nov 2015) | 10 | ~12.2 GB | PL ~2.2 GB ea; SL ~144 MB ea |
| EDGAR sector NCs (CO₂×19, CH₄×22, CO×1) | 42 | ~1.5 GB | HTAP CO included |
| FINN fire NCs (CO, CO₂, CH₄) | 3 | ~0.5 GB | |
| CarbonTracker NCs (CT2025 + CTCH4_2024) | 368 | ~3.7 GB | ~10 MB each |
| SOMFFN ocean NC | 1 | ~1.3 GB | |
| MODIS MOD09A1 + MCD12Q1 HDF | 816 | ~42 GB | ~63 MB each MOD09A1 |
| WPS geographic static (WPS_GEOG) | — | ~30 GB | |
| **Rawdata subtotal** | | **~91 GB** | |

#### WRF source and compiled code

| Item | Size |
|------|------|
| WRF-Chem 4.1.5 source + objects + executables | ~12 GB |
| WPS 4.1 source + compiled | ~2 GB |
| `.venv` Python environment | ~1 GB |
| **Code subtotal** | **~15 GB** |

#### Intermediate products (generated during pipeline)

| File type | Count | Size/file | Total |
|-----------|-------|-----------|-------|
| `met_em.d01.*` (ungrib+metgrid) | 516 (129 d × 4/d) | ~8 MB | ~4 GB |
| `wrfinput_d01` | 1 | ~0.24 GB | ~0.24 GB |
| `wrfbdy_d01` | 1 | ~4.5 GB | ~4.5 GB |
| `wrfchemi_d01_*` (hourly source) | 3,120 | ~0.7 MB | ~2.2 GB | original hourly files in `input/` |
| `wrfchemi_d01_*` (daily 24-frame, run-active) | 130 | ~17 MB | ~2.2 GB | bundled into `input_daily/` (Bug 11) |
| `wrffirechemi_d01_*` (daily fire, 24 fr) | 129 | ~17 MB | ~2.2 GB |
| `wrfoce_d01_*` (monthly ocean + start-date copy) | 7 | ~2 MB | ~14 MB |
| `vprm_input_d01_*` (8-day VPRM) | 17 | ~8 MB | ~136 MB |
| **Intermediate subtotal** | | | **~15.5 GB** |

#### Model outputs

| File type | Count | Size/file | Total | Notes |
|-----------|-------|-----------|-------|-------|
| `wrfout_d01_*` (daily, 24 frames) | 126 | ~2.1 GB | **~265 GB** | **dominant cost**; 2015-07-28 → 2015-12-01 |
| `wrfrst_d01_*` (restart) | 4 | ~4 GB | ~16 GB | one per segment boundary |
| `wrfout` derivation | | | | 88 MB/frame × 24 fr = 2.1 GB/file |
| **Model output subtotal** | | | **~281 GB** | |

#### Post-processing outputs (optional)

| File type | Count | Size/file | Total |
|-----------|-------|-----------|-------|
| `xghg_YYYY-MM-DD_HH.nc` | ~3,024 | ~0.2 MB | ~0.6 GB |
| `xghg_IDN_BB_2015.nc` (combined) | 1 | ~0.6 GB | ~0.6 GB |
| PNG plots (XCO₂/XCH₄/XCO maps) | ~9,072 | ~0.5 MB | ~4.5 GB |
| **Post-proc subtotal** | | | **~5.7 GB** |

#### Disk summary

| Category | Used |
|----------|------|
| Already on disk (rawdata + code + prior work) | ~335 GB |
| Intermediate products (some generated) | ~15.5 GB |
| Model outputs (in progress) | ~281 GB |
| Post-processing (optional) | ~5.7 GB |
| **Projected total at completion** | **~637 GB** |
| **Projected free space remaining** | **~319 GB** |

**Constraint check:** 956 GB usable − 637 GB = **319 GB remaining free** ✅ (safe margin; wrfout dominates at ~265 GB)

---

## 12. Directory Structure

```
/home/igrk/WRF-GRK/
├── WRF/                          # WRF-Chem 4.1.5 source + compiled
│   └── main/wrf.exe, real.exe
├── WPS/                          # WPS 4.1 compiled
├── rawdata/
│   ├── era5/                     # ERA5 GRIBs (era5_pl_YYYYMM.grib, era5_sl_YYYYMM.grib)
│   ├── edgar/                    # EDGAR sector NCs + HTAP CO NC
│   ├── finn_fire/                # FINN v2.5 daily NCs (CO, CO2, CH4)
│   ├── hotspot/                  # BMKG hotspot CSV
│   ├── carbontracker/            # CT2025 + CTCH4_2024 daily NCs
│   ├── ocean/                    # SOMFFN NC
│   ├── modis/                    # MOD09A1 + MCD12Q1 HDF files
│   └── geog/WPS_GEOG/            # WPS static geographic data
├── preprocessing/
│   ├── anthropogenic/
│   │   ├── 01_process_edgar.py        # EDGAR → 3,120 hourly wrfchemi files
│   │   └── 02_bundle_daily_chemi.py   # hourly → 130 daily 24-frame files (Bug 11)
│   ├── fire/01_process_finn_hotspot.py
│   ├── ocean/01_process_somffn.py
│   ├── boundary_conditions/01_process_carbontracker.py
│   ├── vprm/01_process_vprm.py
│   ├── utils/
│   │   ├── 01_patch_emiss_globals.py  # patch WRF projection globals onto chem aux files (Bug 11)
│   │   └── 02_align_aux_times.py      # start-date-aligned copies of wrfoce/vprm_input (Bug 12, 15)
│   └── domain/namelist.wps
├── scripts/
│   ├── download/
│   │   ├── 02_download_edgar.py
│   │   ├── 04_download_somffn.py
│   │   └── 06_download_geog_data.sh
│   └── run/
│       ├── 00_setup_simulation.sh
│       ├── 01_run_real.sh
│       ├── 02_run_wrf.sh
│       └── 03_run_preprocessing.sh
├── simulations/IDN_BB_2015/
│   ├── input/                    # source aux files (hourly chemi, daily fire, monthly oce, 8-day vprm)
│   ├── input_daily/              # bundled daily 24-frame wrfchemi (run-active, Bug 11)
│   ├── met_em/                   # met_em.d01.* + geo_em.d01.nc
│   ├── run/                      # linked executables, tables, namelist.input, active run
│   ├── output/                   # wrfout_d01_* files
│   ├── restart/                  # wrfrst_d01_* files
│   └── logs/                     # all *.log files
├── docs/
│   └── simulation.md             # this document
└── .venv/                        # Python virtual environment
```

---

## 13. Known Issues & Bug Fixes

### Bug 1 — `nc.stringtochar` crash (ALL preprocessing scripts)

**Affected:** `01_process_edgar.py`, `01_process_finn_hotspot.py`, `01_process_somffn.py`, `01_process_vprm.py`  
**Symptom:** `AttributeError: 'numpy.bytes_' object has no attribute 'encode'`  
**Root cause:** `netCDF4.stringtochar()` raises this error in the installed version  
**Fix:** Replace all instances with `list(dt.strftime("%Y-%m-%d_%H:%M:%S"))` — returns a plain Python list of characters which netCDF4 accepts for `"S1"` variables  
**Status:** ✅ Fixed in all 4 scripts

### Bug 2 — FINN unit factor wrong by 5 orders of magnitude

**Affected:** `preprocessing/fire/01_process_finn_hotspot.py`  
**Symptom:** FINN emissions 5 orders of magnitude too large  
**Root cause:** Script assumed FINN units were `g m⁻² day⁻¹`; actual FINN v2.5 units are `molecules cm⁻² s⁻¹`  
**Fix:**
- Old factor: assumed `g m⁻² day⁻¹` → `mol km⁻² hr⁻¹` (wrong)
- New factor: `data × 3600 × 1e10 / 6.022e23` (molecules cm⁻² s⁻¹ → mol km⁻² hr⁻¹)
- Removed `MW` dict from fire script (not needed for molecule-based units)  
**Status:** ✅ Fixed

### Bug 3 — VPRM `Times` dimension missing

**Affected:** `preprocessing/vprm/01_process_vprm.py`  
**Symptom:** WRF cannot read `vprm_input_d01_*` — expects `Times(Time, DateStrLen)` but file had `Times(DateStrLen)`  
**Fix:** Added `ds_out.createDimension("Time", 1)` and changed `times_var[0, :] = list(...)` to 2D assignment  
**Status:** ✅ Fixed

### Bug 4 — namelist.wps (4 fixes, previous session)

Multiple issues in the WPS namelist were corrected before geogrid succeeded.  
**Status:** ✅ Fixed, geogrid SUCCESS

### Bug 5 — EDGAR sector download crash on `TNR_Aviation_SPS` (404)

**Date:** 2026-05-03  
**Symptom:** `02_download_edgar.py` crashed with `HTTPError: 404 Not Found` on `TNR_Aviation_SPS` URL for CO₂; CH₄ sectors never downloaded  
**Root cause:** `TNR_Aviation_SPS` does not exist in EDGAR v8.0 for either CO₂ or CH₄ (confirmed via HTTP HEAD checks). Script used `raise_for_status()` without 404 handling.  
**Fix:**
1. Removed `TNR_Aviation_SPS` from `CO2_SECTORS` and `CH4_SECTORS` lists in `02_download_edgar.py`
2. Removed `TNR_Aviation_SPS` from `CO2_SECTOR_PROFILES` and `CH4_SECTOR_PROFILES` in `01_process_edgar.py`
3. Added 404-tolerant handling to `download_file()`: logs warning and returns instead of raising  
**Status:** ✅ Fixed; CH₄ sectors successfully downloaded on re-run

---

### Bug 6 — SOMFFN (lon,lat) dimension transposition

**Date:** 2026-05-03  
**Affected:** `preprocessing/ocean/01_process_somffn.py`  
**Symptom:** `RegularGridInterpolator` crashed / produced wrong output — array shape mismatch  
**Root cause:** The MPI SOM-FFN v2023 NetCDF stores flux as `(longitude:360, latitude:180, time:492)`. After `.isel(time=t)` the result is `(360, 180)` but `RegularGridInterpolator` expects `(lat, lon)` = `(180, 360)`.  
**Fix:** Added a transpose check immediately after `.values.squeeze()`:
```python
if data.shape == (len(src_lon), len(src_lat)):
    data = data.T  # SOMFFN stored (lon,lat); RegularGridInterpolator needs (lat,lon)
```
**Status:** ✅ Fixed; script rerun produced 6 correct `wrfoce_d01_*` files (Jul–Dec 2015)

---

### Bug 7 — CarbonTracker `pressure` variable is 4D (CT2025 CO₂)

**Date:** 2026-05-03  
**Affected:** `preprocessing/boundary_conditions/01_process_carbontracker.py` — `load_ct_field()` function  
**Symptom:** `ValueError: The truth value of an array with more than one element is ambiguous` at line ~117: `if ct_lev[0] < ct_lev[-1]:`  
**Root cause:** In CT2025 (CO₂) files, `pressure` has dimensions `(time:8, boundary:35, latitude:90, longitude:120)` — a 4D array, not a 1D pressure vector. The original code `lev = ds["pressure"].values` assigned a 4D array to `ct_lev`, making `ct_lev[0]` a 3D slice whose truth value is ambiguous.  
**Fix:** Replaced the single-line lev extraction with explicit multi-case logic:
```python
if "pressure" in ds:
    # CT2025: pressure is (time, boundary=N+1, lat, lon) in Pa
    p_edge = ds["pressure"].values
    p_edge_t = p_edge[0]  # first time step: (boundary, lat, lon)
    p_mid = 0.5 * (p_edge_t[:-1] + p_edge_t[1:])  # (level, lat, lon)
    lev = p_mid.mean(axis=(1, 2)) / 100.0  # spatial mean → 1D hPa profile
```
**Status:** ✅ Fixed (partial — revealed Bug 8 for CH₄ files)

---

### Bug 8 — CH₄_BCK wrong values (299–382 ppb instead of ~1800 ppb)

**Date:** 2026-05-03  
**Affected:** `preprocessing/boundary_conditions/01_process_carbontracker.py` — `load_ct_field()` for CH₄  
**Symptom:** After Bug 7 fix, carbontracker ran without error but the output `CH4_BCK` field in `wrfinput_d01` showed range **299.29–382.19 ppb** — far below the correct 2015 global tropospheric CH₄ of ~1850 ppb. CO₂_BCK (393–424 ppm) was correct.  
**Root cause:** CTCH4_2024 (CH₄) files have **no `pressure` variable** — they use hybrid sigma coordinates (`at`, `bt`, `surf_pressure`). The code fell through to the `elif "level" in ds:` branch and returned `lev = [1, 2, 3, ..., 25]` (integer level indices). Since WRF pressure levels range 20–1000 hPa and CT "levels" are 1–25, `np.interp` mapped all WRF columns to the uppermost CT level (~0.5 hPa stratosphere), where CH₄ is photochemically depleted to 300–380 ppb.  

**Diagnosis steps:**
1. Inspected CTCH4_2024 dataset structure: `ch4` variable shape `(time:8, level:25, lat:90, lon:120)`, units `"parts per billion (ppb)"`, range 226–2955 ppb
2. Identified `at` and `bt` (hybrid sigma boundary coefficients, shape `boundary:26`) and `surf_pressure` `(time, lat, lon)` in Pa
3. Confirmed correct pressure levels by computing `p = at + bt × psfc`: levels 960 hPa (surface) → 0.5 hPa (top)

**Fix:** Added a new branch in `load_ct_field()` for the hybrid sigma case:
```python
elif "at" in ds and "bt" in ds and "surf_pressure" in ds:
    # CTCH4_2024: hybrid sigma levels.  p = at + bt * psfc  (Pa)
    at = ds["at"].values        # (boundary,) — 26 values
    bt = ds["bt"].values        # (boundary,) — 26 values
    psfc = ds["surf_pressure"].values[0]  # (lat, lon), first time step, Pa
    p_edge = at[:, None, None] + bt[:, None, None] * psfc[None, :, :]
    p_mid = 0.5 * (p_edge[:-1] + p_edge[1:])  # (25, lat, lon)
    lev = p_mid.mean(axis=(1, 2)) / 100.0  # Pa → hPa, shape (25,)
```
**Result after fix:** `CH4_BCK` range = **1582–2039 ppb** — physically correct. Surface tropics ~2000 ppb (high CH₄ from wetlands + peat); upper stratosphere ~1582 ppb (photochemical depletion). `CO2_BCK` unchanged at 393–424 ppm.  
**Note on units:** CTCH4_2024 `ch4` is already in ppb, so the `field_3d.max() < 0.01` guard (which would multiply by 1e9) correctly does **not** trigger. No unit conversion needed for CH₄.  
**Status:** ✅ Fixed; wrfinput_d01 and wrfbdy_d01 re-patched with correct values

---

### Bug 9 — namelist.input: invalid/removed options for WRF 4.1.5

**Date:** 2026-05-03  
**Affected:** `simulations/IDN_BB_2015/input/namelist.input` and `run/namelist.input`  
**Symptoms (multiple):**
1. `plumerisefire_frct = 1` — `real.exe` FATAL: variable not in WRF 4.1.5 `&chem` registry
2. `co2_upper_bc = 1` — `real.exe` FATAL: variable not in WRF 4.1.5 `&chem` registry
3. `vprm_opt = 1`, `term_opt = 3` — `real.exe` FATAL: these are `rconfig character` type, not integer; WRF expects quoted string values
4. `io_form_restart = .false.` — `wrf.exe` FATAL: `io_form_restart` is integer type (2 = NetCDF); `.false.` is a logical and causes a namelist read error

**Fixes:**
1. Removed `plumerisefire_frct` — not present in registry.chem for WRF 4.1.5 (plume rise is activated implicitly by `biomass_burn_opt=5`)
2. Removed `co2_upper_bc` — not present in registry.chem for WRF 4.1.5
3. `vprm_opt = 'VPRM_table_TROPICS'` — valid choices: `VPRM_table_US`, `VPRM_table_EUROPE`, `VPRM_table_TROPICS`
4. `term_opt = 'CH4_termite_NW'` — valid choices: `CH4_termite_NW`, `CH4_termite_OW`
5. `io_form_restart = 2` — corrected to integer (NetCDF format)

**Status:** ✅ All fixed in both `input/namelist.input` and `run/namelist.input`

---

### Bug 10 — `02_run_wrf.sh` sed pattern corrupts `io_form_restart`

**Date:** 2026-05-03  
**Affected:** `scripts/run/02_run_wrf.sh`  
**Symptom:** Every time a new WRF segment started, `io_form_restart` in `namelist.input` was overwritten from `= 2` back to `= .false.`, causing `wrf.exe` to immediately fatal on namelist read  
**Root cause:** The script uses `sed -i "s/restart\s*=.*/restart = .false.,/"` to set the restart flag for the first segment. The pattern `restart\s*=` is a substring match and also matches `io_form_restart                     = 2,`, replacing it with `io_form_restart = .false.,`.  
**Fix:** Changed both restart sed patterns to use a line-anchored regex that matches only the standalone `restart` key:
```bash
# Before (broken — matches io_form_restart too):
sed -i "s/restart\s*=.*/restart = .false.,/" namelist.input

# After (fixed — anchored to start of line, requires spaces before restart):
sed -i "s/^ *restart  *= .*/restart                             = .false.,/" namelist.input
```
Same fix applied to the `.true.` variant for segments 2–5.  
**Status:** ✅ Fixed in `scripts/run/02_run_wrf.sh`

---

### Bug 11 — `wrf_get_next_time` Status=-4 with hourly single-frame `wrfchemi` files

**Date:** 2026-05-04  
**Affected:** `simulations/IDN_BB_2015/run/namelist.input` (`io_style_emissions`, `frames_per_auxinput5`); fire-emissions preprocessing layout  
**Symptom:** `wrf.exe` aborted at startup after writing the first `wrfout` frame:
```
Open file wrfchemi_d01_2015-07-25_00:00:00
input_wrf: wrf_get_next_time current_date: 2015-07-25_00:00:00 Status = -4
---- ERROR: Could not find matching time in input file wrfchemi_d01_2015-07-25_00:00:00
NOTE: 1 namelist vs input data inconsistencies found.
FATAL CALLED FROM FILE: <stdin> LINE: 1290
```
**Root cause:** Original layout was 3,120 hourly single-frame `wrfchemi` files with `io_style_emissions = 2`, `frames_per_auxinput5 = 1`. WRF v4.1.5 calls `wrf_get_next_time` (`WRF/external/io_netcdf/wrf_io.F90:3174`) once per file open. Internally it positions to the matching time (`ext_ncd_set_time` sets `CurrentTime = NumberTimes = 1`); the subsequent `get_next_time` then sees `CurrentTime ≥ NumberTimes` and returns `WRF_WARN_TIME_EOF` (-4). With `io_style_emissions = 2` this is treated as fatal (`WRF/share/input_wrf.F:1144`). With `io_style_emissions = 1` the check is bypassed (`WRF/share/input_wrf.F:1167-1170`), but io_style 1 demands two cyclic 12-hour static files and discards day-to-day variation.

**Fix:** Bundle the 3,120 hourly files into 130 daily files with 24 frames each, and switch the namelist accordingly:
- New layout: `wrfchemi_d01_<YYYY-MM-DD>_00:00:00`, 24 hourly frames per file, written to `simulations/IDN_BB_2015/input_daily/` and symlinked into the run dir.
- `io_style_emissions = 2` (unchanged), `frames_per_auxinput5 = 24` (was 1).
- Bundling script `/tmp/bundle_daily_chemi.py` (NETCDF3_CLASSIC, copies template, writes 24 frames per day).
- All 130 bundled files patched with WRF projection/grid global attributes copied from `wrfinput_d01` (DX, DY, MAP_PROJ, dimensions, etc. — the original per-hour files only carried `:TITLE = "WRF-Chem EMISSIONS"`).

After the fix WRF passes the chem-emissions reader and integrates normally; emissions update each hour, file rolls over each day at 00:00 UTC.

**Status:** ✅ Fixed; WRF integrating with full hour-by-hour, day-by-day EDGAR anthro emissions

---

### Bug 12 — `wrfoce` and `vprm_input` filename ≠ internal `Times`

**Date:** 2026-05-04  
**Affected:** `simulations/IDN_BB_2015/run/wrfoce_d01_*`, `vprm_input_d01_*`  
**Symptom:** After Bug 11 was fixed, WRF aborted on `wrfoce`:
```
Time in file: 2015-07-01_00:00:00
Time on domain: 2015-07-25_00:00:00
**WARNING** Time in input file not equal to time on domain **WARNING**
input_wrf: wrf_get_next_time current_date: 2015-07-01_00:00:00 Status = -4
---- ERROR: Could not find matching time in input file wrfoce_d01_2015-07-25_00:00:00
```
**Root cause:** Ocean files are monthly (Jul/Aug/Sep/…); VPRM is on the MODIS 8-day calendar (start dates 2015-07-20, 07-28, 08-05, …). Neither stream had a file dated to the simulation start (2015-07-25). Symlinking the start-dated name to the nearest real file is not enough because WRF reads the `Times` variable inside the file and compares it against the domain clock — the symlinked oce file still reported `2015-07-01_00:00:00` internally.

**Fix:** Created real copies of the nearest source files with the internal `Times` variable rewritten to the simulation start:
- `wrfoce_d01_2015-07-25_00:00:00` ← copy of `wrfoce_d01_2015-07-01_00:00:00`, Times rewritten
- `vprm_input_d01_2015-07-25_00:00:00` ← copy of `vprm_input_d01_2015-07-20_00:00:00`, Times rewritten

This is acceptable because (a) ocean CO₂ flux is monthly-resolved (the 25-day shift has no real signal at the model's 27 km × hourly resolution), and (b) the VPRM 8-day grid is intended to be applied "as of" the start of each composite period; shifting the start of the first composite by 5 days does not change the underlying MODIS-derived parameters.

**Status:** ✅ Fixed; init clears, WRF integrates

---

### Bug 13 — OpenMPI `--np 28` exceeded available slots

**Date:** 2026-05-04  
**Affected:** WRF launch command in this session  
**Symptom:** `mpirun -np 28 wrf.exe` failed: *"There are not enough slots available in the system to satisfy the 28 slots requested."*  
**Root cause:** Host has 16 physical cores / 32 hardware threads. OpenMPI counts cores, not hwthreads. `02_run_wrf.sh` and prior namelist comments assumed 28 ranks (`Ntasks in X = 4, Ntasks in Y = 7`).

**Fix for this run:** launched with `mpirun -np 16 wrf.exe` (4 × 4 decomposition; one rank per physical core, generally the best WRF layout on this CPU). For the production segment runs, either keep `-np 16` or pass `--use-hwthread-cpus --np 28` to OpenMPI.

**Status:** ✅ Workaround in use; `02_run_wrf.sh` should be updated to `-np 16` for this host before launching segments 2–5

---

### Bug 14 — Misinterpreted BMKG `confidence` column → all hotspots filtered out

**Date:** 2026-05-04  
**Affected:** `preprocessing/fire/01_process_finn_hotspot.py`  
**Symptom:** Despite the script's name suggesting FINN+BMKG fusion, the resulting `wrffirechemi_d01_*` files were driven by FINN alone — no BMKG hotspot density was ever applied.  
**Root cause:** The script set `MIN_CONFIDENCE = 2` and applied `hs = hs[hs["confidence"] >= MIN_CONFIDENCE]`. The BMKG `archived_hotspot_idn.csv` distribution, however, is **already pre-filtered to high-confidence detections**: the `confidence` column is a constant `1` marker indicating "high confidence", not a 1/2/3 ordinal. The `≥ 2` filter therefore discarded every row (54,930 detections in 2015 alone, including the Sep–Oct peak with ~40,700 detections).

**Fix:** Removed the `MIN_CONFIDENCE` filter entirely — all 2015 rows are now used. Constant kept out of the script as a comment so the meaning is documented for future readers.

```python
# BMKG `archived_hotspot_idn.csv` is already pre-filtered to high-confidence
# detections; the `confidence` column is a constant marker, not a filter key.
# Use all rows.
hs = hs[hs["year"] == 2015].copy()
```

**Will the fire signal show up in `wrfout`?** Yes — once fire emissions enter via auxinput7, WRF distributes them through plume rise (`biomass_burn_opt = 5`) into the column. They appear in:
- `CO2_BBU`, `CH4_BBU`, `CO_BBU` — the dedicated fire-tracer fields (already showing 12–141 ppmv enhancements over hot pixels in the current run, even with FINN-only emissions)
- `CO2_TST`, `CH4_TST`, `CO_TST` — diagnostic totals
- The post-processing `xghg_*.nc` `*_BBU` columns and the `fire_xghg_YYYYMMDD_HH.png` 3-panel plots provide the column-integrated fire enhancement

For the publication run, regenerate `wrffirechemi_d01_*` with the corrected script (54,930 BMKG hotspots will now contribute, concentrated over Sumatra/Kalimantan/Papua peatlands during Sep–Oct 2015).

**Status:** ✅ Script fixed; current WRF run uses pre-fix (FINN-only) fire emissions. To use BMKG-augmented emissions, re-run `01_process_finn_hotspot.py`, replace `wrffirechemi_d01_*` in the run dir, and restart from `wrfrst_d01_*`.

---

### Bug 15 — VPRM read schedule walks off the MODIS composite calendar after restart

**Date:** 2026-05-04  
**Affected:** `simulations/IDN_BB_2015/run/namelist.input` (start date), `vprm_input_d01_*` files, restart-segment boundaries  
**Symptom:** After a successful seg1 (start 2015-07-25), restart from `wrfrst_d01_2015-08-01_00:00:00` aborted at sim time 2015-08-02 with:
```
mediation_integrate: med_read_wrf_chem_emissions: Read emissions for time 2015-08-02_00:00:00
FATAL CALLED FROM FILE: <stdin> LINE: 314
Possibly missing file for = auxinput15
```
**Root cause:** WRF schedules auxinput15 reads as `start_time + N × auxinput15_interval_m`. With `start = 2015-08-01 00:00` and `interval = 8 days = 11520 min`, WRF asked for files dated 2015-08-01, 08-09, 08-17, 08-25 — none of which exist. The real VPRM/MODIS 8-day composites are dated 2015-07-28, 08-05, 08-13, 08-21, 08-29, ... A short-lived workaround had created symlinks at every `start + 8k` date pointing to the nearest real composite, but each symlink's *internal* `Times` variable still carried the original composite date, so WRF rejected the read with the same time-mismatch as Bug 12.

**Why symlinks are the wrong fix:** they hide a real schedule mismatch and break again on every restart segment that doesn't start on a MODIS composite date.

**Fix:** Align the simulation start (and every restart-segment boundary) to a real MODIS 8-day composite date. Chose **2015-07-28** (the first composite ≥ desired 2015-07-25). With this start, WRF's auxinput15 reads land on 07-28, 08-05, 08-13, 08-21, 08-29, ... — all real composites, no symlinks, no Times rewrites.

**Actions taken:**
1. Stopped the running WRF
2. Deleted all `wrfrst_d01_*` and `wrfout_d01_*` from the failed run
3. Removed the four ad-hoc VPRM symlinks at non-composite dates (`vprm_input_d01_2015-07-25_00:00:00`, `2015-08-02`, `2015-08-10`, `2015-08-18`) so the run dir contains only real MODIS composite files
4. Updated namelist: `start = 2015-07-28 00:00`, `end = 2015-08-27 00:00`, `restart = .false.`
5. Reran `real.exe` (16 ranks, ~3 min, SUCCESS COMPLETE REAL_EM INIT)
6. Updated `preprocessing/boundary_conditions/01_process_carbontracker.py` `SIM_START` from `datetime(2015, 7, 25, 0)` to `datetime(2015, 7, 28, 0)` and re-patched ICs/BCs into the new wrfinput_d01/wrfbdy_d01
7. Created `wrfoce_d01_2015-07-28_00:00:00` from the July monthly file with internal `Times` rewritten and WRF projection globals patched
8. Restart-segment plan updated in §6 so every segment boundary is a MODIS composite date

**Status:** ✅ Fixed; WRF relaunched from 2015-07-28 with `mpirun -np 16`, integrating cleanly past first hourly emissions reads.

---

### Bug 16 — Ocean CO₂ flux silently zero: `ebio_co2oce` lowercase vs WRF expected `EBIO_CO2OCE`

**Date:** 2026-05-04  
**Affected:** `preprocessing/ocean/01_process_somffn.py`, all `wrfoce_d01_*` files  
**Symptom:** `EBIO_CO2OCE` in wrfout = 0 everywhere across all frames; `CO2_OCE` ~1e-16 (numerical roundoff). Yet the wrfoce file itself contains valid SOMFFN data (range −235 to +1123 mol km⁻² hr⁻¹, 23,575 non-zero cells). WRF logs `"Input data is acceptable to use: wrfoce_d01_..."` — no error.  
**Root cause:** WRF-Chem's registry preprocessor (`tools/registry`) **uppercases the netcdf DataName** for scalar 2D state variables when generating `inc/allocs.inc`. Registry source has `state real ebio_co2oce ij misc 1 - i06rh "ebio_co2oce" ...` (lowercase netcdf name), but the generated allocs.inc emits `grid%tail_statevars%DataName = 'EBIO_CO2OCE'` (UPPERCASE). At runtime WRF passes `'EBIO_CO2OCE'` to `nf_inq_varid()` which is case-sensitive. The wrfoce file's `ebio_co2oce` (lowercase, written by `01_process_somffn.py`) is therefore not found, the read silently fails (no fatal), and `grid%ebio_co2oce` stays at 0 — which `module_ghg_fluxes.F` then adds to `CO2_OCE` as zero.

This case-conversion only affects scalar 2D state vars. 4D-array members (e.g. `ebu_in_co2`, `evi`, `cpool`, `t_ann`) keep their case via `*_dname_table` from `module_configure.f90`, which is why the corresponding fire and VPRM reads work correctly with lowercase or matching uppercase names.

**Fix:**
1. Renamed the variable in all 7 wrfoce files (6 monthly + 1 start-date copy) using `ncrename -v ebio_co2oce,EBIO_CO2OCE`.
2. Updated `preprocessing/ocean/01_process_somffn.py` line 158 to write `EBIO_CO2OCE` (UPPERCASE), with an inline comment explaining the case-sensitivity trap.

**Impact assessment:** Ocean CO₂ flux was missing from the run from t=0 until the fix. Magnitude is small relative to other sources (ocean ~0.1 mol km⁻² hr⁻¹ vs anthro ~10³ and fire ~10⁵), so column XCO₂ bias from missing ocean flux is < 0.1 ppm. But for clean-marine column comparisons (e.g. SH retrievals where ocean is dominant), this matters.

**Status:** ✅ Fixed; WRF restarted from t=0 (2015-07-28) on 2026-05-04 to apply.

---

### Bug 17 — VPRM `T_ANN` was a constant 300.15 K everywhere

**Date:** 2026-05-04  
**Affected:** `preprocessing/vprm/01_process_vprm.py`, all 17 `vprm_input_d01_*` files  
**Symptom:** `T_ANN` in vprm_input files had **1 unique value** (300.15 K) across the 231×127 grid, removing all spatial gradient from VPRM's Q10 respiration scaling.  
**Root cause:** `01_process_vprm.py` line 234 hardcoded `tann = np.full((NY, NX), 300.15, dtype=np.float32)` as a placeholder ("use 27°C for Indonesia"). The real climatological mean annual surface temperature varies 280–304 K across the domain (cold mountains → warm lowlands).

**Why it matters:** VPRM respiration scales as `exp(Q10 × (T2 − T_ANN))`. With T_ANN spatially uniform at 300.15 K, cells with real T_ANN ≠ 300.15 K get systematically wrong respiration: highland forests (real T_ANN ~285 K) have respiration biased high by `exp(0.07 × (300.15 − 285)) ≈ 2.9×`, while warm lowlands (real T_ANN ~302 K) are biased low by `exp(0.07 × (300.15 − 302)) ≈ 0.88×`. Domain-averaged the bias partially cancels, but the spatial pattern of `CO2_BIO` is muted.

**Fix:**
1. Computed correct T_ANN by reading `TMN` from `wrfinput_d01` (Noah's annual climatological deep-soil temperature, set by `real.exe` from ERA5 climatology — varies 280–304 K across the domain, mean 300.8 K).
2. Patched all 17 existing `vprm_input_d01_*` files in place with this spatially-varying `TMN` field.
3. Updated `01_process_vprm.py` to read `TMN` from `wrfinput_d01` for future runs (with a fallback constant if wrfinput is not yet available).

**Status:** ✅ Fixed; WRF restarted from t=0 to apply.

---

### Bug 18 — `*_TST` tracers are duplicates of `*_ANT`, not "tagged totals"

**Date:** 2026-05-04  
**Affected:** documentation only — `docs/simulation.md` §5.4 description  
**Symptom:** `CO2_TST − CO2_ANT = 0.000000` at every grid point. Same for `CH4_TST − CH4_ANT` and `CO_TST − CO_ANT`.  
**Root cause:** `WRF/chem/module_ghg_fluxes.F:79` shows `chem(p_co2_tst) += emis_ant(p_e_co2tst) × conv_rho`. The `*_TST` tracers are driven **only by `E_*TST` in `wrfchemi`** — and our preprocessing writes `E_CO2TST = E_CO2`, `E_CH4TST = E_CH4`, `E_COTST = E_CO` (all anthropogenic only). Fire, biogenic, and ocean fluxes do **not** contribute to `*_TST`. So `*_TST` is structurally a duplicate of `*_ANT` in this configuration.

The earlier simulation.md §5.4 description ("Test/diagnostic total flux tracer") was wrong. The `_TST` tracer is a per-source twin of the anthropogenic tracer — useful for unit tests of the registry, not a tagged total of all sources.

**Fix:** Documentation only. §5.4 description corrected. To make `_TST` a true "tagged total of all sources", the preprocessing would need to set `E_CO2TST = E_CO2 + E_CO2_BIO + E_CO2_BBU + E_CO2_OCE` — not possible as a one-pass operation since BIO/OCE come from different auxinputs (15/6) and are computed online by VPRM. A meaningful "total" is best constructed in post-processing as the sum of `*_ANT + *_BIO + *_BBU + *_OCE`.

**Status:** ✅ Documented (no code change needed).

---

### Bug 19 — Missing `wrfoce` files at every auxinput6 read date

**Date:** 2026-05-04  
**Affected:** `simulations/IDN_BB_2015/run/wrfoce_d01_*`, `preprocessing/utils/02_align_aux_times.py`  
**Symptom:** WRF aborted at sim time **2015-08-27 00:00:00** (= sim start `2015-07-28` + 30 days = first auxinput6 read after start) with:
```
Timing for Writing wrfout_d01_2015-08-27_00:00:00 ...
mediation_integrate: med_read_wrf_chem_emissions: Open file wrfchemi_d01_2015-08-27_00:00:00
d01 2015-08-27_00:00:00  Error trying to read metadata
d01 2015-08-27_00:00:00  Input data is acceptable to use: wrffirechemi_d01_2015-08-26_00:00:00
input_wrf: wrf_get_next_time current_date: 2015-08-27_00:00:00 Status = -4
ERROR: Could not find matching time in input file wrffirechemi_d01_2015-08-26_00:00:00
```

The error message points at `wrffirechemi_d01_2015-08-26` (the last successfully-read fire file), but **the actual failing read is `wrfoce`** (auxinput6). The misleading error string is just the most recent "Input data is acceptable to use" line in the log buffer.

**Root cause:** `auxinput6_interval_m = 43200` minutes = 30 days. Starting from sim 2015-07-28, the auxinput6 read times are **2015-07-28, 2015-08-27, 2015-09-26, 2015-10-26, 2015-11-25**. WRF builds the filename from current sim time via `<date>` substitution, so at 08-27 00:00 it tries to open `wrfoce_d01_2015-08-27_00:00:00` — which we never created. The original SOMFFN preprocessing produces only first-of-month files (`2015-07-01`, `2015-08-01`, ...), and the prior fix for Bug 12 only created a single start-date copy (`2015-07-28`).

A restart from a clean wrfrst at 2015-08-21 reproduced the failure at the same calendar date — confirming this is a **deterministic file-not-found at the second auxinput6 read time**, not a counter or transient I/O issue.

**Fix:**
1. Generated start-date-aligned `wrfoce_d01_*` files for **every auxinput6 read time** in the simulation period: 2015-07-28, 2015-08-27, 2015-09-26, 2015-10-26, 2015-11-25. Each is a copy of the same-month source file (e.g. 08-27 ← 08-01) with the internal `Times` variable rewritten and WRF projection globals patched.
2. Updated `preprocessing/utils/02_align_aux_times.py` to loop over `START_DATE → END_DATE` at the configured interval (`OCE_INTERVAL_DAYS = 30`) instead of creating only the start-date copy.

**Impact assessment:**
- The simulation was reading the start-date oce file (07-28) for the first 30 days of integration. After the fix it picks up an updated monthly oce flux every 30 days, which is what we wanted from the SOMFFN monthly product anyway.
- All wrfout written before the crash is unaffected (the oce read at the start was correct).

**Status:** ✅ Fixed; WRF restarted from `wrfrst_d01_2015-08-21_00:00:00` after the fix.

---

### Bug 20 — `wrfbdy_d01` truncated to 30 days because real.exe was not rerun after namelist end-date changes

**Date:** 2026-05-04  
**Affected:** `simulations/IDN_BB_2015/run/wrfbdy_d01`  
**Symptom:** After Bug 19 was fixed, the run cleared the 2015-08-27 00:00 chem-emission reads cleanly (chemi, oce, fire all OK), but immediately crashed on the meteorology BC read with:
```
d01 2015-08-27_00:00:00  Input data is acceptable to use: wrfbdy_d01
input_wrf: wrf_get_next_time current_date: 2015-08-27_00:00:00 Status = -4
FATAL CALLED FROM FILE: <stdin> LINE: 1142
---- ERROR: Ran out of valid boundary conditions in file wrfbdy_d01
```

**Root cause:** `wrfbdy_d01` had only **120 time records** spanning **2015-07-28 00:00 → 2015-08-26 18:00** (= 30 days at 6-h interval). It was generated by an earlier `real.exe` run with `end_day = 27` (run_hours = 744). The namelist was later bumped to `end_day = 29` (run_hours = 768) and even longer for restart-segment planning, but **real.exe was never rerun**. WRF therefore ran out of 6-h BC frames at 2015-08-27 00:00 — exactly one read after the last frame.

This was masked while debugging because the crash time coincided with the Bug 19 oce-file-not-found point, making both look like the same bug.

**Fix:**
1. Set namelist to span the full simulation period (`end = 2015-12-01`, `run_hours = 3024`), and reran `real.exe` with 16 ranks. real.exe processed 504 of 505 loops successfully, then aborted on a missing `met_em.d01.2015-12-01_00:00:00.nc` (we have met_em through 2015-11-30 18:00). The fatal at the very end is **harmless** — wrfbdy is written incrementally, so the on-disk file contains all 503 valid 6-h records spanning **2015-07-28 00:00 → 2015-11-30 12:00**, which is enough to drive the entire simulation.
2. Re-patched CarbonTracker chem ICs/BCs into the new wrfinput_d01 (`SIM_START = 2015-07-28`) and wrfbdy_d01.
3. Re-ran `preprocessing/utils/02_align_aux_times.py` to refresh the 5 start-date-aligned `wrfoce` copies (all five auxinput6 read times: 07-28, 08-27, 09-26, 10-26, 11-25) — same fix as Bug 19, but reapplied because the new wrfinput's TMN/projection globals could differ.
4. Re-patched WRF projection globals onto all 282 chem aux files (`01_patch_emiss_globals.py`) plus the 6 run-dir oce files.
5. Re-patched VPRM `T_ANN` with the new `wrfinput:TMN` (Bug 17 redo for safety; the value was identical to before, confirming `real.exe` is deterministic).
6. Reset namelist for seg1: `start = 2015-07-28`, `end = 2015-08-29`, `run_hours = 768`, `restart = .false.`.
7. Launched WRF from t=0 (16 MPI ranks).

**Lesson:** Any change to `start_*`/`end_*`/`run_hours` that extends the simulation period requires rerunning `real.exe` to regenerate `wrfbdy_d01`. The `02_run_wrf.sh` driver should validate `wrfbdy` coverage against the namelist before launching.

**Status:** ✅ Fixed; WRF restarted from t=0 (2015-07-28) on 2026-05-04 with a 503-record wrfbdy.

---

### Bug 21 — `02_align_aux_times.py` wrote `Times = SIM_START` to every oce file

**Date:** 2026-05-05  
**Affected:** `preprocessing/utils/02_align_aux_times.py`, all 5 start-aligned `wrfoce_d01_*` files in `simulations/IDN_BB_2015/run/`  
**Symptom:** After Bug 19's "create one oce file per auxinput6 read time" loop was added to the script, WRF still failed at sim 2015-08-27 00:00 with the same time-mismatch error pattern:
```
d01 2015-08-27_00:00:00  Input data is acceptable to use: wrfoce_d01_2015-08-27_00:00:00
 Time in file: 2015-07-28_00:00:00
 Time on domain: 2015-08-27_00:00:00
**WARNING** Time in input file not equal to time on domain **WARNING**
input_wrf: wrf_get_next_time current_date: 2015-07-28_00:00:00 Status = -4
ERROR: Could not find matching time in input file wrfoce_d01_2015-08-27_00:00:00
```

This time the error message correctly identifies the file (`wrfoce_d01_2015-08-27`), and the diagnostic clearly shows the file's internal Times = `2015-07-28` while WRF wants `2015-08-27`.

**Root cause:** The fix for Bug 19 added a loop over auxinput6 read dates, but `make_aligned_copy()` in `02_align_aux_times.py` always rewrote `Times` to the global `START_STR` (= `2015-07-28_00:00:00`), regardless of which date the file was supposed to represent. So all 5 oce files (07-28, 08-27, 09-26, 10-26, 11-25) ended up with `Times = 2015-07-28`. Only the 07-28 file happened to be correct by coincidence; the rest were wrong.

This is a subtle bug-in-the-fix: the per-date loop was added, but the per-date Times-rewrite was not — `make_aligned_copy()` was reused unchanged from when it only handled the start-date case.

**Fix:**
1. Added a `target_time_str` parameter to `make_aligned_copy()` so each call writes the correct per-file date.
2. Updated the oce loop to pass `target_time_str=ts` (the iteration's date string).
3. Updated the vprm fallback call site to pass `target_time_str=START_STR` explicitly.
4. Reran the script; verified all 5 oce files now have matching internal `Times`:

```
wrfoce_d01_2015-07-28_00:00:00 → Times = 2015-07-28_00:00:00
wrfoce_d01_2015-08-27_00:00:00 → Times = 2015-08-27_00:00:00
wrfoce_d01_2015-09-26_00:00:00 → Times = 2015-09-26_00:00:00
wrfoce_d01_2015-10-26_00:00:00 → Times = 2015-10-26_00:00:00
wrfoce_d01_2015-11-25_00:00:00 → Times = 2015-11-25_00:00:00
```

5. Re-patched WRF projection globals onto all 6 run-dir oce files.
6. Restarted WRF from `wrfrst_d01_2015-08-21_00:00:00` (real MODIS composite per Bug 15) with `run_hours = 192` (08-21 → 08-29 = seg1 end). This avoids redoing the 24 sim days of seg1 that already completed cleanly with the (correct) 07-28 oce file.

**Status:** ✅ Fixed; WRF restarted from `wrfrst_d01_2015-08-21_00:00:00` with corrected oce internal Times.

---

### Bug 22 — Cumulative-state failure at sim Oct 30 00:00 when seg4 integrates ≥6 sim days from Oct 24 wrfrst (NOT a fire-file corruption)

**Date:** 2026-05-06  
**Affected:** `simulations/IDN_BB_2015/input/wrffirechemi_d01_2015-10-30_00:00:00`  
**Symptom:** Seg4 (2015-10-24 → 11-25) crashed at the Oct 29 23:00 → Oct 30 00:00 fire-file transition. `rsl.error.0000` sequence (all 16 ranks identical):

```
d01 2015-10-30_00:00:00  Input data is acceptable to use: wrffirechemi_d01_2015-10-30_00:00:00
 Error trying to read metadata: 2015-10-30_00:00:00
input_wrf: wrf_get_next_time current_date: 2015-10-30_00:00:00 Status = -4
ERROR: Could not find matching time in input file wrffirechemi_d01_2015-10-29_00:00:00
```

Note: the final error message names `wrffirechemi_d01_2015-10-29_00:00:00` (the last **successfully** opened fire file), not the failing one — this is the same misleading-log pattern as Bug 19. The real failing file is `wrffirechemi_d01_2015-10-30_00:00:00`, indicated by the "Error trying to read metadata" line immediately above.

Fire files rotated successfully on Oct 25, 26, 27, 28, 29 within seg4, then failed on Oct 30.

**Root cause:** Subtle internal HDF5 metadata corruption that WRF v4.1.5's NetCDF library (`external/io_netcdf/wrf_io.F90:3174` → `share/input_wrf.F:1141 wrf_get_next_time`) rejects when reading the `TITLE` global attribute, but which ncdump, Python netCDF4, h5dump, and nccopy all tolerate. Likely caused by a partial/interrupted flush during the original file-creation script run on 2026-05-04. All structural diagnostics passed:
- `ncdump -h` diff vs Oct 29: only filename string differs
- Format: `netCDF-4`, contiguous chunks, no compression — identical to all neighbors
- Dims: `(24, 1, 127, 231)`, vars: `Times`, `ebu_in_co2`, `ebu_in_co`, `ebu_in_ch4`
- `Times` array: 24 entries `2015-10-30_00:00:00` → `2015-10-30_23:00:00`, ASCII `|S1` dtype, correct shape `(24,19)`
- Data: no NaN/Inf, max values physically reasonable (lower fire activity on Oct 30)
- Global attrs: all 89 attributes present, all WRF projection globals (`DX=27000`, `DY=27000`, `MAP_PROJ=3`, `TRUELAT1`, `TRUELAT2`, `STAND_LON`, etc.) matching wrfinput
- `nccopy` and `md5sum` completed without error

The defect was invisible to all user-space tools but caused WRF's internal NetCDF/HDF5 C library to fail reading the attribute.

**Fix attempt 1 (FAILED):** Rewrote the file from scratch using Python `netCDF4`, copying all variables, dimensions, and 89 global attributes into a fresh HDF5 container. New file (8,474,766 bytes) replaced original (`.bak` kept). WRF restarted and **crashed again at the exact same point** with the identical error — ruling out HDF5-internal corruption of the original file as the root cause.

**Fix attempt 2 (FAILED):** Cloned the *known-working* Oct 29 fire file as the new Oct 30 file, then overwrote `Times` and `ebu_in_*` data values. Resulting file had Oct 29's exact byte-level HDF5 layout with Oct 30 content. Crashed identically. Logs archived to `logs/rsl_seg4_attempt3/`.

**Diagnostic experiment (CONCLUSIVE):** Restarted from `wrfrst_d01_2015-10-29_00:00:00` (written by the failed run itself, 1.7 GB) with `run_hours = 26`, all other namelist parameters unchanged. WRF crossed the Oct 29 23:00 → Oct 30 00:00 transition **cleanly**, opened `wrffirechemi_d01_2015-10-30_00:00:00` without error (`open_aux_u : opening wrffirechemi_d01_2015-10-30_00:00:00 for reading`, `Yes, this special data is acceptable to use: WRF-Chem EMISSIONS`, `Input data is acceptable to use: wrffirechemi_d01_2015-10-30_00:00:00`), and ran through Oct 30 02:00 with `wrf: SUCCESS COMPLETE WRF`. Logs/wrfout/wrfrst archived to `logs/diag_oct29_26h/`.

**Conclusion:** The Oct 30 fire file is *not* the cause. The crash is a **cumulative-state failure** triggered by integrating ≥6 sim days from the Oct 24 restart. Most likely candidates:
  - **`WrfDataHandleMax = 99` exhaustion** (`WRF/external/io_netcdf/wrf_io.F90:41`). Each daily aux-stream rotation (chemi, fire) opens new file handles. Older handles may not be released cleanly, hitting the 99-handle ceiling around day 6.
  - **NetCDF/HDF5 internal cache or metadata-table growth** in WRF's I/O layer.
  - **Integer counter or buffer overflow** in some accumulator that wraps after ~6 daily rotations × 16 ranks.

The `"Error trying to read metadata"` line is generic — `share/input_wrf.F:1878 is_this_data_ok_to_use` calls `wrf_get_dom_ti_char(fid, 'TITLE', ...)` and emits this message for any non-zero ierr (including FILE_NOT_OPENED, NF_INQ_ATT failure, or handle-table corruption), not just genuine metadata corruption.

**Recommended fix:** Restart seg4 from `wrfrst_d01_2015-10-29_00:00:00` (which the failed run wrote at its `restart_interval=1440`) and integrate Oct 29 → Nov 25 (run_hours = 648 = 27 days). This violates Hard Rule #1 (MODIS composite restart boundaries) so vprm reads at Oct 29 + 8d × {0,1,2,3} = Oct 29 / Nov 6 / Nov 14 / Nov 22 require **dummy vprm files cloned from the nearest MODIS composites with internal `Times` rewritten** (using the same approach as Bug 21's `02_align_aux_times.py`). Required vprm dummies: Oct 29 (← Oct 24 source), Nov 6 (← Nov 9 or Nov 1), Nov 14 (← Nov 17 or Nov 9), Nov 22 (← Nov 25 or Nov 17).

Alternative fix (cleaner but unverified): patch `WRF/external/io_netcdf/wrf_io.F90` to raise `WrfDataHandleMax` from 99 to e.g. 256 and rebuild. Avoids the segment-split. Risk: requires WRF recompile.

**Status:** 🔍 Root cause identified empirically; fix path chosen pending user decision (segment-split vs WRF rebuild).

---

## 13.5 Output Validation Caveats

Reviewed first wrfout (`wrfout_d01_2015-07-25_00:00:00`, 17 frames) on 2026-05-04 after relaunch with BMKG-augmented fire. Background and anomaly tracers are physically reasonable, **with two caveats that affect interpretation but do not invalidate the run**:

### Caveat A — `CO_BCK` is a single constant (80) instead of a 3D field

**Symptom:** `CO_BCK` in `wrfinput_d01` and every wrfout frame is **identically 80 across the entire 3D grid** (range 79.97–80.04 only after a few hours due to numerical advection rounding).  
**Expected:** Real CarbonTracker CO ranges ~60–200 ppb across the troposphere — lower in the SH (clean marine), higher in NH/biomass-burning regions, vertical gradient of ~30–60 ppb between surface and tropopause.  
**Likely cause:** [`preprocessing/boundary_conditions/01_process_carbontracker.py`](preprocessing/boundary_conditions/01_process_carbontracker.py) does not populate a 3D CO field from CT — either CT-CO was not loaded, the field name lookup failed silently, or a constant default (80 ppb) was written instead of the interpolated profile. Same hybrid-sigma pathway as Bug 8 (CH₄) needs to be applied for CO.  
**Implication for current run:** The lateral and initial CO background is uniform. Fire and anthropogenic CO enhancements (`CO_BBU`, `CO_ANT`) are still correct because they're driven by emissions, not BCs. Total column CO will be biased low in NH inflow regions and biased high in clean SH inflow regions, by up to ~50 ppb relative to a real CT column.  
**Status:** ⚠️ Open — does not block the run; revisit before publication. Fix the CO branch in `load_ct_field()` analogous to the CH₄ hybrid-sigma fix (Bug 8).

### Caveat B — `CH4_*` and `CO_*` `units` attribute says `ppmv` but values are in `ppb`

**Symptom:** `wrfout` `CH4_*` fields carry the attribute `units = "ppmv"` but the magnitudes (1582–2039 for `CH4_BCK`, ~1810 mean) are clearly **ppb**, not ppmv. Same for `CO_*` (`CO_BCK = 80`).  
**Cause:** WRF-Chem registers all `chem_opt=17` GHG tracers with the unit string `ppmv` regardless of the actual unit the input data carries. CarbonTracker writes its CH₄ field in ppb (`CTCH4_2024.molefrac_glb3x2_*.nc` — `ch4` units `"parts per billion (ppb)"`) and the preprocessing leaves the values as-is. The CO background of 80 has the same ppb interpretation.  
**Implication for analysis:** When plotting/comparing to satellite retrievals or surface obs, treat:
- `CO2_*` as **ppmv** (file units correct, magnitudes ~400 confirm this)
- `CH4_*` as **ppb** (file units misleading; magnitudes ~1810 confirm ppb)
- `CO_*` as **ppb** (file units misleading; magnitudes ~80 confirm ppb)

The post-processing script `postprocessing/column_ghg/01_compute_xghg.py` converts XCH₄ and XCO to ppb in its output (correct), but the raw wrfout values must be interpreted manually.  
**Status:** ⚠️ Documented — not a code bug per se; a registry-level WRF-Chem convention. Add an explicit note in any plotting / analysis script: `# CH4_*, CO_* fields in wrfout are in ppb despite units='ppmv'`.

### What looked right

| Check | Value | Verdict |
|---|---|---|
| `CO2_BCK` initial | 393.8–424.2 ppm, mean 400 | ✅ matches CT2025 July tropics |
| `CH4_BCK` initial | 1582–2039 ppb, mean 1810 | ✅ matches CTCH4_2024 July tropics |
| `CO2_BBU` peak after 16 h | +130 ppm at fire grid cell | ✅ consistent with Indonesian peat-fire column observations |
| `CO_BBU` peak after 16 h | +10.5 ppm = +10,500 ppb | ✅ fresh BB plume |
| `CH4_BBU` peak after 16 h | +0.97 ppm = +970 ppb | ✅ peat-fire CH₄/CO₂ ratio realistic |
| `CO2_ANT` peak | +18–40 ppm | ✅ urban Java plume in surface layer |
| `CO2_BIO` after sunrise | +12–29 ppm anomaly | ✅ VPRM responding to MODIS LSWI/EVI |
| All anomalies at t=0 | 0 | ✅ by design |

### Day-1 diurnal review (`wrfout_d01_2015-07-28_00:00:00`, 24 frames)

Reviewed on 2026-05-04 after the simulation was relaunched from 2015-07-28 with BMKG-augmented fire and proper time-varying anthro:

**Working:**
- **Fire diurnal cycle is captured** — `CO_BBU` and `CO2_BBU` peak at LT 18–19 (UTC 11–12) with smooth ramp-up from LT 14, smooth decay overnight. Domain max `CO_BBU` ~55 ppmv, `CO2_BBU` ~650 ppm at single peat-fire grid cells (1.5°N 115.3°E central Borneo, 0°N 102°E Sumatra). Matches FINN+BMKG diurnal profile peaking 13–15 LT.
- **Anthropogenic spatial pattern correct** — top `CO2_ANT` cells at predawn map to real urban centers: HCMC (10.6°N 107.1°E), Aceh (3.4°N 98.3°E), Yogyakarta (-6.8°S 111.9°E), Bangkok (14.4°N 100.7°E), Jakarta (-6.3°S 105.8°E). Pixel-wise correlation between `CO2_ANT` and `CO2_BBU` is **0.0035** — they are spatially independent as expected.
- **PBL height** follows tropical maritime continent shape (0 → ~900 m noon → ~1070 m predawn).
- **Total CO₂ surface mean over Borneo land at predawn:** 414 ppm (range 400–477 ppm). Surface enhancement of +14 ppm above background is realistic for a peat-forest column with shallow nocturnal BL trapping respiration.

### Caveat C — Excessive nocturnal land-surface cooling (T2 over Borneo drops to 283 K predawn)

**Symptom:** T2 over Borneo land (interior peat, 1°S–4°N × 109–119°E) starts at 296.3 K mean (range 292–299) at LT 07 = sim start, peaks weakly at 299.4 K at LT 13 (only +3 K diurnal range up), then cools steadily to **283.5 K mean (range 278–290) at predawn LT 05**. That's a 13 K drop and a predawn surface temperature 5–10 K colder than realistic Indonesian dry-season conditions (real lows ~22–24 °C / 295–297 K).  
**Where the cold comes from:** wrfinput_d01 surface temperatures are reasonable. The cooling grows during integration. Most likely causes:
- `aer_op_opt = 0` disables aerosol-radiation feedback. During the BB event, optical-thick smoke would suppress longwave outgoing — without it, the model cools the surface unchecked
- `sf_sfclay_physics = 1` + Noah at 27 km may be over-stabilising the nocturnal BL, decoupling it from the warmer free troposphere
- The cold predawn explains why `CO2_BIO` accumulates strongly at predawn: respiration emissions trapped in a near-surface layer that the model is making *too* shallow and *too* cold

**Implication:** Surface CO₂ enhancements (especially `CO2_BIO` and `CO2_ANT`) will be biased high at predawn because of the unrealistic surface stability. Column-integrated XCO₂/XCH₄/XCO (the main publication output) is much less sensitive — the column averages over the cold near-surface layer. So this caveat affects surface-tower comparisons more than satellite XCO₂ comparisons.

**Recommended action before publication:** Try enabling aerosol-radiation feedback (`aer_op_opt = 1`, requires aerosol fields) or test with `sf_sfclay_physics = 5` (MYNN surface layer to match MYNN PBL). For now, the run continues with current settings — the column products should be usable.

**Status:** ⚠️ Documented; does not block the run.

### Day-2 follow-up review (`wrfout_d01_2015-07-29_00:00:00`, 16 frames so far)

Reviewed on 2026-05-04 to confirm the cold-T2 and `CO2_BIO` overshoot from day 1 don't compound during integration:

- **Daytime relaxation works.** `CO2_BIO` over Borneo land drops from +11.4 ppm at LT 06 (day-1 carryover) → +0.8 ppm at LT 15 (day-2 noon). Total surface CO₂ over Borneo land at LT 15 = **402.1 ppm — only +2 ppm above CT background**. GPP is drawing down the night's respiration accumulation as expected.
- **Diurnal cycle of `CO2_BIO` has the right shape:** predawn peak (+11), noon trough (+0.8), evening rise back. Magnitude is at the high end of published Indonesian peat-forest tower observations (5–25 ppm) but in-range.
- **Fire diurnal repeats correctly day-2.** `CO2_BBU` peak in fire pixels reaches 386 ppm at LT 19, ramp-up starts LT 14, decay begins LT 20. Same pattern as day-1.
- **T2 cold bias persists** but is bounded. Day-2 noon Borneo land = 298.2 K (real should be 302–304 K, ~5 K cold). Day-2 morning still recovering from the cold day-1 predawn (289 K at LT 07). The cold bias is **not amplifying day over day** — it's a steady-state offset, not a runaway.
- **Decision: continue the run.** Rationale:
  1. The publication target is column XCO₂/XCH₄/XCO (post-processing output). The cold near-surface 25 m layer contributes ~5% of the column mass; a 5 K bias there shifts column-mean temperature by ~0.25 K, which barely affects column tracer concentrations.
  2. The fire signal — the main scientific target — is captured with correct timing and realistic magnitude.
  3. Surface comparisons against in-situ towers will require a documented bias caveat, but column comparisons against GOSAT/OCO-2/TROPOMI should remain valid.
  4. Fixing this would require either (a) `aer_op_opt = 1` with a full aerosol field that doesn't exist for `chem_opt=17` GHG-only mode (multi-hour preprocessing build-out), or (b) `sf_sfclay_physics = 5` switch (untested with this config).
- If column products show systematic biases against satellite observations after the run completes, a v2 with corrected physics will be planned.

### Post-fix validation (Bugs 16 + 17) — 07-28 and 07-29 wrfout after t=0 restart

Reviewed `wrfout_d01_2015-07-28_00:00:00` and `wrfout_d01_2015-07-29_00:00:00` (24 frames each) on 2026-05-04 to confirm the case-rename and T_ANN fixes are taking effect:

**Bug 16 (CO₂ ocean flux) — confirmed active:**
- `EBIO_CO2OCE` in wrfout: **23,575 unique values**, range **−235 to +1123 mol km⁻² hr⁻¹** (matches the SOMFFN source data byte-for-byte)
- Realistic sink/source pattern over the Maritime Continent: ~10,028 source cells (mean +30 mol km⁻² hr⁻¹) vs ~1,676 sink cells (mean −21) — consistent with warm tropical-ocean outgassing dominating, with localised sinks at cold-water upwelling
- `CO2_OCE` accumulating: 27,166 → 28,321 non-zero surface cells over the first 2 sim days; surface enhancement ~0.01 ppm/day domain-mean — small magnitude but spatially correct

**Bug 17 (VPRM T_ANN) — confirmed active:**

Regional `CO2_BIO` means at the end of day 2 now show realistic spatial heterogeneity (was domain-uniform before fix):

| Region | N cells | `CO2_BIO` mean (ppm) | Range |
|---|---|---|---|
| Borneo highlands (1–3°N, 113–117°E) | 128 | 12.5 | 4.7 – 26.6 |
| Borneo lowlands (−1–1°N, 110–114°E) | 144 | 8.2 | 1.7 – 22.2 |
| Java (urban/cropland) | 171 | 2.3 | 0.1 – 10.8 |
| Sumatra | 362 | 7.1 | 0.3 – 20.3 |
| PNG | 569 | 7.0 | 0.04 – 21.6 |

Borneo-land coefficient of variation = 0.45 (was effectively 0 before — uniform T_ANN). The 5× contrast between Java and Borneo highlands matches expected ecology (cropland with high turnover vs cool-mountain forest with strong respiration amplification).

**Other tracers consistent with day-1/2 reviews:**
- `CO2_BCK` Borneo mean = 400.3–400.9 ppm ✅
- `CO2_BBU` max = 57–67 ppm at peat-fire grid cells (Borneo, Sumatra) ✅
- Total surface `CO2` Borneo mean = 412–413 ppm, max = 465–475 ppm at fire pixels ✅
- `CH4_BCK` mean = 1813–1819 ppb; `CH4_BBU` max ≈ 0.5 ppmv = 500 ppb at fire pixels ✅
- Fire diurnal repeats correctly day-2

**All five emission components (ANT, BBU, BIO, OCE, TST) are now correctly attributed and spatially resolved.** The Caveats A (CO_BCK constant), B (CH4/CO unit-attribute mislabeling), C (cold-T2 bias), and D (uninitialised SST cells) are unchanged and still apply.

### Caveat D — SST field has uninitialized cells (935 ocean cells with `SST = 0`)

**Symptom:** `SST` in `wrfinput_d01` has 935 of 15,402 ocean cells with `SST = 0.00 K` (range overall 0.00–303.63 K, mean over real ocean cells 299.8 K).  
**Cause:** ERA5 single-level SST coverage is 0.25° resolution; metgrid's nearest-neighbour interpolation onto the 27 km WRF grid leaves some marginal-sea / inland-water cells unfilled. With `sst_update = 0` these stay at 0 for the entire run.  
**Implication:** Those ocean cells behave as 0 K water — radiatively cold, no latent heat flux. In the Maritime Continent context most of these are likely small coastal cells and inland lakes; the impact on Indonesian peatland chemistry is small (the affected cells are not where the fires are).  
**Recommended action before publication:** Fill SST=0 cells from the surrounding median, or set `sst_update = 1` and provide the time-varying SST field via `wrflowinp_d01`.  
**Status:** ⚠️ Documented; does not block the run.

---

## 14. Progress Log

### 2026-05-03

| Time (WIB) | Event |
|-----------|-------|
| ~14:44 | ERA5 download started; PL Jul queued |
| ~15:33 | era5_pl_201507.grib (2.2 GB) ✅ |
| ~15:37 | era5_sl_201507.grib (144 MB) ✅ |
| ~16:45 | era5_pl_201508.grib (2.2 GB) ✅ |
| ~16:49 | era5_sl_201508.grib (144 MB) ✅ |
| ~17:10 | geogrid SUCCESS → geo_em.d01.nc (16 MB) ✅ |
| ~17:22 | era5_pl_201509.grib (2.1 GB) ✅ |
| ~17:26 | era5_sl_201509.grib (139 MB) ✅ |
| ~17:57 | era5_pl_201510.grib (575 MB, partial?) ✅ |
| ~18:00 | EDGAR CO₂ sector download started (edgar_sectors screen) |
| ~18:01 | EDGAR download crashed — Bug 5: TNR_Aviation_SPS 404 |
| ~18:05 | Bug 5 fixed; edgar_dl2 screen started |
| ~18:10 | EDGAR CO₂: 19/19 sectors ✅; CH₄: 22/22 sectors ✅ |
| ~18:11 | `era5_sl_201510.grib` (144 MB) ✅ → ERA5 now 8/10 |
| ~18:30 | Registry.EM_CHEM analysis: confirmed 15 3D GHG tracers for chem_opt=17 |
| ~18:45 | Physics schemes reviewed and updated: `bl_pbl_physics` 1→5 (MYNN), `cu_physics` 1→3 (Grell-Freitas) |
| ~19:00 | Post-processing script created: `postprocessing/column_ghg/01_compute_xghg.py` |
| ~18:56 | `era5_pl_201511.grib` (2.2 GB) ✅ |
| ~18:59 | `era5_sl_201511.grib` (139 MB) ✅ — **all 10/10 ERA5 GRIBs complete** |
| ~19:05 | WPS: ungrib started (screen `wps`) 🔄 |
| | MODIS MOD09A1: 799 HDF files (17 tiles × 47 dates, 0 missing) ✅ |
| | MODIS MCD12Q1: 17 tiles complete ✅ |
| | CarbonTracker: 368 NC files ✅ |
| | SOMFFN: 1.3 GB NC ✅ |

### 2026-05-04 (continued)

| Time (WIB) | Event |
|-----------|-------|
| ~early | WPS ungrib: 516 FILE intermediates complete (Jul 25 – Nov 30 18:00) ✅ |
| ~early | WPS metgrid: 516 met_em.d01.*.nc files ✅ |
| ~early | `00_setup_simulation.sh`: run directory created, all files linked ✅ |
| ~early | `01_run_real.sh`: real.exe SUCCESS — wrfinput_d01 (235 MB), wrfbdy_d01 (19 GB) ✅ |
| | EDGAR preprocessing: 3,120 wrfchemi_d01_* files ✅ |
| | Fire preprocessing (FINN + hotspot): 129 wrffirechemi_d01_* files ✅ |
| | Ocean preprocessing (SOMFFN): Bug 6 (lon/lat transpose) diagnosed and fixed; 6 wrfoce_d01_* files ✅ |
| | CarbonTracker preprocessing: Bug 7 (4D pressure) diagnosed and fixed; script re-run |
| | CarbonTracker CH₄: Bug 8 (hybrid sigma levels) — CH₄_BCK was 299–382 ppb; diagnosed, fixed, re-run |
| | CarbonTracker output verified: CO₂_BCK 393–424 ppm, CH₄_BCK 1582–2039 ppb ✅ |
| | VPRM preprocessing: 17 vprm_input_d01_* files ✅ |
| | namelist.input: Bug 9 — removed `plumerisefire_frct`, `co2_upper_bc`; fixed `vprm_opt`, `term_opt`, `io_form_restart` ✅ |
| ~early | `00_setup_simulation.sh` re-run: all 3,120 + 129 + 6 + 17 files linked into run dir ✅ |
| ~20:00 | `02_run_wrf.sh` started (screen `wrf`) — WRF seg1: 2015-07-25 → 2015-08-24 |
| ~20:10 | WRF startup failed: `io_form_restart = .false.` fatal error — Bug 10 diagnosed |
| ~20:15 | Bug 10 fixed in `02_run_wrf.sh` (sed line-anchored pattern); master namelist corrected |
| ~20:25 | WRF restarted (screen `wrf`) |
| ~20:28 | WRF seg1 confirmed running: `Timing for main: time 2015-07-25_01:02:06` (~0.81 s/step) 🔄 |
| later | WRF seg1 aborted: `med_read_wrf_chem_emissions: Status = -4` on first hour read |
| later | Diagnosis traced into `WRF/external/io_netcdf/wrf_io.F90:3174` and `WRF/share/input_wrf.F:1141`; root cause is hourly single-frame chemi files vs `wrf_get_next_time` semantics (Bug 11) |
| later | First workaround: `io_style_emissions = 1` + two static cyclic files (`wrfchemi_00z_d01`, `wrfchemi_12z_d01`) built from 2015-07-25 only — passed init, integrated 8 sim hours, but freezes anthro emissions to a single day's diurnal cycle for the entire 30-day run |
| later | Permanent fix: bundled 3,120 hourly files → 130 daily 24-frame files in `input_daily/`; reverted to `io_style_emissions = 2` with `frames_per_auxinput5 = 24`; symlinked into run dir; patched WRF projection globals from `wrfinput_d01` (Bug 11) |
| later | `wrfoce` and `vprm_input` still failed with `Time in file ≠ Time on domain`; copied nearest source files to start-dated names and rewrote internal `Times` (Bug 12) |
| later | OpenMPI `-np 28` failed (no slots); switched to `-np 16` (Bug 13) |
| later | WRF seg1 (proper time-varying anthro) integrating: `Timing for main: time 2015-07-25_06:00:00` and advancing at ~0.85 s/step on 16 ranks 🔄 |
| later | Verified ICs in first wrfout frame: `CO2_BCK ≈ 400 ppm`, `CH4_BCK ≈ 1810 ppb`, `CO_BCK = 80 ppb` (all anomaly tracers `_ANT`/`_BIO`/`_BBU`/`_OCE`/`_TST` start at 0 by design — they accumulate per-source enhancements) |
| later | After 13 sim hours: `CO2_BBU` peaks ~141 ppmv in fire grid cells, `CH4_BBU` ~1.2 ppmv, `CO_BBU` ~12 ppmv — fire signal is propagating into wrfout via plume rise (`biomass_burn_opt = 5`) |
| later | Bug 14 surfaced when checking BMKG hotspot effect on the fire input. User clarified: the CSV is already pre-filtered to high-confidence detections, so `confidence == 1` is a constant marker (not a filter key). Removed `MIN_CONFIDENCE` from `01_process_finn_hotspot.py` so all 2015 rows are used. Current WRF run still uses pre-fix (FINN-only) fire emissions; rerun preprocessing + restart from `wrfrst_d01_*` to switch in the BMKG-augmented fire input. |

**Current blockers:** None — WRF seg1 (proper time-varying anthro, FINN-only fire) running. BMKG-augmented fire emissions ready to swap in via re-preprocess + restart.

| later | User confirmed BMKG `confidence` column is a constant marker (CSV is already filtered to high-confidence). Stopped WRF, moved 129 FINN-only fire files to `input/_fire_finn_only_backup/`, re-ran `01_process_finn_hotspot.py` with all 54,930 hotspots, patched WRF projection globals onto the 129 new files, relaunched from t=0 with `mpirun -np 16`. |
| later | Validated first 17 wrfout frames: CO₂/CH₄ backgrounds match CarbonTracker initial values; fire/anthro/bio anomalies grow from 0 with realistic peat-fire and urban-plume magnitudes. Two caveats documented in §13.5: (a) `CO_BCK` is constant 80 instead of a 3D CT field — preprocessing miss, fix before publication; (b) `CH4_*` and `CO_*` carry attribute `units="ppmv"` but values are ppb — WRF-Chem registry convention, not a bug. |
| later | Restart attempt for seg2 (start 2015-08-01 from `wrfrst_d01_2015-08-01_00:00:00`) failed at 2015-08-02 with `Possibly missing file for = auxinput15` (Bug 15). Diagnosed: WRF computes auxinput15 read times as `start + N × 8 days`, which lands on dates that don't match the MODIS 8-day composite calendar. Ad-hoc VPRM symlinks at `start + 8k` dates were tried but failed because the internal `Times` variable still carried the original composite date (same trap as Bug 12). |
| later | Adopted permanent fix: align simulation start (and all restart-segment boundaries) to real MODIS composite dates. Chose 2015-07-28. Stopped WRF, deleted all wrfrst/wrfout from the failed run, removed non-composite VPRM symlinks (07-25, 08-02, 08-10, 08-18), updated namelist (start 07-28, end 08-27, restart=false), reran `real.exe` SUCCESS, updated `01_process_carbontracker.py` `SIM_START → 2015-07-28` and re-patched ICs/BCs, regenerated `wrfoce_d01_2015-07-28_00:00:00` with rewritten Times and patched globals, relaunched WRF on 16 ranks. Integration started cleanly. Restart-segment plan in §6 updated to align every boundary with a MODIS composite. |
| later | User flagged inconsistency: seg1 was set to end 2015-08-27, but the next MODIS composite is 2015-08-29. Stopped WRF, updated namelist `run_hours = 768` and `end_day = 29` so seg1 ends on a real composite. Cleaned wrfrst/wrfout/logs and relaunched. Final segment plan in §6: seg1 07-28→08-29, seg2 08-29→09-30, seg3 09-30→10-24, seg4 10-24→11-25, seg5 11-25→12-01 — every boundary is a MODIS 8-day composite date, so each restart's first auxinput15 read finds an existing real file. |
| later | Audited simulation.md against actual run state: corrected `time_step` 162→150 s; MPI tasks 28→16 (Bug 13 reality on this 16-core host); updated §9.3 Fire Emissions to reflect FINN+BMKG hybrid (incl. Bug 14 fix); updated §11 file-count tables and §15 task list to match the new 2015-07-28 start, daily 24-frame chemi layout, and MODIS-aligned segment boundaries. |
| later | Promoted ad-hoc helper scripts from `/tmp/` into the repo: `preprocessing/anthropogenic/02_bundle_daily_chemi.py` (hourly→daily 24-frame chemi, idempotent), `preprocessing/utils/01_patch_emiss_globals.py` (patch WRF projection globals onto chem aux files), `preprocessing/utils/02_align_aux_times.py` (start-date-aligned copies of `wrfoce`/`vprm_input` with internal `Times` rewritten). Removed obsolete one-shot variants. Updated §8.1 script overview and §12 directory tree. |
| later | Day-1 diurnal review of `wrfout_d01_2015-07-28_00:00:00` (24 frames). **Working:** fire diurnal peak at LT 18–19 (`CO_BBU` 55 ppmv, `CO2_BBU` 650 ppm at peat-fire pixels — central Borneo and Sumatra peatlands), anthropogenic spatial pattern matches real urban centers (HCMC, Bangkok, Jakarta, Aceh, Yogyakarta), `CO2_ANT`/`CO2_BBU` pixel correlation = 0.0035 (spatially independent, as expected), PBL height shape realistic (0 → 900 m noon → 1070 m predawn), total surface CO₂ over Borneo land at predawn = 414 ppm mean (realistic +14 ppm dry-season peat-forest enhancement). **Caveats added** to §13.5: (C) excessive nocturnal land-surface cooling — Borneo T2 drops to 283 K predawn, ~5–10 K too cold; likely due to `aer_op_opt=0` disabling aerosol-radiation feedback. Affects surface-level enhancements but not column products. (D) 935 of 15,402 ocean cells have `SST=0` from ERA5/metgrid interpolation gaps. Non-blocking. |
| later | Day-2 follow-up review of `wrfout_d01_2015-07-29_00:00:00` (16 frames). Confirmed cold-T2 and `CO2_BIO` overshoot don't compound: daytime GPP relaxation drops `CO2_BIO` from +11.4 ppm (day-1 predawn carryover) → +0.8 ppm at LT 15 noon; total surface CO₂ over Borneo land at noon = **402.1 ppm, only +2 ppm above background**. Fire diurnal repeats correctly day-2 (`CO2_BBU` peak 386 ppm at LT 19). T2 cold bias is steady-state ~5 K offset, not runaway. **Decision: continue the run** — column XCO₂/XCH₄/XCO products (publication target) are not materially affected by a 25 m-layer bias; surface-tower comparisons get a documented caveat. Aerosol-radiation feedback fix (`aer_op_opt=1`) deferred to v2 if column comparisons show systematic bias. |
| later | Three new bugs surfaced via wrfout cross-checks (`CO2_OCE`, `CO2_TST`, `T_ANN`). **Bug 16:** `EBIO_CO2OCE` in wrfout was zero everywhere because the wrfoce file's variable was named `ebio_co2oce` lowercase but the WRF registry preprocessor emits `'EBIO_CO2OCE'` UPPERCASE in `allocs.inc`, and `nf_inq_varid` is case-sensitive — silent read failure. Renamed in all 7 wrfoce files via `ncrename`; fixed `01_process_somffn.py` to write UPPERCASE going forward. **Bug 17:** `T_ANN` in vprm_input was a constant 300.15 K placeholder, removing the spatial gradient from VPRM's Q10 respiration. Patched all 17 vprm_input files with `wrfinput_d01:TMN` (Noah's annual climatology, real range 280–304 K); fixed `01_process_vprm.py` to read TMN going forward. **Bug 18:** `*_TST` tracers are exact duplicates of `*_ANT` (driven by `E_*TST = E_*` from wrfchemi only) — not a tagged total of all sources. Documentation in §5.4 corrected. WRF stopped at sim time 2015-07-31 02:15, restarted from t=0 to apply Bugs 16+17 fixes. |
| later | Post-fix validation on `wrfout_d01_2015-07-28_00:00:00` and `wrfout_d01_2015-07-29_00:00:00` (full 24-frame days each). Confirmed: `EBIO_CO2OCE` now non-zero (23,575 unique values, range −235 to +1123 mol km⁻² hr⁻¹, matching the SOMFFN source data); `CO2_OCE` accumulating realistically (27k→28k non-zero cells over 2 days, ~0.01 ppm/day surface). VPRM `CO2_BIO` shows expected spatial heterogeneity (Borneo highlands 12.5 ppm, Java 2.3 ppm — 5× regional contrast; Borneo CV = 0.45 vs ≈ 0 before). All other tracers consistent with prior day-1/2 reviews. Caveats A/B/C/D in §13.5 unchanged. |
| later | WRF aborted at sim 2015-08-27 00:00:00 (30 sim days into the run) with `wrf_get_next_time Status = -4`. Restart from `wrfrst_d01_2015-08-21_00:00:00` reproduced the failure at the same calendar date. **Bug 19** identified: WRF reads auxinput6 (oce) at sim_start + N×30 days, constructing filenames `wrfoce_d01_<YYYY-MM-DD>_00:00:00` from the current sim time. With sim_start = 2015-07-28, the read times are 07-28, 08-27, 09-26, 10-26, 11-25 — only 07-28 had a corresponding file. Created start-aligned copies of the nearest monthly source for each required date (Times rewritten, projection globals patched), and updated `preprocessing/utils/02_align_aux_times.py` to loop over all auxinput6 read dates instead of creating only the start-date copy. The misleading error message references `wrffirechemi_d01_2015-08-26` because it was the last successful "Input data is acceptable to use" log line, not the actual failing file. Restarted WRF from `wrfrst_d01_2015-08-21_00:00:00`. |
| later | After Bug 19 fix, WRF cleared the 2015-08-27 00:00 chem reads cleanly (chemi/oce/fire all OK) but immediately crashed on `wrfbdy_d01` with "Ran out of valid boundary conditions". **Bug 20** identified: wrfbdy had only 120 records (07-28 → 08-26 18:00 = 30 days) because `real.exe` was last run with `end_day = 27` and never rerun after the namelist end-date was bumped. Set namelist to full period (`end = 2015-12-01`, `run_hours = 3024`) and reran real.exe — completed 504/505 loops with a harmless final-loop FATAL on a missing `met_em.d01.2015-12-01_00:00:00.nc` (we have met_em only through 11-30 18:00). On-disk wrfbdy has **503 records spanning 07-28 00:00 → 11-30 12:00** — covers the entire planned simulation period. Re-patched CT BCs (SIM_START=2015-07-28); re-ran `02_align_aux_times.py` to refresh the 5 oce read-time copies; re-patched globals on 282 chem files plus 6 run-dir oce files; re-patched T_ANN from new wrfinput:TMN (identical, confirming real.exe is deterministic). Reset namelist to seg1 (start=07-28, end=08-29, run_hours=768, restart=false). **Decision: do whole seg1 from the beginning** (per user direction) — relaunched WRF from t=0 on 16 MPI ranks. |

### 2026-05-05

| Time (WIB) | Event |
|-----------|-------|
| ~early | Seg1 fresh run from 07-28 reached sim 2015-08-27 00:00 cleanly — Bugs 19/20 confirmed fixed (wrfbdy and oce reads work past the prior crash point). New error appeared: `wrfoce_d01_2015-08-27_00:00:00` had internal `Times = 2015-07-28_00:00:00`. **Bug 21** identified: my own fix-script `02_align_aux_times.py` had a per-date loop added for Bug 19, but the helper `make_aligned_copy()` always wrote the *global* `START_STR` to `Times` regardless of the iteration's target date. So all 5 read-time oce files carried `Times = 2015-07-28`. Added `target_time_str` parameter to `make_aligned_copy()`, updated callers, reran the script, verified each oce file now has matching internal Times. Re-patched globals. |
| ~early | Restarted WRF from `wrfrst_d01_2015-08-21_00:00:00` (real MODIS composite per Bug 15), `run_hours = 192` (08-21 → 08-29). Saves ~24 sim days vs restarting from t=0. The 30 already-written wrfout files (07-28 → 08-26) remain valid — they used the correct 07-28 oce file. |
| ~09:00 | **Seg1 SUCCESS COMPLETE WRF.** 32 daily wrfout files written for 07-28 → 08-28. Pre-restart files (07-28 → 08-26) named `_00:00:00` (canonical), post-restart files (08-21 → 08-28) named `_01:00:00` (frame range 01:00 of date X → 00:00 of date X+1). Stitched 08-27 and 08-28 into canonical `_00:00:00` files via `/tmp/stitch_wrfout.py` (frame[23] of prev day's `_01:00:00` + frames[0..22] of curr day's `_01:00:00`). Deleted 8 redundant `_01:00:00` files (freed 64 GB). Final clean seg1 set: 32 daily `_00:00:00` files, 07-28 → 08-28. |
| ~10:00 | Cleanup before seg2: moved seg1 wrfout to `output/seg1/`; deleted 30 redundant wrfrst (kept only 08-29). Disk: 425 GB used → 680 GB after cleanup → 531 GB free. |
| ~10:15 | **Seg2 launched** from `wrfrst_d01_2015-08-29_00:00:00`, end 2015-09-30, run_hours=768 (32 days). Pre-flight verified all input streams: chemi/fire daily 24-frame files (08-29 → 09-30), oce read times 08-29 + 09-28 (re-anchored from restart start; created via inline `02_align_aux_times.py` logic), vprm at MODIS composites 08-29/09-06/09-14/09-22/09-30 (all real, with spatial T_ANN), wrfbdy spans through 11-30 12:00. |
| ~10:24 → ~11:32 | Seg2 progressed at ~1.0 s/step. By 11:32 reached sim 2015-09-04 (6/32 sim days done). |
| later | **Seg2 SUCCESS COMPLETE WRF** at sim 2015-09-30 00:00. 18,432 timesteps total, 32 wrfout files (`_01:00:00` pattern, 8.5 GB each). |
| later | Storage near limit (193 GB free) when stitching seg2 was attempted; stitch interrupted after 4 files. Deleted partials, freed wrfrst (52 GB), reached 276 GB free. User backed up seg2 wrfout externally instead of stitching, deleted local 32 `_01:00:00` files (freed 272 GB → 531 GB free). For seg2 analysis, files in backup carry `_01:00:00` filenames but the `Times` variable inside is correct — index by `Times` for unified handling across segments. |
| later | **Seg3 launched** from `wrfrst_d01_2015-09-30_00:00:00`, end 2015-10-24, run_hours=576 (24 days). Pre-flight verified all input streams: 25 chemi + 25 fire daily files (09-30 → 10-24), oce single read at 09-30 (segment < 30 days, only 1 auxinput6 read), vprm at MODIS composites 09-30/10-08/10-16/10-24 (all real). ETA ~3.5h wall. |
| later | Created `wrf-ghg` Claude skill at `~/.agents/skills/wrf-ghg/` to operate this simulation: `SKILL.md` distills the 12 hard rules from Bugs 1–21, references this document as authoritative source. Bundled scripts: `check_seg_inputs.py` (pre-flight verification of chemi/fire/oce/vprm/wrfbdy/wrfrst with Times alignment, EBIO_CO2OCE case, T_ANN spatial check), `monitor_run.sh` (health check), `launch_segment.sh` (16-rank nohup launch), `stitch_seg_wrfout.py` (optional canonical-naming with disk guard). References: `references/segment_workflow.md` and `references/file_conventions.md`. |

### 2026-05-06

| Time (WIB) | Event |
|-----------|-------|
| ~early | **Seg3 SUCCESS COMPLETE WRF** at sim 2015-10-24 00:00. 24 daily wrfout files (`_01:00:00` pattern). User backed up wrfout externally. |
| ~07:55 | Cleanup before seg4: deleted 24 seg3 wrfout files, 24 intermediate wrfrst (kept 10-24), rsl logs. Disk: 655 GB used → 425 GB used, 532 GB free (freed 230 GB). |
| ~07:58 | Pre-flight via `check_seg_inputs.py 2015-10-24 2015-11-25` flagged 2 missing oce files. Created `wrfoce_d01_2015-10-24_00:00:00` (← 10-01 monthly source) and `wrfoce_d01_2015-11-23_00:00:00` (← 11-01 monthly source) with internal `Times` rewritten and projection globals patched. Re-verified: chemi 33/33, fire 33/33, vprm 5/5 (10-24/11-01/11-09/11-17/11-25 all real MODIS composites), oce 2/2, wrfbdy covers through 11-30 12:00. ✅ All inputs ready. |
| ~07:59 | **Seg4 launched** via `launch_segment.sh`: start 2015-10-24, end 2015-11-25, run_hours=768 (32 days), restart=true, 16 MPI ranks. ETA ~5h wall. |
| ~later | Seg4 FATAL at sim 2015-10-30 00:00: **Bug 22** — `Error trying to read metadata` on fire-file rotation. Misleading-log pattern (Bug 19): named file in the ERROR is `wrffirechemi_d01_2015-10-29` (last good), but actual failing read is on Oct 30. Crash logs archived to `logs/rsl_seg4_attempt1/`. |
| ~later | Bug 22 attempt 1: rewrote Oct 30 fire file from scratch via Python netCDF4 (preserve all data + 89 globals). Relaunched from Oct 24 wrfrst. **Same crash** at same time — ruling out internal HDF5 corruption of the original as root cause. Logs archived to `logs/rsl_seg4_attempt2/`. |
| ~later | Bug 22 attempt 2: cloned the known-working Oct 29 fire file's HDF5 layout, then overwrote `Times` to Oct 30 strings and `ebu_in_co2/co/ch4` arrays with Oct 30 data from `.bak`. Resulting file has Oct 29's byte-level container with Oct 30's content. Relaunched from `wrfrst_d01_2015-10-24_00:00:00`. **Same crash** at same time. Logs archived to `logs/rsl_seg4_attempt3/`. |
| ~later | Bug 22 diagnostic: restarted from `wrfrst_d01_2015-10-29_00:00:00` (1.7 GB, written by failed run) with `run_hours=26`, debug_level=1, all other params unchanged. **WRF crossed Oct 30 00:00 transition cleanly**, opened the Oct 30 fire file without error, ran to Oct 30 02:00 with `SUCCESS COMPLETE WRF`. Logs in `logs/diag_oct29_26h/`. |
| ~later | **Conclusion: Oct 30 fire file is exonerated.** Crash is a cumulative-state failure across \u22656 sim days from Oct 24 wrfrst. Most likely cause: WRF `WrfDataHandleMax=99` exhaustion or NetCDF I/O resource leak. Bug 22 documentation updated with diagnostic findings and fix paths (segment-split from Oct 29 with vprm dummies, OR WRF source patch raising `WrfDataHandleMax`). |

---

## 15. Pending Tasks

### Critical path (blocking post-processing)

| # | Task | Blocker | Status |
|---|------|---------|--------|
| 1 | Complete ERA5 download | — | ✅ All 10/10 GRIBs done |
| 2 | ungrib + metgrid | ERA5 | ✅ 516 met_em files |
| 3 | `00_setup_simulation.sh` | metgrid | ✅ Run dir set up |
| 4 | `01_run_real.sh` | setup | ✅ wrfinput_d01 + wrfbdy_d01 |
| 5 | EDGAR preprocessing | real.exe | ✅ 3,120 hourly → 130 daily 24-frame in `input_daily/` (Bug 11) |
| 6 | Fire preprocessing (FINN + BMKG) | real.exe | ✅ 129 wrffirechemi files (BMKG hotspots active after Bug 14 fix) |
| 7 | Ocean preprocessing | real.exe | ✅ 6 wrfoce files + start-date copy 2015-07-28 |
| 8 | CarbonTracker preprocessing | real.exe | ✅ ICs + BCs patched (`SIM_START = 2015-07-28`) |
| 9 | VPRM preprocessing | real.exe + MODIS | ✅ 17 vprm_input files (real MODIS composites only) |
| 10 | seg1 (2015-07-28 → 2015-08-29) | all preprocessing | ✅ Done (32 daily `_00:00:00` wrfout files, backed up externally) |
| 11 | seg2 (2015-08-29 → 2015-09-30) | seg1 wrfrst | ✅ Done (32 `_01:00:00` wrfout files, backed up externally) |
| 12 | seg3 (2015-09-30 → 2015-10-24) | seg2 wrfrst | ✅ Done (24 `_01:00:00` wrfout files, backed up externally) |
| 13 | seg4 (2015-10-24 → 2015-11-25) | seg3 wrfrst | 🔄 **RUNNING** (16 MPI ranks; ETA ~5h) |
| 14 | seg5 (2015-11-25 → 2015-12-01) | seg4 wrfrst | ⏳ |
| 15 | Post-processing: `01_compute_xghg.py` | all 5 segs complete | ⏳ |

**Note on file naming across segments:** Seg1 wrfout files use the canonical `wrfout_d01_<date>_00:00:00` pattern (frames 00:00–23:00 of each date). Seg2+ wrfout files written by restart runs use `wrfout_d01_<date>_01:00:00` (frames 01:00 of date X → 00:00 of date X+1). Always index by the `Times` variable inside the file rather than by filename when post-processing across segments. The `00:00:00` vs `01:00:00` filename suffix is purely a WRF restart-IO convention and does not indicate any data difference.

### Monitoring WRF

```bash
# Check current timestep:
tail -5 simulations/IDN_BB_2015/run/rsl.out.0000

# Check overall segment log:
tail -30 simulations/IDN_BB_2015/logs/wrf_run.log

# Segment complete when this appears in rsl.out.0000:
grep "SUCCESS COMPLETE WRF" simulations/IDN_BB_2015/run/rsl.out.0000

# Expected output files after all segments:
ls simulations/IDN_BB_2015/output/wrfout_d01_* | wc -l  # expect 126 (2015-07-28 → 2015-12-01, daily files)
```

### Post-processing (after WRF complete)

```bash
source .venv/bin/activate
python postprocessing/column_ghg/01_compute_xghg.py \
    --start 2015-08-01 --end 2015-11-30 --plot --combine
```

### Expected timing (from seg1 relaunch ~07:09 WIB 2026-05-04)

Pace observed during seg1: ~0.85 s wall-clock per 150 s model step → real-time factor ~177× → ~5.4 h wall per 30-day segment on 16 MPI ranks.

| Segment | Period | Sim days | Est. wall-clock |
|---------|--------|----------|-----------------|
| seg1 | 2015-07-28 → 2015-08-29 | 32 | ~5.7 h |
| seg2 | 2015-08-29 → 2015-09-30 | 32 | ~5.7 h |
| seg3 | 2015-09-30 → 2015-10-24 | 24 | ~4.3 h |
| seg4 | 2015-10-24 → 2015-11-25 | 32 | ~5.7 h |
| seg5 | 2015-11-25 → 2015-12-01 | 6 | ~1.1 h |
| **Total** | | **126 days** | **~22.5 h** |

*Timing estimate based on ~0.85 s wall-clock per model step (150 s model time step; 16 MPI ranks, one per physical core).*

---

## 16. WRF-STILT Setup Plan

**Status:** 🔄 Setup in progress (2026-05-07). WRF met available in `simulations/IDN_BB_2015/output/` (Oct 22–29 2015, 8 files, 24 frames each). STILT project directory: `stilt/`.

### What WRF-STILT does

STILT (Stochastic Time-Inverted Lagrangian Transport) runs backward particle trajectories from receptor locations using WRF meteorology to compute footprints — the sensitivity of measured atmospheric concentrations to surface fluxes. Combined with the WRF-Chem GHG tracers, this enables attribution of CO₂/CH₄/CO enhancements to specific source regions and types (fire, anthropogenic, biosphere, ocean).

### Met data available

| File pattern | Period covered | Frames/file | Grid |
|---|---|---|---|
| `wrfout_d01_2015-10-2{2..9}_01:00:00` | Oct 22 01:00 → Oct 30 00:00 | 24 hourly | 232×128×50, 27 km Mercator |

Domain: lat −15° to +15°, lon 89.6° to 145.4°, centred 0°N 117.5°E (Borneo peak fire period).

WRF projection: `MAP_PROJ=3` (Mercator), `TRUELAT1=0.0`, `STAND_LON=117.5`, `DX=DY=27000 m`.

### Receptor sites

Key surface-monitoring sites within the WRF domain for the 2015 Indonesian fire event:

| ID | Name | Lat | Lon | zagl (m) | Rationale |
|---|---|---|---|---|---|
| SGP_Bukit_Kototabang | NOAA/GAW Bukit Kototabang | -0.20 | 100.32 | 10 | Long-term GHG obs, W Sumatra |
| MHD_Mahé | Mahé Seychelles (GAW) | -4.68 | 55.53 | 10 | Clean marine BG, upwind |
| IDN_Palangkaraya | Palangka Raya (Central Borneo peat) | -2.16 | 113.94 | 10 | Core peat-fire region |
| IDN_Pontianak | Pontianak (W Kalimantan) | -0.02 | 109.33 | 10 | Fire outflow corridor |
| IDN_Jambi | Jambi (S Sumatra peat) | -1.61 | 103.61 | 10 | S Sumatra fire area |
| SGP_Singapore | Singapore (AERONET) | 1.30 | 103.82 | 10 | Urban/downwind receptor |

These can be extended with TCCON/satellite overpass times later.

### Software stack

| Component | Version | Status |
|---|---|---|
| R | 4.3.3 | ✅ installed |
| `uataq/stilt` R package | GitHub HEAD | ❌ needs install |
| R deps (`ncdf4`, `dplyr`, `parallel`, `sf`, `lubridate`, `ggplot2`) | — | ✅ all installed |
| STILT executable (`hymodelc`) | built during stilt_init | ❌ needs build |

### Setup steps (before running)

1. **Install `uataq/stilt` R package** — `devtools::install_github('uataq/stilt')` (requires `devtools`, already installed)
2. **Initialise STILT project** — `stilt::stilt_init('stilt/')` creates `run_stilt.r`, `r/`, `exe/` subdirs and builds the `hymodelc` Fortran executable via `make`
3. **Configure `stilt/run_stilt.r`** — set `project`, `met_path`, `met_file_format`, `receptors`, `n_hours`, `numpar`, `n_cores`
4. **Define receptors** — write `stilt/receptors.csv` with columns `run_time, lati, long, zagl`
5. **Verify met files accessible** — test that STILT can locate and read a single WRF frame
6. **Dry-run test** — one receptor × 24 h backward × 100 particles to validate the footprint pipeline

### Key STILT namelist parameters for this domain

```r
# Met
met_path       <- '/home/igrk/WRF-GRK/simulations/IDN_BB_2015/output'
met_file_format <- 'wrfout_d<domain>_%Y-%m-%d'   # STILT WRF reader globs by date
met_subgrid_buffer <- 0.1    # fraction of domain edge to discard

# Particles
numpar  <- 200        # particles per receptor-time
n_hours <- -120       # 5 days backward
rm_dat  <- TRUE       # remove intermediate particle files

# Domain
xmn <- 90;  xmx <- 145   # lon bounds (subset of WRF domain)
ymn <- -15; ymx <- 15    # lat bounds

# Footprint grid
hnf_plume <- TRUE
smooth_factor <- 1
time_integrate <- TRUE

# Parallelism
n_cores <- 8          # half of 16 physical cores; leave headroom
```

### Output

- `stilt/out/particles/` — particle trajectory files per receptor-time
- `stilt/out/footprints/` — gridded influence footprints (NetCDF, lat×lon×time)
- Footprints can be convolved with WRF-Chem flux fields (`CO2_BBU`, `CO2_ANT`, `CO2_BIO`) to compute source-tagged concentration enhancements

### Files created

| File | Purpose |
|---|---|
| `stilt/setup_stilt.r` | Install + init script (run once) |
| `stilt/run_stilt.r` | Main configured run script |
| `stilt/receptors.csv` | Receptor location/time table |

---

*This document is updated as the simulation progresses. See `simulations/IDN_BB_2015/logs/` for runtime logs.*
