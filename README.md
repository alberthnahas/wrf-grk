# WRF-GRK — WRF-GHG → STILT → Inversion for the 2015 Indonesia Biomass-Burning Event

Top-down greenhouse-gas modelling pipeline for the **October 2015 El Niño Indonesian fires**,
built at **BMKG** (Indonesian Agency for Meteorology, Climatology and Geophysics). It couples a
WRF-GHG forward simulation, WRF-STILT backward-trajectory footprints, and a Bayesian inversion to
estimate regional fire CO₂ emissions constrained by satellite (OCO-2), reanalysis (CarbonTracker
CT2025), and in-situ (Bukit Kototabang flask) observations.

The 2015 fires were the largest since 1997: an extreme drought (SPI < −2.5) drove intense peat
combustion across Sumatra and Kalimantan, releasing on the order of **0.8–1.0 Pg C**, with peak
daily emissions rivalling the entire EU economy. This repository holds the code, configuration,
and documentation for that study — **not** the large input/output data (see [Data](#data)).

---

## Pipeline

```
   ERA5 · EDGAR · FINN · CarbonTracker · SOMFFN · VPRM        (preprocessing/)
                              │
                              ▼
      WRF-GHG  (WRF-Chem 4.1.5, chem_opt=17)                  (simulations/)
      27 km · 232×128×50 · 2015-07-28 → 2015-12-01
      GHG tracers: CO2/CO/CH4 × {ANT,BIO,OCE,BCK,BBU,TST}
                              │
             ┌────────────────┴─────────────────┐
     meteorology → ARL                   GHG tracer fields
             ▼                                   │  (validation / dormant CO check)
   WRF-STILT footprints (−120 h)                 │
   11 sites · 06 UTC · 0.25° grid       (stilt/) │
             ▼                                   ▼
   Bayesian inversion  ────────────────────────────────────  (inversion/)
   y = CT2025 + BKT flask + OCO-2   →   6 regional CO₂ scaling factors α (posterior)
```

Full scientific rationale, bug history, and run conventions live in [`docs/`](docs/).

---

## Repository layout

Only code, configuration, and documentation are tracked. Everything large or machine-specific is
excluded via [`.gitignore`](.gitignore).

| Path | Contents |
|---|---|
| `preprocessing/` | Input builders: `anthropogenic/` (EDGAR), `fire/` (FINN + hotspots), `ocean/` (SOMFFN), `boundary_conditions/` (CarbonTracker), `vprm/`, `domain/` (WPS) |
| `postprocessing/` | `column_ghg/` — XGHG column computation |
| `scripts/` | Build, download, run, and plotting helpers |
| `stilt/` | WRF→ARL conversion, receptor lists, STILT run scripts (`run_stilt_v3p2_06utc.r`), footprint/quantify plotting; STILT-R framework under `r/` |
| `inversion/` | `run_inversion_v3.py` (Bayesian solver) + OCO-2 downloader and plotting scripts |
| `geojson/` | Indonesia administrative, peat-land, and masking boundaries |
| `docs/` | Scientific and operational documentation (see [Documentation](#documentation)) |
| `AGENTS.md` | Entry notes for AI-assisted work on this repo |

**Not tracked** (regenerate or fetch): `rawdata/`, `simulations/`, `stilt/{met,out,out_06utc}`,
`inversion/data/`, `WRF/`, `WPS/`, `libs/`, `.venv/`.

---

## Data

Input and output data are **not** in the repository (hundreds of GB; the full `wrfout` archive
alone is ~1.1 TB). Sources:

| Dataset | Use | Source |
|---|---|---|
| ERA5 | WRF meteorology / boundary conditions | Copernicus CDS |
| EDGAR | Anthropogenic GHG emissions | EC-JRC |
| FINN v2.5 | Biomass-burning emissions (CO₂/CO/CH₄) | NCAR |
| CarbonTracker CT2025 | Background / pseudo-observations | NOAA GML |
| OCO-2 Lite FP v11.1r | Column XCO₂ constraint | NASA GES DISC (Earthdata login) |
| Bukit Kototabang flask | In-situ CO₂ constraint | NOAA GML |
| MODIS hotspots / MCD | Fire detections, land cover | NASA LAADS/Earthdata |
| SOMFFN | Ocean CO₂ flux | Landschützer et al. |

> Credentials (e.g. NASA Earthdata) are read from `~/.netrc` and are **never** stored in this repo.

---

## Environment

- **WRF-Chem 4.1.5** and **WPS 4.1** (built separately; not vendored here).
- **HYSPLIT/STILT** with the `arw2arl` converter (`stilt/exe/`, binaries not tracked).
- **Python** via the project virtualenv (`.venv/`, not tracked) — `netCDF4`, `numpy`, `scipy`,
  `pandas`, `matplotlib`.
- **R** with the [`uataq/stilt`](https://github.com/uataq/stilt) framework (bundled under `stilt/r/`).

---

## Workflow

1. **Preprocess** inputs → `preprocessing/` (EDGAR, FINN, CarbonTracker, SOMFFN, VPRM, WPS domain).
2. **Run WRF-GHG** in restart segments → hourly `wrfout` (one daily file at 01:00 UTC). See the
   `wrf-ghg` skill and [`docs/simulation.md`](docs/simulation.md).
3. **Convert `wrfout` → ARL** meteorology: `bash stilt/convert_wrf2arl.sh`.
4. **Run STILT** backward footprints (06 UTC receptors): `Rscript stilt/run_stilt_v3p2_06utc.r`.
5. **Run the inversion**: `python3 inversion/run_inversion_v3.py` → posterior α + Tg CO₂ per region
   and diagnostic plots.

Step-by-step operator instructions — including downloads, checks, and the case-study window — are
in [`docs/stilt_inversion_agent_guide.md`](docs/stilt_inversion_agent_guide.md).

---

## Key configuration decisions

- **WRF-GHG tracers** (`chem_opt=17`): CO₂, CO, CH₄ each split into anthropogenic, biogenic, ocean,
  background, biomass-burning (`*_BBU`), and test tracers — enabling source attribution.
- **06 UTC receptors.** STILT receptors are released at 06 UTC (~13:00 WIB, the well-mixed
  afternoon boundary layer), aligned with the ~07 UTC Bukit Kototabang flask and the OCO-2 overpass.
  The earlier 12 UTC (evening, stable-layer) configuration is superseded.
- **Case-study window: Oct 6–27, 2015** — the fire peak with maximum observational overlap (all
  OCO-2 granules Oct 7–26 + the BKT flask Oct 27), staging only ~230 GB of `wrfout`.
- **CO joint constraint disabled** (`USE_CO_I7 = False`): the WRF-`CO_BBU`-vs-model check is
  circular; it is retained behind a flag for future integration of a real satellite CO product.

---

## Documentation

| Document | What it covers |
|---|---|
| [`docs/simulation.md`](docs/simulation.md) | Authoritative WRF-GHG science, run history, and bug log |
| [`docs/plan-stilt_inversion.md`](docs/plan-stilt_inversion.md) | Scientific design of the v3 inversion (network, priors, math, targets) |
| [`docs/stilt_inversion_agent_guide.md`](docs/stilt_inversion_agent_guide.md) | Operator runbook for STILT + inversion (download → run → check) |
| [`docs/stilt-inversion-analysis.md`](docs/stilt-inversion-analysis.md) | Written-up v1/v2 analysis and results |

---

## License

Released under the **GNU General Public License v3.0** — see [`LICENSE`](LICENSE).

## Citation & contact

If you use this work, please cite the repository and the underlying datasets (WRF-Chem, FINN,
CarbonTracker, OCO-2). Maintained by **Alberth Nahas** (BMKG) —
[github.com/alberthnahas/wrf-grk](https://github.com/alberthnahas/wrf-grk).
