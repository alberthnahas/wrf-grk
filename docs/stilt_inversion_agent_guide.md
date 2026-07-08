# STILT & Inversion — Agent Operating Guide

**Audience:** an AI agent picking up the STILT footprint + Bayesian inversion work for the
`IDN_BB_2015` (Indonesia biomass-burning, October 2015) case *after* the WRF-GHG simulation
has produced its `wrfout` files.

**Purpose:** a single, self-contained runbook — download, process, run, check, and record —
so a fresh agent can go from "WRF finished" to "posterior emissions + diagnostics" without
re-deriving the pipeline. This does **not** cover running WRF-GHG itself; for that see the
`wrf-ghg` skill and `docs/simulation.md`.

---

## 0. Read first (do not skip)

Before touching anything, read these in order. They are the authoritative state; this guide
is the procedure.

1. `docs/agent_memory.md` — current live operational state (which segment is running, disk, etc.)
2. `docs/simulation.md` — WRF science/bug history (needed to know which `wrfout` days exist and are valid)
3. `docs/plan-stilt_inversion.md` — the scientific design of the v3 inversion (network, priors, math, targets)
4. `docs/stilt-inversion-analysis.md` — the written-up v1/v2 analysis and results
5. `/home/igrk/.agents/skills/wrf-ghg/SKILL.md` — workflow rules and helper scripts

**Working directory is always `/home/igrk/WRF-GRK`.** Activate the project venv for all Python:
`source .venv/bin/activate`.

**Golden rule:** STILT uses `wrfout` **only for meteorology**. The GHG tracers (`CO2_*`, `CO_*`,
`CH4_*`) are consumed later, by the inversion, not by STILT. Keep these two roles separate in
your head — it prevents most confusion.

---

## 1. Pipeline overview

```
 WRF-GHG wrfout (hourly, 1 daily file @ 01:00 UTC, chem_opt=17)
        │
        ├──[meteorology]──► arw2arl ──► stilt/met/YYYY-MM-DD.arl ──► STILT ──► footprints (*_foot.nc)
        │                                                                              │
        └──[GHG tracers: CO_BBU, CO2_BCK]──────────────────────────┐                   │
                                                                    ▼                   ▼
   External data:  FINN priors, CT2025, OCO-2, BKT flask, MODIS ──► run_inversion_v3.py
                                                                              │
                                                                              ▼
                                              posterior α, Tg CO₂ per region, diagnostics + plots
```

Five stages, executed in order:

| Stage | What | Primary script/tool | Output |
|---|---|---|---|
| 1 | Convert `wrfout` → ARL met | `stilt/convert_wrf2arl.sh` (`arw2arl`) | `stilt/met/*.arl` |
| 2 | Download/prepare external obs | `inversion/get_oco2.sh`, preprocessing scripts | `rawdata/*`, `inversion/data/*` |
| 3 | Run STILT back-trajectories | `stilt/run_stilt_full_oct.r` | `stilt/out/footprints/*_foot.nc` |
| 4 | Sample WRF tracers at receptors | Python (see §5) | receptor CO/CO₂ table |
| 5 | Run Bayesian inversion | `inversion/run_inversion_v3.py` | posterior + `inversion/plots/*` |

---

## 2. Prerequisites — what must be true before you start

- **All required `wrfout` days exist and are valid.** STILT runs a **−120 h (5-day) backward**
  trajectory, so a receptor on day *N* needs met for days *N−5 … N*. For the selected **Oct 6–27**
  window (§2.1) you must stage and convert **Oct 1–27** (~27 daily files). The earliest receptor is
  **Oct 6** (needs Oct 1). Confirm with `docs/agent_memory.md` which segments are complete. A gap on
  *any* single day truncates every trajectory that passes through it — do not start STILT for a
  window whose 5-day lookback is incomplete.
- `wrfout` filenames are `wrfout_d01_YYYY-MM-DD_01:00:00`, one per day, hourly frames
  (`history_interval = 60`, `frames_per_outfile = 24`). Treat the `Times` variable inside the
  file as the time truth.
- Disk: ARL files are ~240 MB/day (~7.5 GB for 31 days). STILT `by-id/` particle output is the
  bulkier consumer. Check free space before launching; do not launch if near the budget.

### 2.1 Selecting the case-study window (why October, and how long)

The full WRF-GHG simulation spans **2015-07-28 → 2015-12-01** (~127 days). At **~8.5 GB/day**
(hourly, 232×128×50, ~15 GHG tracers) that is **~1.1 TB of `wrfout`** — infeasible to keep online
or to run STILT/inversion across in full. The `wrfout` archive is backed up off-box; **only the
chosen window needs to be staged locally.** You do **not** need all 1.1 TB — you need
`[receptor_start − 5 days … receptor_end]` (the −120 h STILT lookback).

**Four constraints pin the window to October 2015:**

1. **Scientific significance (why this event at all).** The 2015 El Niño Indonesian fires were the
   largest since 1997: the "year without rain" (SPI < −2.5) drove extreme peat combustion in
   Sumatra and Kalimantan. Emissions reached **~0.8–1.0 Pg C**, and on peak days Indonesia emitted
   **more CO₂ than the entire EU economy** (Huijnen et al. 2016; Field et al. 2016). This is a
   globally significant carbon anomaly — a defensible, high-impact case study. **October is the
   peak month.**
2. **Satellite constraint (OCO-2).** Usable granules exist only for **Oct 7–26** (GES DISC has none
   over Indonesia for Oct 27–31). This is the column-CO₂ constraint and it hard-caps the useful
   upper end of the window.
3. **Real in-situ constraint (BKT flask).** The single genuine measurement is **Oct 27** (smoke-
   flagged, ~10.4 ppm fire enhancement) — the most valuable observation in the system.
4. **STILT lookback.** −120 h means receptors on day *N* need met back to *N−5*; with met from
   Oct 1 the earliest receptor is **Oct 6**.

**Recommended windows** (all justifiable; **Option A is the one selected for this study**):

| Option | Receptor / inversion window | `wrfout` to stage (recept − 5 d) | Storage | Best for |
|---|---|---|---:|---|
| **✅ A — full obs-overlap (SELECTED)** | **Oct 6–27** (22 d) | Oct 1–27 (~27 d) | **~230 GB** | Robust monthly regional CO₂ budget: uses all 20 OCO-2 granules **and** the BKT flask, spans the full fire peak |
| B — peak fortnight | Oct 14–29 (16 d) | Oct 9–29 (~21 d) | ~180 GB | Focus on the most intense burning + flask; leaner compute, slightly fewer OCO-2 days |
| C — plume / methods demo | Oct 22–29 (8 d) | Oct 17–29 (~13 d) | ~110 GB | Transport/plume case study or method demonstration; too short for a stable budget |

**Chosen: Option A — Oct 6–27, 2015.** It maximizes the independent observational constraint
(satellite + in-situ) on the largest single fire event of the decade while staging only ~230 GB.
This is already set in `run_inversion_v3.py` (`DATES`/`DATES_NOON` = Oct 6–27, and the OCO-2 day
filter capped at Oct 27). **Action for the operator:** stage `wrfout` for **Oct 1–27** (window +
5-day lookback) and convert those days to ARL (§3). The in-progress STILT run may produce
footprints out to Oct 31; that is harmless — the inversion simply uses the Oct 6–27 subset
(22 days × 11 sites = **242 receptor-days**).

> The general method for adapting this to *another* event (different region/period) is in §13.

---

## 3. Stage 1 — Convert `wrfout` to ARL meteorology

STILT/HYSPLIT ingests HYSPLIT-ARL packed met, not NetCDF. `arw2arl` (in `stilt/exe/`) does the
conversion, reading its variable map from `stilt/exe/WRFDATA.CFG`.

### 3.1 Variables extracted (from `WRFDATA.CFG` — do not change without reason)

| Dim | WRF variable | → ARL name | Role |
|---|---|---|---|
| 3D | P | PRES | pressure |
| 3D | T | TEMP | temperature |
| 3D | U, V, W | UWND, VWND, WWND | winds (advection) |
| 3D | W→DIFW | DIFW | vertical mixing |
| 3D | QVAPOR | SPHU | specific humidity |
| 2D | HGT | SHGT | terrain height |
| 2D | PSFC | PRSS | surface pressure |
| 2D | RAIN | TPP1 | precip |
| 2D | PBLH | PBLH | **boundary-layer depth (critical for footprints)** |
| 2D | UST | USTR | friction velocity u\* |
| 2D | SWDOWN | DSWF | shortwave down |
| 2D | HFX, LH | SHTF, LHTF | **surface heat fluxes (drive mixing depth)** |
| 2D | T2, U10, V10 | T02M, U10M, V10M | 2 m/10 m surface fields |

The load-bearing fields for footprint quality are **PBLH, UST, HFX, LH, W**. If any come out
zero/degenerate the footprints are silently wrong — validate them (§7).

### 3.2 Procedure

```bash
cd /home/igrk/WRF-GRK
# 1. Hard-link wrfout to colon-free .nc names (arw2arl dislikes ':' in filenames).
#    Links live in stilt/met/wrfout_links/  e.g. wrfout_d01_2015-10-25.nc
#    (create links for any new days; use hard links, not copies, to save disk)
# 2. Run the batch converter (skips days already converted):
bash stilt/convert_wrf2arl.sh
# Output: stilt/met/YYYY-MM-DD.arl  (+ per-day *_arw2arl.log)
```

`convert_wrf2arl.sh` runs `arw2arl` from `stilt/exe/` (so it finds `WRFDATA.CFG`/`ARLDATA.CFG`),
which writes `ARLDATA.BIN`, then renames it to `stilt/met/<date>.arl`.

### 3.3 Check

```bash
ls -lh stilt/met/*.arl                       # expect ~240 MB each, one per day
grep -iE "ERROR|WARNING" stilt/met/*_arw2arl.log   # should be clean
ls stilt/met/*.arl | wc -l                   # expect 31 for the full run
```

---

## 4. Stage 2 — Download and prepare external observations

The inversion consumes five external datasets. Check `docs/plan-stilt_inversion.md §2.3` for the
authoritative status table; most are already present. Download only what is missing.

| Dataset | Path | How to obtain |
|---|---|---|
| **FINN** fire priors (CO₂, CO, CH₄) | `rawdata/finn_fire/emissions-finnv2.5modvrs_{CO2,CO,CH4}_bb_*_0.1x0.1.nc` | Already downloaded (full-year 2015 daily 0.1°). |
| **CT2025** reanalysis | `rawdata/carbontracker/CT2025.molefrac_glb3x2_2015-10-*.nc` | 31 files present; NOAA CarbonTracker 2025, var `co2` [mol/mol]. |
| **OCO-2** Lite FP v11.1r | `inversion/data/oco2_LtCO2_1510*_B11100Ar_*.nc4` | `bash inversion/get_oco2.sh` (needs Earthdata `~/.netrc`; server has Oct 7–26 only). |
| **BKT flask** (real in-situ) | `inversion/data/co2_bkt_flask.txt` | NOAA GML; Oct 27 2015 value already recorded (409.51 ppm). |
| **MODIS hotspots** | `rawdata/hotspot/archived_hotspot_idn.csv` | Already present (~80k rows). |

### OCO-2 download notes
- Requires a valid NASA Earthdata login in `~/.netrc` + `~/.urs_cookies` (see header of
  `inversion/get_oco2.sh`). If a download returns <1000 bytes it failed (usually auth) — fix
  credentials, don't proceed with a truncated granule.
- Only Oct 7–26 exist on the server; Oct 27–31 have no granules over Indonesia. That is expected,
  not an error — those days are covered by CT2025 in the inversion.

### CT2025 preparation
CT2025 is used as pseudo-observations. If not already processed, the preprocessing lives under
`preprocessing/boundary_conditions/` and the inversion reads the raw `.nc` directly (level 0
~981 hPa, time index 3 ≈ 10:30 UTC). Verify all 31 files are non-truncated.

---

## 5. Stage 3 — Run STILT back-trajectories

### 5.1 Inputs
- **Met:** `stilt/met/YYYY-MM-DD.arl` (Stage 1).
- **Receptors (canonical):** `stilt/receptors_full_oct_06utc.csv` — columns `run_time,lati,long,zagl`.
  11 sites × Oct 6–31 = **286 rows**, all at **06:00 UTC**. `zagl = 10` m (surface in-situ equivalent).
- **Run script (canonical):** `stilt/run_stilt_v3p2_06utc.r`, which reads the 06-UTC receptors and
  writes footprints to **`stilt/out_06utc/`** (kept separate from the older noon run).

**Why 06 UTC — this is the correct receptor time, not an option:**
The single real in-situ constraint, the **Bukit Kototabang (BKT) flask, is sampled at ~07 UTC**.
At ~100 °E (WIB = UTC+7), 06–07 UTC is **≈13:00–14:00 local — the afternoon, well-mixed boundary
layer**, exactly the condition a surface receptor should represent. The earlier **noon-UTC (12 UTC)
run is legacy and superseded**: 12 UTC ≈ 19:00 WIB is evening, a stabilizing/decoupling boundary
layer that does **not** represent well-mixed daytime conditions and does not match the flask
sampling time. **Use the 06-UTC files. Do not use `receptors_full_oct.csv` /
`run_stilt_full_oct.r` / `stilt/out/` for the current inversion** — those are the noon run.

### 5.2 Key parameters (already set in the run script — for reference)
`n_hours=-120`, `numpar=200`, `n_cores=8`, footprint grid `xmn/xmx=90/145`, `ymn/ymx=-15/15`,
`xres=yres=0.25°`, `time_integrate=TRUE`, `met_file_format='%Y-%m-%d.arl'`,
`met_file_tres='24 hours'`. The script **skips simulations whose `by-id/` output already exists**,
so re-runs are cheap and idempotent.

### 5.3 Run

```bash
cd /home/igrk/WRF-GRK/stilt
nohup Rscript run_stilt_v3p2_06utc.r > /tmp/stilt_06utc.log 2>&1 &
# ~3–4 h for the full remaining network on 8 cores. Do NOT pipe the live log through a
# monitoring process that can SIGPIPE the R job — tail the file instead:
tail -f /tmp/stilt_06utc.log
```

### 5.4 Output
- `stilt/out_06utc/footprints/{timestamp}_{lon}_{lat}_{zagl}_foot.nc` — NetCDF, shape `(1,120,220)`,
  variable `foot`, units `ppm / (µmol m⁻² s⁻¹)`.
- `stilt/out_06utc/by-id/{simulation_id}/` — particle trajectories.
- Target count: **286 footprints** (see `plan §4.3`).

### 5.5 Check (mandatory before inversion)

```bash
source /home/igrk/WRF-GRK/.venv/bin/activate
python3 - <<'EOF'
import glob, netCDF4 as nc
fps = sorted(glob.glob('/home/igrk/WRF-GRK/stilt/out_06utc/footprints/*_foot.nc'))
print(f'{len(fps)} footprints found')          # expect 286
bad=[]
for f in fps:
    try:
        with nc.Dataset(f) as ds:
            v = ds['foot'][:]
            if v.max()==0: bad.append(f)        # all-zero = failed transport
    except Exception as e:
        bad.append(f'{f}: {e}')
print(f'{len(bad)} bad footprints:', bad[:5])
EOF
```

Any all-zero or unreadable footprint must be re-run (delete its `by-id/` dir so the script
regenerates it) before proceeding.

---

## 6. Stage 4 — WRF GHG tracers at receptors (optional / dormant)

chem_opt=17 writes these tracers (confirmed present in `wrfout`):

- **CO2:** `CO2_ANT, CO2_BIO, CO2_OCE, CO2_BCK, CO2_TST, CO2_BBU`
- **CO:** `CO_ANT, CO_BCK, CO_BBU, CO_TST`
- **CH4:** `CH4_ANT, CH4_BIO, CH4_BCK, CH4_BBU, CH4_TST`

**Current state — the CO constraint (I7) is DISABLED** (`USE_CO_I7 = False` in
`run_inversion_v3.py`). It was found to be *circular*: WRF `CO_BBU` is produced by forcing FINN
with α=1, then compared against `H_STILT × FINN_CO`, which pulls the posterior toward the WRF
input α rather than truth. The code path is retained behind the flag for future integration of a
**real** satellite CO product (TROPOMI/MOPITT), which would make it a genuine independent
constraint. **Do not re-enable `USE_CO_I7` with WRF `CO_BBU` as the "observation."**

So in the active v3.2 inversion the observation vector is **CT2025 + BKT flask + OCO-2 only** —
no WRF tracer is sampled. If you *do* work in the dormant CO path (or add a real CO obs), the
sampler `load_wrf_tracer_at_site` reads `CO_BBU` at the **lowest model level** (`bottom_top=0`,
≈ surface / 10 m AGL) at **06 UTC** — already fixed to 06 UTC to match receptors/footprints/flask.

**Unit note (for accuracy):** WRF-GHG writes CO₂/CH₄/CO in **ppmv**, and the script keeps the CO
path in ppm throughout (`CO_OBS_SIGMA = 0.05` ppm = 50 ppb). There is **no ppmv↔ppb conversion
bug** — just keep any new obs in the same ppm units, do not mix ppb.

`CO2_BBU` (fire CO₂) and `CO2_BCK` (background CO₂) remain useful for **validation/sanity** — e.g.
checking footprint×FINN enhancement against the WRF forward field — but they are not inputs to the
active solve.

---

## 7. Stage 5 — Run the Bayesian inversion

### 7.1 Inputs consumed by `inversion/run_inversion_v3.py`
- STILT footprints from `stilt/out_06utc/` (§5) — builds the Jacobian **H**.
- FINN priors (§4) — regional prior emissions `E_k`.
- **Observation vector `y` = CT2025 pseudo-obs + BKT flask + OCO-2 grid cells.** (WRF `CO_BBU`
  is **not** used — the CO constraint is disabled, `USE_CO_I7 = False`; see §6.)
- MODIS hotspots — modulate `S_prior` (§5.5 of the plan).
- 6-region definition + prior α/σ — from `plan §6.2` (note: v3.1/v3.2 redrew the Kalimantan/
  Sulawesi region boxes and re-tuned priors relative to the plan — the script is authoritative).

### 7.2 Run

```bash
cd /home/igrk/WRF-GRK
source .venv/bin/activate
python3 inversion/run_inversion_v3.py 2>&1 | tee inversion/logs/inversion_v3.log
```

### 7.3 Outputs
- Posterior α (6 regional scaling factors) + Tg CO₂ per region, in the log.
- `inversion/plots/inversion_fit_v3.png`, `posterior_emissions_v3.png`,
  `inversion_diagnostics_v3.png`, `spatial_prior_vs_posterior_v3.png`.
- Additional figures via `inversion/plot_spatial_v3.py`, `plot_stilt_footprints.py`,
  `plot_literature_comparison.py`.

### 7.4 Validation checks (from `plan §9.3` — the run is not trustworthy until these pass)

1. **χ² reduction:** posterior χ² < prior χ² / 10.
2. **Non-negativity:** all posterior α ≥ 0 (SLSQP-bounded).
3. **Physical ordering:** Kalimantan S+E emission **>>** Sulawesi+East. If Sulawesi > Kalimantan,
   there is a structural problem — stop and diagnose (this was the v2 failure).
4. **BKT fit:** posterior Δc at BKT on Oct 27 within 2σ of the 10.4 ppm observation.
5. **OCO-2 residuals:** mean < 1 ppm (no systematic bias).
6. **Uncertainty reduction:** posterior σ_k < prior σ_k for every region.

If check 3 fails (Sulawesi ≥ Kalimantan), the levers are the region geometry, the OCO-2 IDW
distance cutoff (`OCO2_MAX_DIST_DEG = 6.0`), and the hotspot-modulated priors — **not** the CO
constraint, which is disabled. v3.1/v3.2 already fixed the main v3.0 cause (mis-drawn Kalimantan
peat corridor); re-check `REGIONS`/`HOTSPOT_REGION_BOUNDS` if it recurs.

---

## 8. Consolidated check-gates (run in order; do not skip a gate)

| Gate | Command / check | Pass condition |
|---|---|---|
| G1 met complete | `ls stilt/met/*.arl \| wc -l` | 31 (or all days in your window +5 lookback) |
| G2 met clean | `grep -iE "ERROR" stilt/met/*_arw2arl.log` | no matches |
| G3 footprints present | `ls stilt/out_06utc/footprints/*_foot.nc \| wc -l` | 286 (or expected N) |
| G4 footprints valid | §5.5 script | 0 bad |
| G5 obs present | OCO-2/CT2025/FINN/flask paths | all files, non-truncated |
| G6 time alignment | receptor / footprint / CT2025 / flask all ~06–07 UTC | consistent (see §5.1, §6) |
| G7 inversion sane | §7.4 checks 1–6 | all pass |

---

## 9. Record what you learned — memories

After a meaningful run, persist durable facts to the memory store at
`/home/igrk/.claude/projects/-home-igrk-WRF-GRK/memory/` (one fact per file, with frontmatter),
and add a one-line pointer to `MEMORY.md`. Save memories for:

- **`project`** — the current pipeline state that is *not* derivable from code/git: e.g. "STILT
  full-Oct footprints complete (286/286) as of <date>", "inversion v3 posterior: Kalimantan S+E
  ~X Tg, Sulawesi ~Y Tg", the chosen case-study window (§2.1).
- **`feedback`** — corrections the user gave on how to run things (include **Why** and
  **How to apply**), e.g. a preferred σ_CO, a decision to exclude a site/day.
- **`reference`** — external pointers: OCO-2 GES DISC archive URL, Earthdata account note,
  CT2025 source.

Do **not** memorialize things already in the repo (script contents, file layout, the plan's math)
— point to the doc instead. Before writing, check for an existing memory covering the same fact
and update it rather than duplicating. Link related memories with `[[slug]]`.

Convert relative dates to absolute when saving.

---

## 10. Update the skill

The `wrf-ghg` skill (`/home/igrk/.agents/skills/wrf-ghg/SKILL.md`) already lists STILT/inversion
among its triggers. When you add or change a procedure here:

- Keep **procedure** in this guide (`docs/stilt_inversion_agent_guide.md`) and have the skill
  *point* to it — do not duplicate long content into `SKILL.md`.
- If you create a reusable helper (e.g. a WRF-tracer sampler, a footprint validator), put the
  script in `/home/igrk/.agents/skills/wrf-ghg/scripts/` alongside the existing
  `check_seg_inputs.py`, `monitor_run.sh`, etc., and reference it from the skill.
- Add a reference note under `/home/igrk/.agents/skills/wrf-ghg/references/` for any new
  data-format or workflow gotcha.
- If the STILT/inversion work grows large enough to warrant its own skill, create a
  `stilt-inversion` skill with a `description` covering: convert wrfout→ARL, run STILT footprints,
  sample WRF tracers, run the v3 Bayesian inversion, download OCO-2/CT2025, and interpret
  posterior scaling factors — and have it reference this guide as authoritative.

---

## 11. Known pitfalls (learn from these before repeating them)

| Pitfall | Symptom | Fix |
|---|---|---|
| Incomplete 5-day met lookback | Footprints truncate / go to zero mid-trajectory | Ensure all days *N−5…N* converted before running receptors on day N |
| Colons in wrfout filenames | `arw2arl` fails to open input | Use colon-free hard links in `stilt/met/wrfout_links/` |
| Degenerate met field (PBLH/UST/HFX=0) | Footprints wrong but non-zero | Validate ARL/met fields; re-check the source `wrfout` |
| Re-enabling CO constraint with WRF `CO_BBU` | Posterior pulled toward WRF input α (circular) | Keep `USE_CO_I7 = False`; only enable with a *real* satellite CO obs — see §6 |
| Mixing ppb into the ppm CO path | σ off by 1000× | CO path is all in **ppm** (`CO_OBS_SIGMA = 0.05`); do not introduce ppb |
| IDW "virtual" footprints | Sulawesi > Kalimantan (unphysical) | **Do not interpolate footprints** — use real STILT trajectories only (the entire point of v3) |
| Wrong receptor time (12 UTC) | Evening stable BL, mismatches BKT flask | Use **06 UTC** (`*_06utc`, `out_06utc/`) — already the code default (§5.1) |
| OCO-2 truncated download | granule <1 KB | Fix Earthdata `~/.netrc`, re-download |
| Live-log pipe SIGPIPE | STILT R job killed | Run under `nohup`, `tail -f` the file, never pipe the running job |
| Fire-zone OCO-2 bias | Anomalous XCO₂ >25 ppm from smoke | I8 exclusion already drops cells >25 ppm; keep it |

---

## 12. Quick command reference

```bash
cd /home/igrk/WRF-GRK && source .venv/bin/activate

# Stage 1 — met
bash stilt/convert_wrf2arl.sh
ls stilt/met/*.arl | wc -l

# Stage 2 — obs (only if missing)
bash inversion/get_oco2.sh

# Stage 3 — STILT (06-UTC canonical run)
cd stilt && nohup Rscript run_stilt_v3p2_06utc.r > /tmp/stilt_06utc.log 2>&1 &
tail -f /tmp/stilt_06utc.log
ls out_06utc/footprints/*_foot.nc | wc -l  # expect 286

# Stage 5 — inversion
cd /home/igrk/WRF-GRK
python3 inversion/run_inversion_v3.py 2>&1 | tee inversion/logs/inversion_v3.log
# then review inversion/plots/*_v3.png and validate §7.4
```

---

## 13. Reusing this guide for another case study

This runbook is written for `IDN_BB_2015`, but the pipeline (WRF-GHG → ARL → STILT footprints →
Bayesian inversion) is generic. To adapt it to a new fire episode / region / period, an agent
should change only the case-specific inputs and keep the five-stage procedure intact:

| What changes per case | Where |
|---|---|
| Domain, dates, `wrfout` set | WRF-GHG run (`docs/simulation.md`, `wrf-ghg` skill) — must finish first |
| Receptor network (sites, lat/lon) | `stilt/receptors_*.csv` — pick sites that give real footprint coverage of every emission region (do **not** interpolate) |
| **Receptor time** | Choose the **local afternoon, well-mixed PBL** hour in UTC, and match it to any real in-situ sampling time (for IDN_BB_2015 that is 06 UTC ≈ 13:00 WIB, matching the 07 UTC BKT flask). Re-derive it per case from the site's longitude/time zone — do not copy 06 UTC blindly. |
| Back-trajectory length | `n_hours` — long enough to reach all source regions, short enough to stay in the WRF domain |
| Footprint grid | `xmn/xmx/ymn/ymx/xres/yres` to the new domain |
| Prior emissions + regions | FINN/GFED files, region boxes, prior α/σ in the inversion script |
| Observations | CT2025, satellite (OCO-2/other), in-situ flask, hotspot files for the new period |

**Invariants to preserve in any case study** (these are the lessons, not the parameters):
1. STILT uses `wrfout` for meteorology only; GHG tracers enter at the inversion.
2. Met must cover the full back-trajectory lookback for every receptor day.
3. Receptor time, WRF tracer sampling time, footprint time, and in-situ obs time must **all agree**.
4. Use real STILT footprints for every region — never IDW "virtual" footprints.
5. Keep tracer units consistent (WRF-GHG is ppmv; convert to the obs unit).
6. Validate footprints (no all-zero) and run the §7.4 physical/statistical checks before trusting results.

When you start a new case, copy this file to `docs/stilt_inversion_agent_guide_<CASE>.md`, update
the case-specific rows above, and record the new decisions as memories (§9).

---

*Maintained as the agent-facing runbook for STILT + inversion on IDN_BB_2015. The scientific
rationale lives in `docs/plan-stilt_inversion.md`; the WRF-GHG side lives in `docs/simulation.md`
and the `wrf-ghg` skill. Update this file when the procedure changes.*
