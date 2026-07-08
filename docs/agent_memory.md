# Agent Memory: WRF-GRK / IDN_BB_2015

Last updated: 2026-07-08 Asia/Jakarta.

This file is the repo-local continuity record for future agent sessions. It complements, but does not replace, `docs/simulation.md`. For authoritative scientific history and known bugs, read `docs/simulation.md` first.

## Current Live State

- Case: `simulations/IDN_BB_2015`
- Model: WRF-Chem 4.1.5, GHG tracer mode, `chem_opt=17`
- Domain: 27 km Mercator, 232 x 128 x 49 mass levels in current wrfout files
- Canonical simulation period: `2015-07-28_00:00:00` to `2015-12-01_00:00:00`
- Segment plan:
  - seg1: `2015-07-28 -> 2015-08-29`, `run_hours=768`, cold start
  - seg2: `2015-08-29 -> 2015-09-30`, `run_hours=768`, restart
  - seg3: `2015-09-30 -> 2015-10-24`, `run_hours=576`, restart
  - seg4: `2015-10-24 -> 2015-11-25`, `run_hours=768`, restart
  - seg5: `2015-11-25 -> 2015-12-01`, `run_hours=144`, restart
- Current run at time of this memory update:
  - seg4 is running from `wrfrst_d01_2015-10-24_00:00:00`
  - launch command: `setsid -f mpirun -np 16 ./wrf.exe > wrf.stdout.log 2>&1`
  - host-level process check: 16 `wrf.exe` ranks under `mpirun` PID `4018670`
  - latest checked model time: `2015-10-24_04:17:30`
  - no `FATAL` or `---- ERROR` lines found in current `rsl.out.0000` / `rsl.error.0000`
  - free disk at last check: `526G`
  - active namelist: seg4, `run_hours=768`, `2015-10-24_00:00:00 -> 2015-11-25_00:00:00`, `restart=.true.`
  - note: sandboxed `pgrep` may not see detached WRF ranks; use host/escalated process checks for live runs.

## Rules That Must Not Be Forgotten

- Segment boundaries must be real MODIS 8-day composite dates. This prevents VPRM restart read failures.
- Always check inputs before launch using:

```bash
.venv/bin/python /home/igrk/.agents/skills/wrf-ghg/scripts/check_seg_inputs.py <start> <end>
```

- Use `mpirun -np 16`, not 28 or 32, for WRF on this host.
- Do not launch a segment unless free space is enough. A 32-day segment can require roughly:
  - `~272G` wrfout
  - `~50G` intermediate restarts
  - extra scratch margin
- Prefer `>400G` free before launching a full 32-day segment.
- Back up or transfer segment wrfout before deleting local files.
- Delete intermediate `wrfrst_d01_*` only after the segment is complete and the segment-end restart is confirmed.
- For restart segments, wrfout filenames start at `_01:00:00`. This is normal. Use internal `Times`, not filename, for analysis.
- CH4 and CO tracer values in wrfout are physically in ppb even though WRF unit metadata may say `ppmv`.
- CO background `CO_BCK` is intentionally constant around `80 ppb` because there is no local CarbonTracker CO mole-fraction dataset.

## Seg1 Outcome

- seg1 completed successfully at `2015-08-29_00:00:00`.
- Success line observed:

```text
d01 2015-08-29_00:00:00 wrf: SUCCESS COMPLETE WRF
```

- seg1 runtime estimate: about `6 h 40-45 min` wall-clock.
- seg1 generated:
  - 32 full daily wrfout files from `2015-07-28_00:00:00` through `2015-08-28_00:00:00`
  - 1 final one-frame endpoint file: `wrfout_d01_2015-08-29_00:00:00`
  - final restart: `wrfrst_d01_2015-08-29_00:00:00`
- User backed up seg1 wrfout externally.
- Local seg1 wrfout files were deleted before launching seg2.
- Intermediate seg1 restart files were deleted; only `wrfrst_d01_2015-08-29_00:00:00` was kept for seg2.
- Manifest created before deletion:

```text
simulations/IDN_BB_2015/logs/seg1_wrfout_manifest.csv
```

## Seg2 Preparation And Launch

- Before seg2, local seg1 wrfout and intermediate restarts were removed, freeing disk from about `230G` to about `535G`.
- seg2 namelist was set to:

```text
run_hours = 768
start = 2015-08-29_00:00:00
end = 2015-09-30_00:00:00
restart = .true.
```

- Preflight passed:

```text
chemi: 33/33 OK
fire: 33/33 OK
vprm: 5/5 OK
oce: 2/2 OK (2015-08-29, 2015-09-28)
wrfbdy: covers through 2015-11-30_12:00:00
```

- seg2 startup read:
  - `wrfrst_d01_2015-08-29_00:00:00`
  - `wrfchemi_d01_2015-08-29_00:00:00`
  - `wrffirechemi_d01_2015-08-29_00:00:00`
  - `vprm_input_d01_2015-08-29_00:00:00`
  - `wrfbdy_d01`
- seg2 first wrfout file is `wrfout_d01_2015-08-29_01:00:00`. This is expected for a restart run.

## Seg2 Completion And Cleanup

- seg2 completed successfully at `2015-09-30_00:00:00`.
- Success line observed:

```text
d01 2015-09-30_00:00:00 wrf: SUCCESS COMPLETE WRF
```

- Final seg2 restart exists and is required for seg3:

```text
simulations/IDN_BB_2015/run/wrfrst_d01_2015-09-30_00:00:00
```

- User backed up seg2 wrfout externally before cleanup.
- Verified the last seg2 wrfout filename is normal for restart output:
  - file: `wrfout_d01_2015-09-29_01:00:00`
  - frames: 24
  - first internal `Times`: `2015-09-29_01:00:00`
  - last internal `Times`: `2015-09-30_00:00:00`
- Before cleanup, local free space was only `228G`, not enough for seg3.
- Saved seg2 records:

```text
simulations/IDN_BB_2015/logs/seg2_wrfout_manifest.csv
simulations/IDN_BB_2015/logs/seg2_rsl.out.0000
simulations/IDN_BB_2015/logs/seg2_rsl.error.0000
```

- Deleted local backed-up seg2 wrfout files and intermediate September restarts.
- Kept only:

```text
simulations/IDN_BB_2015/run/wrfrst_d01_2015-09-30_00:00:00
```

- Disk after cleanup: about `514G` free.

## Seg3 Preparation And Launch

- Seg3 preflight passed:

```text
Checking inputs for segment 2015-09-30 -> 2015-10-24
chemi: 25/25 OK
fire: 25/25 OK
vprm: 4/4 OK (MODIS composites in range)
oce: 1/1 OK (auxinput6 read times: ['2015-09-30'])
wrfbdy: covers through 2015-11-30_12:00:00 (segment end 2015-10-24) OK
```

- Seg3 namelist was set to:

```text
run_hours = 576
start = 2015-09-30_00:00:00
end = 2015-10-24_00:00:00
restart = .true.
```

- First attempted background launch with `nohup mpirun ... &` exited before creating `rsl.out.0000`; `wrf.stdout.log` was empty.
- Direct foreground `mpirun -np 16 ./wrf.exe` proved the setup was valid and all 16 ranks started.
- Because the foreground session could not be detached, it was intentionally stopped, and the partial `wrfout_d01_2015-09-30_01:00:00` plus interrupted rsl logs were deleted.
- Clean detached launch succeeded with:

```bash
setsid -f mpirun -np 16 ./wrf.exe > wrf.stdout.log 2>&1
```

- Startup checks showed:
  - restart read: `wrfrst_d01_2015-09-30_00:00:00`
  - accepted inputs at segment start:
    - `wrfchemi_d01_2015-09-30_00:00:00`
    - `wrffirechemi_d01_2015-09-30_00:00:00`
    - `vprm_input_d01_2015-09-30_00:00:00`
    - `wrfbdy_d01`
  - first seg3 wrfout: `wrfout_d01_2015-09-30_01:00:00`
  - latest checked timing line: `2015-09-30_06:10:00`
  - no fatal/error lines.

## Seg3 Completion And Cleanup

- seg3 completed successfully at `2015-10-24_00:00:00`.
- Success line observed:

```text
d01 2015-10-24_00:00:00 wrf: SUCCESS COMPLETE WRF
```

- Final seg3 restart exists and is required for seg4:

```text
simulations/IDN_BB_2015/run/wrfrst_d01_2015-10-24_00:00:00
```

- User backed up seg3 wrfout externally before cleanup.
- Saved seg3 records:

```text
simulations/IDN_BB_2015/logs/seg3_wrfout_manifest.csv
simulations/IDN_BB_2015/logs/seg3_rsl.out.0000
simulations/IDN_BB_2015/logs/seg3_rsl.error.0000
simulations/IDN_BB_2015/logs/seg2_leftover_wrfout_manifest.csv
```

- Deleted local backed-up seg3 wrfout files, two leftover seg2 wrfout files, and intermediate October restarts.
- Kept only:

```text
simulations/IDN_BB_2015/run/wrfrst_d01_2015-10-24_00:00:00
```

- Disk after cleanup: about `528G` free.
- Runtime estimates from saved logs:
  - seg2: `~7.04 h` from summed WRF `Timing for main` lines; manifest timestamps give a lower-bound/near-total of `6 h 48 m 45 s`.
  - seg3: `~5.24 h` from summed WRF `Timing for main` lines; clean-run file timestamps give `5 h 14 m 30 s`.

## Seg4 Preparation And Launch

- Seg4 preflight passed:

```text
Checking inputs for segment 2015-10-24 -> 2015-11-25
chemi: 33/33 OK
fire: 33/33 OK
vprm: 5/5 OK (MODIS composites in range)
oce: 2/2 OK (auxinput6 read times: ['2015-10-24', '2015-11-23'])
wrfbdy: covers through 2015-11-30_12:00:00 (segment end 2015-11-25) OK
```

- Seg4 namelist was set to:

```text
run_hours = 768
start = 2015-10-24_00:00:00
end = 2015-11-25_00:00:00
restart = .true.
```

- Seg4 launched with:

```bash
setsid -f mpirun -np 16 ./wrf.exe > wrf.stdout.log 2>&1
```

- Startup checks showed:
  - restart read: `wrfrst_d01_2015-10-24_00:00:00`
  - accepted inputs at segment start:
    - `wrfchemi_d01_2015-10-24_00:00:00`
    - `wrffirechemi_d01_2015-10-24_00:00:00`
    - `vprm_input_d01_2015-10-24_00:00:00`
    - `wrfbdy_d01`
  - first seg4 wrfout: `wrfout_d01_2015-10-24_01:00:00`
  - latest checked timing line: `2015-10-24_04:17:30`
  - no fatal/error lines.

## CarbonTracker / Boundary Condition Work Done

New script added:

```text
preprocessing/boundary_conditions/02_patch_carbontracker_bdy_timevarying.py
```

Purpose:

- Patch time-varying CarbonTracker CO2 and CH4 lateral boundary conditions into `wrfbdy_d01`.
- Writes all four side arrays:
  - `BXS`, `BXE`, `BYS`, `BYE`
- Computes and writes all four boundary tendency arrays:
  - `BTXS`, `BTXE`, `BTYS`, `BTYE`
- Leaves CO_BCK as constant 80 ppb; no local CT-CO dataset exists.

Important details:

- One corrupt file was found:

```text
rawdata/carbontracker/CTCH4_2024.molefrac_glb3x2_2015-09-26.nc
```

- The time-varying patcher skips unreadable CT files and uses nearest readable day.
- The original `01_process_carbontracker.py` was patched to choose the closest daily CT file instead of the middle file of the month.
- Verified ranges after time-varying patch:
  - CO2_BCK side arrays varied spatially and temporally; BT arrays nonzero.
  - CH4_BCK side arrays varied spatially and temporally; BT arrays nonzero.
  - CO_BCK side arrays were all 80.00; BT arrays all zero.

## Emission/Input Cadence Confirmed

- History output:
  - `history_interval = 60`
  - `frames_per_outfile = 24`
  - wrfout is hourly.
- Fire emissions:
  - `wrffirechemi_d01_<date>`
  - daily files, 24 hourly frames
  - WRF reads hourly (`auxinput7_interval_m=60`)
- Anthropogenic emissions:
  - `wrfchemi_d01_<date>`
  - daily files, 24 hourly frames
  - WRF reads hourly (`auxinput5_interval_m=60`)
- Ocean flux:
  - `wrfoce_d01_<date>`
  - 1-frame files
  - WRF reads every 30 days (`auxinput6_interval_m=43200`)
  - restart segments re-anchor the read dates, so segment-specific `wrfoce` files are required.
- Biogenic/VPRM:
  - `vprm_input_d01_<date>`
  - 1-frame files
  - WRF reads every 8 days (`auxinput15_interval_m=11520`)
  - VPRM input is MODIS-derived:
    - MOD09A1 8-day surface reflectance for EVI/LSWI
    - MCD12Q1 land cover for VPRM vegetation class

## WRF Output Checks Performed

- seg1 GHG fields were checked in wrfout:
  - CO2/CH4 background fields are non-constant and reasonable after CT patch.
  - CO background remains constant near 80 ppb as expected.
  - CO2_BBU, CH4_BBU, CO_BBU fire tracers are active and grow with fire emissions.
- A seg1 hourly surface timeseries was created:

```text
simulations/IDN_BB_2015/plots/wrfout_hourly_trace_20150728_20150824.csv
```

- Summary from complete files `2015-07-28` through `2015-08-24`:
  - CO2 surface mean first/latest complete: about `399.98 -> 403.07 ppm`
  - CO2_BBU surface max peak: about `2180.72 ppm` at `2015-08-15_12:00`
  - CH4_BBU surface max peak: about `18.54 ppb` at `2015-08-15_12:00`
  - CO_BBU surface max peak: about `193.40 ppb` at `2015-08-15_12:00`

## Scripts Added During This Session

### Essential NetCDF Extractor

```text
scripts/extract_wrf_essentials.py
```

Default output:

```text
simulations/IDN_BB_2015/run/wrfessentials_seg1_idnbb2015.nc
```

Purpose:

- Extract important meteo, all GHG tracers, and fields needed for total-column GHG calculation from wrfout.
- Writes compressed NetCDF4 directly using zlib/shuffle/chunking.
- Defaults:
  - `complevel=5`
  - `least_significant_digit=4`
  - refuses to run unless `--min-free-gb` threshold is satisfied

Dry-run found:

```text
33 seg1 files
769 hourly frames
65 variables kept
raw estimate ~123 GiB
compressed estimate roughly 45-90 GiB, plan for ~100 GiB
```

Do not run this during a segment unless disk margin is clearly safe.

### Timeseries Plotter

```text
scripts/plot_wrf_timeseries.py
```

Purpose:

- Read wrfout directly.
- Produce CSV and PNG timeseries for essential meteo and GHG variables.
- Outputs to:

```text
simulations/IDN_BB_2015/plots/wrf_diagnostics/timeseries/
```

Tested successfully for:

```text
CO2_TOTAL CH4_TOTAL CO_TOTAL CO2_BBU CH4_BBU CO_BBU
```

### Spatial Map Plotter

```text
scripts/plot_wrf_spatial_maps.py
```

Purpose:

- Read wrfout directly.
- Create spatial maps using cartopy and repo `geojson/` overlays.
- Uses Indonesia/province/peat overlays:
  - `geojson/indonesia_38prov.geojson`
  - `geojson/Indonesia_peat_lands.json`
  - `geojson/world_without_idn.json`
- Layout inspired by local BMKG hotspot map:
  - broad Indonesia extent
  - gridlines
  - info panel
  - source metadata panel
  - scientific colormaps by variable family
- Outputs to:

```text
simulations/IDN_BB_2015/plots/wrf_diagnostics/maps/
```

Tested successfully with:

```text
CO2_BBU at latest complete seg1 time
```

Example output:

```text
simulations/IDN_BB_2015/plots/wrf_diagnostics/maps/map_CO2_BBU_2015-08-28-230000.png
```

## GeoJSON And Reference Plot Notes

- Repo geojson files:
  - `geojson/Indonesia_peat_lands.json`
  - `geojson/indonesia_38prov.geojson`
  - `geojson/indonesia_kabkota_38prov.geojson`
  - `geojson/world_without_idn.json`
- `geopandas` and `fiona` are not installed in `.venv`.
- `shapely`, `cartopy`, `matplotlib`, `netCDF4`, `numpy`, `pandas` are available.
- Plot scripts therefore use `json` + `shapely` + `cartopy`, not geopandas.
- User reference image exists as:

```text
/home/igrk/Documents/Hotspot/update_hotspot_1.png
```

not `.jpg`.

## Storage/Cleanup History

- Initial large storage users:
  - `simulations/IDN_BB_2015/output` and logs wrfout files were removed for rerun.
  - stale STILT wrfout copies under `stilt/met/wrfout_links` were removed.
- After that cleanup, free space was about `494G`.
- Seg1 local wrfout later occupied about `255.8 GiB`.
- User backed up seg1 wrfout externally.
- Local seg1 wrfout was deleted before seg2 launch.
- Disk after seg2 launch and progress to `2015-09-18`: about `339G` free.
- After seg2 completion, local seg2 wrfout and intermediate restarts were deleted after user backup.
- Disk before seg3 launch: about `514G` free.
- Disk after seg3 startup check: about `512G` free.
- After seg3 completion, local seg3 wrfout, two leftover seg2 wrfout files, and intermediate October restarts were deleted after backup.
- Disk before seg4 launch: about `528G` free.
- Disk after seg4 startup check: about `526G` free.

## Commands For Future Continuation

Check current WRF. For detached seg3, process checks should be run outside the sandbox/with escalation if needed; sandboxed `pgrep` may return zero even while the host run is alive:

```bash
cd /home/igrk/WRF-GRK
pgrep -acu "$USER" wrf.exe || true
grep "Timing for main" simulations/IDN_BB_2015/run/rsl.out.0000 | tail -1
grep -i "FATAL\|---- ERROR" simulations/IDN_BB_2015/run/rsl.out.0000 simulations/IDN_BB_2015/run/rsl.error.0000 2>/dev/null | tail
df -h .
```

On seg4 completion:

```bash
ls -lh simulations/IDN_BB_2015/run/wrfrst_d01_2015-11-25_00:00:00
grep "SUCCESS COMPLETE WRF" simulations/IDN_BB_2015/run/rsl.out.0000
```

Before seg5:

```bash
.venv/bin/python /home/igrk/.agents/skills/wrf-ghg/scripts/check_seg_inputs.py 2015-11-25 2015-12-01
```

Expected seg5 namelist:

```text
run_hours = 144
start = 2015-11-25_00:00:00
end = 2015-12-01_00:00:00
restart = .true.
```

Do not launch seg5 unless seg4 wrfout is backed up/transferred or free disk remains safely high.

## Known Open/Important Caveats

- `CO_BCK` remains constant 80 ppb.
- CH4 and CO wrfout metadata may say `ppmv`, but values are ppb.
- Land-surface T2 cold bias is documented and mostly affects surface comparisons.
- SST has some zero cells from ERA5/metgrid interpolation gaps.
- `*_TST` tracers duplicate anthropogenic tracers, not total emissions.
- Bug 22 from the old run documented a cumulative-state failure in seg4 around Oct 30 when integrating too many days from Oct 24. If it recurs, consult `docs/simulation.md` Bug 22 before choosing a split/rebuild fix.
