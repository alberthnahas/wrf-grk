#!/usr/bin/env Rscript
# WRF-STILT run script for IDN_BB_2015 (Oct 22-29 2015 fire peak)
#
# Met data: ARL format files in stilt/met/ converted from wrfout via arw2arl
# Receptors: 5 Indonesian sites x 6 days (Oct 24-29 noon UTC) = 30 simulations
#
# Usage:
#   cd /home/igrk/WRF-GRK/stilt
#   Rscript run_stilt.r
#
# Outputs:
#   out/particles/   particle trajectory files
#   out/footprints/  gridded influence footprints (NetCDF)

# ── Paths (must be set before source('r/dependencies.r')) ─────────────────────
stilt_wd  <- '/home/igrk/WRF-GRK/stilt'
output_wd <- file.path(stilt_wd, 'out')
lib.loc   <- .libPaths()[1]

# ── Parallel simulation settings ──────────────────────────────────────────────
n_cores            <- 8      # 8 of 16 cores; leaves headroom for other tasks
n_nodes            <- 1
processes_per_node <- n_cores
slurm              <- FALSE
slurm_options      <- list()

# ── Receptors ─────────────────────────────────────────────────────────────────
# Read from receptors.csv: run_time, lati, long, zagl
receptors <- read.csv(file.path(stilt_wd, 'receptors.csv'),
                      stringsAsFactors = FALSE)
receptors$run_time <- as.POSIXct(receptors$run_time, tz = 'UTC')

# ── Footprint grid ────────────────────────────────────────────────────────────
hnf_plume      <- TRUE
projection     <- '+proj=longlat'
smooth_factor  <- 1
time_integrate <- TRUE
xmn  <-  90.0;  xmx <- 145.0   # lon bounds (WRF domain)
ymn  <- -15.0;  ymx <-  15.0   # lat bounds (WRF domain)
xres <- 0.25;   yres <- 0.25   # footprint resolution (degrees)

# ── Meteorological data ───────────────────────────────────────────────────────
# ARL files named YYYY-MM-DD.arl, one per day, in stilt/met/
met_path        <- file.path(stilt_wd, 'met')
met_file_format <- '%Y-%m-%d.arl'
met_file_tres   <- '24 hours'
met_subgrid_buffer <- 0.1
met_subgrid_enable <- FALSE
met_subgrid_levels <- NA
n_met_min          <- 1

# ── Model control ─────────────────────────────────────────────────────────────
n_hours       <- -120     # backward 5 days
numpar        <- 200
rm_dat        <- TRUE
run_foot      <- TRUE
run_trajec    <- TRUE
simulation_id <- NA
timeout       <- 3600
varsiwant     <- c('time', 'indx', 'long', 'lati', 'zagl', 'foot', 'mlht',
                   'dens', 'samt', 'sigw', 'tlgr')

# ── Transport and dispersion settings ─────────────────────────────────────────
capemin     <- -1
cmass       <- 0
conage      <- 48
cpack       <- 1
delt        <- 1
dxf         <- 1
dyf         <- 1
dzf         <- 0.01
efile       <- ''
emisshrs    <- 0.01
frhmax      <- 3
frhs        <- 1
frme        <- 0.1
frmr        <- 0
frts        <- 0.1
frvs        <- 0.01
hscale      <- 10800
ichem       <- 8
idsp        <- 2
initd       <- 0
k10m        <- 1
kagl        <- 1
kbls        <- 1
kblt        <- 5
kdef        <- 0
khinp       <- 0
khmax       <- 9999
kmix0       <- 150
kmixd       <- 3
kmsl        <- 0
kpuff       <- 0
krand       <- 4
krnd        <- 6
kspl        <- 1
kwet        <- 1
kzmix       <- 0
maxdim      <- 1
maxpar      <- numpar
mgmin       <- 10
mhrs        <- 9999
nbptyp      <- 1
ncycl       <- 0
ndump       <- 0
ninit       <- 1
nstr        <- 0
nturb       <- 0
nver        <- 0
outdt       <- 0
p10f        <- 1
pinbc       <- ''
pinpf       <- ''
poutf       <- ''
qcycle      <- 0
rhb         <- 80
rht         <- 60
splitf      <- 1
tkerd       <- 0.18
tkern       <- 0.18
tlfrac      <- 0.1
tout        <- 0
tratio      <- 0.75
tvmix       <- 1
veght       <- 0.5
vscale      <- 200
vscaleu     <- 200
vscales     <- -1
wbbh        <- 0
wbwf        <- 0
wbwr        <- 0
wvert       <- FALSE
w_option    <- 0
zicontroltf <- 0
ziscale     <- rep(list(rep(1, 24)), nrow(receptors))
z_top       <- 25000

# ── Transport error settings ──────────────────────────────────────────────────
horcoruverr <- NA
siguverr    <- NA
tluverr     <- NA
zcoruverr   <- NA
horcorzierr <- NA
sigzierr    <- NA
tlzierr     <- NA

# ── User hooks ────────────────────────────────────────────────────────────────
before_trajec    <- function() {output}
before_footprint <- function() {output}

# ── Source STILT functions from r/ ────────────────────────────────────────────
setwd(stilt_wd)
source('r/dependencies.r')

# ── Create output directories ─────────────────────────────────────────────────
system(paste0('rm -r ', output_wd, '/footprints'), ignore.stderr = TRUE)
if (run_trajec) {
  system(paste0('rm -r ', output_wd, '/by-id'),    ignore.stderr = TRUE)
  system(paste0('rm -r ', output_wd, '/met'),      ignore.stderr = TRUE)
  system(paste0('rm -r ', output_wd, '/particles'), ignore.stderr = TRUE)
}
for (d in c('by-id', 'particles', 'footprints')) {
  d <- file.path(output_wd, d)
  if (!file.exists(d)) dir.create(d, recursive = TRUE)
}

# ── Run trajectory simulations ─────────────────────────────────────────────────
stilt_apply(FUN = simulation_step,
            simulation_id      = simulation_id,
            slurm              = slurm,
            slurm_options      = slurm_options,
            n_cores            = n_cores,
            n_nodes            = n_nodes,
            processes_per_node = processes_per_node,
            before_footprint   = list(before_footprint),
            before_trajec      = list(before_trajec),
            lib.loc            = lib.loc,
            capemin     = capemin,
            cmass       = cmass,
            conage      = conage,
            cpack       = cpack,
            delt        = delt,
            dxf         = dxf,
            dyf         = dyf,
            dzf         = dzf,
            efile       = efile,
            emisshrs    = emisshrs,
            frhmax      = frhmax,
            frhs        = frhs,
            frme        = frme,
            frmr        = frmr,
            frts        = frts,
            frvs        = frvs,
            hnf_plume   = hnf_plume,
            hscale      = hscale,
            ichem       = ichem,
            idsp        = idsp,
            initd       = initd,
            k10m        = k10m,
            kagl        = kagl,
            kbls        = kbls,
            kblt        = kblt,
            kdef        = kdef,
            khinp       = khinp,
            khmax       = khmax,
            kmix0       = kmix0,
            kmixd       = kmixd,
            kmsl        = kmsl,
            kpuff       = kpuff,
            krand       = krand,
            krnd        = krnd,
            kspl        = kspl,
            kwet        = kwet,
            kzmix       = kzmix,
            maxdim      = maxdim,
            maxpar      = maxpar,
            met_file_format    = met_file_format,
            met_path           = met_path,
            met_subgrid_buffer = met_subgrid_buffer,
            met_subgrid_enable = met_subgrid_enable,
            met_subgrid_levels = met_subgrid_levels,
            met_file_tres      = met_file_tres,
            mgmin       = mgmin,
            mhrs        = mhrs,
            n_hours     = n_hours,
            n_met_min   = n_met_min,
            nbptyp      = nbptyp,
            ncycl       = ncycl,
            ndump       = ndump,
            ninit       = ninit,
            nstr        = nstr,
            nturb       = nturb,
            numpar      = numpar,
            nver        = nver,
            outdt       = outdt,
            output_wd   = output_wd,
            p10f        = p10f,
            pinbc       = pinbc,
            pinpf       = pinpf,
            poutf       = poutf,
            projection  = projection,
            qcycle      = qcycle,
            r_run_time  = receptors$run_time,
            r_lati      = receptors$lati,
            r_long      = receptors$long,
            r_zagl      = receptors$zagl,
            rhb         = rhb,
            rht         = rht,
            rm_dat      = rm_dat,
            run_foot    = run_foot,
            run_trajec  = run_trajec,
            smooth_factor = smooth_factor,
            splitf      = splitf,
            stilt_wd    = stilt_wd,
            time_integrate = time_integrate,
            timeout     = timeout,
            tkerd       = tkerd,
            tkern       = tkern,
            tlfrac      = tlfrac,
            tout        = tout,
            tratio      = tratio,
            tvmix       = tvmix,
            varsiwant   = list(varsiwant),
            veght       = veght,
            vscale      = vscale,
            vscaleu     = vscaleu,
            vscales     = vscales,
            w_option    = w_option,
            wbbh        = wbbh,
            wbwf        = wbwf,
            wbwr        = wbwr,
            wvert       = wvert,
            xmn  = xmn,  xmx  = xmx,
            ymn  = ymn,  ymx  = ymx,
            xres = xres, yres = yres,
            z_top          = z_top,
            zicontroltf    = zicontroltf,
            ziscale        = ziscale,
            horcoruverr    = horcoruverr,
            siguverr       = siguverr,
            tluverr        = tluverr,
            zcoruverr      = zcoruverr,
            horcorzierr    = horcorzierr,
            sigzierr       = sigzierr,
            tlzierr        = tlzierr)
