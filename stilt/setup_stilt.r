#!/usr/bin/env Rscript
# WRF-STILT setup script for IDN_BB_2015
# Run once to install the uataq/stilt package and initialise the project.
#
# Usage:
#   cd /home/igrk/WRF-GRK
#   Rscript stilt/setup_stilt.r
#
# After this completes successfully, the stilt/ directory will contain:
#   exe/hymodelc     - compiled STILT Fortran executable
#   r/               - STILT R source functions
#   run_stilt.r      - pre-configured run script (edit before running)

project_dir <- '/home/igrk/WRF-GRK/stilt'

# ── 1. Install uataq/stilt from GitHub ────────────────────────────────────────
if (!requireNamespace('uataq', quietly = TRUE)) {
  message('Installing uataq/uataq ...')
  if (!requireNamespace('devtools', quietly = TRUE))
    install.packages('devtools', repos = 'https://cloud.r-project.org')
  devtools::install_github('uataq/uataq')
} else {
  message('uataq already installed.')
}

library(uataq)

# ── 2. Initialise STILT project ───────────────────────────────────────────────
# stilt_init() clones uataq/stilt into project_dir and compiles permute.f90.
message('Initialising STILT project in: ', project_dir)
stilt_init(project_dir)

message('\nSetup complete.')
message('Next step: configure stilt/run_stilt.r, then run:')
message('  Rscript stilt/run_stilt.r')
