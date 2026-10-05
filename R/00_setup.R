# =============================================================================
# 00_setup.R — run this BEFORE the workshop
# Introduction to meta-analysis in R · Galilean School
# =============================================================================

# Packages we need --------------------------------------------------------------
pkgs <- c(
  "metafor",      # the main meta-analysis package (Viechtbauer, 2010)
  "dplyr",        # data wrangling
  "clubSandwich"  # (optional) cluster-robust tests for the bonus exercise
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

# Check ---------------------------------------------------------------------------
library(metafor)
library(dplyr)

packageVersion("metafor")  # >= 4.0 is needed (vcalc() is used in the bonus part)

# Working directory ---------------------------------------------------------------
# Open the project file 'gali-meta-analysis.Rproj' in RStudio, OR set the working
# directory to the workshop folder (the one that contains data/ and R/), e.g.:
# setwd("~/Desktop/gali")

file.exists("data/psilodep.csv")  # should print TRUE

# If you see TRUE and no errors above, you are ready!
