# =============================================================================
# 02_psilodep_live.R — LIVE-CODING SCRIPT (fill in the ___ together)
# Introduction to meta-analysis in R · Galilean School
#
# Is psilocybin-assisted therapy more effective than control conditions in
# reducing depressive symptoms?
#
# Data:     data-depression-psiloctr (Metapsy / Sypres Collaboration)
#           https://github.com/metapsy-project/data-depression-psiloctr
# Pipeline: https://sypres.io/docs/datasets/psilodep-meta-analysis/
# Book:     Harrer et al. (2021) Doing Meta-Analysis with R — https://doing-meta.guide/
#
# Solutions: R/02_psilodep_solutions.R
# =============================================================================

library(metafor)
library(dplyr)

# -----------------------------------------------------------------------------
# STEP 1 — Import & explore the data
# -----------------------------------------------------------------------------
# The file uses ";" as separator and "," as decimal mark
dat <- ___("data/psilodep.csv", na.strings = "NA")

dim(dat)
n_distinct(dat$study)
table(dat$study)

# Why so many rows per study? Explore these columns:
table(dat$outcome_type)
table(dat$instrument)
table(dat$condition_arm2)

dat |>
  filter(study == "Raison 2023") |>
  select(study, outcome_type, instrument, time_weeks,
         primary_instrument, primary_timepoint)

# -----------------------------------------------------------------------------
# STEP 2 — Select ONE effect size per study
# -----------------------------------------------------------------------------
# Keep: primary instrument, primary time point, means & SDs (msd / imsd),
# pre-crossover data, Goodwin 2022 25 mg vs 1 mg only,
# and exclude Carhart-Harris 2021 and Krempien 2023.
dat_main <- dat |>
  filter(
    primary_instrument == ___,
    primary_timepoint == ___,
    outcome_type %in% c(___, ___),
    is.na(post_crossover) | post_crossover == 0,
    !study %in% c("Carhart-Harris 2021", "Krempien 2023"),
    !(multi_arm1 %in% "10 mg" | multi_arm2 %in% "10 mg")
  )

nrow(dat_main)
n_distinct(dat_main$study)   # should be equal to nrow()!

# -----------------------------------------------------------------------------
# STEP 3 — Compute effect sizes (Hedges' g) with escalc()
# -----------------------------------------------------------------------------
dat_main <- escalc(
  measure = ___,
  m1i = mean_arm1, sd1i = sd_arm1, n1i = n_arm1,
  m2i = ___,       sd2i = ___,     n2i = ___,
  data = dat_main,
  slab = study
)

dat_main |> select(study, yi, vi)

# Sanity check against the pre-computed Metapsy values (.g)
round(cbind(ours = dat_main$yi, metapsy = dat_main$.g), 3)

# -----------------------------------------------------------------------------
# STEP 4 — Fit the random-effects model
# -----------------------------------------------------------------------------
res <- rma(___, ___, data = dat_main, method = ___, test = ___)
res

confint(res)    # CIs for tau^2 and I^2
predict(res)    # prediction interval

# -----------------------------------------------------------------------------
# STEP 5 — Forest plot
# -----------------------------------------------------------------------------
forest(___,
       header = c("Study", "Hedges' g [95% CI]"),
       xlab = "Hedges' g (negative favours psilocybin)",
       addpred = TRUE,
       order = "obs",
       shade = TRUE)

# -----------------------------------------------------------------------------
# STEP 6 — Sensitivity analyses
# -----------------------------------------------------------------------------
# 6a) Equal-effects model
res_ee <- rma(yi, vi, data = dat_main, method = ___)
res_ee

# 6b) Influence diagnostics
inf <- ___(res)
inf
plot(inf)

# 6c) Leave-one-out
___(res)

# 6d) Without the most influential study
res_noout <- rma(yi, vi, data = dat_main, method = "REML", test = "knha",
                 subset = study != ___)
res_noout

# -----------------------------------------------------------------------------
# STEP 7 — Subgroup analysis & meta-regression
# -----------------------------------------------------------------------------
table(dat_main$rating)

res_rating <- rma(yi, vi, mods = ~ ___, data = dat_main,
                  method = "REML", test = "knha")
res_rating

# Pooled effect per subgroup
rma(yi, vi, mods = ~ 0 + rating, data = dat_main, method = "REML", test = "knha")

# Continuous moderator: publication year (centred at 2016)
res_year <- rma(yi, vi, mods = ~ I(year - 2016), data = dat_main,
                method = "REML", test = "knha")
res_year
regplot(res_year, xlab = "Years since 2016", ylab = "Hedges' g")

# -----------------------------------------------------------------------------
# STEP 8 — Small-study effects / publication bias
# -----------------------------------------------------------------------------
___(res, xlab = "Hedges' g")

funnel(res, level = c(90, 95, 99), shade = c("white", "gray55", "gray75"),
       refline = 0, legend = TRUE)

___(res, model = "lm")         # Egger's regression test

trimfill(rma(yi, vi, data = dat_main, method = "REML"))

# =============================================================================
# YOUR TURN
# =============================================================================

# Exercise 1 — Re-run the main model WITHOUT studies at high risk of bias
#              (column 'rob'). Does the conclusion change?


# Exercise 2 — Do parallel and crossover trials (column 'design') differ?


# Exercise 3 — Response rates as a Risk Ratio.
#   - same selection as Step 2, but outcome_type == "response"
#   - escalc(measure = "RR", ai = event_arm1, n1i = totaln_arm1,
#                            ci = event_arm2, n2i = totaln_arm2, ...)
#   - pooled results are on the LOG scale: use predict(res, transf = exp)


# BONUS — Use ALL post-dose time points (time_days > 0) and fit a three-level
#         model with rma.mv(..., random = ~ 1 | study / es_id).
#         See https://doing-meta.guide/multilevel-ma.html
