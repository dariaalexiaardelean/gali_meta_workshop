# =============================================================================
# 02_psilodep_solutions.R — FULL SOLUTIONS
# Introduction to meta-analysis in R · Galilean School
#
# Is psilocybin-assisted therapy more effective than control conditions in
# reducing depressive symptoms?
#
# Data:     data-depression-psiloctr (Metapsy / Sypres Collaboration)
#           https://github.com/metapsy-project/data-depression-psiloctr
# Pipeline: follows https://sypres.io/docs/datasets/psilodep-meta-analysis/
#           (they use metapsyTools; here we reproduce it step by step with metafor)
# Book:     Harrer et al. (2021) Doing Meta-Analysis with R — https://doing-meta.guide/
# =============================================================================

library(metafor)
library(dplyr)

# -----------------------------------------------------------------------------
# STEP 1 — Import & explore the data
# -----------------------------------------------------------------------------
# The file uses ";" as separator and "," as decimal mark -> read.csv2()
dat <- read.csv2("data/psilodep.csv", na.strings = "NA")

dim(dat)                       # 221 rows, 74 columns
n_distinct(dat$study)          # ...but only 15 studies!
table(dat$study)

# Why so many rows? One row = one comparison x instrument x time point x outcome
table(dat$outcome_type)        # msd, imsd, change, response, remission
table(dat$instrument)          # MADRS, HAM-D, BDI, QIDS-SR, ...
table(dat$condition_arm2)      # control types: placebo, niacin, waitlist, low-dose psilocybin...

# Look at one study in detail
dat |>
  filter(study == "Raison 2023") |>
  select(study, outcome_type, instrument, time_weeks,
         primary_instrument, primary_timepoint)

# -----------------------------------------------------------------------------
# STEP 2 — Select ONE effect size per study (Sypres inclusion rules)
# -----------------------------------------------------------------------------
# Standard meta-analysis assumes INDEPENDENT effect sizes. We therefore keep:
#   - the primary instrument at the primary time point
#   - continuous outcomes reported as means & SDs (msd / imsd)
#   - pre-crossover data only (crossover trials)
#   - Goodwin 2022: only 25 mg vs 1 mg (drop comparisons involving the 10 mg arm)
#   - exclude Carhart-Harris 2021 and Krempien 2023 (as in the Sypres pipeline)
dat_main <- dat |>
  filter(
    primary_instrument == 1,
    primary_timepoint == 1,
    outcome_type %in% c("msd", "imsd"),
    is.na(post_crossover) | post_crossover == 0,
    !study %in% c("Carhart-Harris 2021", "Krempien 2023"),
    !(multi_arm1 %in% "10 mg" | multi_arm2 %in% "10 mg")
  )

nrow(dat_main)                 # 12 rows ...
n_distinct(dat_main$study)     # ... from 12 studies: one per study

dat_main |>
  select(study, condition_arm2, instrument,
         n_arm1, mean_arm1, sd_arm1, n_arm2, mean_arm2, sd_arm2)

# -----------------------------------------------------------------------------
# STEP 3 — Compute effect sizes with escalc()
# -----------------------------------------------------------------------------
# measure = "SMD" -> bias-corrected standardised mean difference = Hedges' g
# Arm 1 = psilocybin, arm 2 = control. Lower scores = fewer symptoms, so
# NEGATIVE g = psilocybin better.
dat_main <- escalc(
  measure = "SMD",
  m1i = mean_arm1, sd1i = sd_arm1, n1i = n_arm1,
  m2i = mean_arm2, sd2i = sd_arm2, n2i = n_arm2,
  data = dat_main,
  slab = study                 # study labels for the plots
)

# yi = effect size (g), vi = sampling variance
dat_main |> select(study, yi, vi)

# Sanity check: compare with the effect sizes pre-computed by Metapsy (.g)
round(cbind(ours = dat_main$yi, metapsy = dat_main$.g), 3)

# -----------------------------------------------------------------------------
# STEP 4 — Fit the random-effects model
# -----------------------------------------------------------------------------
# method = "REML" -> estimator for tau^2
# test = "knha"   -> Knapp-Hartung adjustment (t-distribution, better with few studies)
res <- rma(yi, vi, data = dat_main, method = "REML", test = "knha")
res

# How to read it:
# - estimate (mu): the AVERAGE true effect across studies, with SE, CI and p-value
# - tau^2 / tau:   between-study variance / SD of the TRUE effects
# - I^2:           % of the observed variability that is due to true heterogeneity
# - Q-test:        is there more variability than expected by chance alone?

# Confidence intervals for the heterogeneity statistics
confint(res)

# Prediction interval: range of TRUE effects to expect in a NEW, comparable study
predict(res)

# -----------------------------------------------------------------------------
# STEP 5 — Forest plot
# -----------------------------------------------------------------------------
forest(res,
       header = c("Study", "Hedges' g [95% CI]"),
       xlab = "Hedges' g (negative favours psilocybin)",
       addpred = TRUE,          # show the prediction interval
       order = "obs",           # sort by effect size
       shade = TRUE)

# -----------------------------------------------------------------------------
# STEP 6 — Sensitivity analyses
# -----------------------------------------------------------------------------
# 6a) Equal-effects model (assumes ONE true effect for all studies)
res_ee <- rma(yi, vi, data = dat_main, method = "EE")
res_ee

# 6b) Influence diagnostics: is any study driving the result?
inf <- influence(res)
inf
plot(inf)

# 6c) Leave-one-out analysis
leave1out(res)

# 6d) Re-fit without the most influential study (Davis 2021)
res_nodavis <- rma(yi, vi, data = dat_main, method = "REML", test = "knha",
                   subset = study != "Davis 2021")
res_nodavis

# 6e) Use the pre-computed Metapsy effect sizes -> reproduces Sypres exactly
#     (g = -0.90, tau^2 = 0.117). Small differences from our g come from the
#     formula used for the sampling variance.
rma(.g, sei = .g_se, data = dat_main, method = "REML", test = "knha")

# -----------------------------------------------------------------------------
# STEP 7 — Subgroup analysis & meta-regression
# -----------------------------------------------------------------------------
# Does the effect differ between clinician-rated and self-report instruments?
table(dat_main$rating)

# Subgroup analysis = meta-regression with a categorical moderator
res_rating <- rma(yi, vi, mods = ~ rating, data = dat_main,
                  method = "REML", test = "knha")
res_rating
# -> 'Test of Moderators' tells whether the subgroups differ

# Pooled effect in each subgroup (no intercept parametrisation)
rma(yi, vi, mods = ~ 0 + rating, data = dat_main, method = "REML", test = "knha")

# Continuous moderator: publication year (centred, so the intercept is the
# predicted effect for a study published in 2016, the earliest in our data)
range(dat_main$year)
res_year <- rma(yi, vi, mods = ~ I(year - 2016), data = dat_main,
                method = "REML", test = "knha")
res_year
regplot(res_year, xlab = "Years since 2016", ylab = "Hedges' g")

# -----------------------------------------------------------------------------
# STEP 8 — Small-study effects / publication bias
# -----------------------------------------------------------------------------
funnel(res, xlab = "Hedges' g")

# Contour-enhanced funnel plot (centred at 0: shaded = regions of significance)
funnel(res, level = c(90, 95, 99), shade = c("white", "gray55", "gray75"),
       refline = 0, legend = TRUE)

# Egger's regression test (classic version, as used by Sypres: p ~ 0.21)
regtest(res, model = "lm")

# Trim-and-fill (needs a model without Knapp-Hartung)
trimfill(rma(yi, vi, data = dat_main, method = "REML"))

# NB: with k = 12 studies these tests have very little power!

# =============================================================================
# EXERCISES
# =============================================================================

# -----------------------------------------------------------------------------
# Exercise 1 — Risk of bias
# Re-run the main model excluding studies at HIGH risk of bias (column 'rob').
# -----------------------------------------------------------------------------
table(dat_main$rob)
res_rob <- rma(yi, vi, data = dat_main, method = "REML", test = "knha",
               subset = rob != "High")
res_rob

# -----------------------------------------------------------------------------
# Exercise 2 — Study design
# Do parallel-group and crossover trials ('design') give different results?
# -----------------------------------------------------------------------------
table(dat_main$design)
res_design <- rma(yi, vi, mods = ~ design, data = dat_main,
                  method = "REML", test = "knha")
res_design

# -----------------------------------------------------------------------------
# Exercise 3 — A dichotomous outcome: response rates (Risk Ratio)
# Same selection rules, but outcome_type == "response".
# Columns: event_arm1, totaln_arm1, event_arm2, totaln_arm2.
# Hint: escalc(measure = "RR", ai = , n1i = , ci = , n2i = ) — log scale!
# -----------------------------------------------------------------------------
dat_resp <- dat |>
  filter(
    outcome_type == "response",
    primary_instrument == 1,
    primary_timepoint == 1,
    is.na(post_crossover) | post_crossover == 0,
    !study %in% c("Carhart-Harris 2021", "Krempien 2023"),
    !(multi_arm1 %in% "10 mg" | multi_arm2 %in% "10 mg")
  )

dat_resp <- escalc(measure = "RR",
                   ai = event_arm1, n1i = totaln_arm1,
                   ci = event_arm2, n2i = totaln_arm2,
                   data = dat_resp, slab = study)

res_resp <- rma(yi, vi, data = dat_resp, method = "REML", test = "knha")
res_resp
predict(res_resp, transf = exp)   # back-transform: RR ~ 2.8
forest(res_resp, atransf = exp, refline = 0,
       xlab = "Risk Ratio (log scale; > 1 favours psilocybin)")

# Try the same for remission (outcome_type == "remission"): RR ~ 4.2

# -----------------------------------------------------------------------------
# BONUS — Using ALL time points: a three-level (CHE) model
# Studies contribute several effect sizes (different weeks) -> NOT independent.
# A multilevel model nests effect sizes within studies; vcalc() adds an
# assumed correlation (rho = 0.6) between sampling errors of the same study.
# See: https://doing-meta.guide/multilevel-ma.html
# -----------------------------------------------------------------------------
dat_time <- dat |>
  filter(
    outcome_type %in% c("msd", "imsd"),
    primary_instrument == 1,
    time_days > 0,                               # all post-dose time points
    is.na(post_crossover) | post_crossover == 0,
    !study %in% c("Carhart-Harris 2021", "Krempien 2023"),
    !(multi_arm1 %in% "10 mg" | multi_arm2 %in% "10 mg")
  )

dat_time <- escalc(measure = "SMD",
                   m1i = mean_arm1, sd1i = sd_arm1, n1i = n_arm1,
                   m2i = mean_arm2, sd2i = sd_arm2, n2i = n_arm2,
                   data = dat_time, slab = study)
dat_time$es_id <- seq_len(nrow(dat_time))       # one id per effect size
table(dat_time$study)                           # 39 effect sizes, 12 studies

V <- vcalc(vi, cluster = study, obs = es_id, rho = 0.6, data = dat_time)

res_che <- rma.mv(yi, V, random = ~ 1 | study / es_id,
                  data = dat_time, test = "t")
res_che
robust(res_che, cluster = study, clubSandwich = TRUE)  # cluster-robust inference

# Does the effect fade over time?
res_time <- rma.mv(yi, V, mods = ~ time_weeks, random = ~ 1 | study / es_id,
                   data = dat_time, test = "t")
res_time
