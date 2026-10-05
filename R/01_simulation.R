# =============================================================================
# 01_simulation.R — building intuition by SIMULATING meta-analyses
# Introduction to meta-analysis in R · Galilean School
#
# Adapted from the "metasimulation" workshop by Filippo Gambarota & Gianmarco
# Altoè (University of Padova): https://stat-teaching.github.io/metasimulation/
# The sim_studies() function below is a simplified version of theirs.
# =============================================================================

library(metafor)
set.seed(2026)

# -----------------------------------------------------------------------------
# 1. One study = one noisy estimate of the true effect
# -----------------------------------------------------------------------------
# Simulate k two-group studies (treatment vs control), each with n participants
# per group. The outcome is standardised (SD = 1), so the mean difference is
# already on the "Cohen's d" scale.
#   es   = true (average) effect
#   tau2 = between-study variance (0 = equal-effects world)

sim_studies <- function(k, es, tau2 = 0, n) {
  n <- rep_len(n, k)
  delta <- rnorm(k, 0, sqrt(tau2))              # study-specific deviations
  yi <- vi <- numeric(k)
  for (i in 1:k) {
    ctr <- rnorm(n[i], 0, 1)                    # control group
    trt <- rnorm(n[i], es + delta[i], 1)        # treatment group
    yi[i] <- mean(trt) - mean(ctr)              # observed effect
    vi[i] <- var(trt) / n[i] + var(ctr) / n[i]  # its sampling variance
  }
  data.frame(study = paste("Study", 1:k), n = n, yi = yi, vi = vi)
}

# Small vs large studies, same true effect (0.5) ------------------------------
small <- sim_studies(k = 10, es = 0.5, n = 20)
large <- sim_studies(k = 10, es = 0.5, n = 300)

par(mfrow = c(1, 2))
forest(small$yi, small$vi, slab = small$study, refline = 0.5, header = "n = 20 per group", xlim = c(-3, 4), at = seq(-1, 2, 0.5))
forest(large$yi, large$vi, slab = large$study, refline = 0.5, header = "n = 300 per group", xlim = c(-3, 4), at = seq(-1, 2, 0.5))
par(mfrow = c(1, 1))
# -> Same true effect, but small studies scatter a lot more (sampling error).

# -----------------------------------------------------------------------------
# 2. Combining studies: inverse-variance weighting
# -----------------------------------------------------------------------------
dat <- sim_studies(k = 15, es = 0.5, n = rpois(15, 30) + 10)

mean(dat$yi)                       # naive mean: every study counts the same
wi <- 1 / dat$vi                   # precise studies get MORE weight
sum(wi * dat$yi) / sum(wi)         # weighted mean = equal-effects estimate

rma(yi, vi, data = dat, method = "EE")   # same number, plus SE, CI, Q, I^2

# -----------------------------------------------------------------------------
# 3. Equal-effects vs random-effects
# -----------------------------------------------------------------------------
# Now the true effects differ across studies: mu = 0.5, tau = 0.3
dat_re <- sim_studies(k = 30, es = 0.5, tau2 = 0.3^2, n = rpois(30, 30) + 10)

fit_ee <- rma(yi, vi, data = dat_re, method = "EE")
fit_re <- rma(yi, vi, data = dat_re, method = "REML")

fit_ee
fit_re
# -> Similar point estimate, but the RE model has a WIDER CI because it also
#    accounts for the between-study variance (tau^2).

# The weights become more balanced under RE:
round(cbind(EE = weights(fit_ee), RE = weights(fit_re)), 1)

# Prediction interval: where we expect the TRUE effect of a NEW study to be
predict(fit_re)

# -----------------------------------------------------------------------------
# 4. Your turn
# -----------------------------------------------------------------------------
# a) Re-run section 3 with tau2 = 0. What happens to tau^2, I^2 and the PI?
# b) Re-run with tau2 = 0.5^2. What happens?
# c) Re-run with tau2 = 0.3^2 but n = 1000 per group. Does tau^2 change?
#    Does I^2 change? (Hint: I^2 is NOT an absolute measure of heterogeneity!)
