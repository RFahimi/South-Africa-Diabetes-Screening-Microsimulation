# =============================================================================
#  COST-EFFECTIVENESS OF FAST-AND-FRUGAL-TREE (FFT) AND MACHINE-LEARNING (ML)
#  BASED TYPE 2 DIABETES SCREENING IN SOUTH AFRICA
#
#  Annotated analysis script (requires separately obtained input files).
#
#  Manuscript : Cost-Effectiveness of Fast-and-Frugal Tree-Based and Machine
#               Learning-Based Strategies for Type 2 Diabetes Screening in
#               South Africa (Frontiers in Clinical Diabetes and Healthcare)
#  Written    : 2026-07-19
#
#  ---------------------------------------------------------------------------
#  WHAT THIS SCRIPT DOES
#  ---------------------------------------------------------------------------
#  Runs the entire economic evaluation end to end:
#     1. Vectorised individual-level microsimulation engine (with common
#        random numbers)
#     2. Probabilistic parameter distributions
#     3. Internal self-tests (engine invariants + sampler moments)
#     4. Base case, efficient frontier, and first-order stability check
#     5. One-way deterministic sensitivity analysis (tornado) + scenarios
#     6. Probabilistic sensitivity analysis (CE plane, CEAC/CEAF, credible
#        intervals, probability cost-effective)
#     7. Expected value of perfect information (EVPI)
#     8. Supplementary Table S1 (PSA distributions)
#     9. Figure 1 (cost-effectiveness plane)
#
#  It supersedes and consolidates:
#     R/engine_vectorized.R, R/distributions.R, R/run_basecase.R, R/run_dsa.R,
#     R/run_psa.R, R/run_evpi.R, R/build_s1.R, tests/test_engine.R,
#     tests/test_distributions.R, Figure1_CE_plane_corrected_final.R
#  (The earlier per-person script Objective5_microsimulation_final_revision.R
#   is superseded entirely - see the note on common random numbers below.)
#
#  ---------------------------------------------------------------------------
#  THE KEY METHODOLOGICAL POINT: COMMON RANDOM NUMBERS (CRN)
#  ---------------------------------------------------------------------------
#  Incremental QALY differences between these strategies are small (~0.03-0.04
#  per person). In the original implementation each strategy was simulated on
#  an INDEPENDENT random-number stream, so the between-strategy Monte-Carlo
#  noise (SD ~0.033 across seeds) was as large as the signal being measured.
#  The identity of the "efficient" strategy was therefore unstable: repeating
#  the analysis under different seeds returned FFT screening, FFT case-finding
#  or ML targeting in different runs.
#
#  This script evaluates EVERY strategy on the SAME simulated individuals using
#  the SAME pre-generated random draws (common random numbers - a standard
#  variance-reduction technique). A transition fires when a FIXED uniform draw
#  falls below a PARAMETER-DEPENDENT threshold, so when a strategy (or a PSA
#  draw) changes a probability, only the threshold moves - not the noise.
#
#  Effect: the incremental-QALY SD across seeds falls from ~0.033 to ~0.003,
#  and FFT CASE-FINDING is the sole strategy on the efficient frontier in
#  30/30 random seeds. This is why the conclusions differ from the original
#  submission. Absolute costs and QALYs are essentially unchanged; only the
#  noise-sensitive increments moved.
#
#  ---------------------------------------------------------------------------
#  HOW TO RUN
#  ---------------------------------------------------------------------------
#     Rscript Objective5_full_analysis_v2.R
#  from the project root (the directory containing df_t2dm.csv).
#  Requires: R >= 4.0, tidyverse (dplyr, ggplot2, purrr, readr, tidyr), scales.
#  Runtime : ~2 minutes (the 1,000-draw PSA dominates).
#  All randomness is seeded, so repeated runs are bit-identical.
#
#  INPUTS  : df_t2dm.csv, ml_predictions_sa_fixed.csv, fft_predictions_sa.csv
#  OUTPUTS : output/*.csv, output/figures/fig4-fig7.{png,pdf},
#            Figure 1_corrected.{png,pdf},
#            output/manuscript_edits/supplementary_table_S1_distributions.{csv,md}
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

dir.create("output",                     showWarnings = FALSE)
dir.create("output/figures",             showWarnings = FALSE, recursive = TRUE)
dir.create("output/manuscript_edits",    showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()


# =============================================================================
# SECTION 1.  MICROSIMULATION ENGINE
# =============================================================================
# A four-state, annual-cycle, individual-level model. Everything is vectorised
# over the cohort (one pass per YEAR, not per PERSON), which is what makes a
# 1,000-draw PSA feasible: a single cohort-run is ~0.02 s rather than minutes.

# ---- Health states (absorbing: DEAD) ----------------------------------------
STATE_NO_DM   <- 0L   # no diabetes
STATE_UNDX_DM <- 1L   # undiagnosed diabetes
STATE_DX_DM   <- 2L   # diagnosed diabetes (assumed to be receiving care)
STATE_DEAD    <- 3L

STATE_NAMES <- c("No diabetes", "Undiagnosed diabetes", "Diagnosed diabetes", "Death")
STRATEGIES  <- c("baseline", "ml_targeted", "fft_casefinding",
                 "fft_screening", "fft_ruleout")

# ---- Base-case parameters ---------------------------------------------------
# SINGLE SOURCE OF TRUTH. The DSA perturbs one field at a time; the PSA draws
# many fields jointly. Every value is a calibrated model assumption unless
# noted - this is why the uncertainty analysis matters so much (see Section 5-7)
# and why the manuscript frames these as assumptions, not empirical estimates.
make_base_params <- function() {
  list(
    time_horizon      = 20L,    # years
    discount_rate     = 0.03,   # per annum, applied to costs and QALYs
    incidence_target  = 0.012,  # calibrated mean annual incidence (SADHS/IDF)
    baseline_dx_share = 0.55,   # share of prevalent diabetes already diagnosed

    # Annual probability of being detected. 'baseline' applies to everyone not
    # targeted by the active strategy. These are TARGETING-INTENSITY parameters,
    # NOT the sensitivity/specificity of the FFT rules or the ML classifier.
    p_detect_baseline      = 0.08,
    p_detect_ml_high       = 0.25,
    p_detect_casefind_high = 0.30,   # highest yield -> drives the main result
    p_detect_screen_high   = 0.22,
    p_detect_ruleout_low   = 0.02,   # rule-out REDUCES detection for low-risk

    # Assumed relative risk of diabetes onset among targeted individuals.
    # Unvalidated: the PSA lets these reach 1.0 (no effect) in the upper tail.
    rr_onset_ml       = 0.80,
    rr_onset_casefind = 0.82,
    rr_onset_screen   = 0.85,

    # Health-state utilities. The undiagnosed -> diagnosed GAP (0.78 -> 0.84) is
    # the mechanism through which detection generates QALYs, so it is a primary
    # driver of every result in this analysis.
    u_nodm = 0.95, u_undx = 0.78, u_dx = 0.84, u_dead = 0.00,

    # Annual health-state costs (2021 USD, payer perspective).
    c_nodm = 0, c_undx = 40, c_dx = 90, c_dead = 0,

    # Annual per-person programme costs by strategy.
    pc_baseline = 0, pc_ml = 12, pc_casefind = 10, pc_screen = 8, pc_ruleout = 5,

    # All-cause mortality: logistic in age, with additive log-odds excess for
    # undiagnosed and diagnosed diabetes.
    mort_intercept = -9, mort_slope = 0.075,
    mort_exc_undx  = 0.25, mort_exc_dx = 0.15
  )
}

# ---- Cohort construction ----------------------------------------------------
# Loads the South African DHS extract and joins the pre-generated FFT flags and
# ML risk scores by individual id.
#
# NOTE ON ml_risk: these are IN-SAMPLE predictions from a model fitted to this
# same survey, so the ML arm's targeting is optimistic relative to prospective
# deployment (disclosed in the manuscript Limitations). The score is strongly
# zero-inflated (60.3% exact zeros), which is why modelled incident diabetes
# concentrates in the minority with non-zero predicted risk.
prepare_cohort <- function(data_path = ".") {
  df <- read_csv(file.path(data_path, "df_t2dm.csv"), show_col_types = FALSE) %>%
    filter(countryname == "South Africa") %>%
    mutate(id = row_number()) %>%
    select(-countryname)

  ml_pred  <- read_csv(file.path(data_path, "ml_predictions_sa_fixed.csv"),
                       show_col_types = FALSE)
  fft_pred <- read_csv(file.path(data_path, "fft_predictions_sa.csv"),
                       show_col_types = FALSE)

  df <- df %>%
    left_join(ml_pred,  by = "id") %>%
    left_join(fft_pred, by = "id") %>%
    mutate(
      has_diabetes = as.integer(has_diabetes),
      # Age categories are mapped to representative ages for the mortality model.
      age_num = case_when(
        age_cate == "Young Adult"              ~ 27,
        age_cate == "Adult"                    ~ 42,
        age_cate == "Middle-Aged /Older Adult" ~ 60,
        TRUE                                   ~ 42
      )
    )

  stopifnot(!any(is.na(df$ml_risk)),      !any(is.na(df$fft_casefind)),
            !any(is.na(df$fft_screen)),   !any(is.na(df$fft_ruleout)))

  # Individual onset risk is expressed RELATIVE to the mean among people without
  # diabetes, so it can later be rescaled to hit any incidence target.
  mean_ml_non_dm <- mean(df$ml_risk[df$has_diabetes == 0L], na.rm = TRUE)
  df$ml_ratio    <- df$ml_risk / mean_ml_non_dm

  # "High risk" under the ML strategy = top quintile of predicted risk.
  df$ml_highrisk <- as.integer(df$ml_risk >= quantile(df$ml_risk, 0.80, na.rm = TRUE))
  df
}

# ---- Common random numbers --------------------------------------------------
# Generated ONCE and reused across all strategies AND all PSA draws. This is the
# heart of the method (see header). max_horizon must cover the longest horizon
# tested in the DSA (30 years) so horizon scenarios reuse the same draws.
make_rng <- function(n, max_horizon = 30L, seed = 1234L) {
  set.seed(seed)
  list(
    U_init  = runif(n),                                        # initial dx status
    U_death = matrix(runif(n * max_horizon), nrow = n),        # annual mortality
    U_trans = matrix(runif(n * max_horizon), nrow = n)         # onset OR detection
  )
}

# ---- Single-strategy cohort simulation --------------------------------------
# Returns one row per person. Semantics deliberately mirror the original
# per-person model: accrual happens at the START of the cycle (so a person who
# dies during a year still accrues that year), mortality is resolved before
# disease transitions, and onset/detection are mutually exclusive within a cycle.
run_cohort <- function(df, strategy, params, rng) {
  n <- nrow(df); H <- params$time_horizon; r <- params$discount_rate

  # -- Per-person annual onset probability -----------------------------------
  # Scaled to the incidence target, then bounded to keep probabilities sane.
  onset_raw <- pmin(pmax(df$ml_ratio * params$incidence_target, 0.0001), 0.20)

  # Targeted strategies are assumed to reduce onset among those they target.
  rr <- rep(1, n)
  if (strategy == "ml_targeted")     rr[df$ml_highrisk == 1L]  <- params$rr_onset_ml
  if (strategy == "fft_casefinding") rr[df$fft_casefind == 1L] <- params$rr_onset_casefind
  if (strategy == "fft_screening")   rr[df$fft_screen == 1L]   <- params$rr_onset_screen
  p_onset_vec <- pmin(pmax(onset_raw * rr, 0.0001), 0.20)

  # -- Per-person annual detection probability -------------------------------
  # Everyone starts at the background rate; the strategy raises it for those it
  # targets (or LOWERS it, for rule-out - which is why rule-out can do harm).
  p_detect_vec <- rep(params$p_detect_baseline, n)
  if (strategy == "ml_targeted")     p_detect_vec[df$ml_highrisk == 1L]  <- params$p_detect_ml_high
  if (strategy == "fft_casefinding") p_detect_vec[df$fft_casefind == 1L] <- params$p_detect_casefind_high
  if (strategy == "fft_screening")   p_detect_vec[df$fft_screen == 1L]   <- params$p_detect_screen_high
  if (strategy == "fft_ruleout")     p_detect_vec[df$fft_ruleout == 1L]  <- params$p_detect_ruleout_low

  program_cost <- switch(strategy,
    baseline        = params$pc_baseline,
    ml_targeted     = params$pc_ml,
    fft_casefinding = params$pc_casefind,
    fft_screening   = params$pc_screen,
    fft_ruleout     = params$pc_ruleout)

  # Indexed by state + 1 (states are 0-based).
  state_cost <- c(params$c_nodm, params$c_undx, params$c_dx, params$c_dead)
  state_util <- c(params$u_nodm, params$u_undx, params$u_dx, params$u_dead)

  # -- Initial state (uses the shared CRN draw, so identical across strategies)
  state <- integer(n)                                  # defaults to STATE_NO_DM
  diab  <- df$has_diabetes == 1L
  dx0   <- diab & (rng$U_init < params$baseline_dx_share)
  state[dx0]          <- STATE_DX_DM
  state[diab & !dx0]  <- STATE_UNDX_DM

  age <- df$age_num
  cost <- numeric(n); qaly <- numeric(n); ly <- numeric(n)
  years_undx <- numeric(n); incident_dm <- integer(n)
  detected <- as.integer(state == STATE_DX_DM)   # prevalent diagnosed count once

  # -- Annual cycle ----------------------------------------------------------
  for (t in 0:(H - 1L)) {
    dfac  <- 1 / ((1 + r)^t)                     # discount factor for this year
    alive <- state != STATE_DEAD
    s1    <- state + 1L

    # (a) Accrue costs/QALYs/life-years for everyone alive at cycle start.
    cost[alive] <- cost[alive] + (state_cost[s1[alive]] + program_cost) * dfac
    qaly[alive] <- qaly[alive] + state_util[s1[alive]] * dfac
    ly[alive]   <- ly[alive] + dfac
    is_undx <- alive & (state == STATE_UNDX_DM)
    years_undx[is_undx] <- years_undx[is_undx] + 1

    # (b) Mortality (vectorised logistic; clamped to a plausible range).
    lp <- params$mort_intercept + params$mort_slope * age +
          (state == STATE_UNDX_DM) * params$mort_exc_undx +
          (state == STATE_DX_DM)   * params$mort_exc_dx
    p_death <- pmin(pmax(plogis(lp), 0.0001), 0.40)
    died <- alive & (rng$U_death[, t + 1L] < p_death)
    state[died] <- STATE_DEAD

    # (c) Disease transitions among survivors.
    #     IMPORTANT: eligibility is fixed BEFORE any state update, so a person
    #     who develops diabetes this cycle cannot also be detected this cycle.
    surv      <- alive & !died
    was_no_dm <- surv & (state == STATE_NO_DM)
    was_undx  <- surv & (state == STATE_UNDX_DM)

    onset  <- was_no_dm & (rng$U_trans[, t + 1L] < p_onset_vec)
    detect <- was_undx  & (rng$U_trans[, t + 1L] < p_detect_vec)

    state[onset]        <- STATE_UNDX_DM
    incident_dm[onset]  <- incident_dm[onset] + 1L
    state[detect]       <- STATE_DX_DM
    detected[detect]    <- detected[detect] + 1L

    age[alive] <- age[alive] + 1
  }

  tibble(id = df$id, strategy = strategy,
         total_cost = cost, total_qaly = qaly, total_ly = ly,
         years_undx = years_undx, incident_dm = incident_dm,
         detected_cases = detected, final_state = state)
}

# ---- Run all strategies on the SAME rng (this is the CRN pairing) ------------
simulate_all <- function(df, params, rng, strategies = STRATEGIES) {
  map_dfr(strategies, function(s) {
    res <- run_cohort(df, s, params, rng)
    tibble(strategy            = s,
           mean_cost           = mean(res$total_cost),
           mean_qaly           = mean(res$total_qaly),
           mean_ly             = mean(res$total_ly),
           mean_years_undx     = mean(res$years_undx),
           incident_dm_rate    = mean(res$incident_dm),
           detected_cases_mean = mean(res$detected_cases),
           death_rate          = mean(res$final_state == STATE_DEAD))
  })
}

# ---- Incremental results vs the no-screening comparator ---------------------
icer_vs_baseline <- function(summary_tbl) {
  base_row <- summary_tbl %>% filter(strategy == "baseline")
  summary_tbl %>%
    filter(strategy != "baseline") %>%
    mutate(incr_cost = mean_cost - base_row$mean_cost,
           incr_qaly = mean_qaly - base_row$mean_qaly,
           ICER      = incr_cost / incr_qaly) %>%
    select(strategy, mean_cost, mean_qaly, incr_cost, incr_qaly, ICER)
}

# Net monetary benefit. Used instead of the ICER wherever incremental QALYs can
# cross zero (ICERs are undefined/uninterpretable there - e.g. dominated arms).
nmb <- function(mean_qaly, mean_cost, lambda) lambda * mean_qaly - mean_cost

# ---- Efficient frontier -----------------------------------------------------
# Walks the cost-ordered options, repeatedly taking the lowest-ICER extension.
# This removes BOTH strictly dominated options (more cost, fewer QALYs) and
# EXTENDEDLY dominated options (beaten by a mixture of two others) - the latter
# is what excludes FFT screening in the base case.
frontier <- function(sm) {
  pts <- sm %>% select(strategy, mean_cost, mean_qaly) %>%
    arrange(mean_cost, desc(mean_qaly))
  current <- pts %>% slice(1); eff <- current$strategy
  repeat {
    cand <- pts %>% filter(mean_cost > current$mean_cost,
                           mean_qaly > current$mean_qaly)
    if (nrow(cand) == 0) break
    cand <- cand %>% mutate(icer = (mean_cost - current$mean_cost) /
                                   (mean_qaly - current$mean_qaly))
    nxt <- cand %>% arrange(icer) %>% slice(1)
    eff <- c(eff, nxt$strategy); current <- nxt
  }
  eff
}


# =============================================================================
# SECTION 2.  PROBABILISTIC SENSITIVITY ANALYSIS: DISTRIBUTIONS
# =============================================================================
# Samplers are implemented in base R (no extra dependencies).
#
# IMPORTANT CAVEAT, stated plainly in the manuscript: most parameters here are
# calibrated ASSUMPTIONS, not empirically estimated quantities. The PSA
# therefore propagates ASSUMED uncertainty. It cannot manufacture empirical
# confidence that the inputs do not have - which is exactly why the EVPI in
# Section 7 is the more decision-relevant output.

# Beta by method of moments (for probabilities and utilities, bounded 0-1).
beta_from_mean_sd <- function(n, mean, sd) {
  vmax <- mean * (1 - mean)
  v <- min(sd^2, 0.99 * vmax)                  # guard against infeasible SDs
  k <- mean * (1 - mean) / v - 1
  rbeta(n, mean * k, (1 - mean) * k)
}

# Gamma by method of moments (for costs and rates: positive, right-skewed).
gamma_from_mean_sd <- function(n, mean, sd) {
  rgamma(n, shape = (mean / sd)^2, rate = mean / sd^2)
}

# Lognormal specified by its median and an upper quantile. Used for the onset
# relative risks so that RR = 1 (i.e. NO preventive effect) sits at the upper
# 2.5% tail - an honest representation of an unvalidated assumed effect.
lnorm_from_median_upper <- function(n, median, upper, p = 0.975) {
  sigma <- (log(upper) - log(median)) / qnorm(p)
  rlnorm(n, meanlog = log(median), sdlog = sigma)
}

# ---- The PSA specification --------------------------------------------------
# One row per parameter. Doubles as the source for Supplementary Table S1, so
# the table can never drift out of step with what was actually sampled.
psa_spec <- function() {
  tribble(
    ~param,                    ~base,   ~dist,     ~arg1,   ~arg2,  ~basis,
    "p_detect_baseline",       0.08,    "beta",    0.08,    0.024,  "Assumption; CV 30% (unvalidated)",
    "p_detect_ml_high",        0.25,    "beta",    0.25,    0.075,  "Assumption; CV 30% (unvalidated)",
    "p_detect_casefind_high",  0.30,    "beta",    0.30,    0.090,  "Assumption; CV 30% (unvalidated)",
    "p_detect_screen_high",    0.22,    "beta",    0.22,    0.066,  "Assumption; CV 30% (unvalidated)",
    "p_detect_ruleout_low",    0.02,    "beta",    0.02,    0.006,  "Assumption; CV 30% (unvalidated)",
    "baseline_dx_share",       0.55,    "beta",    0.55,    0.0825, "Assumption; CV 15%",
    "incidence_target",        0.012,   "gamma",   0.012,   0.0018, "Calibration anchor SADHS/IDF; CV 15%",
    "u_nodm",                  0.95,    "beta",    0.95,    0.02,   "HRQoL literature; SE 0.02",
    "u_undx",                  0.78,    "beta",    0.78,    0.04,   "HRQoL literature; SE 0.04",
    "u_gain",                  0.06,    "gamma",   0.06,    0.02,   "Diagnosed-utility gain (>=0 preserves ordering)",
    "c_undx",                  40,      "gamma",   40,      10,     "SA costing literature; CV 25%",
    "c_dx",                    90,      "gamma",   90,      22.5,   "SA costing literature; CV 25%",
    "pc_ml",                   12,      "gamma",   12,      2.4,    "Assumption; CV 20%",
    "pc_casefind",             10,      "gamma",   10,      2.0,    "Assumption; CV 20%",
    "pc_screen",               8,       "gamma",   8,       1.6,    "Assumption; CV 20%",
    "pc_ruleout",              5,       "gamma",   5,       1.0,    "Assumption; CV 20%",
    "rr_onset_ml",             0.80,    "lnorm",   0.80,    1.00,   "Lognormal; RR=1 at upper 2.5% (unvalidated effect)",
    "rr_onset_casefind",       0.82,    "lnorm",   0.82,    1.00,   "Lognormal; RR=1 at upper 2.5% (unvalidated effect)",
    "rr_onset_screen",         0.85,    "lnorm",   0.85,    1.00,   "Lognormal; RR=1 at upper 2.5% (unvalidated effect)",
    "mort_intercept",         -9.0,     "norm",   -9.0,     0.20,   "Log-odds; SD 0.20",
    "mort_exc_undx",           0.25,    "norm",    0.25,    0.05,   "Log-odds excess; SD 0.05",
    "mort_exc_dx",             0.15,    "norm",    0.15,    0.05,   "Log-odds excess; SD 0.05"
  )
}

# ---- Draw n parameter sets --------------------------------------------------
# Two deliberate modelling choices, both stated in the manuscript:
#   (1) UTILITY ORDERING IS PRESERVED. We sample the diagnosed-utility GAIN as a
#       non-negative quantity, so diagnosis never reduces quality of life in any
#       draw. This treats "diagnosis helps" as a structural assumption; the
#       DSA stress-tests it separately with a no-utility-gain scenario.
#   (2) Onset RRs are capped at 1 (a "preventive" effect cannot be harmful).
draw_psa_params <- function(n_draws, seed = 20260719L) {
  set.seed(seed)
  s <- psa_spec()
  cols <- lapply(seq_len(nrow(s)), function(i) {
    r <- s[i, ]
    switch(r$dist,
      beta  = beta_from_mean_sd(n_draws, r$arg1, r$arg2),
      gamma = gamma_from_mean_sd(n_draws, r$arg1, r$arg2),
      lnorm = lnorm_from_median_upper(n_draws, r$arg1, r$arg2),
      norm  = rnorm(n_draws, r$arg1, r$arg2))
  })
  names(cols) <- s$param
  D <- as_tibble(cols) %>%
    mutate(
      u_dx              = pmin(u_undx + u_gain, u_nodm),   # ordering preserved
      rr_onset_ml       = pmin(rr_onset_ml, 1),
      rr_onset_casefind = pmin(rr_onset_casefind, 1),
      rr_onset_screen   = pmin(rr_onset_screen, 1)
    )
  D$draw <- seq_len(n_draws)
  D
}

# Merge one draw into a full parameter list.
# NOTE: discount rate and time horizon are ANALYTIC CONVENTIONS, not uncertain
# quantities, so they are held fixed here and varied deterministically instead
# (standard ISPOR/CHEERS practice).
params_from_draw <- function(base, draw_row) {
  p <- base
  for (nm in names(draw_row)) if (nm %in% names(p)) p[[nm]] <- draw_row[[nm]]
  p$time_horizon  <- base$time_horizon
  p$discount_rate <- base$discount_rate
  p
}


# =============================================================================
# SECTION 3.  SELF-TESTS
# =============================================================================
# Fast invariant checks. They run every time so a silent regression in the
# engine or the samplers cannot propagate into published numbers.

cat("\n=== [1/8] Loading cohort and running self-tests ===\n")
df   <- prepare_cohort(".")
rng  <- make_rng(nrow(df), seed = 1234L)
base_params <- make_base_params()
cat(sprintf("Cohort: n = %s individuals\n", comma(nrow(df))))

# (a) With onset switched off, nobody should develop diabetes. A handful can
#     still appear because the onset probability is clamped at a 0.0001 floor.
p_test <- base_params; p_test$incidence_target <- 1e-9
stopifnot(sum(run_cohort(df, "baseline", p_test, rng)$incident_dm) < 0.01 * nrow(df))

# (b) Discounting and mortality must pull life-years below the raw horizon.
res_test <- run_cohort(df, "baseline", base_params, rng)
stopifnot(all(res_test$total_ly <= base_params$time_horizon + 1e-8),
          mean(res_test$total_ly) < base_params$time_horizon)

# (c) CRN determinism: identical parameters + identical draws => identical output.
stopifnot(identical(run_cohort(df, "fft_screening", base_params, rng)$total_qaly,
                    run_cohort(df, "fft_screening", base_params, rng)$total_qaly))

# (d) NMB arithmetic.
stopifnot(nmb(1, 100, 6000) == 5900)

# (e) Sampler moments and the preserved utility ordering.
set.seed(1)
.b <- beta_from_mean_sd(2e5, 0.25, 0.075)
stopifnot(abs(mean(.b) - 0.25) < 0.01, abs(sd(.b) - 0.075) < 0.01, all(.b >= 0 & .b <= 1))
.g <- gamma_from_mean_sd(2e5, 40, 10)
stopifnot(abs(mean(.g) - 40) < 0.5, abs(sd(.g) - 10) < 0.5, all(.g > 0))
.l <- lnorm_from_median_upper(2e5, 0.80, 1.00)
stopifnot(abs(median(.l) - 0.80) < 0.01, abs(quantile(.l, 0.975) - 1.00) < 0.03)
.D <- draw_psa_params(5000)
stopifnot(all(.D$u_dx >= .D$u_undx), all(.D$u_nodm >= .D$u_dx),
          all(.D$rr_onset_screen <= 1))
rm(p_test, res_test, .b, .g, .l, .D)
cat("All self-tests passed.\n")


# =============================================================================
# SECTION 4.  BASE CASE, EFFICIENT FRONTIER, AND STABILITY
# =============================================================================

cat("\n=== [2/8] Base case ===\n")
summary_tbl <- simulate_all(df, base_params, rng)
icer_tbl    <- icer_vs_baseline(summary_tbl)

# Final health-state distribution (feeds manuscript Table 7 and Figure 2).
state_dist <- map_dfr(STRATEGIES, function(s) {
  res <- run_cohort(df, s, base_params, rng)
  tibble(strategy = s, final_state = res$final_state) %>%
    count(strategy, final_state) %>%
    group_by(strategy) %>% mutate(prop = n / sum(n)) %>% ungroup()
}) %>% mutate(state = STATE_NAMES[final_state + 1L])

write_csv(summary_tbl, "output/basecase_summary.csv")
write_csv(icer_tbl,    "output/basecase_icer.csv")
write_csv(state_dist,  "output/basecase_state_distribution.csv")

print(summary_tbl); print(icer_tbl)

# ---- Conclusion gate --------------------------------------------------------
# Asserts the CRN-corrected findings. If these ever fail, STOP - do not adjust
# the engine to restore the old result; re-examine the diagnostic instead.
fr <- frontier(summary_tbl); active_eff <- setdiff(fr, "baseline")
cf  <- icer_tbl %>% filter(strategy == "fft_casefinding")
scr <- icer_tbl %>% filter(strategy == "fft_screening")
ro  <- icer_tbl %>% filter(strategy == "fft_ruleout")
mlt <- icer_tbl %>% filter(strategy == "ml_targeted")
gate <- c(
  casefinding_sole_efficient       = identical(active_eff, "fft_casefinding"),
  casefinding_positive_qaly        = cf$incr_qaly > 0,
  casefinding_icer_plausible       = cf$ICER > 2500 && cf$ICER < 5000,
  screening_positive_not_efficient = scr$incr_qaly > 0 && !("fft_screening" %in% fr),
  ml_positive_but_dominated        = mlt$incr_qaly > 0 && !("ml_targeted" %in% fr),
  ruleout_harmful                  = ro$incr_qaly < 0
)
cat("\nEfficient frontier:", paste(fr, collapse = " -> "), "\n"); print(gate)
if (!all(gate)) stop("CRN conclusion gate FAILED - investigate before proceeding.")
cat("Gate PASSED: FFT case-finding is the sole efficient active strategy.\n")

# ---- First-order (stochastic) stability -------------------------------------
# Re-runs the base case under 30 independent CRN streams. This is the evidence
# that the conclusion is a signal and not seed luck.
cat("\n=== [3/8] Stability across 30 seeds ===\n")
stab <- map_dfr(1:30, function(sd_i) {
  rng_s <- make_rng(nrow(df), seed = sd_i)
  sm_s  <- simulate_all(df, base_params, rng_s)
  icer_vs_baseline(sm_s) %>%
    filter(strategy == "fft_casefinding") %>%
    transmute(seed = sd_i, incr_cost, incr_qaly, ICER,
              sole_efficient = identical(setdiff(frontier(sm_s), "baseline"),
                                         "fft_casefinding"))
})
write_csv(stab, "output/basecase_stability.csv")
cat(sprintf("Case-finding sole-efficient in %d/%d seeds; ICER range %.0f-%.0f\n",
            sum(stab$sole_efficient), nrow(stab), min(stab$ICER), max(stab$ICER)))


# =============================================================================
# SECTION 5.  ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS (TORNADO) + SCENARIOS
# =============================================================================
# Outcome is the incremental NET MONETARY BENEFIT of the efficient strategy at
# 1x GDP per capita, NOT the ICER: NMB stays interpretable when incremental
# QALYs cross zero, where an ICER would flip sign meaninglessly.

cat("\n=== [4/8] One-way DSA ===\n")
LAMBDA  <- 6000                  # willingness-to-pay, ~1x GDP per capita (USD)
PRIMARY <- "fft_casefinding"     # the efficient strategy under CRN

cf_nmb <- function(params, lambda = LAMBDA) {
  s <- icer_vs_baseline(simulate_all(df, params, rng)) %>% filter(strategy == PRIMARY)
  nmb(s$incr_qaly, s$incr_cost, lambda)
}
base_nmb <- cf_nmb(base_params)

# Case-finding vs baseline depends only on case-finding's own parameters plus
# the shared ones, so screening/ML/rule-out parameters are not varied here.
ranges <- tribble(
  ~param,                    ~low,    ~high,   ~label,
  "p_detect_casefind_high",  0.21,    0.39,    "FFT case-finding detection prob (0.30 +/-30%)",
  "p_detect_baseline",       0.056,   0.104,   "Background detection (0.08 +/-30%)",
  "rr_onset_casefind",       0.65,    1.00,    "Onset RR, case-finding (0.82; 0.65-1.00)",
  "u_undx",                  0.72,    0.84,    "Undiagnosed utility (0.78; 0.72-0.84)",
  "u_dx",                    0.80,    0.90,    "Diagnosed utility (0.84; 0.80-0.90)",
  "c_dx",                    67.5,    112.5,   "Diagnosed-state cost (90 +/-25%)",
  "c_undx",                  30,      50,      "Undiagnosed-state cost (40 +/-25%)",
  "pc_casefind",             8,       12,      "FFT case-finding programme cost (10 +/-20%)",
  "baseline_dx_share",       0.44,    0.66,    "Baseline diagnosed share (0.55 +/-20%)",
  "incidence_target",        0.009,   0.015,   "Diabetes incidence target (0.012 +/-25%)",
  "mort_exc_undx",           0.125,   0.375,   "Excess mortality, undiagnosed (0.25 +/-50%)",
  "discount_rate",           0.00,    0.06,    "Discount rate (0.03; 0-6%)",
  "time_horizon",            10,      30,      "Time horizon (20; 10-30 yr)"
)

tornado <- map_dfr(seq_len(nrow(ranges)), function(i) {
  r  <- ranges[i, ]
  lo <- base_params; lo[[r$param]] <- if (r$param == "time_horizon") as.integer(r$low)  else r$low
  hi <- base_params; hi[[r$param]] <- if (r$param == "time_horizon") as.integer(r$high) else r$high
  tibble(param = r$param, label = r$label,
         nmb_low = cf_nmb(lo), nmb_high = cf_nmb(hi))
})

# Structural stress test of the one assumption the PSA deliberately holds fixed:
# that diagnosis improves quality of life at all.
noreward <- base_params; noreward$u_dx <- base_params$u_undx
tornado <- bind_rows(tornado, tibble(
  param = "no_utility_gain", label = "No utility gain from diagnosis (u_dx = u_undx)",
  nmb_low = cf_nmb(noreward), nmb_high = base_nmb))

# ---- Scenario table (manuscript Table 6) ------------------------------------
scen_run <- function(mods) {
  p <- base_params; for (nm in names(mods)) p[[nm]] <- mods[[nm]]
  icer_vs_baseline(simulate_all(df, p, rng)) %>%
    filter(strategy == PRIMARY) %>% transmute(incr_cost, incr_qaly, ICER)
}
scenarios <- list(
  "Base case"                    = list(),
  "Programme costs +20%"         = list(pc_ml=14.4, pc_casefind=12, pc_screen=9.6, pc_ruleout=6),
  "Programme costs -20%"         = list(pc_ml=9.6,  pc_casefind=8,  pc_screen=6.4, pc_ruleout=4),
  "Active detection +20%"        = list(p_detect_ml_high=0.30, p_detect_casefind_high=0.36,
                                        p_detect_screen_high=0.264, p_detect_ruleout_low=0.024),
  "Active detection -20%"        = list(p_detect_ml_high=0.20, p_detect_casefind_high=0.24,
                                        p_detect_screen_high=0.176, p_detect_ruleout_low=0.016),
  "No onset reduction"           = list(rr_onset_ml=1, rr_onset_casefind=1, rr_onset_screen=1),
  "Lower diagnosis utility gain" = list(u_dx=0.81),
  "Discount rate 5%"             = list(discount_rate=0.05),
  "Time horizon 10 years"        = list(time_horizon=10L),
  "Time horizon 30 years"        = list(time_horizon=30L)
)
scen_tbl <- imap_dfr(scenarios, ~ scen_run(.x) %>% mutate(scenario = .y)) %>%
  select(scenario, incr_cost, incr_qaly, ICER)
write_csv(scen_tbl, "output/dsa_scenarios.csv")
print(scen_tbl)

# ---- Figure 4: tornado ------------------------------------------------------
tornado <- tornado %>%
  mutate(lo = pmin(nmb_low, nmb_high), hi = pmax(nmb_low, nmb_high),
         span = hi - lo) %>%
  arrange(span)                                   # widest bar plotted at the top
write_csv(tornado %>% mutate(base_nmb = base_nmb), "output/dsa_oneway.csv")

lvl <- tornado$label
p_tornado <- ggplot(tornado) +
  geom_vline(xintercept = base_nmb, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey60") +      # cost-effectiveness boundary
  geom_rect(aes(xmin = lo, xmax = hi,
                ymin = match(label, lvl) - 0.4, ymax = match(label, lvl) + 0.4),
            fill = "#4C72B0", alpha = 0.85) +
  scale_y_continuous(breaks = seq_along(lvl), labels = lvl) +
  labs(x = sprintf("Incremental net monetary benefit at WTP $%s/QALY (US$)",
                   format(LAMBDA, big.mark = ",")),
       y = NULL,
       title = "One-way sensitivity analysis: FFT case-finding vs baseline",
       subtitle = "Dashed line = base-case NMB; solid grey = zero (cost-effectiveness boundary)") +
  theme_minimal(base_size = 11)
ggsave("output/figures/fig4_tornado.png", p_tornado, width = 10.5, height = 6, dpi = 300)
ggsave("output/figures/fig4_tornado.pdf", p_tornado, width = 10.5, height = 6)

flip <- tornado %>% filter(lo < 0, hi > 0) %>% pull(label)
cat(sprintf("Base-case NMB at $%s/QALY: %.0f\n", format(LAMBDA, big.mark=","), base_nmb))
cat("Parameters that reverse the conclusion:",
    if (length(flip)) paste(flip, collapse = "; ") else "none", "\n")


# =============================================================================
# SECTION 6.  PROBABILISTIC SENSITIVITY ANALYSIS
# =============================================================================
# The SAME rng is reused for every draw, so between-draw variation reflects
# PARAMETER uncertainty rather than simulation noise. Raise N_DRAWS to check
# convergence (1,000 -> 2,000 moves P(optimal) by <2 percentage points).

cat("\n=== [5/8] PSA (this is the slow step) ===\n")
N_DRAWS <- 1000L
D <- draw_psa_params(N_DRAWS, seed = 20260719L)

psa <- map_dfr(seq_len(N_DRAWS), function(d_i) {
  p  <- params_from_draw(base_params, as.list(D[d_i, ]))
  sm <- simulate_all(df, p, rng)
  bl <- sm %>% filter(strategy == "baseline")
  sm %>% mutate(draw      = d_i,
                incr_cost = mean_cost - bl$mean_cost,
                incr_qaly = mean_qaly - bl$mean_qaly)
})
write_csv(psa, "output/psa_draws_results.csv")

active <- c("fft_screening", "fft_casefinding", "fft_ruleout", "ml_targeted")

# ---- Incremental estimates with 95% credible intervals ----------------------
incr_summary <- psa %>% filter(strategy %in% active) %>%
  group_by(strategy) %>%
  summarise(mean_incr_cost = mean(incr_cost),
            cost_lo = quantile(incr_cost, 0.025), cost_hi = quantile(incr_cost, 0.975),
            mean_incr_qaly = mean(incr_qaly),
            qaly_lo = quantile(incr_qaly, 0.025), qaly_hi = quantile(incr_qaly, 0.975),
            .groups = "drop")
write_csv(incr_summary, "output/psa_incremental_summary.csv")

# ---- CEAC: probability each strategy is OPTIMAL (max NMB) at each WTP -------
lambdas  <- seq(0, 30000, by = 500)
nmb_wide <- psa %>% select(draw, strategy, mean_cost, mean_qaly)

ceac <- map_dfr(lambdas, function(l) {
  nmb_wide %>% mutate(nmb = mean_qaly * l - mean_cost) %>%
    group_by(draw) %>% slice_max(nmb, n = 1, with_ties = FALSE) %>% ungroup() %>%
    count(strategy) %>% mutate(lambda = l, prob = n / N_DRAWS)
})
ceac_full <- expand_grid(lambda = lambdas, strategy = unique(psa$strategy)) %>%
  left_join(ceac %>% select(lambda, strategy, prob), by = c("lambda", "strategy")) %>%
  mutate(prob = replace_na(prob, 0))
write_csv(ceac_full, "output/psa_ceac.csv")

# CEAF: which strategy has the highest EXPECTED NMB at each WTP (the decision
# you would actually take, as distinct from how confident you can be in it).
ceaf <- map_dfr(lambdas, function(l) {
  psa %>% mutate(nmb = mean_qaly * l - mean_cost) %>%
    group_by(strategy) %>% summarise(enmb = mean(nmb), .groups = "drop") %>%
    slice_max(enmb, n = 1, with_ties = FALSE) %>% mutate(lambda = l)
})
write_csv(ceaf, "output/psa_ceaf.csv")

# ---- Probability table at the two policy thresholds -------------------------
prob_at <- function(l) {
  opt <- nmb_wide %>% mutate(nmb = mean_qaly * l - mean_cost) %>%
    group_by(draw) %>% slice_max(nmb, n = 1, with_ties = FALSE) %>% ungroup() %>%
    count(strategy) %>% mutate(p_optimal = n / N_DRAWS) %>% select(strategy, p_optimal)
  ce <- psa %>% filter(strategy %in% active) %>% group_by(strategy) %>%
    summarise(p_ce_vs_baseline = mean(incr_qaly * l - incr_cost > 0), .groups = "drop")
  full_join(ce, opt, by = "strategy") %>%
    mutate(lambda = l, p_optimal = replace_na(p_optimal, 0))
}
prob_tbl <- bind_rows(prob_at(6000), prob_at(21000))
write_csv(prob_tbl, "output/psa_probability_table.csv")
print(incr_summary); print(prob_tbl)

# ---- Figure 5: cost-effectiveness plane -------------------------------------
p_ce_plane <- ggplot(psa %>% filter(strategy %in% active),
                     aes(incr_qaly, incr_cost, colour = strategy)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 6000,  intercept = 0, linetype = "dashed") +   # 1x GDP
  geom_abline(slope = 21000, intercept = 0, linetype = "dotted") +   # 3x GDP
  geom_point(alpha = 0.25, size = 0.7) +
  labs(x = "Incremental QALYs vs no active screening",
       y = "Incremental cost vs no active screening (US$)", colour = "Strategy",
       title = "Probabilistic cost-effectiveness plane (1,000 Monte Carlo draws)",
       subtitle = "Dashed = $6,000/QALY (1x GDP); dotted = $21,000/QALY (3x GDP)") +
  theme_minimal(base_size = 11)
ggsave("output/figures/fig5_ce_plane.png", p_ce_plane, width = 8, height = 6.5, dpi = 300)
ggsave("output/figures/fig5_ce_plane.pdf", p_ce_plane, width = 8, height = 6.5)

# ---- Figure 6: CEAC ---------------------------------------------------------
p_ceac <- ggplot(ceac_full, aes(lambda, prob, colour = strategy)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 6000,  linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 21000, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(labels = dollar) +
  labs(x = "Willingness to pay (US$/QALY)", y = "Probability strategy is optimal",
       colour = "Strategy", title = "Cost-effectiveness acceptability curves",
       subtitle = "Vertical lines: $6,000 (1x GDP) and $21,000 (3x GDP) per QALY") +
  theme_minimal(base_size = 11)
ggsave("output/figures/fig6_ceac.png", p_ceac, width = 8, height = 6, dpi = 300)
ggsave("output/figures/fig6_ceac.pdf", p_ceac, width = 8, height = 6)

# ---- Sanity check: the PSA cloud should be centred on the base case ---------
det <- icer_vs_baseline(simulate_all(df, base_params, rng)) %>% filter(strategy == PRIMARY)
cen <- incr_summary %>% filter(strategy == PRIMARY)
cat(sprintf("Sanity: PSA mean incr QALY %.4f vs base case %.4f\n",
            cen$mean_incr_qaly, det$incr_qaly))
if (abs(cen$mean_incr_qaly - det$incr_qaly) > 0.02)
  warning("PSA centroid drifts from the base case by >0.02 QALY - check distributions.")


# =============================================================================
# SECTION 7.  EXPECTED VALUE OF PERFECT INFORMATION
# =============================================================================
# EVPI = E[max NMB] - max E[NMB]: the expected loss from having to decide under
# current uncertainty. It monetises what resolving that uncertainty is worth.
#
# NOTE ON SHAPE: EVPI here RISES across the plausible WTP range rather than
# peaking and falling. That is a real feature, not a bug - the QALY-maximising
# strategy keeps changing across draws (case-finding 595, ML 365, screening 40
# of 1,000), so decision uncertainty never resolves.

cat("\n=== [6/8] EVPI ===\n")
ELIGIBLE_POP <- 35e6    # approx. SA adults 20+; REPLACE with a cited estimate
POP_FACTOR   <- ELIGIBLE_POP        # one-time prevalent-population decision

evpi <- map_dfr(lambdas, function(l) {
  w <- psa %>% mutate(nmb = mean_qaly * l - mean_cost)
  enmb_perfect <- w %>% group_by(draw) %>% summarise(m = max(nmb), .groups = "drop") %>%
    pull(m) %>% mean()                                   # E[max NMB]
  enmb_current <- w %>% group_by(strategy) %>%
    summarise(e = mean(nmb), .groups = "drop") %>% pull(e) %>% max()   # max E[NMB]
  tibble(lambda = l,
         evpi_per_person = enmb_perfect - enmb_current,
         evpi_population = (enmb_perfect - enmb_current) * POP_FACTOR)
})
write_csv(evpi, "output/evpi_curve.csv")

p_evpi <- ggplot(evpi, aes(lambda, evpi_population / 1e6)) +
  geom_line(linewidth = 0.9, colour = "#C44E52") +
  geom_vline(xintercept = 6000,  linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 21000, linetype = "dotted", colour = "grey50") +
  scale_x_continuous(labels = dollar) +
  labs(x = "Willingness to pay (US$/QALY)", y = "Population EVPI (US$ millions)",
       title = "Expected value of perfect information",
       subtitle = sprintf("Scaled to %s eligible adults (one-time decision)",
                          comma(ELIGIBLE_POP))) +
  theme_minimal(base_size = 11)
ggsave("output/figures/fig7_evpi.png", p_evpi, width = 8, height = 6, dpi = 300)
ggsave("output/figures/fig7_evpi.pdf", p_evpi, width = 8, height = 6)

cat(sprintf("EVPI at $6,000/QALY: $%.2f per person (population $%s)\n",
            evpi$evpi_per_person[evpi$lambda == 6000],
            comma(round(evpi$evpi_population[evpi$lambda == 6000]))))


# =============================================================================
# SECTION 8.  SUPPLEMENTARY TABLE S1
# =============================================================================
# Built directly from psa_spec() so the published table cannot disagree with
# what was actually sampled.

cat("\n=== [7/8] Supplementary Table S1 ===\n")
s1 <- psa_spec() %>%
  mutate(distribution = recode(dist, beta = "Beta", gamma = "Gamma",
                               lnorm = "Lognormal", norm = "Normal"),
         dispersion = case_when(
           dist %in% c("beta", "gamma", "norm") ~ sprintf("mean %s, SD %s", arg1, arg2),
           dist == "lnorm" ~ sprintf("median %s, upper 97.5%% = %s", arg1, arg2))) %>%
  select(Parameter = param, `Base value` = base, Distribution = distribution,
         Dispersion = dispersion, `Basis / source` = basis)
write_csv(s1, "output/manuscript_edits/supplementary_table_S1_distributions.csv")

# ---- Table footnote ---------------------------------------------------------
# Written to its own file so that the manuscript builder (integrate_manuscript.py)
# reads exactly this text: the note in the Word document can never drift from
# the note generated alongside the table.
s1_footnote <- c(
  paste(
    "Abbreviations: CV = coefficient of variation (the standard deviation",
    "expressed as a percentage of the mean); SD = standard deviation;",
    "SE = standard error; RR = relative risk; HRQoL = health-related quality of",
    "life; IDF = International Diabetes Federation; SADHS = South Africa",
    "Demographic and Health Survey; SA = South Africa."),
  paste(
    "Notes: Distributions were chosen to respect each parameter's natural range —",
    "Beta for probabilities and utilities (bounded 0–1), Gamma for costs and rates",
    "(positive and right-skewed), lognormal for relative risks, and Normal for",
    "parameters expressed on the log-odds scale. Beta and Gamma distributions were",
    "fitted by method of moments to the stated mean and standard deviation. Where",
    "dispersion is reported as a CV, the standard deviation equals the CV multiplied",
    "by the base value (for example, a base value of 0.30 with CV 30% corresponds to",
    "an SD of 0.09); where a parameter is not naturally expressed as a proportion of",
    "its mean, dispersion is given directly as an SD or SE. For the onset relative",
    "risks, the lognormal distribution was specified by its median and upper 97.5th",
    "percentile, set so that RR = 1 (no preventive effect) lies at the upper tail.",
    "u_gain denotes the utility gain conferred by diagnosis; it is drawn as a",
    "non-negative quantity and the diagnosed-state utility is obtained as",
    "u_dx = min(u_undx + u_gain, u_nodm), so the ordering u_nodm ≥ u_dx ≥ u_undx",
    "holds in every draw. The discount rate and time horizon are analytic conventions",
    "rather than uncertain quantities and were varied deterministically (Table 6 and",
    "Figure 4) rather than sampled here. The analysis used 1,000 Monte Carlo draws",
    "under common random numbers. Parameter names correspond to the variables in the",
    "analysis code. Several parameters are calibrated model assumptions rather than",
    "empirically estimated quantities; the distributions therefore represent assumed",
    "rather than observed uncertainty, and the probabilistic results should be",
    "interpreted accordingly.")
)
writeLines(s1_footnote, "output/manuscript_edits/supplementary_table_S1_footnote.txt")

writeLines(
  c("| Parameter | Base value | Distribution | Dispersion | Basis / source |",
    "|---|---|---|---|---|",
    apply(s1, 1, function(r) sprintf("| %s | %s | %s | %s | %s |",
                                     r[1], r[2], r[3], r[4], r[5])),
    "", s1_footnote),
  "output/manuscript_edits/supplementary_table_S1_distributions.md")


# =============================================================================
# SECTION 9.  FIGURE 1 - COST-EFFECTIVENESS PLANE (BASE CASE)
# =============================================================================
# Deterministic companion to the probabilistic Figure 5. Points are taken
# straight from icer_tbl so the figure cannot drift from the results table.

cat("\n=== [8/8] Figure 1 ===\n")
WTP_LOW  <- 6000     # ~1x GDP per capita
WTP_HIGH <- 21000    # ~3x GDP per capita

# Arithmetic guard on the threshold lines (they are drawn as cost = WTP * QALY).
stopifnot(isTRUE(all.equal(WTP_LOW  * 0.04, 240)),
          isTRUE(all.equal(WTP_HIGH * 0.04, 840)))

fig1_labels <- c(fft_ruleout = "FFT Rule-out", ml_targeted = "ML Targeted",
                 fft_casefinding = "FFT Case-finding", fft_screening = "FFT Screening")
ce_data <- icer_tbl %>%
  transmute(strategy = unname(fig1_labels[strategy]), incr_qaly, incr_cost)
plot_data <- bind_rows(tibble(strategy = "Baseline", incr_qaly = 0, incr_cost = 0),
                       ce_data) %>%
  mutate(strategy = factor(strategy,
    levels = c("Baseline", "FFT Rule-out", "ML Targeted",
               "FFT Case-finding", "FFT Screening")))

label_data <- plot_data %>%
  mutate(label_dy = case_when(strategy == "ML Targeted"      ~  80,
                              strategy == "FFT Case-finding" ~ -80,
                              strategy == "FFT Screening"    ~ -80,
                              strategy == "FFT Rule-out"     ~  80,
                              TRUE                           ~ -80),
         label_dx = case_when(strategy == "FFT Screening" ~ -0.004,
                              strategy == "Baseline"      ~  0.005,
                              TRUE                        ~  0.000),
         label_y = incr_cost + label_dy, label_x = incr_qaly + label_dx)

x_limits <- c(-0.02, 0.055); y_limits <- c(-300, 900)
wtp_lines <- tibble(x = seq(x_limits[1], x_limits[2], length.out = 2)) %>%
  crossing(wtp = c(WTP_LOW, WTP_HIGH)) %>%
  mutate(y = wtp * x,
         wtp_label = factor(wtp, levels = c(WTP_LOW, WTP_HIGH),
                            labels = c(paste0("WTP = ", dollar(WTP_LOW),  "/QALY"),
                                       paste0("WTP = ", dollar(WTP_HIGH), "/QALY"))))
shaded_region <- tibble(x = c(0, x_limits[2], x_limits[2]),
                        y = c(0, 0, WTP_HIGH * x_limits[2]))
strategy_colors <- c("Baseline" = "grey45", "FFT Rule-out" = "#D62728",
                     "ML Targeted" = "#F4A460", "FFT Case-finding" = "#1F5FA6",
                     "FFT Screening" = "#2CA02C")

p_fig1 <- ggplot() +
  geom_polygon(data = shaded_region, aes(x, y), fill = "#DCE6F0", alpha = 0.6) +
  geom_line(data = wtp_lines, aes(x, y, linetype = wtp_label),
            color = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_point(data = plot_data, aes(incr_qaly, incr_cost, color = strategy), size = 4) +
  geom_text(data = label_data, aes(label_x, label_y, label = strategy, color = strategy),
            size = 4.2, show.legend = FALSE) +
  scale_color_manual(values = strategy_colors, guide = "none") +
  scale_linetype_manual(name = NULL, values = c("dashed", "dotted")) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0),
                     labels = function(v) paste0("$", comma(v))) +
  labs(title = "Cost-Effectiveness Plane", x = "Incremental QALYs",
       y = "Incremental Cost (USD)") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(color = "black", size = 13),
        panel.grid = element_blank(),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 12),
        plot.margin = margin(10, 14, 10, 10))

ggsave("Figure 1_corrected.png", p_fig1, width = 9, height = 9 * 1519 / 1800,
       dpi = 600, bg = "white")
ggsave("Figure 1_corrected.pdf", p_fig1, width = 9, height = 9 * 1519 / 1800,
       device = cairo_pdf)


# =============================================================================
#  DONE
# =============================================================================
cat("\n=============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Elapsed: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("\nHeadline results\n")
cat(sprintf("  Efficient frontier      : %s\n", paste(fr, collapse = " -> ")))
cat(sprintf("  FFT case-finding ICER   : $%.0f/QALY\n", cf$ICER))
cat(sprintf("  Stability               : sole-efficient in %d/30 seeds\n",
            sum(stab$sole_efficient)))
cat(sprintf("  P(optimal) at $6,000    : %.1f%%\n",
            100 * prob_tbl$p_optimal[prob_tbl$lambda == 6000 &
                                     prob_tbl$strategy == PRIMARY]))
cat(sprintf("  P(cost-effective) $6,000: %.1f%%\n",
            100 * prob_tbl$p_ce_vs_baseline[prob_tbl$lambda == 6000 &
                                            prob_tbl$strategy == PRIMARY]))
cat("\nNOTE: Figures 2 and 3 are produced by make_figures_2_3.py (matplotlib),\n")
cat("      which reads output/basecase_state_distribution.csv and\n")
cat("      output/basecase_summary.csv written above.\n")
cat("=============================================================\n")
