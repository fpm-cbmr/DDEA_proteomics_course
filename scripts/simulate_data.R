# =============================================================================
# simulate_data.R
# -----------------------------------------------------------------------------
# Build a synthetic single-fiber proteomics dataset for the DDEA proteomics
# course. The structure (8 subjects, 2 timepoints, ~95 fibers each, 2373
# proteins) and per-gene mean/SD priors are taken from the real REX dataset,
# but every value is regenerated and subject IDs are anonymized.
#
# Outputs (relative to project root):
#   data/single_fiber_proteomics.rds   list(intensities, metadata)
#   results/model_results.rds          named list (one entry per gene)
#                                       each entry: list(model_results = <tbl>)
#                                       matching the real GAM-results schema
#
# Notes on the model_results schema (matches the real REX file):
#   gene
#   Estimate_(Intercept)   Std. Error_(Intercept)   t value_(Intercept)   Pr(>|t|)_(Intercept)
#   Estimate_time1         Std. Error_time1         t value_time1         Pr(>|t|)_time1
#   p_val_fiber_type       <- p-value for association of intensity with PC1
#   p_val_interaction      <- p-value for the time x PC1 interaction
#
# The real pipeline fits a GAM with s(PC1, by = time). Here we fit lightweight
# linear surrogates per protein on the synthetic data:
#   - time effect from lm(subject_mean ~ time)
#   - PC1 effect from lm(intensity ~ PC1) on fiber-level data
#   - interaction from lm(intensity ~ time * PC1)
# This is enough for course visualisations because the synthetic data are
# generated WITH those exact effects injected, so every "significantly
# associated with PC1" protein really does vary with PC1 in the matrix.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ---- 0. Config --------------------------------------------------------------

REAL_PATH <- "C:/Users/jns822/Desktop/Scripts/REX_single_fiber/data/proteome/data_model_proteome_T1_T3.rds"
OUT_DATA  <- "data/single_fiber_proteomics.rds"
OUT_MODEL <- "results/model_results.rds"

set.seed(20260527)

N_SUBJECTS      <- 8
TIMES           <- c("PRE", "POST")
FIBERS_PER_CELL <- 95
NA_RATE         <- 0.05

# Effect fractions & magnitudes (chosen so simulated p-value densities roughly
# match the real REX results: ~94% fiber-type, ~38% interaction).
FRAC_DA              <- 0.18  # proteins with non-zero PRE vs POST main effect
TIME_EFFECT_SD       <- 0.6
FIBERTYPE_FRAC       <- 0.90  # proteins loaded on PC1
FIBERTYPE_LOADING_SD <- 0.04
FRAC_INTERACTION     <- 0.30  # proteins with PC1 x time interaction
INTERACTION_COEF_SD  <- 0.020

SUBJECT_RE_SD        <- 0.20
RESID_NOISE          <- 0.30

# ---- 1. Per-gene priors from the real data ----------------------------------

cat("Loading real data (priors only)...\n")
real <- readRDS(REAL_PATH)

gene_prior <- real %>%
  group_by(genes) %>%
  summarise(
    mean_int = mean(intensities, na.rm = TRUE),
    sd_int   = sd(intensities,   na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  filter(!is.na(mean_int), !is.na(sd_int)) %>%
  arrange(desc(mean_int))

rm(real); gc(verbose = FALSE)

genes   <- gene_prior$genes
n_genes <- length(genes)
cat(sprintf("Priors built for %d proteins.\n", n_genes))

# ---- 2. Design --------------------------------------------------------------

subject_ids <- sprintf("S%02d", seq_len(N_SUBJECTS))

metadata <- expand_grid(
    subject_id = subject_ids,
    time       = TIMES,
    fiber_num  = seq_len(FIBERS_PER_CELL)
  ) %>%
  mutate(fiber_id = sprintf("%s_%s_%03d", subject_id, time, fiber_num)) %>%
  select(fiber_id, subject_id, time, fiber_num)

n_fibers <- nrow(metadata)

# Latent fiber-type axis (becomes PC1). Mixture of slow / hybrid / fast.
mix <- sample(c(-1, 0, 1), n_fibers, replace = TRUE, prob = c(0.45, 0.10, 0.45))
metadata$PC1_latent <- mix * 20 + rnorm(n_fibers, 0, 6)

# Time indicator: 0 / 1
time_indicator <- ifelse(metadata$time == "POST", 1, 0)

# Subject-level global shifts (small, on top of per-protein subject RE)
subject_global_shift <- setNames(rnorm(N_SUBJECTS, 0, 0.10), subject_ids)
fiber_subject_idx <- match(metadata$subject_id, subject_ids)

# Per-protein truths (the values we INJECT and that the fitted model should recover)
is_fibertype   <- runif(n_genes) < FIBERTYPE_FRAC
ft_loading     <- ifelse(is_fibertype, rnorm(n_genes, 0, FIBERTYPE_LOADING_SD), 0)

is_da          <- runif(n_genes) < FRAC_DA
true_log2fc    <- ifelse(is_da, rnorm(n_genes, 0, TIME_EFFECT_SD), 0)

is_interaction <- runif(n_genes) < FRAC_INTERACTION
true_interact  <- ifelse(is_interaction, rnorm(n_genes, 0, INTERACTION_COEF_SD), 0)

cat(sprintf("Truth: %d DA (time), %d fiber-type, %d interaction.\n",
            sum(is_da), sum(is_fibertype), sum(is_interaction)))

# ---- 3. Generate intensities ------------------------------------------------

cat("Simulating intensity matrix...\n")
intensities <- matrix(
  NA_real_,
  nrow = n_genes, ncol = n_fibers,
  dimnames = list(genes, metadata$fiber_id)
)

pc1_vec <- metadata$PC1_latent

for (i in seq_len(n_genes)) {
  mu     <- gene_prior$mean_int[i]
  sigma  <- gene_prior$sd_int[i]
  s_re   <- rnorm(N_SUBJECTS, 0, SUBJECT_RE_SD)

  baseline <- mu +
              s_re[fiber_subject_idx] +
              subject_global_shift[fiber_subject_idx] +
              ft_loading[i]    * pc1_vec +
              true_log2fc[i]   * time_indicator +
              true_interact[i] * pc1_vec * time_indicator

  noise_sd <- sqrt(sigma^2 + RESID_NOISE^2)
  intensities[i, ] <- baseline + rnorm(n_fibers, 0, noise_sd)
}

# MNAR-style missingness — low intensities more likely to be NA
cat("Injecting MNAR missingness...\n")
flat   <- as.vector(intensities)
ranks  <- rank(flat, ties.method = "first")
prob_na <- NA_RATE * 2 * (1 - ranks / length(ranks))
prob_na <- pmin(pmax(prob_na, 0), 1)
flat[runif(length(flat)) < prob_na] <- NA_real_
intensities <- matrix(flat, nrow = n_genes, ncol = n_fibers,
                      dimnames = dimnames(intensities))
cat(sprintf("Realized NA fraction: %.3f\n", mean(is.na(intensities))))

# Compute realized PC1 from the synthetic matrix; store in metadata.
cat("Computing PC1 of synthetic matrix...\n")
imp <- intensities
row_means <- rowMeans(imp, na.rm = TRUE)
na_idx <- which(is.na(imp), arr.ind = TRUE)
imp[na_idx] <- row_means[na_idx[, "row"]]
pc <- prcomp(t(imp), center = TRUE, scale. = FALSE, rank. = 1)
metadata$PC1 <- as.numeric(pc$x[, 1])
metadata$PC1_latent <- NULL

# Use the REALIZED PC1 (not the latent) for the per-protein PC1 / interaction tests,
# so the model sees what a user analysing the data would see.
PC1_real <- metadata$PC1

# ---- 4. Per-protein model_results -------------------------------------------

cat("Fitting per-protein models for model_results...\n")
t_factor <- factor(metadata$time, levels = c("PRE", "POST"))
# `time1` in the real schema corresponds to the second level vs. the first.

t0 <- Sys.time()
model_results_list <- vector("list", n_genes)
names(model_results_list) <- genes

# Pre-build the long-form helpers used inside the loop
subj_time_grid <- metadata %>% select(fiber_id, subject_id, time)

for (i in seq_len(n_genes)) {
  y_fiber <- intensities[i, ]
  ok      <- !is.na(y_fiber)

  # ---- (a) time main effect from subject means ----
  subj_means <- tapply(y_fiber[ok], list(metadata$subject_id[ok], metadata$time[ok]),
                       mean, na.rm = TRUE)
  # subj_means is 8 x 2 matrix with columns "POST","PRE" (alphabetical) — be explicit:
  pre  <- subj_means[, "PRE"]
  post <- subj_means[, "POST"]
  diffs <- post - pre
  n_pair <- sum(!is.na(diffs))
  est_time  <- mean(diffs, na.rm = TRUE)
  sd_time   <- sd(diffs, na.rm = TRUE)
  se_time   <- sd_time / sqrt(n_pair)
  t_time    <- est_time / se_time
  p_time    <- 2 * pt(-abs(t_time), df = n_pair - 1)

  intercept_est <- mean(c(pre, post), na.rm = TRUE)  # grand mean
  intercept_se  <- sd(c(pre, post), na.rm = TRUE) / sqrt(length(c(pre, post)))
  intercept_t   <- intercept_est / intercept_se
  intercept_p   <- 2 * pt(-abs(intercept_t),
                          df = length(c(pre, post)) - 1)

  # ---- (b) PC1 main effect (proxy for GAM s(PC1)) ----
  fit_pc <- tryCatch(
    .lm.fit(cbind(1, PC1_real[ok]), y_fiber[ok]),
    error = function(e) NULL
  )
  if (!is.null(fit_pc)) {
    rss     <- sum(fit_pc$residuals^2)
    n_obs   <- length(fit_pc$residuals)
    df_res  <- n_obs - 2
    sigma2  <- rss / df_res
    xtx_inv <- solve(crossprod(cbind(1, PC1_real[ok])))
    coef_pc <- fit_pc$coefficients[2]
    se_pc   <- sqrt(sigma2 * xtx_inv[2, 2])
    t_pc    <- coef_pc / se_pc
    p_fiber_type <- 2 * pt(-abs(t_pc), df = df_res)
  } else {
    p_fiber_type <- NA_real_
  }

  # ---- (c) PC1 x time interaction ----
  X <- cbind(
    1,
    time_indicator[ok],
    PC1_real[ok],
    time_indicator[ok] * PC1_real[ok]
  )
  fit_ix <- tryCatch(.lm.fit(X, y_fiber[ok]), error = function(e) NULL)
  if (!is.null(fit_ix)) {
    rss    <- sum(fit_ix$residuals^2)
    df_res <- length(fit_ix$residuals) - 4
    sigma2 <- rss / df_res
    xtx_inv <- solve(crossprod(X))
    coef_ix <- fit_ix$coefficients[4]
    se_ix   <- sqrt(sigma2 * xtx_inv[4, 4])
    t_ix    <- coef_ix / se_ix
    p_inter <- 2 * pt(-abs(t_ix), df = df_res)
  } else {
    p_inter <- NA_real_
  }

  mr <- tibble(
    gene                       = genes[i],
    `Estimate_(Intercept)`     = intercept_est,
    `Std. Error_(Intercept)`   = intercept_se,
    `t value_(Intercept)`      = intercept_t,
    `Pr(>|t|)_(Intercept)`     = intercept_p,
    Estimate_time1             = est_time,
    `Std. Error_time1`         = se_time,
    `t value_time1`            = t_time,
    `Pr(>|t|)_time1`           = p_time,
    p_val_fiber_type           = p_fiber_type,
    p_val_interaction          = p_inter
  )

  model_results_list[[i]] <- list(model_results = mr)

  if (i %% 500 == 0) {
    cat(sprintf("  fitted %d / %d (%.1f s elapsed)\n",
                i, n_genes, as.numeric(Sys.time() - t0, units = "secs")))
  }
}
cat(sprintf("Per-protein fitting done in %.1f s.\n",
            as.numeric(Sys.time() - t0, units = "secs")))

# ---- 5. Save ----------------------------------------------------------------

dir.create(dirname(OUT_DATA),  showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_MODEL), showWarnings = FALSE, recursive = TRUE)

saveRDS(list(intensities = intensities, metadata = metadata), OUT_DATA)
saveRDS(model_results_list, OUT_MODEL)

cat("\nWrote:\n  ", OUT_DATA, "\n  ", OUT_MODEL, "\n")

# ---- 6. Brief sanity summary ------------------------------------------------

flat_mr <- bind_rows(lapply(model_results_list, `[[`, "model_results"))
cat("\n=== synthetic model_results summary ===\n")
cat("rows:", nrow(flat_mr), "\n")
cat("time main effect significant (Pr(>|t|)_time1 < 0.05): ",
    sum(flat_mr$`Pr(>|t|)_time1` < 0.05, na.rm = TRUE), "\n", sep = "")
cat("p_val_fiber_type < 0.05: ",
    sum(flat_mr$p_val_fiber_type   < 0.05, na.rm = TRUE), "\n", sep = "")
cat("p_val_interaction < 0.05: ",
    sum(flat_mr$p_val_interaction  < 0.05, na.rm = TRUE), "\n", sep = "")

cat("Done.\n")
