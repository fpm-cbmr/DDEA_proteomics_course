# data_simulation_v1.R
# Option B — PCA reconstruction with calibrated residual noise
#
# PURPOSE:
#   Generate a public-safe synthetic dataset that preserves the dominant
#   biological structure of the real single-fiber proteomics data (fiber-type
#   separation, time effect, subject variation, protein-protein correlations)
#   while ensuring no individual fiber's real intensities are present in the
#   output.
#
# HOW TO RUN:
#   Source this script from the project root after sourcing or running the
#   tutorial up to and including the PCA chunk, so that pca_res, imp,
#   intensities, and metadata are already in your R session.
#   Alternatively, set LOAD_DATA <- TRUE below to load fresh from disk.
#
# OUTPUT:
#   data/single_fiber_proteomics_simulated.rds
#   Same list structure as the real data: list(intensities, metadata)
#   Ready to swap in as data/single_fiber_proteomics.rds for public release.
#
# TUNING:
#   The only knob is K (number of PCs retained). Increase K to preserve more
#   fine-grained structure; decrease K for stronger privacy. The diagnostic
#   plots printed at the end help you decide.

library(tidyverse)

# ── Optional: load data fresh from disk ───────────────────────────────────────
LOAD_DATA <- TRUE   # set TRUE if running independently of the tutorial session

if (LOAD_DATA) {
  fiber_data  <- readRDS("data/single_fiber_proteomics.rds")
  intensities <- fiber_data$intensities
  metadata    <- fiber_data$metadata

  imp <- intensities |>
    PhosR::tImpute() |>
    as.data.frame()

  set.seed(42)

  pca_res <- prcomp(t(imp), scale. = TRUE)

  metadata <- metadata |>
    left_join(
      as_tibble(pca_res$x, rownames = "fiber_id") |> select(fiber_id, PC1, PC2),
      by = "fiber_id"
    )
}

# ── Sanity check that required objects exist ──────────────────────────────────
stopifnot(
  exists("pca_res"),
  exists("imp"),
  exists("intensities"),
  exists("metadata")
)

# =============================================================================
# STEP 0: Tuning parameter
# =============================================================================
# K = number of PCs to retain for reconstruction.
# These PCs carry the biology the tutorial depends on. The remaining PCs
# are replaced by calibrated synthetic noise, severing the link to individual
# fiber profiles.
#
# Guidelines:
#   K = 30-50  — good balance; captures fiber type, time, subject structure
#   K < 20     — too few; PCA scatter and UMAP may look noisier than real data
#   K > 100    — individual-level information leaks back in; weaker privacy

K <- 50

# =============================================================================
# STEP 1: Reconstruct in PCA-scaled space  (fibers × proteins)
# =============================================================================
# prcomp(t(imp), scale. = TRUE) centers and scales by protein before computing
# PCs, so reconstruction must happen in that same scaled space.
#
#   pca_res$x        — scores matrix    (fibers   × all PCs)
#   pca_res$rotation — loadings matrix  (proteins × all PCs)
#   pca_res$center   — per-protein means (used to unscale at the end)
#   pca_res$scale    — per-protein SDs   (used to unscale at the end)

reconstructed_scaled <- pca_res$x[, 1:K] %*% t(pca_res$rotation[, 1:K])
# dim: fibers × proteins, in centered-and-scaled space

var_captured <- 100 * sum(pca_res$sdev[1:K]^2) / sum(pca_res$sdev^2)
cat(sprintf("K = %d PCs capture %.1f%% of total variance.\n", K, var_captured))

# =============================================================================
# STEP 2: Estimate per-protein residual noise from the discarded PCs
# =============================================================================
# The residuals = real (scaled) data − K-PC reconstruction.
# Their SD per protein tells us how much variance the discarded PCs carry
# for that protein. We use this to calibrate the synthetic noise so the
# simulated data has the same total variance structure as the real data.

data_scaled      <- scale(t(imp), center = pca_res$center, scale = pca_res$scale)
residuals_scaled <- data_scaled - reconstructed_scaled
residual_sd      <- apply(residuals_scaled, 2, sd, na.rm = TRUE)

cat(sprintf(
  "Residual SD: median = %.4f, range = [%.4f, %.4f]\n",
  median(residual_sd), min(residual_sd), max(residual_sd)
))

# =============================================================================
# STEP 3: Generate synthetic noise
# =============================================================================
# Draw iid N(0, residual_sd_j) noise for each protein j independently.
# This preserves the per-protein noise magnitude while breaking any
# individual-level covariance structure present in the tail PCs.

set.seed(42)
n_fibers   <- nrow(pca_res$x)
n_proteins <- nrow(pca_res$rotation)

noise_scaled <- matrix(
  rnorm(n_fibers * n_proteins),
  nrow = n_fibers,
  ncol = n_proteins
) * matrix(residual_sd, nrow = n_fibers, ncol = n_proteins, byrow = TRUE)

# =============================================================================
# STEP 4: Combine and unscale back to original intensity space
# =============================================================================

sim_scaled            <- reconstructed_scaled + noise_scaled
sim_fibers_x_proteins <- t(t(sim_scaled) * pca_res$scale + pca_res$center)

# Transpose to proteins × fibers — same layout as `intensities`
sim_intensities <- t(sim_fibers_x_proteins)
rownames(sim_intensities) <- rownames(intensities)
colnames(sim_intensities) <- colnames(intensities)

cat(sprintf(
  "Simulated intensity matrix: %d proteins × %d fibers\n",
  nrow(sim_intensities), ncol(sim_intensities)
))

# =============================================================================
# STEP 5: Build simulated metadata
# =============================================================================
# Preserve the group structure (time, fiber_type, fiber_id) but replace real
# subject IDs with synthetic donor labels so no participant can be identified.

real_subjects <- unique(metadata$subject_id)
fake_subjects <- paste0("donor_", seq_along(real_subjects))
subject_map   <- setNames(fake_subjects, real_subjects)

sim_metadata <- metadata |>
  mutate(subject_id = unname(subject_map[subject_id]))

cat("Subject ID mapping (real → simulated):\n")
print(data.frame(real = real_subjects, simulated = fake_subjects))

# =============================================================================
# STEP 6: Package in the same list structure as the real data
# =============================================================================

sim_fiber_data <- list(
  intensities = sim_intensities,
  metadata    = sim_metadata
)

cat("\nsim_fiber_data structure:\n")
str(sim_fiber_data, max.level = 1)

# =============================================================================
# STEP 7: Diagnostic plots
# =============================================================================
# Compare simulated vs real on four dimensions. All should look similar.
# If not, adjust K upward and re-run.

# 1. PCA scatter: should reproduce fiber-type separation and time gradient
sim_pca    <- prcomp(t(sim_intensities), scale. = TRUE)
sim_scores <- as_tibble(sim_pca$x, rownames = "fiber_id") |>
  left_join(sim_metadata |> select(fiber_id, time, subject_id), by = "fiber_id")

p_real <- ggplot(metadata, aes(PC1, PC2, colour = time)) +
  geom_point(size = 0.6, alpha = 0.6) +
  theme_minimal() +
  labs(title = "Real data — PC1 vs PC2")

p_sim <- ggplot(sim_scores, aes(PC1, PC2, colour = time)) +
  geom_point(size = 0.6, alpha = 0.6) +
  theme_minimal() +
  labs(title = sprintf("Simulated (K = %d) — PC1 vs PC2", K))

gridExtra::grid.arrange(p_real, p_sim, ncol = 2)

# 2. Variance explained: bars should track closely for PC1–10
var_real <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
var_sim  <- (sim_pca$sdev^2) / sum(sim_pca$sdev^2)

tibble(
  PC      = rep(1:20, 2),
  var_exp = c(var_real[1:20], var_sim[1:20]),
  source  = rep(c("Real", "Simulated"), each = 20)
) |>
  ggplot(aes(PC, var_exp, fill = source)) +
  geom_col(position = "dodge") +
  scale_x_continuous(breaks = 1:20) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Variance explained per PC: real vs simulated",
       y = "% variance", fill = NULL) |>
  print()

# 3. Marginal intensity distributions: density curves should overlap well
set.seed(7)
check_proteins <- rownames(intensities)[sample(nrow(intensities), 6)]

bind_rows(
  as_tibble(t(intensities[check_proteins, ])) |> mutate(source = "Real"),
  as_tibble(t(sim_intensities[check_proteins, ])) |> mutate(source = "Simulated")
) |>
  pivot_longer(-source, names_to = "protein", values_to = "intensity") |>
  ggplot(aes(intensity, colour = source)) +
  geom_density() +
  facet_wrap(~ protein, scales = "free") +
  theme_minimal() +
  labs(title = "Intensity distributions: real vs simulated", colour = NULL) |>
  print()

# 4. Protein-protein correlations: scatter should hug the y = x diagonal
#    Target: R² > 0.9 across the upper triangle
set.seed(1)
sample_prots <- sample(rownames(intensities), 100)
cor_real     <- cor(t(intensities[sample_prots, ]),     use = "pairwise.complete.obs")
cor_sim      <- cor(t(sim_intensities[sample_prots, ]), use = "pairwise.complete.obs")

cor_df <- tibble(
  real      = cor_real[upper.tri(cor_real)],
  simulated = cor_sim[upper.tri(cor_sim)]
)
r2 <- cor(cor_df$real, cor_df$simulated)^2
cat(sprintf("\nProtein-protein correlation R² (real vs simulated): %.3f\n", r2))
cat(if (r2 > 0.9) "✓ R² > 0.9 — structure well preserved.\n"
    else "⚠ R² ≤ 0.9 — consider increasing K.\n")

cor_df |>
  ggplot(aes(real, simulated)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_abline(colour = "firebrick", linetype = "dashed") +
  annotate("text", x = -0.5, y = 0.9, label = sprintf("R² = %.3f", r2), size = 4) +
  theme_minimal() +
  labs(title = "Protein-protein correlations: real vs simulated",
       x = "Real correlation", y = "Simulated correlation") |>
  print()

# =============================================================================
# STEP 8: Save (only after inspecting the diagnostics above)
# =============================================================================
# Acceptance criteria before saving:
#   ✓ PCA scatter shows same fiber-type separation as real data
#   ✓ Variance-explained bars track closely for PC1–10
#   ✓ Density curves overlap well for most proteins
#   ✓ Correlation R² > 0.9
#
# The output file uses a "_simulated" suffix so it does NOT overwrite the real
# data. Rename it manually (or update the path in the .qmd) once satisfied.

sim_fiber_data$metadata <- sim_fiber_data$metadata |>
  dplyr::select(!c(PC1, PC2))

out_path <- "data/single_fiber_proteomics_simulated.rds"
saveRDS(sim_fiber_data, out_path)
cat(sprintf("\nSaved: %s\n", out_path))
cat("Inspect the diagnostics above before replacing the real data file.\n")
cat("To use in the tutorial: rename to data/single_fiber_proteomics.rds\n")
