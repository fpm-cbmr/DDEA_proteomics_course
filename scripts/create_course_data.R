# =============================================================================
# create_course_data.R
# -----------------------------------------------------------------------------
# Build the course-ready dataset and model-results files from the REX
# single-fiber proteome. The only transformations applied to the real data are:
#   - drop the precomputed PC1 column (PC1 is recomputed in the Quarto doc)
#   - pivot the long-form tibble into a (proteins x fibers) matrix
#   - assemble a per-fiber metadata table
#   - anonymise subject identifiers (REX## -> S##) in every place they appear
#
# Outputs (relative to project root):
#   data/single_fiber_proteomics.rds   list(intensities, metadata)
#   results/model_results.rds          named list, one entry per gene
#                                       (passthrough of the real GAM results;
#                                        no subject identifiers exist inside)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
})

# ---- 0. Paths ---------------------------------------------------------------

REAL_INTENS  <- "C:/Users/jns822/Desktop/Scripts/REX_single_fiber/data/proteome/data_model_proteome_T1_T3.rds"
REAL_MODEL   <- "C:/Users/jns822/Desktop/Scripts/REX_single_fiber/data/proteome/proteome_gam_interaction_results_T1_T3.rds"
OUT_DATA     <- "data/single_fiber_proteomics.rds"
OUT_MODEL    <- "results/model_results.rds"

# ---- 1. Load real long-form intensities, drop PC1 ---------------------------

cat("Loading real long-form intensities...\n")
long <- readRDS(REAL_INTENS) |>
  select(-PC1)   # PC1 is recomputed downstream in the Quarto doc

stopifnot(all(c("genes", "sample_id", "intensities",
                "subject_id", "time") %in% names(long)))

cat(sprintf("  rows: %d   genes: %d   sample_ids: %d   subjects: %d\n",
            nrow(long),
            dplyr::n_distinct(long$genes),
            dplyr::n_distinct(long$sample_id),
            dplyr::n_distinct(long$subject_id)))

# ---- 2. Build subject anonymisation map (REX## -> S##) ----------------------

real_subjects <- sort(unique(long$subject_id))
subj_map <- setNames(sprintf("S%02d", seq_along(real_subjects)), real_subjects)

cat("\nSubject anonymisation map:\n")
print(subj_map)

# Apply the mapping inside sample_id strings too. Real sample_id format is
# "<REX##>_<PRE|POST>_<fiber_num>", so we only touch the leading subject token.
anonymise_sample_id <- function(s) {
  parts <- str_split_fixed(s, "_", 3)
  paste(subj_map[parts[, 1]], parts[, 2], parts[, 3], sep = "_")
}

long_anon <- long |>
  mutate(
    subject_id = unname(subj_map[subject_id]),
    sample_id  = anonymise_sample_id(sample_id)
  )

# Sanity: no REX identifier survives anywhere
stopifnot(!any(grepl("REX", long_anon$sample_id)))
stopifnot(!any(grepl("REX", long_anon$subject_id)))

# ---- 3. Pivot to (proteins x fibers) intensity matrix -----------------------

cat("\nPivoting long form to (proteins x fibers) matrix...\n")
genes  <- sort(unique(long_anon$genes))
fibers <- sort(unique(long_anon$sample_id))

intensities <- matrix(
  NA_real_,
  nrow = length(genes), ncol = length(fibers),
  dimnames = list(genes, fibers)
)
intensities[cbind(match(long_anon$genes,     genes),
                  match(long_anon$sample_id, fibers))] <- long_anon$intensities

cat(sprintf("  intensities matrix: %d x %d  (NA fraction: %.3f)\n",
            nrow(intensities), ncol(intensities), mean(is.na(intensities))))

# ---- 4. Per-fiber metadata --------------------------------------------------

metadata <- long_anon |>
  distinct(sample_id, subject_id, time) |>
  rename(fiber_id = sample_id) |>
  mutate(
    fiber_num = as.integer(str_split_fixed(fiber_id, "_", 3)[, 3]),
    # T1 -> pre, T3 -> post (matches the PRE/POST tokens already inside sample_id)
    time      = factor(recode(time, T1 = "pre", T3 = "post"),
                       levels = c("pre", "post"))
  ) |>
  select(fiber_id, subject_id, time, fiber_num) |>
  arrange(subject_id, time, fiber_num)

# Re-order the columns of `intensities` to match metadata$fiber_id row order
intensities <- intensities[, metadata$fiber_id]

stopifnot(identical(colnames(intensities), metadata$fiber_id))
stopifnot(nrow(metadata) == ncol(intensities))

cat(sprintf("  metadata: %d fibers x %d columns (%s)\n",
            nrow(metadata), ncol(metadata),
            paste(colnames(metadata), collapse = ", ")))
cat("  fibers per (subject, time):\n")
print(metadata |> count(subject_id, time) |> tidyr::pivot_wider(
  names_from = time, values_from = n
))

# ---- 5. Load real model results (passthrough) -------------------------------

cat("\nLoading real model results...\n")
model_results <- readRDS(REAL_MODEL)
cat(sprintf("  %d genes in model_results\n", length(model_results)))

# Sanity: nothing inside references REX subject IDs (the structure is per-gene
# stats; subject info would only appear in raw fits, not in the summary tibbles).
flat_check <- bind_rows(lapply(model_results, `[[`, "model_results"))
stopifnot(!any(grepl("REX", unlist(flat_check), useBytes = TRUE)))

# ---- 6. Save ----------------------------------------------------------------

dir.create(dirname(OUT_DATA),  showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_MODEL), showWarnings = FALSE, recursive = TRUE)

fiber_data <- list(
  intensities = intensities,
  metadata    = metadata
)

saveRDS(fiber_data,    OUT_DATA)
saveRDS(model_results, OUT_MODEL)

cat("\nWrote:\n  ", OUT_DATA, "\n  ", OUT_MODEL, "\n")
cat("Done.\n")
