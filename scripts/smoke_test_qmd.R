# Mirrors the early Quarto chunks against the freshly exported real data to
# confirm the .qmd will render. Runs: load -> impute -> prcomp -> attach PC1
# to metadata -> flatten model_results -> grab a top hit and inspect.
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(tidyr) })

# ---- load ----
fiber_data  <- readRDS("data/single_fiber_proteomics.rds")
intensities <- fiber_data$intensities
metadata    <- fiber_data$metadata

cat("intensities:", dim(intensities), "\n")
cat("metadata cols:", paste(colnames(metadata), collapse = ", "), "\n")
cat("time levels:"); print(unique(metadata$time))
cat("subjects:"); print(unique(metadata$subject_id))
stopifnot(!"PC1" %in% colnames(metadata))

# ---- mean-impute ----
imp <- intensities
rm_ <- rowMeans(imp, na.rm = TRUE)
ix  <- which(is.na(imp), arr.ind = TRUE)
imp[ix] <- rm_[ix[, "row"]]

# ---- prcomp ----
pca_res <- prcomp(t(imp), scale. = TRUE)
cat("variance explained PC1-5: ",
    paste(round((pca_res$sdev^2 / sum(pca_res$sdev^2))[1:5] * 100, 1), collapse = " "),
    " %\n", sep = "")

# ---- attach PC1/PC2 to metadata ----
metadata <- metadata %>%
  left_join(as_tibble(pca_res$x, rownames = "fiber_id") %>%
              select(fiber_id, PC1, PC2),
            by = "fiber_id")
cat("metadata cols after attach:", paste(colnames(metadata), collapse = ", "), "\n")

# ---- flatten model_results ----
model_list <- readRDS("results/model_results.rds")
da_results <- bind_rows(lapply(model_list, `[[`, "model_results")) %>%
  mutate(
    adj_p_time        = p.adjust(`Pr(>|t|)_time1`,  method = "BH"),
    adj_p_fiber_type  = p.adjust(p_val_fiber_type,  method = "BH"),
    adj_p_interaction = p.adjust(p_val_interaction, method = "BH")
  )

cat("\nDA hits at adj_p_time < 0.05:        ", sum(da_results$adj_p_time < 0.05, na.rm = TRUE), "\n")
cat("Fiber-type hits at adj_p_fiber_type < 0.05:", sum(da_results$adj_p_fiber_type < 0.05, na.rm = TRUE), "\n")
cat("Interaction hits at adj_p_interaction < 0.05:", sum(da_results$adj_p_interaction < 0.05, na.rm = TRUE), "\n")

# ---- linkage spot-check: a top fiber-type hit vs metadata$PC1 ----
top_ft <- da_results %>% arrange(p_val_fiber_type) %>% slice(1) %>% pull(gene)
vals <- intensities[top_ft, ]
ok   <- !is.na(vals)
cat(sprintf("\nTop fiber-type gene: %s  | cor(intensity, PC1) = %.3f\n",
            top_ft, cor(vals[ok], metadata$PC1[ok])))

# ---- top time hit: PRE vs POST means ----
top_time <- da_results %>% arrange(`Pr(>|t|)_time1`) %>% slice(1) %>% pull(gene)
vals <- intensities[top_time, ]
m <- aggregate(vals, by = list(time = metadata$time), FUN = mean, na.rm = TRUE)
cat(sprintf("\nTop time gene: %s  | pre mean = %.2f  post mean = %.2f  delta = %.2f\n",
            top_time, m$x[m$time == "pre"], m$x[m$time == "post"],
            m$x[m$time == "post"] - m$x[m$time == "pre"]))
