# Detailed inspection — counts, distributions, design.
suppressPackageStartupMessages(library(dplyr))

path <- "C:/Users/jns822/Desktop/Scripts/REX_single_fiber/data/proteome/data_model_proteome_T1_T3.rds"
d <- readRDS(path)

cat("\n=== unique counts ===\n")
cat("subjects   :", dplyr::n_distinct(d$subject_id), "\n")
cat("sample_ids :", dplyr::n_distinct(d$sample_id), "\n")
cat("genes      :", dplyr::n_distinct(d$genes), "\n")
cat("times      :"); print(table(d$time))

cat("\n=== sample_id parsing ===\n")
parts <- do.call(rbind, strsplit(unique(d$sample_id), "_"))
cat("first 5 parsed sample_ids:\n"); print(head(parts, 5))
cat("unique tokens in 2nd field:\n"); print(table(parts[, 2]))

cat("\n=== fibers per subject ===\n")
print(summary(table(unique(d[, c("subject_id", "sample_id")])$subject_id)))

cat("\n=== fibers per subject x time ===\n")
fst <- d %>% distinct(subject_id, sample_id, time) %>%
  count(subject_id, time, name = "n_fibers")
print(summary(fst$n_fibers))

cat("\n=== intensities ===\n")
cat("NA fraction in long form:", mean(is.na(d$intensities)), "\n")
cat("quantiles:\n"); print(quantile(d$intensities, c(0, .01, .25, .5, .75, .99, 1), na.rm = TRUE))

cat("\n=== per-gene mean & sd (first 10 by mean) ===\n")
g <- d %>% group_by(genes) %>%
  summarise(mean_int = mean(intensities, na.rm = TRUE),
            sd_int   = sd(intensities, na.rm = TRUE),
            n_obs    = sum(!is.na(intensities)),
            .groups  = "drop") %>%
  arrange(desc(mean_int))
print(head(g, 10))
cat("\n=== per-gene n_obs distribution ===\n")
print(summary(g$n_obs))
cat("\n=== per-gene sd distribution ===\n")
print(summary(g$sd_int))

cat("\n=== PC1 ===\n")
cat("quantiles:\n"); print(quantile(d$PC1, c(0, .25, .5, .75, 1), na.rm = TRUE))
cat("PC1 is per-fiber or per-row? (range of unique PC1 per sample_id)\n")
pc1_per_fiber <- d %>% group_by(sample_id) %>%
  summarise(n_unique_PC1 = dplyr::n_distinct(PC1), .groups = "drop")
print(table(pc1_per_fiber$n_unique_PC1))
