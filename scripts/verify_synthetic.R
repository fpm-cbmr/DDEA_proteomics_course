# Quick sanity check on the simulated outputs.
fd <- readRDS("data/single_fiber_proteomics.rds")
da <- readRDS("results/da_results.rds")

cat("\n=== fiber_data list ===\n")
cat("names:"); print(names(fd))
cat("intensities dim:", dim(fd$intensities), "\n")
cat("metadata  rows:", nrow(fd$metadata), "  cols:", paste(colnames(fd$metadata), collapse = ", "), "\n")

cat("\n=== anonymity check ===\n")
cat("unique subject_ids:\n"); print(unique(fd$metadata$subject_id))
cat("any sample id contains 'REX'? ", any(grepl("REX", fd$metadata$fiber_id)), "\n", sep = "")

cat("\n=== design check ===\n")
print(table(fd$metadata$subject_id, fd$metadata$time))

cat("\n=== intensity distribution ===\n")
cat("NA fraction:", mean(is.na(fd$intensities)), "\n")
print(quantile(fd$intensities, c(0, .01, .25, .5, .75, .99, 1), na.rm = TRUE))

cat("\n=== PC1 in metadata ===\n")
print(summary(fd$metadata$PC1))

cat("\n=== top DA hits ===\n")
print(head(da, 12))

cat("\n=== DA result distribution ===\n")
cat("n significant adj_p<0.05:", sum(da$adj_p_value < 0.05, na.rm = TRUE), "\n")
cat("n with |log2FC|>0.5 and adj_p<0.05:",
    sum(abs(da$log2FC) > 0.5 & da$adj_p_value < 0.05, na.rm = TRUE), "\n")
