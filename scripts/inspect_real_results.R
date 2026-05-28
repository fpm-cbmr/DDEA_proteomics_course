path <- "C:/Users/jns822/Desktop/Scripts/REX_single_fiber/data/proteome/proteome_gam_interaction_results_T1_T3.rds"
obj  <- readRDS(path)

cat("\n=== top-level ===\n")
cat("class:"); print(class(obj))
cat("length:", length(obj), "\n")
cat("first 10 names:\n"); print(head(names(obj), 10))

cat("\n=== first entry structure (depth 2) ===\n")
str(obj[[1]], max.level = 2)

cat("\n=== names of entries inside one gene ===\n")
print(names(obj[[1]]))

cat("\n=== model_results for first 3 genes ===\n")
for (g in head(names(obj), 3)) {
  cat(sprintf("\n--- %s ---\n", g))
  mr <- obj[[g]]$model_results
  cat("class:"); print(class(mr))
  if (is.data.frame(mr)) {
    cat("dim:", dim(mr), "\n")
    cat("colnames:\n"); print(colnames(mr))
    print(mr)
  } else {
    str(mr, max.level = 2)
  }
}

cat("\n=== aggregate column types across all genes ===\n")
# Pick one model_results to learn the schema
mr1 <- obj[[1]]$model_results
if (is.data.frame(mr1)) {
  cat("column types:\n"); print(sapply(mr1, class))
}

cat("\n=== distribution of estimate_time across all genes (if present) ===\n")
all_cols <- names(obj[[1]]$model_results)
cat("Checking column:", all_cols, "\n")
if ("estimate_time" %in% all_cols) {
  vals <- sapply(obj, function(x) x$model_results$estimate_time[1])
  cat("estimate_time summary:\n"); print(summary(vals))
}
if ("p_val_fiber_type" %in% all_cols) {
  vals <- sapply(obj, function(x) x$model_results$p_val_fiber_type[1])
  cat("p_val_fiber_type summary:\n"); print(summary(vals))
  cat("n significant <0.05:", sum(vals < 0.05, na.rm = TRUE), "\n")
}
if ("p_val_interaction" %in% all_cols) {
  vals <- sapply(obj, function(x) x$model_results$p_val_interaction[1])
  cat("p_val_interaction summary:\n"); print(summary(vals))
  cat("n significant <0.05:", sum(vals < 0.05, na.rm = TRUE), "\n")
}
