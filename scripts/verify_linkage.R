# Confirms that the model_results faithfully describe the intensity matrix:
#   - top fiber-type hits show a clear intensity vs PC1 correlation
#   - top time hits show a PRE vs POST mean shift
#   - top interaction hits show a PC1-slope difference between PRE and POST
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

fd  <- readRDS("data/single_fiber_proteomics.rds")
mr  <- readRDS("results/model_results.rds")

flat <- bind_rows(lapply(mr, `[[`, "model_results"))
cat("flat results rows:", nrow(flat), "\n")
cat("columns:\n"); print(colnames(flat))

I <- fd$intensities
M <- fd$metadata

show_correlation <- function(gene) {
  vals <- I[gene, ]
  ok   <- !is.na(vals)
  cat(sprintf("  cor(intensity, PC1) = %.3f  (n=%d)\n",
              cor(vals[ok], M$PC1[ok]), sum(ok)))
}

show_time_shift <- function(gene) {
  vals <- I[gene, ]
  d <- data.frame(intensity = vals, time = M$time)
  m <- d %>% group_by(time) %>%
    summarise(mean_int = mean(intensity, na.rm = TRUE), .groups = "drop")
  cat(sprintf("  mean(PRE)=%.3f  mean(POST)=%.3f  delta=%.3f\n",
              m$mean_int[m$time == "PRE"], m$mean_int[m$time == "POST"],
              m$mean_int[m$time == "POST"] - m$mean_int[m$time == "PRE"]))
}

show_interaction <- function(gene) {
  vals <- I[gene, ]; ok <- !is.na(vals)
  pre  <- M$time == "PRE"  & ok
  post <- M$time == "POST" & ok
  s_pre  <- coef(lm(vals[pre]  ~ M$PC1[pre]))[2]
  s_post <- coef(lm(vals[post] ~ M$PC1[post]))[2]
  cat(sprintf("  PC1 slope PRE=%.4f  POST=%.4f  diff=%.4f\n",
              s_pre, s_post, s_post - s_pre))
}

cat("\n=== top 5 fiber-type hits ===\n")
top_ft <- flat %>% arrange(p_val_fiber_type) %>% head(5)
for (g in top_ft$gene) {
  cat(g, "  p_val_fiber_type =", signif(flat$p_val_fiber_type[flat$gene == g], 3), "\n")
  show_correlation(g)
}

cat("\n=== bottom 5 fiber-type hits (should have low |cor|) ===\n")
bot_ft <- flat %>% arrange(desc(p_val_fiber_type)) %>% head(5)
for (g in bot_ft$gene) {
  cat(g, "  p_val_fiber_type =", signif(flat$p_val_fiber_type[flat$gene == g], 3), "\n")
  show_correlation(g)
}

cat("\n=== top 5 time hits ===\n")
top_t <- flat %>% arrange(`Pr(>|t|)_time1`) %>% head(5)
for (g in top_t$gene) {
  cat(g, "  Estimate_time1 =", signif(flat$Estimate_time1[flat$gene == g], 3),
      "  p =", signif(flat$`Pr(>|t|)_time1`[flat$gene == g], 3), "\n")
  show_time_shift(g)
}

cat("\n=== top 5 interaction hits ===\n")
top_i <- flat %>% arrange(p_val_interaction) %>% head(5)
for (g in top_i$gene) {
  cat(g, "  p_val_interaction =", signif(flat$p_val_interaction[flat$gene == g], 3), "\n")
  show_interaction(g)
}
