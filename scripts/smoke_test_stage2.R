# Walk through Stage 2 — confirms model_results loads and gives sensible
# numbers for exercise answers.
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

fd <- readRDS("data/single_fiber_proteomics.rds")
intensities <- fd$intensities
metadata    <- fd$metadata

model_list <- readRDS("results/model_results.rds")
cat("model_list length:", length(model_list), "\n")
cat("entries per gene: "); print(names(model_list[[1]]))
cat("model_results columns:\n"); print(colnames(model_list[[1]]$model_results))

da_results <- bind_rows(lapply(model_list, `[[`, "model_results")) %>%
  mutate(
    adj_p_time        = p.adjust(`Pr(>|t|)_time1`,  method = "BH"),
    adj_p_fiber_type  = p.adjust(p_val_fiber_type,  method = "BH"),
    adj_p_interaction = p.adjust(p_val_interaction, method = "BH")
  )

cat("\n=== Hit counts ===\n")
print(da_results %>%
  summarise(
    n_total           = dplyr::n(),
    n_time            = sum(adj_p_time        < 0.05, na.rm = TRUE),
    n_fiber_type      = sum(adj_p_fiber_type  < 0.05, na.rm = TRUE),
    n_interaction     = sum(adj_p_interaction < 0.05, na.rm = TRUE),
    n_time_up         = sum(adj_p_time < 0.05 & Estimate_time1 > 0, na.rm = TRUE),
    n_time_down       = sum(adj_p_time < 0.05 & Estimate_time1 < 0, na.rm = TRUE)
  ))

cat("\n=== Top 5 time hits ===\n")
print(da_results %>% arrange(`Pr(>|t|)_time1`) %>%
  select(gene, Estimate_time1, `Pr(>|t|)_time1`, adj_p_time) %>% head(5))

cat("\n=== Top 5 fiber_type hits ===\n")
print(da_results %>% arrange(p_val_fiber_type) %>%
  select(gene, p_val_fiber_type, adj_p_fiber_type) %>% head(5))

cat("\n=== Top 5 interaction hits ===\n")
print(da_results %>% arrange(p_val_interaction) %>%
  select(gene, p_val_interaction, adj_p_interaction) %>% head(5))

cat("\n=== A NON-significant time protein ===\n")
print(da_results %>% filter(adj_p_time > 0.5) %>% slice(1) %>%
  select(gene, Estimate_time1, adj_p_time))

cat("\n=== A NON-significant fiber_type protein ===\n")
print(da_results %>% filter(adj_p_fiber_type > 0.5) %>% slice(1) %>%
  select(gene, p_val_fiber_type, adj_p_fiber_type))

# Verify a top time hit really shows a pre/post shift in the matrix
top_time <- da_results %>% arrange(`Pr(>|t|)_time1`) %>% slice(1) %>% pull(gene)
vals <- intensities[top_time, ]
m_pre  <- mean(vals[metadata$time == "pre"],  na.rm = TRUE)
m_post <- mean(vals[metadata$time == "post"], na.rm = TRUE)
cat(sprintf("\nTop time hit %s: pre mean = %.2f, post mean = %.2f, delta = %.2f\n",
            top_time, m_pre, m_post, m_post - m_pre))
