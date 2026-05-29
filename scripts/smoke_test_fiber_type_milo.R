# Walk through:
#  - prcomp PCA -> inject into Seurat -> RunUMAP
#  - FindClusters -> MYH7/MYH2 -> fiber_type label
#  - Build Milo with sample = subject x time x fiber_type
#  - testNhoods for main time effect AND the fiber_type:time interaction
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

fd <- readRDS("data/single_fiber_proteomics.rds")
intensities <- fd$intensities
metadata    <- fd$metadata

stopifnot("MYH7" %in% rownames(intensities), "MYH2" %in% rownames(intensities))

# --- PCA (Section 1.1) ---
imp <- intensities
rm_ <- rowMeans(imp, na.rm = TRUE)
ix  <- which(is.na(imp), arr.ind = TRUE)
imp[ix] <- rm_[ix[, "row"]]
pca_res <- prcomp(t(imp), scale. = TRUE)

# --- Seurat + UMAP (Section 1.3) with injected prcomp PCs ---
md_seu <- metadata %>% column_to_rownames("fiber_id")
seu <- Seurat::CreateSeuratObject(counts = imp, meta.data = md_seu)
Seurat::VariableFeatures(seu) <- rownames(seu)
seu[["RNA"]]$data <- seu[["RNA"]]$counts

pc_emb <- pca_res$x[, 1:30]
pc_ldg <- pca_res$rotation[, 1:30]
colnames(pc_emb) <- paste0("PC_", seq_len(ncol(pc_emb)))
colnames(pc_ldg) <- paste0("PC_", seq_len(ncol(pc_ldg)))
seu[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pc_emb[colnames(seu), ],
  loadings   = pc_ldg,
  assay      = Seurat::DefaultAssay(seu),
  key        = "PC_"
)

set.seed(42)
seu <- Seurat::RunUMAP(seu, dims = 1:20, verbose = FALSE)
seu <- Seurat::FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- Seurat::FindClusters(seu, resolution = 0.5, verbose = FALSE)

cat("\n=== Clusters ===\n")
print(table(seu$seurat_clusters))

# --- Fiber-type assignment (Section 1.4) ---
cluster_myh <- tibble(
    fiber_id = colnames(seu),
    cluster  = as.character(seu$seurat_clusters),
    MYH7     = intensities["MYH7", colnames(seu)],
    MYH2     = intensities["MYH2", colnames(seu)]
  ) %>%
  group_by(cluster) %>%
  summarise(MYH7 = mean(MYH7, na.rm = TRUE),
            MYH2 = mean(MYH2, na.rm = TRUE),
            n    = dplyr::n(),
            .groups = "drop") %>%
  mutate(fiber_type = ifelse(MYH2 > MYH7, "fast", "slow"))

cat("\n=== Cluster -> fiber_type ===\n"); print(cluster_myh)

fiber_type_map <- setNames(cluster_myh$fiber_type, cluster_myh$cluster)
seu$fiber_type <- unname(fiber_type_map[as.character(seu$seurat_clusters)])
metadata <- metadata %>%
  left_join(tibble(fiber_id = colnames(seu), fiber_type = seu$fiber_type),
            by = "fiber_id")

cat("\n=== fiber_type distribution ===\n")
print(table(metadata$fiber_type, metadata$time))

# --- MILO with subject x time x fiber_type as the sample ---
md <- metadata %>%
  mutate(sample_id = paste(subject_id, time, fiber_type, sep = "_")) %>%
  as.data.frame()
rownames(md) <- md$fiber_id

imp_mat <- as.matrix(imp)[, md$fiber_id]
sce <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(logcounts = imp_mat),
  colData = md
)
SingleCellExperiment::reducedDim(sce, "PCA")  <- pca_res$x[colnames(sce), 1:30]
SingleCellExperiment::reducedDim(sce, "UMAP") <- Seurat::Embeddings(seu, "umap")[colnames(sce), ]

milo_obj <- miloR::Milo(sce)
milo_obj <- miloR::buildGraph(milo_obj, k = 30, d = 20, reduced.dim = "PCA")
milo_obj <- miloR::makeNhoods(milo_obj, prop = 0.1, k = 30, d = 20,
                              refined = TRUE, reduced_dims = "PCA")
milo_obj <- miloR::countCells(milo_obj, meta.data = md, sample = "sample_id")
milo_obj <- miloR::buildNhoodGraph(milo_obj)

design <- md %>%
  distinct(sample_id, subject_id, time, fiber_type) %>%
  as.data.frame()
rownames(design) <- design$sample_id
design$time       <- factor(design$time,       levels = c("pre",  "post"))
design$fiber_type <- factor(design$fiber_type, levels = c("fast", "slow"))

cat("\n=== design ===\n"); print(design)

cat("\n=== Running main time effect ===\n")
da_time <- miloR::testNhoods(milo_obj,
                             design    = ~ subject_id + fiber_type + time,
                             design.df = design)
cat("DA neighbourhoods (SpatialFDR < 0.1):",
    sum(da_time$SpatialFDR < 0.1, na.rm = TRUE), "/", nrow(da_time), "\n")

cat("\n=== Running stratified contrasts per fiber type ===\n")
design$group <- factor(
  paste(design$fiber_type, design$time, sep = "_"),
  levels = c("fast_pre", "fast_post", "slow_pre", "slow_post")
)

da_fast <- miloR::testNhoods(
  milo_obj,
  design          = ~ 0 + group + subject_id,
  design.df       = design,
  model.contrasts = "groupfast_post - groupfast_pre"
)

da_slow <- miloR::testNhoods(
  milo_obj,
  design          = ~ 0 + group + subject_id,
  design.df       = design,
  model.contrasts = "groupslow_post - groupslow_pre"
)

da_counts <- tibble(
  fiber_type = c("fast", "slow"),
  n_total    = c(nrow(da_fast), nrow(da_slow)),
  n_sig      = c(sum(da_fast$SpatialFDR < 0.1, na.rm = TRUE),
                 sum(da_slow$SpatialFDR < 0.1, na.rm = TRUE)),
  n_up       = c(sum(da_fast$SpatialFDR < 0.1 & da_fast$logFC > 0, na.rm = TRUE),
                 sum(da_slow$SpatialFDR < 0.1 & da_slow$logFC > 0, na.rm = TRUE)),
  n_down     = c(sum(da_fast$SpatialFDR < 0.1 & da_fast$logFC < 0, na.rm = TRUE),
                 sum(da_slow$SpatialFDR < 0.1 & da_slow$logFC < 0, na.rm = TRUE))
)
cat("\nPer-fiber-type significant nhood counts:\n")
print(da_counts)

cat("\nAll steps completed cleanly.\n")
