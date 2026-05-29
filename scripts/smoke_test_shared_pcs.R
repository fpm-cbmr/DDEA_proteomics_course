# Mirrors the new section 1.3 -> 1.4 chain on the real data:
#  - prcomp PCA once
#  - inject into Seurat, RunUMAP with seed 42
#  - inject same PCs + Seurat's UMAP into SCE
#  - confirm the two reducedDims/Embeddings line up

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr)
})

fd <- readRDS("data/single_fiber_proteomics.rds")
intensities <- fd$intensities
metadata    <- fd$metadata

# --- mean-impute and prcomp (Section 1.1) ---
imp <- intensities
rm_ <- rowMeans(imp, na.rm = TRUE)
ix  <- which(is.na(imp), arr.ind = TRUE)
imp[ix] <- rm_[ix[, "row"]]

pca_res <- prcomp(t(imp), scale. = TRUE)

# --- Seurat: inject prcomp PCs, then RunUMAP (Section 1.3) ---
md <- metadata %>% column_to_rownames("fiber_id")
seu <- Seurat::CreateSeuratObject(counts = imp, meta.data = md)
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

seurat_umap <- Seurat::Embeddings(seu, "umap")
cat("Seurat UMAP dim:", dim(seurat_umap), "\n")
cat("Seurat UMAP first 3 rows:\n"); print(head(seurat_umap, 3))

# --- MILO: inject prcomp PCs + reuse Seurat UMAP (Section 1.4) ---
md_milo <- metadata %>%
  mutate(sample_id = paste(subject_id, time, sep = "_")) %>%
  as.data.frame()
rownames(md_milo) <- md_milo$fiber_id

imp_mat <- as.matrix(imp)[, md_milo$fiber_id]

sce <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(logcounts = imp_mat),
  colData = md_milo
)

SingleCellExperiment::reducedDim(sce, "PCA")  <- pca_res$x[colnames(sce), 1:30]
SingleCellExperiment::reducedDim(sce, "UMAP") <-
  Seurat::Embeddings(seu, "umap")[colnames(sce), ]

sce_umap <- SingleCellExperiment::reducedDim(sce, "UMAP")
cat("\nSCE UMAP dim:", dim(sce_umap), "\n")
cat("SCE UMAP first 3 rows:\n"); print(head(sce_umap, 3))

cat("\nMax abs diff between Seurat and SCE UMAP coords: ",
    max(abs(seurat_umap[rownames(sce_umap), ] - sce_umap)), "\n", sep = "")

# Quick build + count + buildNhoodGraph to confirm MILO chain still works
milo_obj <- miloR::Milo(sce)
milo_obj <- miloR::buildGraph(milo_obj, k = 30, d = 20, reduced.dim = "PCA")
milo_obj <- miloR::makeNhoods(milo_obj, prop = 0.1, k = 30, d = 20,
                              refined = TRUE, reduced_dims = "PCA")
milo_obj <- miloR::countCells(milo_obj, meta.data = md_milo, sample = "sample_id")
milo_obj <- miloR::buildNhoodGraph(milo_obj)

cat("\nMILO neighbourhoods:", ncol(miloR::nhoods(milo_obj)), "\n")
cat("All steps completed cleanly.\n")
