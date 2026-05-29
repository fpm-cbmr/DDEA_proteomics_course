# =============================================================================
# qmd_to_colab.R
# -----------------------------------------------------------------------------
# Convert a Quarto .qmd into a Google Colab-ready R notebook (.ipynb).
#
# Splits the .qmd into alternating markdown cells and R code cells; sets the
# kernelspec to "ir" (Colab's R kernel); and prepends a setup cell that
# installs every package referenced anywhere in the document and clones the
# course data from the user's GitHub fork.
#
# Usage:
#   Rscript scripts/qmd_to_colab.R \
#     [input.qmd] [output.ipynb] [github_user/repo]
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(stringr)
})

args        <- commandArgs(trailingOnly = TRUE)
qmd_path    <- if (length(args) >= 1) args[1] else "proteomics_data_visualization.qmd"
out_path    <- if (length(args) >= 2) args[2] else "proteomics_data_visualization.ipynb"
gh_repo     <- if (length(args) >= 3) args[3] else "fpm_cbmr/DDEA_proteomics_course"

stopifnot(file.exists(qmd_path))

# ---- 1. Parse .qmd into a sequence of (type, source) cells ----------------

lines <- readLines(qmd_path, warn = FALSE)

# Strip the YAML front-matter (between leading --- and the next ---)
if (length(lines) > 0 && lines[1] == "---") {
  close_idx <- which(lines == "---")[2]
  lines <- lines[seq.int(close_idx + 1, length(lines))]
}

cells   <- list()
buf_md  <- character()
buf_r   <- character()
in_code <- FALSE

flush_md <- function() {
  if (length(buf_md) == 0) return()
  src <- paste(buf_md, collapse = "\n")
  # Skip cells that are only whitespace / horizontal rules
  if (nchar(trimws(gsub("-", "", src))) == 0) {
    buf_md <<- character()
    return()
  }
  cells[[length(cells) + 1]] <<- list(
    cell_type = "markdown",
    metadata  = setNames(list(), character()),
    source    = strsplit(paste0(src, "\n"), "(?<=\n)", perl = TRUE)[[1]]
  )
  buf_md <<- character()
}

flush_r <- function() {
  if (length(buf_r) == 0) return()
  # Drop Quarto chunk-options lines (`#| ...`) — Colab doesn't read them.
  body <- buf_r[!grepl("^#\\|", buf_r)]
  src  <- paste(body, collapse = "\n")
  if (nchar(trimws(src)) == 0) {
    buf_r <<- character()
    return()
  }
  cells[[length(cells) + 1]] <<- list(
    cell_type       = "code",
    metadata        = setNames(list(), character()),
    execution_count = NULL,
    outputs         = list(),
    source          = strsplit(paste0(src, "\n"), "(?<=\n)", perl = TRUE)[[1]]
  )
  buf_r <<- character()
}

for (ln in lines) {
  if (!in_code && grepl("^```\\{r", ln)) {
    flush_md()
    in_code <- TRUE
    next
  }
  if (in_code && grepl("^```\\s*$", ln)) {
    flush_r()
    in_code <- FALSE
    next
  }
  if (in_code) buf_r <- c(buf_r, ln) else buf_md <- c(buf_md, ln)
}
flush_md(); flush_r()

cat(sprintf("Parsed %d cells (%d markdown, %d code).\n",
            length(cells),
            sum(sapply(cells, function(c) c$cell_type == "markdown")),
            sum(sapply(cells, function(c) c$cell_type == "code"))))

# ---- 2. Discover packages used in the document -----------------------------
# Anything that looks like `pkg::fn(` or `library(pkg)` is treated as a dep.

all_r_src <- unlist(lapply(cells[sapply(cells, function(c) c$cell_type == "code")],
                            function(c) c$source))
text <- paste(all_r_src, collapse = "\n")

ns_pkgs  <- unique(regmatches(text, gregexpr("\\b[A-Za-z][A-Za-z0-9.]*(?=::)",
                                              text, perl = TRUE))[[1]])
lib_pkgs <- unique(regmatches(text,
                              gregexpr("(?<=library\\()[A-Za-z][A-Za-z0-9.]*",
                                       text, perl = TRUE))[[1]])

pkgs <- sort(unique(c(ns_pkgs, lib_pkgs)))
cat("Packages discovered:\n"); print(pkgs)

# Manually split into CRAN vs Bioconductor (override list)
bioc_pkgs <- intersect(pkgs, c(
  "clusterProfiler", "org.Hs.eg.db", "limma", "miloR",
  "SingleCellExperiment", "PhosR", "EnhancedVolcano",
  "ComplexHeatmap", "fgsea", "msigdbr", "scater", "BiocParallel"
))
cran_pkgs <- setdiff(pkgs, bioc_pkgs)
# Drop base / always-available packages from the CRAN list
cran_pkgs <- setdiff(cran_pkgs, c("base", "stats", "utils", "grDevices",
                                   "graphics", "methods", "ggplot2"))

# ---- 3. Build a Colab setup cell -------------------------------------------

setup_lines <- c(
  "# === Google Colab setup (run this cell first) ===",
  "# Installs every package the rest of the notebook needs, then clones",
  "# the course repo so the data and scripts are available.",
  "",
  "# 1) CRAN packages",
  sprintf("cran_pkgs <- c(%s)",
          paste(sprintf('\"%s\"', cran_pkgs), collapse = ", ")),
  "missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))",
  "if (length(missing_cran)) install.packages(missing_cran, Ncpus = 4)",
  "",
  "# 2) Bioconductor packages",
  "if (!requireNamespace(\"BiocManager\", quietly = TRUE))",
  "  install.packages(\"BiocManager\")",
  sprintf("bioc_pkgs <- c(%s)",
          paste(sprintf('\"%s\"', bioc_pkgs), collapse = ", ")),
  "missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))",
  "if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)",
  "",
  "# 3) Clone the course repo so data/ and results/ are available",
  sprintf("repo_dir <- \"%s\"", basename(gh_repo)),
  "if (!dir.exists(repo_dir)) {",
  sprintf("  system(\"git clone --depth 1 https://github.com/%s.git\")",
          gh_repo),
  "}",
  "setwd(repo_dir)",
  "",
  "cat(\"Working directory:\", getwd(), \"\\n\")",
  "cat(\"Contents of data/:\\n\"); print(list.files(\"data\"))",
  "cat(\"Contents of results/:\\n\"); print(list.files(\"results\"))"
)

setup_cell <- list(
  cell_type       = "code",
  metadata        = setNames(list(), character()),
  execution_count = NULL,
  outputs         = list(),
  source          = strsplit(paste0(paste(setup_lines, collapse = "\n"), "\n"),
                             "(?<=\n)", perl = TRUE)[[1]]
)

intro_md <- list(
  cell_type = "markdown",
  metadata  = setNames(list(), character()),
  source    = strsplit(paste0(
    "# Proteomics Data Visualization — Google Colab\n",
    "\n",
    "**DDEA Proteomics Course — Hands-On Tutorial**\n",
    "\n",
    "This is the Google Colab version of `proteomics_data_visualization.qmd`. ",
    "Run the setup cell below once at the start of the session to install ",
    "all dependencies and pull the course data.\n",
    "\n",
    "**Make sure your runtime is set to R**: ",
    "Runtime → Change runtime type → Runtime type: R.\n"
  ), "(?<=\n)", perl = TRUE)[[1]]
)

cells <- c(list(intro_md, setup_cell), cells)
cat(sprintf("Final notebook has %d cells.\n", length(cells)))

# ---- 4. Assemble the .ipynb JSON -------------------------------------------

nb <- list(
  cells = cells,
  metadata = list(
    kernelspec = list(
      display_name = "R",
      language     = "R",
      name         = "ir"
    ),
    language_info = list(
      codemirror_mode = "r",
      file_extension  = ".r",
      mimetype        = "text/x-r-source",
      name            = "R",
      pygments_lexer  = "r",
      version         = paste(R.version$major, R.version$minor, sep = ".")
    ),
    colab = list(
      provenance = list(),
      toc_visible = TRUE
    )
  ),
  nbformat       = 4L,
  nbformat_minor = 4L
)

writeLines(
  jsonlite::toJSON(nb, auto_unbox = TRUE, pretty = 2, null = "null"),
  out_path
)
cat(sprintf("\nWrote: %s\n", out_path))
