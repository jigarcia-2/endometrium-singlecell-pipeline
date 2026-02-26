############################################################
# Endometrial Organoid scRNA-seq (EEOs) - SeuratDisk Pipeline
# Based on García-Alonso et al., 2021
#
# Validated environment:
# - Windows, R 4.3.3 (Rtools43)
# - Seurat 5.2.1 / SeuratObject 5.0.2
# - SeuratDisk 0.0.0.9021 (GitHub)
#
# Note (Windows): Seurat plots require ggplot2 >= 3.5.2 & < 4.0.0
# (see scr/00_fix_windows_packages.R if needed)
############################################################

## 0) Reproducibility and global options ----
set.seed(1234)
options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(pkgType = "binary")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

## 1) Project paths (fixed) ----
base_dir <- getwd()

data_raw_dir       <- file.path(base_dir, "data", "raw")
fig_dir            <- file.path(base_dir, "outputs", "figures")
table_dir          <- file.path(base_dir, "outputs", "tables")
obj_dir            <- file.path(base_dir, "outputs", "objects")

dir.create(data_raw_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,            recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,          recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir,            recursive = TRUE, showWarnings = FALSE)

## 2) Required packages ----
# Install required packages if not already available.
# ggplot2 version is validated separately (see section 3b).

cran_pkgs <- c(
  "Seurat", "SeuratObject", "dplyr", "patchwork",
  "hdf5r")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, type = "binary")
  }
}

# SeuratDisk from GitHub (only if missing)
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", type = "binary")
  }
  remotes::install_github("mojaveazure/seurat-disk", upgrade = "never")
}

## 3)3a- Load libraries ----
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(ggplot2)
library(hdf5r)

##3b Version check: ggplot2 compatibility (Windows safeguard) ----
# On Windows (R 4.3.3), Seurat plotting may be unstable with ggplot2 >= 4.0.0.
# This check prevents silent crashes by enforcing a supported version range.

ggv <- packageVersion("ggplot2")
message("ggplot2:      ", as.character(ggv))

if (ggv < "3.5.2" || ggv >= "4.0.0") {
  stop(
    "Unsupported ggplot2 version detected (", as.character(ggv), ").\n",
    "Required: ggplot2 >= 3.5.2 and < 4.0.0.\n",
    "If needed, run: source('scr/00_fix_windows_packages.R') in a clean R session."
  )
}

## 4) Record package versions ----
message("=== Versions ===")
message("R:            ", getRversion())
message("Seurat:       ", as.character(packageVersion("Seurat")))
message("SeuratObject: ", as.character(packageVersion("SeuratObject")))
message("SeuratDisk:   ", as.character(packageVersion("SeuratDisk")))
message("ggplot2:      ", as.character(packageVersion("ggplot2")))
message("hdf5r:        ", as.character(packageVersion("hdf5r")))

# Record session information (reproducibility) ----
session_file <- file.path(base_dir, "outputs", "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = session_file)
message("sessionInfo saved to: ", session_file)

## 5) Input Data Validation ----
h5ad_file <- file.path(data_raw_dir, "endometrium_organoid.h5ad")
if (!file.exists(h5ad_file)) {
  stop("File not found:\n", h5ad_file,
       "\n\nPlease move 'endometrium_organoid.h5ad' into the 'data/raw/' directory and rerun the pipeline.")
}
message("Input dataset confirmed: ", h5ad_file)

# Define expected output path (.h5seurat will be created next to the .h5ad file)
created_h5seurat <- sub("\\.h5ad$", ".h5seurat", h5ad_file)
message("Expected .h5seurat path: ", created_h5seurat)

## 6) Convert .h5ad -> .h5seurat (if needed) ----
# Convert only when the .h5seurat file is missing.
# This keeps the pipeline re-runnable without trying to recreate the same file.

if (!file.exists(created_h5seurat)) {
  message("Converting .h5ad -> .h5seurat (first run)...")
  SeuratDisk::Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)
} else {
  message("Using existing .h5seurat (skipping conversion): ", created_h5seurat)
}

# Post-conversion validation
if (!file.exists(created_h5seurat)) {
  stop(
    "Conversion unsuccessful: expected .h5seurat file not detected.\n",
    "Potential causes: SeuratDisk/hdf5r incompatibility, a locked destination file, ",
    "or a corrupted .h5ad input."
  )
}

## 7) Object Initialization ----
# Load the .h5seurat file into memory as a Seurat object.
# UpdateSeuratObject() ensures structural compatibility with the
# current Seurat version in use.

message("Loading Seurat object into R environment...")
gc()  # Perform garbage collection to optimize memory allocation

pbmc <- SeuratDisk::LoadH5Seurat(file = created_h5seurat)
pbmc <- SeuratObject::UpdateSeuratObject(pbmc)

message("Seurat object successfully initialized.")
message("Dimensions — Cells: ", ncol(pbmc), 
        " | Features (genes): ", nrow(pbmc))

# Optional checkpoint: serialize raw loaded object for rapid reloading
saveRDS(pbmc, file = file.path(obj_dir, "seurat_loaded_raw.rds"))
message("Checkpoint saved to: outputs/objects/seurat_loaded_raw.rds")

## 8) Experimental filtering: retain control (non-inhibitor) cells ----

# Inspect distribution of experimental conditions
table(pbmc$Inhibitor)

# Validate presence of required metadata column
if (!"Inhibitor" %in% colnames(pbmc@meta.data)) {
  stop("Required metadata column 'Inhibitor' not found in Seurat object.")
}

# Subset cells based on primary experimental condition
pbmc_ctrl <- subset(pbmc, subset = Inhibitor == "Ctrl")

message("Cells retained (Inhibitor == 'Ctrl'): ", ncol(pbmc_ctrl))


## 9) QC plot (Control cells) ----
Idents(pbmc_ctrl) <- "Ctrl"

vln_qc <- VlnPlot(
  pbmc_ctrl,
  features = required_qc,
  ncol = 3,
  pt.size = 0,
  cols = "tomato"
) &
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )

vln_qc <- vln_qc +
  patchwork::plot_annotation(
    title = "Quality control metrics",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 20,
        face = "bold",
        hjust = 0.5   # centrar
      )
    )
  )

ggsave(
  filename = file.path(fig_dir, "QC_violinplot.png"),
  plot = vln_qc,
  width = 12, height = 4, dpi = 300
)

message("QC summaries:")
print(summary(pbmc_ctrl@meta.data[[qc_n_genes]]))
print(summary(pbmc_ctrl@meta.data[[qc_n_counts]]))
print(summary(pbmc_ctrl@meta.data[[qc_percent_mito]]))


## QC scatter plots (control cells) ----
# Standardize QC column names for plotting (Seurat tutorials expect these)
# - nCount_RNA: total RNA counts (UMIs) per cell
# - nFeature_RNA: number of detected genes per cell
# - percent.mt: fraction of mitochondrial transcripts per cell

if (!"nCount_RNA" %in% colnames(pbmc_ctrl@meta.data) && !is.na(qc_n_counts)) {
  pbmc_ctrl$nCount_RNA <- pbmc_ctrl@meta.data[[qc_n_counts]]
}
if (!"nFeature_RNA" %in% colnames(pbmc_ctrl@meta.data) && !is.na(qc_n_genes)) {
  pbmc_ctrl$nFeature_RNA <- pbmc_ctrl@meta.data[[qc_n_genes]]
}
if (!"percent.mt" %in% colnames(pbmc_ctrl@meta.data) && !is.na(qc_percent_mito)) {
  pbmc_ctrl$percent.mt <- pbmc_ctrl@meta.data[[qc_percent_mito]]
}

# Clean identity (avoids showing "AnnData" as legend group)
Idents(pbmc_ctrl) <- "Ctrl"

# Scatter 1: counts vs mitochondrial fraction
p_sc1 <- FeatureScatter(pbmc_ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggplot2::ggtitle("QC: RNA counts vs mitochondrial fraction") +
  ggplot2::labs(x = "nCount_RNA (total counts / UMIs)", y = "percent.mt (mitochondrial percentage)") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    legend.position = "none"
  )

# Scatter 2: counts vs detected genes
p_sc2 <- FeatureScatter(pbmc_ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggplot2::ggtitle("QC: RNA counts vs detected genes") +
  ggplot2::labs(x = "nCount_RNA (total counts / UMIs)", y = "nFeature_RNA (detected genes)") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    legend.position = "none"
  )

# Save each 
ggsave(file.path(fig_dir, "QC_scatter_nCount_vs_percentMT.png"),
       p_sc1, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "QC_scatter_nCount_vs_nFeature.png"),
       p_sc2, width = 7, height = 6, dpi = 300)

# Save both together in one figure 
p_scatter_combined <- p_sc1 + p_sc2 + patchwork::plot_layout(ncol = 2)

ggsave(file.path(fig_dir, "QC_scatter_combined.png"),
       p_scatter_combined, width = 14, height = 6, dpi = 300)

message("QC scatter plots saved (individual + combined).")

## 10) Normalization, HVG selection, scaling, and PCA (control cells) ----
# This section performs the standard Seurat preprocessing workflow on the
# control (Ctrl) subset:
# 1) Library-size normalization (LogNormalize)
# 2) Identification of highly variable genes (HVGs) using VST
# 3) Feature-wise scaling (z-scoring)
# 4) Principal component analysis (PCA) for downstream clustering/UMAP

# 10.1 Normalize counts (log-normalization with a fixed scale factor)
pbmc_ctrl <- NormalizeData(
  pbmc_ctrl,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# 10.2 Identify highly variable features (HVGs)
pbmc_ctrl <- FindVariableFeatures(
  pbmc_ctrl,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# 10.3 Scale highly variable genes prior to PCA
# Z-score scaling is applied to HVGs to standardize gene expression
# and ensure comparability across cells.
pbmc_ctrl <- ScaleData(
  pbmc_ctrl,
  features = VariableFeatures(pbmc_ctrl),
  verbose = FALSE
)

# 10.4 Run PCA on HVGs to capture the major axes of biological variance
# while reducing technical noise.
pbmc_ctrl <- RunPCA(
  pbmc_ctrl,
  features = VariableFeatures(pbmc_ctrl),
  verbose = FALSE
)

## PCA dimensionality selection ----
# The number of PCs retained is chosen using the elbow criterion and documented
# in the README to keep the workflow transparent and reproducible.
DIMS_USE <- 1:9
RESOLUTION <- 0.5

elbow_plot <- ElbowPlot(pbmc_ctrl) +
  ggplot2::ggtitle("Elbow Plot for Principal Component Selection") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(fig_dir, "ElbowPlot.png"),
  plot = elbow_plot,
  width = 6,
  height = 4,
  dpi = 300,
  bg = "white"
)

## 11) Graph-based clustering and UMAP embedding ----
# 11.1 Construct nearest-neighbor graph in PCA space
pbmc_ctrl <- FindNeighbors(pbmc_ctrl, dims = DIMS_USE, verbose = FALSE)

# 11.2 Perform graph-based clustering

pbmc_ctrl <- FindClusters(pbmc_ctrl, resolution = RESOLUTION, verbose = FALSE)

# 11.3 Compute UMAP embedding for visualization
pbmc_ctrl <- RunUMAP(pbmc_ctrl, dims = DIMS_USE, verbose = FALSE)
# 11.4 Visualize clusters in UMAP space
umap_clusters <- DimPlot(
  pbmc_ctrl,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  ggtitle(paste0(
    "UMAP of control cells | PCs: 1–", max(DIMS_USE),
    " | Resolution: ", RESOLUTION
  ))

ggsave(
  file.path(fig_dir, "UMAP_seurat_clusters.png"),
  umap_clusters,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)
## 12) Differential expression analysis: marker genes per cluster ----
# Identify cluster-specific marker genes using Seurat's Wilcoxon rank-sum test.
# Only positively enriched markers are retained (only.pos = TRUE),
# and genes must be expressed in at least 30% of cells within a cluster (min.pct = 0.3).
# This step provides biologically interpretable signatures for downstream
# annotation and functional characterization of epithelial cell states.

all_clusters <- FindAllMarkers(pbmc_ctrl, min.pct = 0.3, only.pos = TRUE, verbose = FALSE)
out_csv <- file.path(table_dir, "all_clusters_control.csv")
write.csv(all_clusters, file = out_csv, row.names = FALSE)
message("Markers CSV saved to: ", out_csv)

## 13) Heatmap of top marker genes per cluster ----
# To facilitate biological interpretation, we select the top 10 positively
# enriched marker genes for each cluster based on the highest average log2
# fold change (avg_log2FC).
#
# This approach prioritizes genes with the strongest differential expression
# signal within each cluster, highlighting transcriptional programs that
# define cluster identity.
#
# The resulting heatmap visualizes scaled expression patterns of these
# marker genes across all control cells, enabling qualitative assessment
# of cluster-specific gene signatures.

top10_markers <- all_clusters %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)

p_heat <- DoHeatmap(
  pbmc_ctrl,
  features = unique(top10_markers$gene),
  size = 3
) +
  ggtitle("Top marker genes per cluster")

ggsave(
  file.path(fig_dir, "Heatmap_top10_markers_per_cluster.png"),
  p_heat,
  width = 12,
  height = 10,
  dpi = 300
)
## 14) UMAP visualization with biological cell type annotation ----
# To assess the biological coherence of the unsupervised clustering,
# we visualize the UMAP embedding colored by the provided
# "Celltype_without_inhibitors" annotation.
#
# This allows comparison between computationally derived clusters
# and previously curated biological cell states, facilitating
# biological interpretation of cluster identities.

umap_celltypes <- DimPlot(
  pbmc_ctrl,
  reduction = "umap",
  group.by = "Celltype_without_inhibitors",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  ggtitle("UMAP projection colored by annotated cell states")

ggsave(
  file.path(fig_dir, "UMAP_Celltype_without_inhibitors.png"),
  umap_celltypes,
  width = 9,
  height = 7,
  dpi = 300
)

#UMAP colored by hormonal condition
p_hormones <- DimPlot(
  pbmc_ctrl,
  reduction = "umap",
  group.by = "Hormones",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  ggplot2::ggtitle("UMAP — hormonal condition") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Combine side-by-side and save
p_compare <- umap_celltypes | p_hormones
ggsave(
  filename = file.path(fig_dir, "UMAP_cellstate_vs_hormones.png"),
  plot = p_compare,
  width = 14, height = 7, dpi = 300,
  bg = "white"
)
## 15) Marker visualization: stem/progenitor and proliferation candidates ----
# Objective: map canonical epithelial/stromal markers and proliferation-associated genes
# onto the UMAP embedding to contextualize their expression across clusters/cell states.
#
# Markers:
# - Epithelial / stem-related: CDH2, EPCAM
# - Stromal (expected low in epithelial organoids): VIM, SUSD2, PDGFRB
# - Proliferation / cycling programs: S100P, CDH1, NEDD9

genes_to_plot <- c("CDH2", "EPCAM", "VIM", "SUSD2", "PDGFRB", "S100P", "CDH1", "NEDD9")

# Safety check: keep only genes that exist in the dataset features
genes_present <- genes_to_plot[genes_to_plot %in% rownames(pbmc_ctrl)]
genes_missing <- setdiff(genes_to_plot, genes_present)

if (length(genes_missing) > 0) {
  warning("Markers not found in this object (skipping): ", paste(genes_missing, collapse = ", "))
}

## 15.1 FeaturePlots: spatial expression on UMAP ----
# For each selected marker gene, generate a UMAP-based FeaturePlot.
# Expression values are visualized as a continuous gradient across cells.
# The cutoffs (q5–q95) reduce the impact of extreme outliers and
# improve visual interpretability of gene expression patterns.

for (g in genes_present) {
  p <- FeaturePlot(
    pbmc_ctrl,
    features = g,
    max.cutoff = "q95",
    min.cutoff = "q5"
  ) +
    ggtitle(paste0("FeaturePlot - ", g))
  
  ggsave(
    filename = file.path(fig_dir, paste0("FeaturePlot_", g, ".png")),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

## 15.2 DotPlot: marker summary across clusters ----
# Generate a DotPlot summarizing marker expression across Seurat clusters.
# Dot size represents the percentage of cells expressing the gene.
# Dot color reflects scaled average expression within each cluster.
# This provides a compact comparative overview of marker distribution.

if (length(genes_present) > 0) {
  p_dot <- DotPlot(
    pbmc_ctrl,
    features = genes_present,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("DotPlot of selected stem/proliferation markers by cluster") +
    theme_classic() +   # ← fondo blanco limpio
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(fig_dir, "DotPlot_selected_markers.png"),
    plot = p_dot,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"   # ← asegura fondo blanco en el PNG
  )
}

## 16) Save final processed object + session info ----
saveRDS(pbmc_ctrl, file = file.path(obj_dir, "seurat_final_processed.rds"))
message(" Final object saved: outputs/objects/seurat_final_processed.rds")

writeLines(capture.output(sessionInfo()), con = file.path(table_dir, "sessionInfo.txt"))
message("pipeline completed successfully.")
