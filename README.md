# Endometrial Organoid scRNA-seq Analysis  
## Stem/Progenitor-Associated Programs in Human Endometrial Epithelium

Single-cell RNA-seq analysis of human endometrial epithelial organoids to characterize cellular heterogeneity and explore transcriptional programs linked to epithelial plasticity, proliferation, and tissue regeneration.

This project reconstructs and modernizes my original Master’s thesis into a fully reproducible **R + Seurat pipeline**, implemented using best practices in single-cell transcriptomics and high-dimensional data analysis.

**Dataset scale:** ~118,672 cells × 27,254 genes (sparse matrix)

This project demonstrates how clinically relevant biological questions 
can be translated into scalable, reproducible computational workflows.

---

## Clinical & Scientific Motivation

The human endometrium regenerates cyclically, suggesting the presence of resident epithelial stem/progenitor populations.

Understanding these programs is relevant for:

- Fertility & implantation biology  
- Endometrial regeneration  
- Reproductive disorders  
- Early tumorigenic mechanisms  
- Translational precision medicine  

This project evaluates whether canonical epithelial markers and proliferation-associated genes show convergent spatial patterns across epithelial subpopulations.

---

## Experimental Design Overview

- Human endometrial epithelial organoids
  - Timepoints: Day0, Day2, Day6
  - Hormonal conditions: Control, Estrogen (E), Estrogen + Progesterone (E+P)
  - Pathway inhibition: NOTCH (DBZ), WNT (XAV)
- Single-cell RNA sequencing
- Dimensionality reduction (PCA, UMAP)
- Graph-based clustering (Seurat)
- Differential expression analysis
- Marker-driven biological interpretation

---

## Methods Overview

Pipeline implemented entirely in **R (Seurat 5)**:

1. Import `.h5ad`
2. Convert to `.h5seurat`
3. Subset control condition
4. Quality control inspection
5. Log-normalization (scale factor 10,000)
6. Highly variable gene selection (VST)
7. Scaling
8. PCA
9. Graph-based clustering
10. UMAP embedding
11. Differential expression (`FindAllMarkers`, Wilcoxon test)
12. Marker visualization (FeaturePlot, DotPlot, Heatmap)
13. Export tables, figures, and session metadata

---

## Differential Expression Strategy

Cluster-specific markers were identified using:

- `FindAllMarkers()`  
- Wilcoxon rank-sum test  
- Positive markers only  
- `min.pct = 0.3`

All results are programmatically exported as `.csv` files for traceability.

---

## Functional Annotation (Interpretative Layer)

Top markers were manually reviewed using:

- UniProt  
- Human Protein Atlas  
- GeneCards  
- Selected peer-reviewed literature  

Genes were grouped into high-level themes (proliferation, differentiation, ciliation, secretion, immune modulation, tissue remodeling).

This interpretative step is intentionally separated from the computational pipeline to preserve reproducibility.

---

# Key Results

## Clustering

**12 epithelial subpopulations** identified via UMAP.

Functional programs include:

- Proliferation  
- Differentiation  
- Implantation-related activity  
- Immune signaling  
- Tissue repair  

---

## Marker Expression Insights

- **CDH2** enriched in angiogenesis/tissue repair clusters  
- **EPCAM** broadly expressed, higher in implantation and cell-cycle clusters  
- **SUSD2** and **PDGFRB** minimally expressed (consistent with epithelial model)  
- **VIM** broadly expressed (epithelial–mesenchymal plasticity)  

Proliferation-associated genes:

- **S100P** enriched in proliferative clusters  
- **CDH1** broadly expressed across epithelial states  
- **NEDD9** enriched in tissue remodeling clusters  

---

## Interpretation

Overlap between canonical epithelial markers (**CDH2, EPCAM**) and proliferation genes (**CDH1, NEDD9**) suggests a potential association between epithelial plasticity and proliferative programs.

Findings are hypothesis-generating and require functional in vivo validation.

---

# Repository Structure

```text
.
endometrium-singlecell-pipeline/
│
├── .gitignore  
├── README.md
├── .gitignore
│
├── data/
│ └── raw/ # External dataset (not versioned)
│
├── src/
│ ├── organoids-script.R
│ └── 00_fix_windows_packages.R
│
├── outputs/
│ ├── figures/ # Generated plots
│ ├── tables/ # DEG results + sessionInfo.txt
│ └── objects/ # Saved Seurat .rds objects
```
---

## Reproducibility & Environment

Validated on:

- R 4.3.3 (Windows UCRT)
- Seurat 5.2.1
- SeuratDisk 0.0.0.9021
- ggplot2 3.5.2 (pinned for stability)

### Windows Compatibility Patch

`00_fix_windows_packages.R` pins a stable version of ggplot2 
to prevent known Seurat plotting issues on Windows (R 4.3.3).

This script is optional and only required on Windows systems.

Full session metadata available in:

```text
.
outputs/tables/sessionInfo.txt
```
---
## Data Availability

Raw data files are not included in this repository due to size restrictions.

Please download the original `.h5ad` file from the publication source and place it in:

data/raw/

Expected filename:
endometrium_organoid.h5ad

---

# Technical Skills Demonstrated

- Single-cell RNA-seq analysis
- Graph-based clustering
- Differential gene expression
- High-dimensional sparse matrix handling
- Cross-format data conversion (.h5ad → .h5seurat)
- Reproducible research engineering
- Dependency/version control
- Biological interpretation of transcriptomic data

---

## Author

Jimena Taciana Garcia, MD  
Bioinformatics | Data Science | Reproductive Genomics