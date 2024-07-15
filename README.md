# HW_6B.R Analysis

This R script performs an analysis on PBMC and genotype data using Seurat and other R packages. Below is a summary of the key steps and analyses performed:

### Library Imports and Setup
- Set the working directory using `setwd()`.
- Imported essential libraries such as `dplyr`, `Seurat`, `patchwork`, and `readr`.

### Data Loading and Initialization
- Loaded the PBMC dataset using `Read10X()` and created a Seurat object with `CreateSeuratObject()`.
- Examined a few genes in the first thirty cells and compared the object size of dense and sparse matrices.

### Quality Control (QC)
- Added columns to the object metadata for QC stats, specifically the percentage of mitochondrial genes.
- Visualized QC metrics with a violin plot using `VlnPlot()`.
- Created scatter plots to visualize feature-feature relationships using `FeatureScatter()`.

### Data Preprocessing
- Filtered the data based on QC metrics.
- Normalized the data using `NormalizeData()`.
- Identified variable features with `FindVariableFeatures()` and visualized them.

### Dimensionality Reduction
- Scaled the data and performed PCA with `RunPCA()`.
- Visualized PCA results with `DimPlot()` and `DimHeatmap()`.
- Conducted JackStraw analysis and created an ElbowPlot to determine significant PCs.

### Clustering and UMAP
- Found neighbors and clusters using `FindNeighbors()` and `FindClusters()`.
- Ran UMAP for visualization and plotted the clusters.

### Additional Analysis
- Loaded genotype data using `read_csv()`.
- Filtered clinical data and calculated the mean number of cigarettes per day for a specific demographic.
- Created a histogram of cigarettes per day.
- Performed PCA on genotype data and visualized it with `plot_ly()`.

## Key Findings
- Successfully processed and visualized PBMC data through QC, normalization, and clustering.
- Identified significant PCs and clusters using PCA and UMAP.
- Calculated demographic-specific statistics from clinical data.
- Visualized PCA results of genotype data to understand population diversity.

This summary encapsulates the major components and findings of the analysis. For detailed code and results, please refer to the specific sections in the script.
``` &#8203;:citation[oaicite:0]{index=0}&#8203;

