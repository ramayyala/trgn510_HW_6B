setwd('G:/My Drive/Academics/TRGN 510/HW_6B')
library(dplyr)
library(Seurat)
library(patchwork)
# Q1-Q3

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


#Q4
library(readr)
kg_genotypes <- read_csv("trgn510/kg.genotypes.csv")
View(kg_genotypes)

#Q5
library(dplyr)

trgn.clinical<- read.csv('trgn510/trgn599.clinical.tsv',sep="\t")
test <- filter(trgn.clinical, gender =="female" & race=="white") 
mean(test$cigarettes_per_day,na.rm=TRUE)

#Q6
hist(trgn.clinical$cigarettes_per_day)


#Q7
library(dplyr)
library(plotly)
genotypes <- read.csv('trgn510/kg.genotypes.csv',header = TRUE, sep = ",", 
                      quote = "\"", dec = ".", fill = TRUE, row.names = 1)
population_details <- read.csv('trgn510/kg.poplist.csv',header = TRUE, sep = ",", 
                               quote = "\"", dec = ".", fill = TRUE)
PopList <- data.frame(population_details[,c(4,3,2)])
mini<-genotypes[c(1:10000),]
#Takes 60 minutes to do 
pca_anal <- prcomp(mini)
PCs <- data.frame(Sample=row.names(pca_anal$rotation),
                  pca_anal$rotation[,1:6],row.names = NULL) 
people_PCA_diversity<-inner_join(PopList, PCs, b, by = "Sample")
plot_ly(people_PCA_diversity, x = ~PC1, y = ~PC2, z = ~PC4, color = ~Cohort)

