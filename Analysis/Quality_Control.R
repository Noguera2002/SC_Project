# Load libraries
library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(BiocSingular)
library(scDblFinder)
library(SingleCellExperiment)
library(scater)
library(org.Hs.eg.db)
library(Matrix)
library(scales)
library(AnnotationDbi)
library(org.Hs.eg.db)

### 1. Load the RNA Seurat object 
HCC_data <- readRDS("Data/scRNA.rds") 
head(HCC_data@meta.data) # Investigate metadata


### 2. Quality Control

# Mitochondrial genes. They indicate the presence of dead cells.
# Ribosomal genes. They are highly expressed in all cells. Bad for normalization.
# Thankfully the Seurat metadata comes with the percentage of mitochondrial genes.

# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(HCC_data),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Identify ribosomal genes
ribo_genes <- names(gene_symbols[grep("^RP[SL]", gene_symbols)])

# Calculate percent ribosomal content
HCC_data <- PercentageFeatureSet(
  HCC_data,
  features = ribo_genes,
  col.name = "percent_ribo"
)

# Visualize QC metrics as a violin plot
VlnPlot(HCC_data, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 3)

# Removal of low-quality cells
HCC_data.filtered <- subset(
  HCC_data,
  subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &
    percent_mt < 2 &   # I believe we can go higher, like 5
    percent_ribo < 20   # I believe we can go lower, like 2 or 5
)

# Visualize QC metrics after filtering
VlnPlot(HCC_data.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 3)

# Removal of low-quality genes
C <- HCC_data.filtered@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

# Replace Ensembl IDs with Gene Symbols for boxplot visualization
gene_names <- gene_symbols[rownames(C)] # Map Ensembl IDs to gene symbols
rownames(C) <- ifelse(is.na(gene_names), rownames(C), gene_names) # Replace or keep Ensembl IDs

# Visualize the most expressed genes as a boxplot
par(mar = c(4, 8, 2, 1))
boxplot(
  as.matrix(t(C[most_expressed, ])),
  cex = 0.1,
  las = 1,
  xlab = "% total count per cell",
  col = (scales::hue_pal())(20)[20:1],  # Color palette using scales::hue_pal()
  horizontal = TRUE,
  names = rownames(C[most_expressed, ])  # Gene names as labels
)
# No need to filter genes

# Compare the number of cells and genes before and after filtering
dim(HCC_data)
dim(HCC_data.filtered)



### 3. Normalization
HCC_data.filtered <- NormalizeData(HCC_data.filtered, normalization.method = "LogNormalize", scale.factor = 10000)



### 4. Feature selection
HCC_data.filtered <- FindVariableFeatures(HCC_data.filtered, selection.method = "vst", nfeatures = 2000)



### 5. Scaling
HCC_data.filtered <- ScaleData(HCC_data.filtered) # No need to adjust based on ribo (plot showed it was not necessary) 



### 6. Visualization
# Linear dimension reduction
HCC_data.filtered <- RunPCA(HCC_data.filtered, features = VariableFeatures(object = HCC_data.filtered))

# Elbow plot
ElbowPlot(HCC_data.filtered) 

# Lets see how many PCs we should use since elbow is not perfect elbow. Cover 80% of the variance
# Calculate variance explained by each PC

#Extract the PCA results (standard deviations of the PCs)
stdev <- HCC_data.filtered[["pca"]]@stdev

# Calculate the variance explained by each principal component
variance <- stdev^2 / sum(stdev^2)

# Calculate the cumulative variance explained
cumulative_var <- cumsum(variance)

# Find the number of components responsible for 80% of the variability
which(cumulative_var >= 0.75)[1]

# We need 16 dimensions

#Cluster cells
HCC_data.filtered <- FindNeighbors(HCC_data.filtered, dims = 1:16)
HCC_data.filtered <- FindClusters(HCC_data.filtered, resolution = c(0.3, 0.5, 0.7))

#Non linear dimentional reduction
HCC_data.filtered <- RunUMAP(HCC_data.filtered, dims = 1:16)
p1<- DimPlot(HCC_data.filtered, reduction = "umap")

DimPlot(HCC_data.filtered, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(HCC_data.filtered, group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(HCC_data.filtered, group.by = "RNA_snn_res.0.7", label = TRUE)


### 7. Doublet detection
#Convert your seurat object to a sce object
library(SingleCellExperiment)
library(scDblFinder)

# 1. Convert to SCE object

sce <- as.SingleCellExperiment(HCC_data.filtered)

# 2. Create a "sample" variable combining age group and donor ID to account for batch effects
sce$sample <- paste0(HCC_data.filtered$age_group, "_", HCC_data.filtered$donor_id)

# 3. Run scDblFinder for doublet detection
set.seed(10010101)  # For reproducibility
sce <- scDblFinder(sce, samples = "sample")  # Use "sample" to handle batches

# 4. Inspect the doublet classification results
doublet_table <- table(sce$scDblFinder.class)
doublet_percentage <- sum(sce$scDblFinder.class == "doublet") / length(sce$scDblFinder.class) * 100

# Print doublet results
print(doublet_table)
cat(sprintf("Doublet percentage: %.2f%%\n", doublet_percentage))

# 5. Add doublet classifications back to Seurat object
HCC_data.filtered$scDblFinder.class <- sce$scDblFinder.class
HCC_data.filtered$scDblFinder.score <- sce$scDblFinder.score

# 6. Visualize doublet scores on UMAP
# Add UMAP coordinates to SCE for visualization (if not already present)
if (is.null(reducedDim(sce, "UMAP"))) {
  reducedDim(sce, "UMAP") <- HCC_data.filtered@reductions$umap@cell.embeddings
}

# Plot UMAP with doublet scores
plotUMAP(sce, colour_by = "scDblFinder.score")

# 7. Visualize clusters and doublets in Seurat
# Plot Seurat clusters and highlight doublets
DimPlot(HCC_data.filtered, group.by = "seurat_clusters", label = TRUE)
FeaturePlot(HCC_data.filtered, features = "scDblFinder.score")


# very low percentage of doublets, so we can proceed with the analysis. No need to remove them.

# Save the filtered Seurat object
saveRDS(HCC_data.filtered, "Data/scRNA_filtered.rds")
