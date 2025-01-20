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
library(SingleCellExperiment)
library(scDblFinder)
library(grid)
library(biomaRt)

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
vln_plot <- VlnPlot(HCC_data, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 3)
ggsave("Results/QC_violin_RNA.png", plot = vln_plot, width = 10, height = 7)
vln_plot

# Removal of low-quality cells
HCC_data.filtered <- subset(
  HCC_data,
  subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &
    percent_mt < 2 &   
    percent_ribo < 20   
)

# Visualize QC metrics after filtering
vln_plot_2 <- VlnPlot(HCC_data.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), ncol = 3)
ggsave("Results/QC_violin_RNA_2.png", plot = vln_plot_2, width = 10, height = 7)
vln_plot_2

# Add an extra filtration step to keep genes that are present in at least 10 cells
HCC_data.filtered <- HCC_data.filtered[rowSums(HCC_data.filtered@assays$RNA@counts > 0) >= 10, ]

# Removal of low-quality genes 
C <- HCC_data.filtered@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

# Replace Ensembl IDs with Gene Symbols for boxplot visualization
gene_names <- gene_symbols[rownames(C)] # Map Ensembl IDs to gene symbols
rownames(C) <- ifelse(is.na(gene_names), rownames(C), gene_names) # Replace or keep Ensembl IDs

# Visualize the most expressed genes as a boxplot
# Set up the file path where the boxplot will be saved
png("Results/most_expressed_boxplot.png", width = 800, height = 600)
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
dev.off()
# Show the plot in RStudio after saving it
par(mar = c(4, 8, 2, 1))  # Ensure the plot area is set
boxplot(
  as.matrix(t(C[most_expressed, ])),
  cex = 0.1,
  las = 1,
  xlab = "% total count per cell",
  col = (scales::hue_pal())(20)[20:1],  # Color palette using scales::hue_pal()
  horizontal = TRUE,
  names = rownames(C[most_expressed, ])  # Gene names as labels
)


# We will filter the MALAT1 gene since it is over-expressed and it probably will not provide significant biological information (long non-coding RNA) 
# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Query for the Ensembl gene ID of MALAT1
malat1_ensg <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                     filters = "external_gene_name", 
                     values = "MALAT1", 
                     mart = ensembl)
# Extract Ensembl gene ID for MALAT1
malat1_ensg_id <- malat1_ensg$ensembl_gene_id[1]  # assuming the first result is MALAT1
# Remove MALAT1 from the counts matrix directly in HCC_data.filtered
HCC_data.filtered@assays$RNA@counts <- HCC_data.filtered@assays$RNA@counts[!rownames(HCC_data.filtered@assays$RNA@counts) %in% malat1_ensg_id, ]
# Also remove it from the data slot (for normalization) in HCC_data.filtered
HCC_data.filtered@assays$RNA@data <- HCC_data.filtered@assays$RNA@data[!rownames(HCC_data.filtered@assays$RNA@data) %in% malat1_ensg_id, ]


# Check how the most expressed look now
C2 <- HCC_data.filtered@assays$RNA@counts
C2 <- Matrix::t(Matrix::t(C2)/Matrix::colSums(C2)) * 100
most_expressed_2 <- order(apply(C2, 1, median), decreasing = T)[20:1]

# Replace Ensembl IDs with Gene Symbols for boxplot visualization
gene_names <- gene_symbols[rownames(C2)] # Map Ensembl IDs to gene symbols
rownames(C2) <- ifelse(is.na(gene_names), rownames(C2), gene_names) # Replace or keep Ensembl IDs

par(mar = c(4, 8, 2, 1))  # Ensure the plot area is set
boxplot(
  as.matrix(t(C2[most_expressed_2, ])),
  cex = 0.1,
  las = 1,
  xlab = "% total count per cell",
  col = (scales::hue_pal())(20)[20:1],  # Color palette using scales::hue_pal()
  horizontal = TRUE,
  names = rownames(C2[most_expressed_2, ])  # Gene names as labels
)


# Compare the number of cells and genes before and after filtering
dim(HCC_data)
dim(HCC_data.filtered)


### 3. Normalization
# Remove the corresponding feature metadata (for MALAT1)
HCC_data.filtered@assays$RNA@meta.features <- HCC_data.filtered@assays$RNA@meta.features[!rownames(HCC_data.filtered@assays$RNA@meta.features) %in% malat1_ensg_id, ]
HCC_data.filtered <- NormalizeData(HCC_data.filtered, normalization.method = "LogNormalize", scale.factor = 10000)


### 4. Feature selection
HCC_data.filtered <- FindVariableFeatures(HCC_data.filtered, selection.method = "vst", nfeatures = 2000)


### 5. Scaling
HCC_data.filtered <- ScaleData(HCC_data.filtered) # No need to adjust based on ribo (plot showed it was not necessary) 


### 6. Visualization
# Linear dimension reduction
HCC_data.filtered <- RunPCA(HCC_data.filtered, features = VariableFeatures(object = HCC_data.filtered))

# Elbow plot
elbow_plot <- ElbowPlot(HCC_data.filtered)
ggsave("Results/elbow_plot.png", plot = elbow_plot)
elbow_plot

# Lets see how many PCs we should use since elbow is not perfect elbow. Cover 80% of the variance
# Calculate variance explained by each PC

#Extract the PCA results (standard deviations of the PCs)
stdev <- HCC_data.filtered[["pca"]]@stdev

# Calculate the variance explained by each principal component
variance <- stdev^2 / sum(stdev^2)

# Calculate the cumulative variance explained
cumulative_var <- cumsum(variance)

# Find the number of components responsible for 80% of the variability
which(cumulative_var >= 0.80)[1]

# We need 16 dimensions

#Cluster cells
HCC_data.filtered <- FindNeighbors(HCC_data.filtered, dims = 1:16)
HCC_data.filtered <- FindClusters(HCC_data.filtered, resolution = c(0.3, 0.5, 0.7))

#Non linear dimentional reduction
HCC_data.filtered <- RunUMAP(HCC_data.filtered, dims = 1:16)
p1<- DimPlot(HCC_data.filtered, reduction = "umap")

plot_0.3 <- DimPlot(HCC_data.filtered, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap")
ggsave("Results/dimplot_res_0.3.png", plot = plot_0.3, width = 8, height = 6)
plot_0.3

plot_0.5 <- DimPlot(HCC_data.filtered, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "umap")
ggsave("Results/dimplot_res_0.5.png", plot = plot_0.5, width = 8, height = 6)
plot_0.5

plot_0.7 <- DimPlot(HCC_data.filtered, group.by = "RNA_snn_res.0.7", label = TRUE, reduction = "umap")
ggsave("Results/dimplot_res_0.7.png", plot = plot_0.7, width = 8, height = 6)
plot_0.7


# ### 7. Doublet detection
# #Convert your seurat object to a sce object
# # 1. Convert to SCE object
# sce <- as.SingleCellExperiment(HCC_data.filtered)
# 
# # 2. Create a "sample" variable combining age group and donor ID to account for batch effects
# sce$sample <- paste0(HCC_data.filtered$age_group, "_", HCC_data.filtered$donor_id)
# 
# # 3. Run scDblFinder for doublet detection
# set.seed(10010101)  # For reproducibility
# sce <- scDblFinder(sce, samples = "sample")  # Use "sample" to handle batches
# 
# # 4. Inspect the doublet classification results
# doublet_table <- table(sce$scDblFinder.class)
# doublet_percentage <- sum(sce$scDblFinder.class == "doublet") / length(sce$scDblFinder.class) * 100
# 
# # Print doublet results
# print(doublet_table)
# cat(sprintf("Doublet percentage: %.2f%%\n", doublet_percentage))
# 
# # 5. Add doublet classifications back to Seurat object
# HCC_data.filtered$scDblFinder.class <- sce$scDblFinder.class
# HCC_data.filtered$scDblFinder.score <- sce$scDblFinder.score
# 
# # 6. Visualize doublet scores on UMAP
# # Add UMAP coordinates to SCE for visualization (if not already present)
# if (is.null(reducedDim(sce, "UMAP"))) {
#   reducedDim(sce, "UMAP") <- HCC_data.filtered@reductions$umap@cell.embeddings
# }
# 
# # Plot UMAP with doublet scores
# scDblFinder_score_plot <- plotUMAP(sce, colour_by = "scDblFinder.score")
# ggsave("Results/scDblFinder_score.png", plot = scDblFinder_score_plot, width = 8, height = 6)
# scDblFinder_score_plot
# 
# # 7. Visualize clusters and doublets in Seurat
# # Plot Seurat clusters and highlight doublets
# DimPlot(HCC_data.filtered, group.by = "seurat_clusters", label = TRUE)
# scDblFinder_score_plot2 <- FeaturePlot(HCC_data.filtered, features = "scDblFinder.score")
# ggsave("Results/scDblFinder_score2.png", plot = scDblFinder_score_plot2, width = 8, height = 6)
# scDblFinder_score_plot2
# 
# # 8. Remove doublets
# # Remove doublets from the Seurat object
# HCC_data.filtered_no_doublets <- subset(HCC_data.filtered, subset = scDblFinder.class != "doublet")
# # Check the number of cells before and after removal
# cat(sprintf("Cells before removal: %d\n", ncol(HCC_data.filtered)))
# cat(sprintf("Cells after removal: %d\n", ncol(HCC_data.filtered_no_doublets)))
# # Visualize after doublet removal
# scDblFinder_score_plot_no_doublets <- FeaturePlot(HCC_data.filtered_no_doublets, features = "scDblFinder.score")
# ggsave("Results/scDblFinder_score_no_doublets.png", plot = scDblFinder_score_plot_no_doublets, width = 8, height = 6)
# scDblFinder_score_plot_no_doublets
# 
# dim(HCC_data.filtered_no_doublets)
# 
# # Save the filtered Seurat object
# saveRDS(HCC_data.filtered_no_doublets, "Data/scRNA_filtered_no_doublets.rds")






### 7. Doublet detection
# 1. Convert your Seurat object to a SingleCellExperiment (SCE) object
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

# 6. **Add UMAP coordinates to SCE for visualization (if not already present)**
# Explicitly assign UMAP coordinates from Seurat to SCE
if (is.null(reducedDim(sce, "UMAP"))) {
  reducedDim(sce, "UMAP") <- HCC_data.filtered@reductions$umap@cell.embeddings
}
# Remove JOINT_WNN_UMAP if it's not needed
if ("JOINT_WNN_UMAP" %in% names(reducedDims(sce))) {
  reducedDim(sce, "JOINT_WNN_UMAP") <- NULL
}

# Plot UMAP with doublet scores
scDblFinder_score_plot <- plotUMAP(sce, colour_by = "scDblFinder.score")
ggsave("Results/scDblFinder_score.png", plot = scDblFinder_score_plot, width = 8, height = 6)
scDblFinder_score_plot

# 7. Visualize clusters and doublets in Seurat
# Plot Seurat clusters and highlight doublets
scDblFinder_score_plot2 <- FeaturePlot(HCC_data.filtered, features = "scDblFinder.score", reduction = "umap")
ggsave("Results/scDblFinder_score2.png", plot = scDblFinder_score_plot2, width = 8, height = 6)
scDblFinder_score_plot2

# 8. Remove doublets
# Remove doublets from the Seurat object
HCC_data.filtered_no_doublets <- subset(HCC_data.filtered, subset = scDblFinder.class != "doublet")
# Check the number of cells before and after removal
cat(sprintf("Cells before removal: %d\n", ncol(HCC_data.filtered)))
cat(sprintf("Cells after removal: %d\n", ncol(HCC_data.filtered_no_doublets)))
# Visualize after doublet removal
scDblFinder_score_plot_no_doublets <- FeaturePlot(HCC_data.filtered_no_doublets, features = "scDblFinder.score", reduction = "umap")
ggsave("Results/scDblFinder_score_no_doublets.png", plot = scDblFinder_score_plot_no_doublets, width = 8, height = 6)
scDblFinder_score_plot_no_doublets

dim(HCC_data.filtered_no_doublets)

# Save the filtered Seurat object
saveRDS(HCC_data.filtered_no_doublets, "Data/scRNA_filtered_no_doublets.rds")