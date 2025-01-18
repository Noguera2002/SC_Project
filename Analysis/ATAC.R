# Install required packages if not already installed
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("Signac", quietly = TRUE)) {
  install.packages("Signac")
}
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  install.packages("hdf5r")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
  BiocManager::install("JASPAR2020")
}
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}

library(scDblFinder)
library(SingleCellExperiment)
library(scater)
library(hdf5r)
library(Seurat)
library(Signac)
library(ggplot2)
library(JASPAR2020)
library(Matrix)

### 1. Load the ATAC Seurat object 
atac_data<- readRDS("Data/ATAC.rds")
head(atac_data@meta.data)

### 2. Quality control
# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(atac_data),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Identify ribosomal genes
ribo_genes <- names(gene_symbols[grep("^RP[SL]", gene_symbols)])

# Calculate percent ribosomal content
atac_data <- PercentageFeatureSet(
  atac_data,
  features = ribo_genes,
  col.name = "percent_ribo"
)

# Dropped the  RNA counts so that it doesn't affect the analysis on the ATAC (variable features etc)
atac_data@meta.data <- atac_data@meta.data[, !(colnames(atac_data@meta.data) %in% c("nCount_RNA", "nFeature_RNA"))]

data.filtered <- subset(atac_data, subset = nCount_ATAC > 200  & nCount_ATAC < 100000 & nucleosome_signal<3 &percent_mt < 5)


# Cells that are present in less than 10 samples are excluded

# Filter genes expressed in at least 10 cells

# Efficiently calculate the number of cells expressing each gene in the sparse matrix
data.filtered <- data.filtered[rowSums(data.filtered@assays$RNA@counts > 0) >= 10, ]


VlnPlot(data.filtered, features = c("nFeature_ATAC", "nCount_ATAC","TSS_percentile","nucleosome_signal", "percent_mt", "percent_ribo"), ncol = 5)

dim(data)
dim(data.filtered)


#This is how we perform normalization in scATAC data
data.filtered <- RunTFIDF(data.filtered, method = 3)
data.filtered <- FindTopFeatures(data.filtered, min.cutoff="q75") # use the top 25%
data.filtered


# Singular Value Decomposition (dimensionality reduction)
data.filtered <- RunSVD(data.filtered)

# We can keep the following but they do say how many dimensions use  == 10
#Estimate the cummulative Variance explained to determine the dimensions for the UMAP

# Extract singular values from the SVD reduction
#svd_values <- data.filtered@reductions$lsi@stdev
# Compute variance explained
#variance_explained <- (svd_values^2) / sum(svd_values^2)
# Compute cumulative variance explained
#cumulative_variance <- cumsum(variance_explained)
# View cumulative variance explained
#cumulative_variance

# Create a data frame for plotting
#variance_df <- data.frame(
 # Component = 1:length(cumulative_variance),
#  CumulativeVariance = cumulative_variance
#)

# Plot using ggplot2
#ggplot(variance_df, aes(x = Component, y = CumulativeVariance)) +
 # geom_line() +
#  geom_point() +
 # labs(
 #   title = "Cumulative Variance Explained by SVD Components",
 #   x = "Component",
 #   y = "Cumulative Variance Explained"
 # ) +
 # theme_minimal()


# We exclude the first dimension as this is typically correlated with sequencing depth

data.filtered <- RunUMAP(data.filtered, reduction = "lsi", dims = 2:10, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
p1 <- DimPlot(data.filtered, reduction = "umap") 
p1

#Doublet identification



data.filtered$sample <- paste0(data.filtered$donor_id, "_", data.filtered$age_group)
sce <- as.SingleCellExperiment(data.filtered)
#Store the highly variable features as top.var (required to use with scDblFinder)
top.var <- VariableFeatures(data.filtered)


set.seed(123)
sce.dbl <- scDblFinder(sce, samples="sample", clusters=colLabels(sce))
#Plot the doublets

plotUMAP(sce.dbl, colour_by="scDblFinder.score")
plotUMAP(sce.dbl, colour_by="scDblFinder.class")

table(sce.dbl$scDblFinder.class)

# Keep only singlets
sce.filtered <- sce.dbl[, sce.dbl$scDblFinder.class == "singlet"]
plotUMAP(sce.filtered, colour_by="scDblFinder.score")

# Convert filtered SCE to Seurat object
data.filtered <- as.Seurat(sce.filtered)

# Proceed with UMAP and plotting
data.filtered <- RunUMAP(data.filtered, reduction = "LSI", dims = 2:10, 
                         reduction.name = "umap.atac", reduction.key = "umapatac_")

# Plot the UMAP
p2 <- DimPlot(data.filtered, group.by = "age_group",reduction = "UMAP")  # Use "umap.atac" for the custom reduction name
p2
DimPlot(data.filtered, reduction = "UMAP", split.by = "age_group")


saveRDS(data.filtered, "ATAC_filtered.rds")


---------------------------------------------
  
# Cell annotation
---------------------------------------------

data.filtered<- readRDS("ATAC_filtered.rds")

#split the data by the age group

data.filtered[["RNA"]] <- split(data.filtered[["RNA"]], f = data.filtered$age_group)

data.filtered[["RNA"]] <- JoinLayers(data.filtered[["RNA"]])

#This is how we perform normalization in scATAC data
DefaultAssay(data.filtered) <- "RNA"
data.filtered[["RNA"]] <- SelectAssayLayer(data.filtered[["RNA"]], layer = "data")

data.filtered <- RunTFIDF(data.filtered, method = 3)

# Singular Value Decomposition (dimensionality reduction) instead of PCA
data.filtered <- RunSVD(data.filtered)












#Create the geneActivity matrix
# estimate the transcriptional activity of each gene by quantifying ATAC-seq counts in the 2 kb-upstream region and gene body, using the GeneActivity() function in the Signac package. 

#gene.activities <- GeneActivity(data.filtered, features = VariableFeatures(pbmc.rna))



