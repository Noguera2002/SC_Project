# Install necessary packages (if not already installed)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("EnsDb.Hsapiens.v75")  # For GRCh37/hg19

# Load required libraries
library(SingleR)
library(celldex)
library(org.Hs.eg.db)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
library(tidyr)
library(multtest)
library(metap)
library(tibble)
library(AnnotationDbi)  # For mapping IDs
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)
library(clustree)

# Step 1: Read and visualize the data
scRNA <- readRDS("../Data/scRNA_integrated.rds")

# Visualize clustering resolutions using clustree
clustree_plot <- clustree(scRNA, prefix = "RNA_snn_res.")

# Save the plot to a file
ggsave("clustree_plot.png", plot = clustree_plot, width = 10, height = 8, dpi = 300)

Idents(scRNA) <- "RNA_snn_res.0.3"

# Step 2: Load the reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Optional: Inspect the reference dataset
ref_metadata <- as.data.frame(colData(ref))
# print(ref_metadata)  # Uncomment to explore the reference labels

# Step 3: Prepare scRNA.counts
# Extract counts matrix from Seurat object
scRNA.counts <- GetAssayData(scRNA, slot = "counts")

# Convert ENSG identifiers in scRNA.counts to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = rownames(scRNA.counts),
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# Replace rownames with gene symbols
rownames(scRNA.counts) <- gene_symbols

# Remove rows with NA (no mapping to gene symbols)
scRNA.counts <- scRNA.counts[!is.na(rownames(scRNA.counts)), ]

# Step 4: Standardize gene names (case sensitivity)
rownames(scRNA.counts) <- toupper(rownames(scRNA.counts))
rownames(ref) <- toupper(rownames(ref))

# Step 5: Find common genes
common_genes <- intersect(rownames(scRNA.counts), rownames(ref))

# Subset both datasets to the common genes
scRNA.counts <- scRNA.counts[common_genes, ]
ref <- ref[common_genes, ]

# Step 6: Run SingleR
pred <- SingleR(test = scRNA.counts, 
                ref = ref, 
                labels = ref$label.main,
                de.method = "wilcox")

# Output predictions
print(pred)

saveRDS(pred, "cell_type_pred_singleR.rds")


pred <- readRDS("cell_type_pred_singleR.rds")
# Add SingleR predictions to metadata
scRNA@meta.data$singleR.labels <- pred$pruned.labels[match(rownames(scRNA@meta.data), rownames(pred))]


# Verify the label distribution
table(scRNA@meta.data$singleR.labels)

# Plot UMAP with SingleR labels, including "Unknown"
DimPlot(scRNA, group.by = "singleR.labels")

# Diagnostics for SingleR annotation
# 1. View scores
print(pred$scores)

# 2. Plot score heatmap
plotScoreHeatmap(pred)

# 3. Plot delta distribution
delta_scores <- plotDeltaDistribution(pred)
ggsave("delta_scores.png", plot = delta_scores, width = 10, height = 8, dpi = 300)

saveRDS(scRNA, "annotation.rds")


----------------------------------

###Manual annotation ##




set.seed(123)
library(Seurat)
library(harmony)

integrate.data <- readRDS("../Data/scRNA_integrated.rds")

#Identify the conserved markers
DimPlot(integrate.data, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(integrate.data, reduction = 'umap', group.by = "age_group")

# Ensure the Seurat object is integrated and Idents are properly set
library(Seurat)

# Get the list of all clusters
Idents(integrate.data) <- "seurat_clusters"
clusters <- unique(Idents(integrate.data))

# Create an empty list to store the conserved markers for each cluster
conserved_markers_list <- list()

# Loop through each cluster and find conserved markers
for (cluster in clusters) {
  conserved_markers <- FindConservedMarkers(
    integrate.data,
    ident.1 = cluster,
    grouping.var = "age_group",  # Replace with the variable of interest
    only.pos = TRUE,
    min.pct = 0.25,
    min.diff.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  # Store the results in the list
  conserved_markers_list[[paste0("Cluster_", cluster)]] <- conserved_markers
}

# To view the results for a specific cluster, e.g., Cluster 1:
print(conserved_markers_list$Cluster_1)



# Set the active identity class to the integrated clustering results
Idents(integrate.data) <- integrate.data$RNA_snn_res.0.3
celltype <- rep(NA, length = ncol(integrate.data))

# Map clusters to cell types
celltype[which(Idents(integrate.data) %in% c(11))] <- 'RG'
celltype[which(Idents(integrate.data) %in% c('5_3'))] <- 'IPC'
celltype[which(Idents(integrate.data) %in% c(8))] <- 'IN-MGE'
celltype[which(Idents(integrate.data) %in% c(10))] <- 'IN-CGE'
celltype[which(Idents(integrate.data) %in% c(3))] <- 'IN-fetal'
celltype[which(Idents(integrate.data) %in% c(0,14))] <- 'EN-fetal-late'
celltype[which(Idents(integrate.data) %in% c(6,12))] <- 'EN'
celltype[which(Idents(integrate.data) %in% c('5_0','5_1','5_2','5_4'))] <- 'EN-fetal-early'

celltype[which(Idents(integrate.data) %in% c(4,9,15))] <- 'Astrocytes'
celltype[which(Idents(integrate.data) %in% c(1,21,22))] <- 'Oligodendrocytes'
celltype[which(Idents(integrate.data) %in% c(2,19))] <- 'OPC'
celltype[which(Idents(integrate.data) %in% c(7,13,17,23))] <- 'Microglia'
celltype[which(Idents(integrate.data) %in% c(16))] <- 'Endothelial'
celltype[which(Idents(integrate.data) %in% c(18))] <- 'Pericytes'
celltype[which(Idents(integrate.data) %in% c(20))] <- 'VSMC'



# Convert cell type annotations to a factor
celltype <- factor(celltype, 
                   levels = c('RG', 'IPC', 'EN-fetal-early', 'EN-fetal-late', 'EN', 
                              'IN-fetal', 'IN-MGE', 'IN-CGE',
                              'OPC', 'Astrocytes', 'Oligodendrocytes', 'Microglia', 
                              'Endothelial', 'Pericytes', 'VSMC'), 
                   ordered = TRUE)


integrate.data$celltype.manual <- celltype


# Check for unmapped clusters
if (any(is.na(celltype))) {
  warning("Some clusters are not assigned a cell type.")
}

# Visualize with DimPlot
DimPlot(integrate.data, group.by = "celltype.manual", reduction = "umap", label = TRUE) 

# Optionally save the results to a file
saveRDS(conserved_markers_list, file = "conserved_markers_all_clusters.rds")

