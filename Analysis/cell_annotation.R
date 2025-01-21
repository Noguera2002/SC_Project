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
scRNA <- readRDS("./annotation_2.rds")


DimPlot(scRNA, reduction = "umap_harmony", group.by = "RNA_snn_res.0.2", label = TRUE)
# Visualize clustering resolutions using clustree
clustree_plot <- clustree(scRNA, prefix = "RNA_snn_res.")

# Save the plot to a file
ggsave("clustree_plot.png", plot = clustree_plot, width = 10, height = 8, dpi = 300)

Idents(scRNA) <- "RNA_snn_res.0.2"

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
singleRplot <- DimPlot(scRNA, group.by = "singleR.labels")

singleRplot_umap <- DimPlot(scRNA, reduction = "umap", group.by = "singleR.labels")
ggsave("SingleR_umap.png", plot = singleRplot_umap, width = 10, height = 8, dpi = 300)

# Diagnostics for SingleR annotation
# 1. View scores
print(pred$scores)

# 2. Plot score heatmap
plotScoreHeatmap(pred)

# 3. Plot delta distribution
delta_scores <- plotDeltaDistribution(pred)
delta_scores
ggsave("delta_scores.png", plot = delta_scores, width = 10, height = 8, dpi = 300)

saveRDS(scRNA, "annotation.rds")





annotation <- readRDS("Analysis/annotation_2.rds")

cluster7_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 7,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster6_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 6,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)


cluster5_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 5,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster4_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 4,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster3_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 3,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster2_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 2,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster1_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 1,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)



cluster0_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 0,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster8_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 8,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)



cluster9_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 9,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster10_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 10,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster11_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 11,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster12_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 12,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)

cluster13_conserved_markers <- FindConservedMarkers(annotation,
                                                   ident.1 = 13,
                                                   grouping.var = "age_group",
                                                   only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster14_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 14,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)


cluster15_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 15,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)
cluster16_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 16,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)

cluster17_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 17,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)
cluster18_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 18,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)


----------------------------------

###Manual annotation ##




set.seed(123)
library(Seurat)
library(harmony)

integrate.data <- readRDS("./Data/scRNA_integrated_no_doublets.rds")

#Identify the conserved markers
DimPlot(integrate.data, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)
DimPlot(integrate.data, reduction = 'umap', group.by = "age_group")

# Ensure the Seurat object is integrated and Idents are properly set
library(Seurat)

# Get the list of all clusters

# Set the active identity class to the integrated clustering results
# Set the active identity class to the integrated clustering results
Idents(integrate.data) <- integrate.data$RNA_snn_res.0.2

# Create a mapping of cluster IDs to cell type names
integrate.data <- RenameIdents(object = integrate.data,
  "0" = "Astrocyte",
  "1" = "Endothelial cell",
  "2" = "Oligodendrocyte precursor cell",
  "3" = "Oligodendrocyte",
  "4" = "Glutamatergic neuron",
  "5" = "Microglial cell",
  "6" = "Pyramidal cell", #
  "7" = "Gabaergic neuron",#
  "8" = "Neuroblast",
  "9" = "CD4+, alpha-beta T cell",
  "10" = "Pyramidal Cell",
  "11" = "Oligodendrocyte precursor cell",
  "12" = "NA",
  "13" = "Endothelial cell",
  "14" = "Pericyte cell",
  "15" = "Oligodendrocyte precursor cell",
  "16" = "Microglial cell",
  "17" = "Fibroblast",
  "18" = "Astrocyte"
)

# Assign cell types based on cluster IDs
integrate.data$cell_type_manual <- Idents(integrate.data)

# visualize data

clusters <- DimPlot(integrate.data, reduction = 'umap_harmony', label = TRUE)
celltype <- DimPlot(integrate.data, reduction = 'umap_harmony', group.by = 'cell_type_manual', label = TRUE)
ggsave(filename = 'cell_annotation.png', plot = celltype)





# Optionally save the results to a file
saveRDS(integrate.data, file ="annotated_cells.rds")


-----------------------------------------------------------------------------------

  
##Mapping of the genes

# Load the necessary library
library(org.Hs.eg.db)

# List of ENSG IDs
ensg_ids <- row.names(cluster7_conserved_markers)

# Map ENSG IDs to gene symbols
gene_names <- AnnotationDbi::select(org.Hs.eg.db, keys = ensg_ids, keytype = "ENSEMBL", columns = "SYMBOL")

# Print the results
print(gene_names)
gene_names_only <- gene_names$SYMBOL
gene_names_only# Filter out the None values and structure them one under the other


