### Data Integration

# Load required libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
library(tidyr)
library(multtest)
library(metap)
library(tibble)
library(future)

# Set global options
set.seed(123)
plan("multisession", workers = 4)  # Parallel processing
options(future.globals.maxSize = 8000 * 1024^2)  # Increase memory allocation

# Load preprocessed Seurat object
scRNA <- readRDS("Data/scRNA_filtered_no_doublets.rds")

# Visualize data
Idents(scRNA) <- "RNA_snn_res.0.2"
umap_age_group <- DimPlot(scRNA, reduction = 'umap', group.by = 'age_group')
ggsave("Results/dimplot_age_group.png", plot = umap_age_group, width = 8, height = 6)

umap_split_age_group <- DimPlot(scRNA, reduction = 'umap', split.by = 'age_group')
ggsave("Results/dimplot_split_age_group.png", plot = umap_split_age_group, width = 8, height = 6)

### Preprocess Data for Integration

# Split data based on age group
scRNA[["RNA"]] <- split(scRNA[["RNA"]], f = scRNA$age_group)

# Normalize, find variable features, scale, and run PCA
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA)

# Determine number of PCs
stdev <- scRNA[["pca"]]@stdev
variance <- stdev^2 / sum(stdev^2)
cumulative_var <- cumsum(variance)
num_pcs <- which(cumulative_var >= 0.80)[1]  # PCs for 80% variance
cat(sprintf("Number of PCs for 80%% variance: %d\n", num_pcs))  # Should print 16

### Integrate Data Using Harmony and CCA

# CCA Integration
scRNA <- IntegrateLayers(
  object = scRNA,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "CCA_Integration",
  verbose = FALSE
)

# Save intermediate results
saveRDS(scRNA, "Data/scRNA_integrated_cca_no_doublets2.rds")

# Harmony Integration
scRNA <- IntegrateLayers(
  object = scRNA,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "Harmony_Integration",
  verbose = FALSE
)

# Reset future plan
plan("sequential")

# Re-join split layers
scRNA[["RNA"]] <- JoinLayers(scRNA[["RNA"]])

# Save integrated data
saveRDS(scRNA, "Data/scRNA_integrated_harmony_no_doublets2.rds")

# Run UMAP on integrated reductions
scRNA <- RunUMAP(scRNA, dims = 1:num_pcs, reduction = "CCA_Integration", reduction.name = "umap_cca")
scRNA <- RunUMAP(scRNA, dims = 1:num_pcs, reduction = "Harmony_Integration", reduction.name = "umap_harmony")

# Visualize UMAPs
p1 <- DimPlot(scRNA, reduction = 'umap', group.by = 'age_group') + ggtitle("UMAP Raw Data")
ggsave("Results/umap_raw_data2.png", plot = p1, width = 8, height = 6)

p2 <- DimPlot(scRNA, reduction = "umap_cca", group.by = "age_group") + ggtitle("UMAP CCA")
ggsave("Results/umap_cca2.png", plot = p2, width = 8, height = 6)

p3 <- DimPlot(scRNA, reduction = "umap_harmony", group.by = "age_group") + ggtitle("UMAP Harmony")
ggsave("Results/umap_harmony2.png", plot = p3, width = 8, height = 6)

# Combine UMAP plots
combined_umap_plot <- wrap_plots(p1, p2, p3, nrow = 1) + plot_layout(guides = "collect")
ggsave("Results/combined_umap_plots2.png", plot = combined_umap_plot, width = 24, height = 8)
combined_umap_plot

### Find Neighbors and Clusters for Harmony and CCA Reductions

# Harmony clustering
scRNA.harmony <- FindNeighbors(scRNA, dims = 1:num_pcs, reduction = "Harmony_Integration")
scRNA.harmony <- FindClusters(scRNA.harmony, resolution = 0.2)

# CCA clustering
scRNA.cca <- FindNeighbors(scRNA, dims = 1:num_pcs, reduction = "CCA_Integration")
scRNA.cca <- FindClusters(scRNA.cca, resolution = 0.2)



#################
# Harmony clustering
scRNA.harmony <- FindNeighbors(scRNA, 
                               dims = 1:num_pcs, 
                               reduction = "Harmony_Integration", 
                               graph.name = c("Harmony_nn", "Harmony_snn"))
scRNA.harmony <- FindClusters(scRNA.harmony, 
                              resolution = 0.2, 
                              graph.name = "Harmony_snn")

# CCA clustering
scRNA.cca <- FindNeighbors(scRNA, 
                           dims = 1:num_pcs, 
                           reduction = "CCA_Integration", 
                           graph.name = c("CCA_nn", "CCA_snn"))
scRNA.cca <- FindClusters(scRNA.cca, 
                          resolution = 0.2, 
                          graph.name = "CCA_snn")
#################

# Save clustering results
saveRDS(scRNA.harmony, "Data/scRNA_harmony_clustered.rds")
saveRDS(scRNA.cca, "Data/scRNA_cca_clustered.rds")

### Visualize Clustering Results

# Assign identities for Harmony and CCA
Idents(scRNA.harmony) <- "Harmony_snn_res.0.2"
Idents(scRNA.cca) <- "CCA_snn_res.0.2"

# Plot Harmony clustering
DimPlot(scRNA.harmony, reduction = "umap_harmony", split.by = "age_group")
#ggsave("Results/umap_harmony_clusters2.png", width = 8, height = 6)

# Plot CCA clustering
DimPlot(scRNA.cca, reduction = "umap_cca", split.by = "age_group")
#ggsave("Results/umap_cca_clusters2.png", width = 8, height = 6)



# CELL ANNOTATION PART

# Step 2: Load the reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Optional: Inspect the reference dataset
ref_metadata <- as.data.frame(colData(ref))
# print(ref_metadata)  # Uncomment to explore the reference labels

# Step 3: Prepare scRNA.harmony.counts
# Extract counts matrix from Seurat object
scRNA.harmony.counts <- GetAssayData(scRNA.harmony, slot = "counts")

# Convert ENSG identifiers in scRNA.harmony.counts to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = rownames(scRNA.harmony.counts),
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# Replace rownames with gene symbols
rownames(scRNA.harmony.counts) <- gene_symbols

# Remove rows with NA (no mapping to gene symbols)
scRNA.harmony.counts <- scRNA.harmony.counts[!is.na(rownames(scRNA.harmony.counts)), ]

# Step 4: Standardize gene names (case sensitivity)
rownames(scRNA.harmony.counts) <- toupper(rownames(scRNA.harmony.counts))
rownames(ref) <- toupper(rownames(ref))

# Step 5: Find common genes
common_genes <- intersect(rownames(scRNA.harmony.counts), rownames(ref))

# Subset both datasets to the common genes
scRNA.harmony.counts <- scRNA.harmony.counts[common_genes, ]
ref <- ref[common_genes, ]

# Step 6: Run SingleR
pred <- SingleR(test = scRNA.harmony.counts, 
                ref = ref, 
                labels = ref$label.main,
                de.method = "wilcox")

# Output predictions
print(pred)

saveRDS(pred, "Analysis/cell_type_pred_singleR_2.rds")


pred <- readRDS("Analysis/cell_type_pred_singleR_2.rds")
# Add SingleR predictions to metadata
scRNA.harmony@meta.data$singleR.labels <- pred$pruned.labels[match(rownames(scRNA.harmony@meta.data), rownames(pred))]


# Verify the label distribution
table(scRNA.harmony@meta.data$singleR.labels)

# Plot UMAP with SingleR labels, including "Unknown"
singleRplot <- DimPlot(scRNA.harmony, group.by = "singleR.labels")

singleRplot_umap <- DimPlot(scRNA.harmony, reduction = "umap_harmony", group.by = "singleR.labels")
ggsave("Results/SingleR_umap.png", plot = singleRplot_umap, width = 10, height = 8, dpi = 300)
singleRplot_umap

# Diagnostics for SingleR annotation
# 1. View scores
print(pred$scores)

# 2. Plot score heatmap
plotScoreHeatmap(pred)

# 3. Plot delta distribution
delta_scores <- plotDeltaDistribution(pred)
delta_scores
ggsave("Results/delta_scores.png", plot = delta_scores, width = 10, height = 8, dpi = 300)

saveRDS(scRNA.harmony, "Analysis/annotation_2.rds")

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

cluster19_conserved_markers <- FindConservedMarkers(annotation,
                                                    ident.1 = 19,
                                                    grouping.var = "age_group",
                                                    only.pos = TRUE,min.pct = 0.25,  min.diff.pct = 0.25,
                                                    logfc.threshold = 0.25)


### AFTER FINDING CONSERVED GENES ### CHANGE THE ONES THAT ARE HERE ### ASK ANTONI ###

# Create a mapping of cluster IDs to cell type names
annotation <- RenameIdents(object = annotation,
                               "0" = "Neurons",
                               "1" = "Glioblasts",
                               "2" = "Erythroid_progenitor_cells",
                               "3" = "Erythroid_progenitor_cells",
                               "4" = "Ciliated_cells",
                               "5" = "CD4+, alpha-beta T cell",
                               "6" = "Stem_cells",
                               "7" = "Neuronal_stem_Cells",
                               "8" = "ZBTB32+ B cell",
                               "9" = "Astrocytes",
                               "10" = "Myeloblasts",
                               "11" = "Glioblasts",
                               "12" = "Microglia",
                               "13" = "Endothelial cells",
                               "14" = "NA",
                               "15" = "Fibroblasts",
                               "16" = "Multi-ciliated epithelial cell",
                               "17" = "CD8+ T cell",
                               "18" = "CD8+ T cell"
)

# Assign cell types based on cluster IDs
annotation$cell_type_manual <- Idents(annotation)

# visualize data

clusters <- DimPlot(annotation, reduction = 'umap_harmony', label = TRUE)
celltype <- DimPlot(annotation, reduction = 'umap_harmony', group.by = 'cell_type_manual', label = TRUE)
ggsave(filename = 'cell_annotation.png', plot = celltype)

# Optionally save the results to a file
saveRDS(annotation, file ="Analysis/annotated_cells.rds")