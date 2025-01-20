### Data Integration 

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

set.seed(123)

# Perform integration of cells across conditions to enable meaningful cross-sample comparisons
# Perform clustering and marker identification workflow with integrated data
# Identify cell types based on known gene marker expression
# Identify cell types using a previously annotated dataset as reference

scRNA <- readRDS("Data/scRNA_filtered_no_doublets.rds") 

# 0 Set the Idents to resolution 0.2
Idents(scRNA) <- "RNA_snn_res.0.2"

# 1 Lets visualize the data
umap_age_group <- DimPlot(scRNA, reduction = 'umap', group.by = 'age_group') # group by age_group
ggsave("Results/dimplot_age_group.png", plot = umap_age_group, width = 8, height = 6)
umap_age_group
umap_split_age_group <- DimPlot(scRNA, reduction = 'umap', split.by = 'age_group')
ggsave("Results/dimplot_split_age_group.png", plot = umap_split_age_group, width = 8, height = 6)
umap_split_age_group

### 2. Preprocess the Data for Integration

# Lets split the data based on age_group
scRNA[["RNA"]] <- split(scRNA[["RNA"]], f = scRNA$age_group)

# Normalize, find variable features, scale, and run PCA
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA)

# Choose the number of PCs to use for clustering. (Different than before we split)

#Extract the PCA results (standard deviations of the PCs)
stdev <- scRNA[["pca"]]@stdev

# Calculate the variance explained by each principal component
variance <- stdev^2 / sum(stdev^2)

# Calculate the cumulative variance explained
cumulative_var <- cumsum(variance)

# Find the number of components responsible for 80% of the variability
which(cumulative_var >= 0.80)[1]

# Number of PCs to use for clustering is 16


### 3. Integrate the Data Using Harmony and CCA

library(future)
plan("multisession", workers =4)
options(future.globals.maxSize = 8000 * 1024^2)

# a Anchor-Based CCA Integration
scRNA <- IntegrateLayers(
  object = scRNA,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "CCA_Integration",
  verbose = FALSE
)

#Save the data
saveRDS(scRNA, "Data/scRNA_integrated_no_doublets.rds")

# b Harmony Integration
scRNA <- IntegrateLayers(
  object = scRNA,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "Harmony_Integration",
  verbose = FALSE
)

# Set the future plan back to "sequential" after running itegration
plan("sequential")

#Re-join the split layers
scRNA[["RNA"]] <- JoinLayers(scRNA[["RNA"]])

#Save the data
saveRDS(scRNA, "Data/scRNA_integrated_no_doublets.rds")

# Run UMAP on the CCA and Harmony Integrated Reductions
scRNA <- RunUMAP(scRNA, dims = 1:16, reduction = "CCA_Integration", reduction.name = "umap_cca")
scRNA <- RunUMAP(scRNA, dims = 1:16, reduction = "Harmony_Integration", reduction.name = "umap_harmony")

# Create UMAP Plots
p1 <- DimPlot(scRNA, reduction = 'umap', group.by = 'age_group') + ggtitle("UMAP Raw Data")
ggsave("Results/umap_raw_data.png", plot = p1, width = 8, height = 6)
p2 <- DimPlot(scRNA, reduction = "umap_cca", group.by = "age_group") + ggtitle("UMAP CCA")
ggsave("Results/umap_cca.png", plot = p2, width = 8, height = 6)
p3 <- DimPlot(scRNA, reduction = "umap_harmony", group.by = "age_group") + ggtitle("UMAP Harmony")
ggsave("Results/umap_harmony.png", plot = p3, width = 8, height = 6)

# Combine the Plots in a Single Row
combined_umap_plot <- wrap_plots(p1, p2, p3, nrow = 1) + plot_layout(guides = "collect")
ggsave("Results/combined_umap_plots.png", plot = combined_umap_plot, width = 24, height = 8)
combined_umap_plot

#Save the data
saveRDS(scRNA, "Data/scRNA_integrated_no_doublets.rds")
