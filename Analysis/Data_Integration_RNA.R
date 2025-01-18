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

# Perform integration of cells across conditions to enable meaningful cross-sample comparisons
# Perform clustering and marker identification workflow with integrated data
# Identify cell types based on known gene marker expression
# Identify cell types using a previously annotated dataset as reference

scRNA <- readRDS("Data/scRNA_filtered.rds") 

# 1 Lets visualize the data

DimPlot(scRNA, reduction = 'umap', group.by = 'age_group') # group by age_group
DimPlot(scRNA, reduction = 'umap', split.by = 'age_group')



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
which(cumulative_var >= 0.75)[1]

# Number of PCs to use for clustering is 16



### 3. Integrate the Data Using Harmony and CCA

# library(future)
# plan("multisession", workers = 4)
# options(future.globals.maxSize = 8000 * 1024^2)

# a Anchor-Based CCA Integration
scRNA <- IntegrateLayers(
  object = scRNA,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "CCA_Integration",
  verbose = FALSE
)

# b Harmony Integration
scRNA <- IntegrateLayers(
  object = scRNA,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "Harmony_Integration",
  verbose = FALSE
)

# Set the future plan back to "sequential" after running itegration
# plan("sequential")

#Re-join the split layers
scRNA[["RNA"]] <- JoinLayers(scRNA[["RNA"]])

#Save the data
saveRDS(scRNA, "Data/scRNA_integrated.rds")

# Run UMAP on the CCA and Harmony Integrated Reductions
scRNA <- RunUMAP(scRNA, dims = 1:16, reduction = "CCA_Integration", reduction.name = "umap_cca")
scRNA <- RunUMAP(scRNA, dims = 1:16, reduction = "Harmony_Integration", reduction.name = "umap_harmony")

# Create UMAP Plots
p1 <- DimPlot(scRNA, reduction = 'umap', group.by = 'age_group') + ggtitle("UMAP Raw Data")
p2 <- DimPlot(scRNA, reduction = "umap_cca", group.by = "age_group") + ggtitle("UMAP CCA")
p3 <- DimPlot(scRNA, reduction = "umap_harmony", group.by = "age_group") + ggtitle("UMAP Harmony")

# Combine the Plots in a Single Row
wrap_plots(p1, p2, p3, nrow = 1) + plot_layout(guides = "collect")























