BiocManager::install("slingshot")
BiocManager::install("tradeSeq")
BiocManager::install("TSCAN")

suppressPackageStartupMessages({
  library(slingshot)
  library(Seurat)
  library(SingleCellExperiment)
  library(tradeSeq)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
  
})

annotated_data <- readRDS("./Analysis/annotation.rds")



#Question: How do progenitor populations like RG (Radial Glia) or IPC (Intermediate Progenitor Cells)
#transition into mature neuronal and glial cell types (e.g., EN, IN-MGE, IN-CGE, Astrocytes, Oligodendrocytes) across age groups?

HSC_CD34 <- subset(annotated_data, subset = singleR.labels == "Embryonic_stem_cells")

#Remove layers
HSC_CD34@assays$RNA@layers$scale.data.childhood <- NULL
HSC_CD34@assays$RNA@layers$scale.data.adulthood <- NULL
HSC_CD34@assays$RNA@layers$scale.data.infancy <- NULL
HSC_CD34@assays$RNA@layers$scale.data.adolescence <- NULL
HSC_CD34@assays$RNA@layers$scale.data$early_fetal <- NULL
HSC_CD34@assays$RNA@layers$scale.data$late_fetal <- NULL
HSC_CD34@assays$RNA@layers$scale.data <- NULL


#remove reductions
HSC_CD34@reductions <- list()

HSC_CD34 <- FindVariableFeatures(HSC_CD34)
HSC_CD34 <- ScaleData(HSC_CD34)
HSC_CD34 <- RunPCA(HSC_CD34, npcs = 20)

ElbowPlot(HSC_CD34)

# We decided to use 12 PCs for following analysis
HSC_CD34 <- FindNeighbors(HSC_CD34, dims = 1:20)


HSC_CD34 <- FindClusters(HSC_CD34, resolution = c(0.3, 0.5, 0.7))
HSC_CD34 <- RunUMAP(HSC_CD34, dims = 1:20)


DimPlot(HSC_CD34, reduction = "umap", group.by = "RNA_snn_res.0.5", label = TRUE)

Idents(HSC_CD34) <- "RNA_snn_res.0.5"


pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

dimred <- HSC_CD34@reductions$umap@cell.embeddings
clustering <- HSC_CD34$RNA_snn_res.0.5

var_features<- VariableFeatures(HSC_CD34)
counts<- as.matrix(HSC_CD34@assays$RNA$counts[var_features, ])


set.seed(123)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "1")#define where to start the trajectories

lineages



# Plot the lineages
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.7, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(lineages), lwd=2, type = 'lineages', col = c("black"))


# Find differentially expressed genes between the lineages


# We will use the 200 most variable genes instead of the 200 previously calculated for demonstration purposes and time constraint. 
RG <- FindVariableFeatures(RG, nfeatures = 200)
var_features<- VariableFeatures(RG)
# Subset the count matrix with these variable features
counts <- as.matrix(RG@assays$RNA$counts[var_features, ])

curves <- getCurves(SlingshotDataSet(lineages), approx_points = 500, thresh = 0.01, stretch = 0.8, allow.breaks = TRUE, shrink = 0.99)
curves


library(ggplot2)
set.seed(42)

# Two ways to provide the required input to fitGAM.
# Approach 1
sce <- fitGAM(counts = as.matrix(counts), sds = curves)

plotGeneCount(curves, counts, clusters = clustering, models = sce)



different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)

# helper function for plotting a gene's differential expression
plot_differential_expression <- function(feature_id) {
  patchwork::wrap_plots(
    plotGeneCount(curves, counts, clusters = clustering, models = sce, gene = feature_id) + theme(legend.position = "none"),
    plotSmoothers(sce, counts, gene = feature_id)
  )
}

# GZMB, which stands for Granzyme B, is primarily expressed by cytotoxic T lymphocytes, so it follows the biology that it would be expressed in more differentiated T cells.
print(feature_id)


plot_differential_expression(feature_id)



library(scater)
RG$clusters <- Idents(RG)
sce <- as.SingleCellExperiment(RG)
by.cluster <- aggregateAcrossCells(sce, ids=RG$RNA_snn_res.0.5)
centroids <- reducedDim(by.cluster, "UMAP")
# Set clusters=NULL as we have already aggregated above.
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst


line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")

plotUMAP(sce, colour_by="clusters") + 
  geom_line(data=line.data, mapping=aes(x=umap_1, y=umap_2, group=edge))



map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="UMAP", clusters = sce$clusters)

tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)

#Visualize the orderings
common.pseudo <- averagePseudotime(tscan.pseudo) 
plotUMAP(sce, colour_by=I(common.pseudo), text_by="clusters", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=umap_1, y=umap_2, group=edge))


saveRDS(annotated_data, "trajectory.rds")

