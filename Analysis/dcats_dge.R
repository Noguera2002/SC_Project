# Load required libraries
suppressPackageStartupMessages({
  library(DCATS)
  library(Seurat)
  library(SingleCellExperiment)
  library(clusterProfiler)
  library(dplyr)
  library(enrichplot)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(tidyverse)
  library(DESeq2)
})

# Load the dataset
data_dge <- readRDS("annotated_cells.rds")

# View metadata to identify columns for conditions and cell types
head(data_dge@meta.data)
unique(data_dge@meta.data$age_group)
table(data_dge@meta.data$age_group)


# Create a unique identifier for each cell
data_dge$id <- paste0(data_dge$age_group, "_", data_dge$donor_id)
unique(data_dge@meta.data$donor_id)

# Convert Seurat object to SingleCellExperiment object
sce <- as.SingleCellExperiment(data_dge)

# Estimate the similarity matrix using the SNN graph and cell type annotations
knn_mat <- knn_simMat(
  data_dge@graphs$RNA_snn,
  data_dge@meta.data$cell_type_manual
)

# Create a contingency table of sample IDs and cell types
count_mat <- table(data_dge$id, data_dge$cell_type_manual)

# Create the design matrix for age groups
sim_design <- data.frame(condition = data_dge$age_group)
sim_design <- data.frame(
  condition = data_dge@meta.data$age_group[match(rownames(count_mat), data_dge$id)]
)
rownames(sim_design) <- rownames(count_mat)

# Map age groups to numeric values
numeric_mapping <- c(
  "early fetal" = 1,
  "late fetal" = 2,
  "infancy" = 3,
  "childhood" = 4,
  "adolescence" = 5,
  "adulthood" = 6
)
sim_design$condition <- as.numeric(factor(sim_design$condition, levels = names(numeric_mapping)))

# Run DCATS analysis
results <- dcats_GLM(
  count_mat,
  sim_design,
  knn_mat
)
print(results)

---------------------------------------------------------------------------------

  
  
  
  

  # Load the RDS file
data_dge <- readRDS("./annotated_cells.rds")

# Ensure meta.data is a data.frame
meta.data <- data_dge@meta.data

# Step 1: Calculate abundance and apply closure to 100
abundance <- meta.data %>%
  group_by(age_group, cell_type_manual) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(age_group) %>%
  mutate(relative_abundance_closure = (count / sum(count, na.rm = TRUE)) * 100) %>%
  ungroup()

# Step 2: Calculate the geometric mean (GM) of the closed abundance
# Function to calculate geometric mean
gm <- function(x) exp(mean(log(x[x > 0])))  # Exclude zeros to avoid log issues

# Step 3: Apply GM-based recalculation
abundance <- abundance %>%
  group_by(age_group) %>%
  mutate(gm_value = gm(relative_abundance_closure),
         final_relative_abundance = (relative_abundance_closure / gm_value) * 100) %>%
  ungroup()

# Print the results
print(abundance)


library(ggplot2)

# Stacked barplot for relative_abundance_closure
ggplot(abundance, aes(x = age_group, y = relative_abundance_closure, fill = cell_type_manual)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Stacked Barplot of Relative Abundance (Closure to 100%)",
    x = "Age Group",
    y = "Relative Abundance (%)",
    fill = "Cell Type"
  ) +
  theme_minimal()

ggsave(filename = 'compositional_plot.png', plot = celltype)

#___________________________________________________-


data_dge <- readRDS("./annotated_cells.rds")
unique(data_dge@meta.data$donor_id)
unique(data_dge@meta.data$age_group)

# Extract sample names from counts
# Create metadata from the column names of counts
# Define the age_groups and their new classes
#data_dge@meta.data$age_class <- ifelse(data_dge@meta.data$age_group %in% 
                                         #c("early fetal", "late fetal", "infancy"), 
                                       #"young", 
                                       #"old")

data_dge@meta.data


library(DESeq2)

# Aggregate counts to sample level
counts <- AggregateExpression(
  data_dge, 
  group.by = c("cell_type_manual", "donor_id"),
  assays = "RNA",
  return.seurat = FALSE
)
counts <- counts$RNA

counts


# Transpose and process the count matrix
counts.t <- t(counts)
counts.t <- as.data.frame(counts.t)

# Split the data frame by condition

# get values where to split
splitRows <- gsub('_.*', '', rownames(counts.t))

# split data.frame
cts.split <- split.data.frame(counts.t,
                              f = factor(splitRows))

# Modify and transpose each subset
cts.split.modified <- lapply(cts.split, function(x) {
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})

# Extract a specific cell type
counts_neuroblast <- cts.split.modified$`Oligodendrocyte precursor cell`

# Generate sample-level metadata

# Generate sample-level metadata
colData <- data.frame(samples = colnames(counts_neuroblast))

# Classify samples into "young" and "old"
colData <- colData %>%
  mutate(condition = ifelse(grepl("EaFet|LaFet|Inf|Chil", samples), "young", "old")) %>%
  mutate(condition = factor(condition, levels = c("young", "old"))) %>% # Convert to factor
  column_to_rownames(var = "samples") # Set samples as rownames


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_neuroblast,
  colData = colData,
  design = ~ condition
)

dds
# Filter genes with fewer than 10 reads
# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Check the coefficients for the comparison
resultsNames(dds)

vst_data <- vst(dds, blind = TRUE)
pca_data <- assay(vst_data)
pca_results <- prcomp(t(pca_data), scale. = TRUE)
summary(pca_results)


library(ggplot2)
pca_df <- as.data.frame(pca_results$x)
pca_df$condition <- colData$condition

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of RNA-Seq Data", x = "PC1", y = "PC2") +
  theme_minimal()


# Generate results object
res <- results(dds, name = "condition_old_vs_young")

summary(res)


# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 


# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


## Order results by log fold change
top20_sig_genes <- sig_res %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::pull(gene) %>%
  head(n=20)



# Step 8: Visualize the results (MA plot)
plotMA(res, main = "DESeq2: MA Plot", ylim = c(-5, 5))

# Step 9: Volcano plot
library(ggplot2)
ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("gray", "red"))
