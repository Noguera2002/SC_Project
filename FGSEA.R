# Load required libraries
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(GSEABase)
library(fgsea)
library(ggplot2)

# 1. Ensure column names in `res` are unique
colnames(res) <- make.unique(colnames(res))  # Ensure unique column names

# 2. Handle existing 'genes' column
if ("genes" %in% colnames(res)) {
  res <- res %>%
    rename(existing_genes = genes)  # Rename the existing 'genes' column
}


# Add row names as 'genes'
res <- res %>%
  rownames_to_column(var = "genes") %>%  # Add row names to a new 'genes' column
  mutate(genes = as.character(genes))   # Ensure 'genes' is a character vector

# 2. Map Ensembl Gene IDs to Entrez IDs
entrez_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res$genes),  # Unique gene identifiers
  columns = c("ENSEMBL", "ENTREZID"),
  keytype = "ENSEMBL"
) %>%
  as_tibble() %>%
  filter(!is.na(ENSEMBL) & !is.na(ENTREZID))  # Remove missing mappings

# Join with `res` and filter out missing data
res <- res %>%
  inner_join(entrez_mapping, by = c("genes" = "ENSEMBL")) %>%
  filter(!is.na(pvalue) & !is.na(log2FoldChange))  # Remove rows with missing stats

# 3. Create ranking metric for GSEA
res <- res %>%
  mutate(stat_sig = -log10(pvalue) * sign(log2FoldChange))  # Combine p-value and direction
rankData <- res$stat_sig
names(rankData) <- res$ENTREZID  # Use Entrez IDs as names
rankData <- sort(rankData, decreasing = TRUE)  # Sort for fgsea input

# 4. Load gene sets
gmt_file <- "Data/brain_geneset.gmt"  # Replace with your actual path
if (!file.exists(gmt_file)) stop("Gene set file not found: ", gmt_file)
genesets <- getGmt(gmt_file)  # Load .gmt file
message("Gene sets loaded successfully.")

# 5. Perform GSEA using fgseaMultilevel
fgseaRes <- fgseaMultilevel(
  pathways = genesets,
  stats = rankData,
  minSize = 15,
  maxSize = 500,
  scoreType = "std"  # Use standard scoring for balanced analysis
)

# 6. Process results
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))  # Sort by Normalized Enrichment Score (NES)

# 7. Visualize top 10 upregulated and bottom 10 downregulated pathways
top_10 <- fgseaResTidy %>%
  slice_max(order_by = NES, n = 10)

bottom_10 <- fgseaResTidy %>%
  slice_min(order_by = NES, n = 10)

top_bottom_10 <- bind_rows(top_10, bottom_10)

# Plot the results
ggplot(top_bottom_10, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj < 0.05), show.legend = TRUE) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray")) +
  labs(
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    title = "Top 10 and Bottom 10 Pathways (GSEA)",
    fill = "Significant (padj < 0.05)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )
