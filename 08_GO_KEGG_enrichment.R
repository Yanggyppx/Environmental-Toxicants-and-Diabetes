# 08_GO_KEGG_enrichment.R
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

#-------------------------------------------------------------------------------
# 1. Load Intersection Genes
#-------------------------------------------------------------------------------

genes_df <- readxl::read_xlsx("results/intersection_genes.xlsx")
gene_symbols <- genes_df$Gene

cat("Number of genes for enrichment:", length(gene_symbols), "\n")

#-------------------------------------------------------------------------------
# 2. Convert Gene Symbols to Entrez IDs
#-------------------------------------------------------------------------------

gene_mapping <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

entrez_ids <- unique(gene_mapping$ENTREZID)
cat("Mapped Entrez IDs:", length(entrez_ids), "\n")

#-------------------------------------------------------------------------------
# 3. GO Enrichment Analysis
#-------------------------------------------------------------------------------

GO_results <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

go_df <- as.data.frame(GO_results)
go_sig <- go_df[go_df$p.adjust < 0.05, ]

cat("Significant GO terms:", nrow(go_sig), "\n")

# Select top 5 per category
go_top <- go_sig %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5)

#-------------------------------------------------------------------------------
# 4. Visualize GO Results
#-------------------------------------------------------------------------------

go_top$GeneRatio <- sapply(go_top$GeneRatio, function(x) eval(parse(text = x)))

p_go <- ggplot(go_top, aes(x = GeneRatio, y = reorder(Description, -log10(pvalue)))) +
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_color_gradient(low = "#137CBD", high = "#D53E51") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  labs(x = "Gene Ratio", y = "", size = "Gene Count", color = "-log10(p-value)") +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("figures/GO_enrichment.pdf", p_go, width = 9, height = 6)

#-------------------------------------------------------------------------------
# 5. KEGG Enrichment Analysis
#-------------------------------------------------------------------------------

KEGG_results <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

kegg_df <- as.data.frame(KEGG_results)
kegg_sig <- kegg_df[kegg_df$p.adjust < 0.05, ]

cat("Significant KEGG pathways:", nrow(kegg_sig), "\n")

# Visualize top 10 KEGG pathways
kegg_top <- head(kegg_sig, 10)
kegg_top$GeneRatio <- sapply(kegg_top$GeneRatio, function(x) eval(parse(text = x)))

p_kegg <- ggplot(kegg_top, aes(x = reorder(Description, GeneRatio), y = GeneRatio)) +
  geom_col(aes(fill = -log10(pvalue)), width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "#A0BCC2", high = "#D8A7B1") +
  labs(x = "", y = "Gene Ratio", fill = "-log10(p)") +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text.y = element_text(size = 10)
  )

ggsave("figures/KEGG_enrichment.pdf", p_kegg, width = 8, height = 5)

#-------------------------------------------------------------------------------
# 6. Save Results
#-------------------------------------------------------------------------------

writexl::write_xlsx(go_sig, "results/GO_enrichment_results.xlsx")
writexl::write_xlsx(kegg_sig, "results/KEGG_enrichment_results.xlsx")

cat("\nEnrichment analysis completed!\n")