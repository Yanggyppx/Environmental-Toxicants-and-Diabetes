# 07_Venn_diagram.R
library(ggVennDiagram)
library(ggvenn)
library(VennDiagram)
library(ggplot2)
library(grid)

#-------------------------------------------------------------------------------
# 1. Load Gene Lists
#-------------------------------------------------------------------------------

cat("===== Loading Gene Lists =====\n")

# Load CTD genes (from exposure analysis)
ctd_genes <- readRDS("results/CTD_gene_lists.rds")
exposure_genes <- ctd_genes$Combined

# Load disease-associated genes (from GEO differential expression)
# Adjust path based on your data
deg_data <- read.table("data/GEO/GSE277813_DEG.csv", 
                       header = TRUE, sep = "\t", row.names = 1)
disease_genes <- rownames(deg_data)

cat("Exposure-related genes (CTD):", length(exposure_genes), "\n")
cat("Disease-related genes (GEO):", length(disease_genes), "\n")

#-------------------------------------------------------------------------------
# 2. Two-Set Venn Diagram (Exposure vs Disease)
#-------------------------------------------------------------------------------

cat("\n===== Creating Two-Set Venn Diagram =====\n")

gene_list_2 <- list(
  "Exposure" = exposure_genes,
  "Diabetes" = disease_genes
)

# Find intersection
intersection_genes <- intersect(exposure_genes, disease_genes)
cat("Intersection genes:", length(intersection_genes), "\n")

# Method 1: ggvenn (simple and clean)
p_venn2 <- ggvenn(
  gene_list_2,
  show_percentage = TRUE,
  fill_color = c("#137CBD", "#D53E51"),
  stroke_size = 0.8,
  text_size = 5,
  set_name_size = 5
) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("figures/Venn_two_sets.pdf", p_venn2, width = 6, height = 5)

#-------------------------------------------------------------------------------
# 3. Four-Set Venn Diagram (Individual Chemicals)
#-------------------------------------------------------------------------------

cat("\n===== Creating Four-Set Venn Diagram =====\n")

gene_list_4 <- list(
  Antimony = ctd_genes$Antimony,
  Glycidamide = ctd_genes$Glycidamide,
  EthyleneOxide = ctd_genes$EthyleneOxide,
  Uranium = ctd_genes$Uranium
)

p_venn4 <- ggVennDiagram(
  gene_list_4,
  category.names = c("Sb", "Glycidamide", "EthyleneOxide", "U"),
  set_color = "black",
  set_size = 4,
  label = "both",
  label_alpha = 0,
  label_size = 4,
  edge_size = 1
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "") +
  theme(
    legend.position = "none",
    text = element_text(family = "serif")
  )

ggsave("figures/Venn_four_chemicals.pdf", p_venn4, width = 8, height = 7)

#-------------------------------------------------------------------------------
# 4. High-Quality VennDiagram Package Version
#-------------------------------------------------------------------------------

cat("\n===== Creating Publication-Quality Venn Diagram =====\n")

venn_plot <- venn.diagram(
  x = gene_list_2,
  category.names = c("Exposure-related\nGenes", "Diabetes-related\nGenes"),
  filename = NULL,
  
  # Colors

fill = c("#E63946", "#457B9D"),
  alpha = 0.5,
  
  # Borders
  lwd = 2,
  col = c("#C1121F", "#1D3557"),
  
  # Labels
  cex = 1.8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Category labels
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  
  # Other settings
  margin = 0.1,
  scaled = FALSE
)

pdf("figures/Venn_publication.pdf", width = 8, height = 6)
grid.draw(venn_plot)
dev.off()

#-------------------------------------------------------------------------------
# 5. Save Intersection Genes
#-------------------------------------------------------------------------------

cat("\n===== Saving Results =====\n")

# Save intersection genes
intersection_df <- data.frame(Gene = intersection_genes)
writexl::write_xlsx(intersection_df, "results/intersection_genes.xlsx")
write.csv(intersection_df, "results/intersection_genes.csv", row.names = FALSE)

# Create summary statistics
venn_summary <- data.frame(
  Category = c("Exposure only", "Disease only", "Intersection", "Total unique"),
  Count = c(
    length(setdiff(exposure_genes, disease_genes)),
    length(setdiff(disease_genes, exposure_genes)),
    length(intersection_genes),
    length(union(exposure_genes, disease_genes))
  )
)

write.csv(venn_summary, "results/Venn_summary.csv", row.names = FALSE)

cat("\nVenn Diagram Summary:\n")
print(venn_summary)

cat("\n===== Venn Diagram Analysis Completed =====\n")
cat("Intersection genes saved:", length(intersection_genes), "genes\n")