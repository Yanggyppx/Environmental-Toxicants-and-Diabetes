# 06_CTD_gene_extraction.R
library(readxl)
library(dplyr)

#-------------------------------------------------------------------------------
# 1. Set Working Directory and Load CTD Data
#-------------------------------------------------------------------------------

# Note: CTD data should be downloaded from https://ctdbase.org/
# Search each chemical and download the "Gene Interactions" file

cat("===== Loading CTD Gene Interaction Data =====\n\n")

#-------------------------------------------------------------------------------
# 2. Load and Process Each Chemical
#-------------------------------------------------------------------------------

# Antimony (Sb)
cat("Processing: Antimony (Sb)... ")
sb_data <- read_excel("data/CTD/Sb-ctd.xlsx")
sb_data <- sb_data[sb_data$`Interaction Count` >= 1, ]
sb_genes <- unique(sb_data$`Gene Symbol`)
cat(length(sb_genes), "genes\n")

# Uranium (U)
cat("Processing: Uranium (U)... ")
u_data <- read_excel("data/CTD/Uranium-ctd.xlsx")
u_data <- u_data[u_data$`Interaction Count` >= 1, ]
u_genes <- unique(u_data$`Gene Symbol`)
cat(length(u_genes), "genes\n")

# Glycidamide
cat("Processing: Glycidamide... ")
g_data <- read_excel("data/CTD/Glycidamide-ctd.xlsx")
g_data <- g_data[g_data$`Interaction Count` >= 1, ]
g_genes <- unique(g_data$`Gene Symbol`)
cat(length(g_genes), "genes\n")

# Ethylene Oxide
cat("Processing: Ethylene Oxide... ")
e_data <- read_excel("data/CTD/EthyleneOxide-ctd.xlsx")
e_data <- e_data[e_data$`Interaction Count` >= 1, ]
e_genes <- unique(e_data$`Gene Symbol`)
cat(length(e_genes), "genes\n")

#-------------------------------------------------------------------------------
# 3. Combine All Genes
#-------------------------------------------------------------------------------

cat("\n===== Combining Gene Lists =====\n")

# Union of all genes
all_genes <- unique(c(sb_genes, u_genes, g_genes, e_genes))
cat("Total unique genes:", length(all_genes), "\n")

# Create summary table
gene_summary <- data.frame(
  Chemical = c("Antimony", "Uranium", "Glycidamide", "Ethylene Oxide", "Combined"),
  Gene_Count = c(length(sb_genes), length(u_genes), length(g_genes), 
                 length(e_genes), length(all_genes))
)

print(gene_summary)

#-------------------------------------------------------------------------------
# 4. Create Gene-Chemical Mapping Matrix
#-------------------------------------------------------------------------------

gene_chemical_matrix <- data.frame(
  Gene = all_genes,
  Antimony = all_genes %in% sb_genes,
  Uranium = all_genes %in% u_genes,
  Glycidamide = all_genes %in% g_genes,
  EthyleneOxide = all_genes %in% e_genes
)

# Count how many chemicals each gene is associated with
gene_chemical_matrix$Chemical_Count <- rowSums(gene_chemical_matrix[, 2:5])

# Sort by chemical count
gene_chemical_matrix <- gene_chemical_matrix %>% 
  arrange(desc(Chemical_Count))

cat("\nGenes associated with multiple chemicals:\n")
print(head(gene_chemical_matrix[gene_chemical_matrix$Chemical_Count > 1, ], 20))

#-------------------------------------------------------------------------------
# 5. Save Results
#-------------------------------------------------------------------------------

# Save individual gene lists
gene_lists <- list(
  Antimony = sb_genes,
  Uranium = u_genes,
  Glycidamide = g_genes,
  EthyleneOxide = e_genes,
  Combined = all_genes
)

saveRDS(gene_lists, "results/CTD_gene_lists.rds")

# Save combined gene list
write.csv(data.frame(Gene = all_genes), 
          "results/CTD_combined_genes.csv", row.names = FALSE)

# Save gene-chemical matrix
write.csv(gene_chemical_matrix, 
          "results/CTD_gene_chemical_matrix.csv", row.names = FALSE)

# Save summary
write.csv(gene_summary, "results/CTD_summary.csv", row.names = FALSE)

cat("\n===== CTD Gene Extraction Completed =====\n")
cat("Combined gene list saved with", length(all_genes), "genes\n")