#12_GEO_analysis.R
#convertFPKMtoTPM
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(stringr)
library(tidyverse)

exp_FPKM=read.table("GSE277813_99_expression_matrix.txt",sep="\t",header=T)
FPKMToTPM<-function(fpkm)
{
	exp(log(fpkm)-log(sum(fpkm))+log(1e6))
}
exprSetTPM<-apply(exp_FPKM,2,FPKMToTPM)
colSums(exprSetTPM)
k=round(exp_TPM[,1],1)==round(exprSetTPM[,1],1);table(k)
summary(as.vector(exprSetTPM))
expr_norm<- normalizeBetweenArrays(exprSetTPM, method ="quantile")
sample_names <- colnames(expr_norm)
group<- ifelse(grepl("Control", sample_names),"Control",
				ifelse(grepl("T2DM", sample_names),"T2DM","Others"))
group<- factor(group, levels = c("Control","T2DM")) 
table(group)
design<- model.matrix(~0+ group)
colnames(design) <- c("Control","T2DM") 
contrast_matrix <- makeContrasts(T2DM_vs_Control= T2DM - Control,levels = design)
fit <- lmFit(expr_norm, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
deg_T2DM_vs_Control <- topTable(fit2,
						   coef ="T2DM_vs_Control",
						   number = Inf,
						   adjust.method ="BH",
						   sort.by ="p")
add_deg_labels <- function(deg_df, logfc_threshold =1.2, padj_threshold =0.05) {
	deg_df%>%
	mutate(
		label = case_when(
			FC >= fc_threshold & adj.P.Val < padj_threshold ~ "up",
			FC <= -fc_threshold & adj.P.Val < padj_threshold~ "down",
			TRUE~ "nosig")
		)
	}
deg_T2DM_vs_Control <- add_deg_labels(deg_T2DM_vs_Control)
write.csv(deg_T2DM_vs_Control,"T2DM_vs_Control_DEGs.csv", row.names = TRUE)

deg_data=deg_T2DM_vs_Control

comparisons = list(c("Control", "T2DM"))
deg_data$significance <- "Not Significant"
deg_data$significance[deg_data$padj < 0.05 & deg_data$FoldChange > 1.2] <- "Upregulated"
deg_data$significance[deg_data$padj < 0.05 & deg_data$FoldChange < -1.2] <- "Downregulated"
top_genes <- deg_data[deg_data$padj < 0.05 & abs(deg_data$log2FoldChange) > 2, ]
if(nrow(top_genes) > 10) {
    top_genes <- top_genes[order(top_genes$padj), ][1:10, ]
}
#Volcano_Plots
p <- ggplot(deg_data, aes(x = log2FoldChange, y = -log10(P-value), color = significance)) +
     geom_point(alpha = 0.6, size = 1) +
     scale_color_manual(values = c("Downregulated" = "#6fa9b7", 
                                 "Not Significant" = "grey", 
                                 "Upregulated" = "#e89db3")) +
     labs(title = title,
         x = "log2 Fold Change",
         y = "-log10(P-value)") +
     theme_classic() +
     theme(legend.position = "bottom") +
     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
     geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50")
     if(nrow(top_genes) > 0) {
    	p <- p + geom_text(data = top_genes, 
    		aes(label = gene_id), 
                size = 3, vjust = 1, hjust = 0.5, check_overlap = TRUE)
            }

ggsave(p,"Volcano_Plots.pdf", width = 7, height = 7)
#heatmap_Plots
deg_data
sig_res_df <- deg_data[which(deg_data$adj.P.Val< 0.05 & abs(deg_data$FoldChange) >1.2), ]
interested_genes <-row.names(sig_res_df)
heatmap_matrix <- expr_norm[interested_genes, ]
annotation_col<-data.frame(Group=group)
rownames(annotation_col)<-colnames(expr_norm)
pdf("heatmap.pdf",height=12,width=7)
pheatmap(heatmap_matrix,
		 scale ="row", 
		 annotation_col = annotation_col,
		 cluster_rows = TRUE,
		 cluster_cols = FALSE,
		 show_rownames = FALSE,
		 show_colnames = FALSE,
		 fontsize_row = 8,
		 fontsize_col = 10, 
		 main =" Control vs Dibetes (Blood) ",
		 border_color ="NA", 
		 color = colorRampPalette(c("#6fa9b7","white","#e89db3"))(100))
dev.off()
