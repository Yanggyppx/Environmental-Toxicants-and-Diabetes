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

#External validation
library(tidyverse)
library(edgeR)
library(DESeq2)
library(FactoMineR)
library(org.Mm.eg.db)
library(stringr)
library(stringi)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(patchwork)
library(ggplotify)
library(GenomicFeatures)
rm(list=ls())
counts<-read.table("GSE203624_Raw_gene_counts_matrix.txt",sep="\t",header=T)
counts<-counts[c("Geneid","D137","D138","D139","D140","D141","D142","D143","D144","T73","T74","T75","T76","T77","T78","T79","T80")]
matrix0 <- order(rowMeans(counts[,-1]),decreasing = T)
expr_ordered <- counts[matrix0,]
keep=!duplicated(expr_ordered$Geneid)
expr_max=expr_ordered[keep,]
counts=expr_max[c("D137","D138","D139","D140","D141","D142","D143","D144","T73","T74","T75","T76","T77","T78","T79","T80")]
rownames(counts)=expr_max$Geneid
head(counts)
txdb <- makeTxDbFromGFF("Rattus_norvegicus.mRatBN7.2.113.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by="gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens[1:10]
gene_length <- sapply(exons_gene_lens,function(x){x})
id_length <- as.data.frame(gene_length)
head(id_length)
dim(id_length)
library("biomaRt")
listMarts()
ensembl <- useMart("ensembl", dataset='rnorvegicus_gene_ensembl',host="https://dec2021.archive.ensembl.org/")
ensembl
my_gene <- rownames(id_length)
rat_symbol<-getBM(
  attributes=c('ensembl_gene_id', 'external_gene_name'),
  filters="ensembl_gene_id",
  values = my_gene,
  mart = ensembl)
head(rat_symbol)
dim(rat_symbol)
rownames(rat_symbol)=rat_symbol$ensembl_gene_id
rat_symbol=rat_symbol[sort(rownames(rat_symbol)),]
id_length$ensembl_gene_id=rownames(id_length)
id_symbol_length=merge(rat_symbol,id_length,by="ensembl_gene_id")
colnames(id_symbol_length)=c("ensembl","Geneid","Length")
geneid_efflen <- subset(id_symbol_length,select = c("Geneid","Length"))
colnames(geneid_efflen) <- c("geneid","efflen")
geneid_efflen=unique(geneid_efflen)
geneid_efflen <- geneid_efflen[!duplicated(geneid_efflen$geneid),]
counts.filters=counts[rownames(counts)%in%geneid_efflen$geneid,]
counts.filters=na.omit(counts.filters)
geneid_efflen=geneid_efflen[geneid_efflen$geneid%in%rownames(counts.filters),]
efflen <- geneid_efflen[match(rownames(counts.filters),geneid_efflen$geneid) , "efflen"]
counts2TPM2 <- function(count=count, efflength=efflen){
	RPK <- count/(efflength/1000) 
	PMSC_rpk <- sum(RPK)/1e6
	RPK/PMSC_rpk}
tpm_raw <- apply(counts.filters, 2, counts2TPM2, efflength = efflen)
head(tpm_raw)
colSums(tpm_raw)
dat <- log2(tpm_raw+1)
head(dat)
######Deseq2
> head(counts.filters)
counts.filters=counts.filters[c("T73","T74","T75","T76","T77","T78","T79","T80","D137","D138","D139","D140","D141","D142","D143","D144")]
Group<-factor(c(rep("Ta",8),rep("DU",8)),levels=c("Ta","DU"))
library(DESeq2)
exp_Count=counts.filters
colData<-data.frame(row.names=colnames(exp_Count),group_list=Group)
dds<-DESeqDataSetFromMatrix(countData=exp_Count,
							colData=colData,
							design=~group_list)
dds$group_list
dds<-DESeq(dds)
##results extracts a result table from a DESeq analysis
res<-results(dds,contrast=c("group_list","DU","Ta"))
resOrdered<-res[order(res$padj),]
head(resOrdered)
DEG<-as.data.frame(resOrdered)
DEG_DESeq2<-na.omit(DEG)
logFC<-1
Pvalue<-0.05
DEG_DESeq2$significant <- ifelse(DEG_DESeq2$padj<0.05& abs(DEG_DESeq2$log2FoldChange) >=1,
						ifelse(DEG_DESeq2$log2FoldChange >=1,"Up","Down"),"Not significant")
library(ggplot2)
DEG=DEG_DESeq2
mytheme<-theme_bw()+
		 theme(legend.key=element_rect(fill='transparent'),
		       legend.background=element_rect(fill='transparent'),
		       legend.position=c(0.15,0.9),
		       legend.title=element_blank(),
		       legend.text=element_text(size=20,margin=margin(t=6)),
		       axis.text.x=element_text(hjust=0.5,size=18),
		       axis.text.y=element_text(size=18),
		       axis.title.x=element_text(size=20),
		       axis.title.y=element_text(size=20),
		       axis.line=element_line(size=1),
		       plot.title=element_text(size=24,hjust=0.5))
Pvalue=0.05
log2FC=1
DEG$Group=DEG$significant
DEG$Group=factor(DEG$Group,levels=c("Up","Not significant","Down"))
volcano=ggplot(DEG,aes(x=log2FoldChange,y=-log10(padj),color=Group))+
			   geom_point(alpha=0.75,size=4)+
			   labs(x=bquote(~Log[2]~"(foldchange)"),
			        y=bquote(~-Log[10]~italic("Adjust P-value")),
			        title="Test")+      scale_colour_manual(name="",values=alpha(c("#EB4232","#d8d8d8","#2DB2EB"),0.7))+
			        scale_x_continuous(limits=c(-8,8),breaks=seq(-8,8,by=2))+		        scale_y_continuous(expand=expansion(add=c(0,0)),limits=c(0,40),breaks=seq(0,40,by=10))+geom_hline(yintercept=c(-log10(Pvalue)),size=0.7,color="black",lty="dashed")+
geom_vline(xintercept=c(-log2FC,log2FC),size=0.7,color="black",lty="dashed")+
			        mytheme
volcano
ggsave(volcano,file="volcano.pdf")
DEG$Symbol<-rownames(DEG)
Up<-filter(DEG,Group=='Up')%>%distinct(Symbol,.keep_all=T)%>%top_n(5,-log10(padj))
Down<-filter(DEG,Group=='Down')%>%distinct(Symbol,.keep_all=T)%>%top_n(5,-log10(padj))
head(Up);head(Down)
library(ggrepel)
##geom_text_repel()
p2<-volcano+
	geom_point(data=Up,aes(x=log2FoldChange,y=-log10(padj)),
	color='#EB4232',size=7.5,alpha=0.2)+
	geom_text_repel(data=Up,aes(x=log2FoldChange,y=-log10(pvalue),label=Symbol),
	seed=23456,color='black',show.legend=FALSE,
	min.segment.length=0,	segment.linetype=1,
	force=2,	force_pull=2,	size=8,	box.padding=unit(2,"lines"),
	point.padding=unit(1,"lines"),	max.overlaps=Inf)
p2
ggsave(p2,file="volcano.upgene.pdf")
p3<-volcano+
	geom_point(data=Down,aes(x=log2FoldChange,y=-log10(padj)),
	color='#EB4232',size=7.5,alpha=0.2)+
	geom_text_repel(data=Up,aes(x=log2FoldChange,y=-log10(pvalue),label=Symbol),
	seed=23456,color='black',show.legend=FALSE,
	min.segment.length=0,	segment.linetype=1,
	force=2,	force_pull=2,	size=8,	box.padding=unit(2,"lines"),
	point.padding=unit(1,"lines"),	max.overlaps=Inf)
p3
ggsave(p3,file="volcano.downgene.pdf")
library(org.Mm.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(ggplot2)
head(deg)
go_ythdf2<-deg$Symbol
go_ythdf2_id_trance <- bitr(go_ythdf2,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Mm.eg.db",drop = T)
f <- as.data.frame(go_ythdf2_id_trance[,2])
colnames(f)[1] <- ''
colnames(f) <- "V1"
EG2Ensembl=toTable(org.Mm.egENSEMBL)
f=f$V1  
geneLists=data.frame(ensembl_id=f)
results=merge(geneLists,EG2Ensembl,by='ensembl_id',all.x=T)
id=na.omit(results$gene_id)
All <- enrichGO(OrgDb="org.Mm.eg.db", gene = id, ont = "ALL",pvalue=0.05) 
GO_readable <- setReadable(All, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
GO_enrich = data.frame(GO_readable) 
head(GO_enrich)
GO_enrich$xais=-log10(GO_enrich$pvalue)
write.csv(GO_enrich,'GSE203624.DEG.GO_enrich.csv')
KEGG <- enrichKEGG(gene= id, organism  = 'mmu',pvalue=1) 
KEGG  <- setReadable(KEGG , OrgDb = org.Mm.eg.db, keyType="ENTREZID")
KEGG_enrich  = data.frame(KEGG) 
write.csv(KEGG_enrich,'GSE203624.DEG.KEGG_enrich.csv')
#plot:
table=data.frame(Description=c('extracellular matrix organization',
	                            'humoral immune response',
                                'cellular respiration',
                            	'negative regulation of lipid storage',
                                'negative regulation of apoptotic signaling pathway',
                                'blood coagulation',
                                'ERK1 and ERK2 cascade',
                                'response to metal ion'),
value=c(4.89024939912296,4.67063963887676,3.94812416880096,3.58639043853712,3.32546188525697,3.1705553261584,3.1385152908998,3.09692874515406))
table$Description=factor(table$Description,levels=c('response to metal ion',
													'ERK1 and ERK2 cascade',
												    'blood coagulation',
												    'negative regulation of apoptotic signaling pathway',
												    'negative regulation of lipid storage',
												    'cellular respiration',
												     'humoral immune response',
												     'extracellular matrix organization'))
g=ggplot(table,aes(x=value,y=Description))+geom_col(fill="#2978b2")+theme_bw() + theme(panel.grid=element_blank())
ggsave('GSE203624.DEG.GO.ggplot.pdf',g,height=7,width=9)
C：intersect.gene
data.tpm=readRDS("GSE203624.depleted_uranium.log2tpm.rds")
library(stringr)
target=data.tpm[rownames(data.tpm)%in%target.gene$Gene_formatted,]
write.table(target,"GSE203624.depleted_uranium.Rat.intersectDEG.csv",sep='\t',col.names=T,row.names=F,quote = FALSE,na='')
library(tidyverse)
library(paletteer)
library(aplot)
library(ggplot2)
gene_order <- c(
  "MAP1LC3A", "TGIF2", "WASL", "MAPK12", "CEBPA", "SPATC1L", "HMGA1", 
  "GSTP1", "CHCHD5", "SULT1A1", "IFRD2", "MPG", "GTPBP6", "EPOR", 
  "ALDH3B2", "BAX", "NR1H2", "CYP27A1", "TNFRSF25", "SELENBP1", 
  "PRDX2", "ST6GALNAC4", "RNF123", "MCOLN1", "MT2A", "ALPL", "SULT1A4", "OSM",
  "GSPT1", "CDKN3", "EIF1", "MAPK1", "PPP2CA", "PDCD6IP", "SRRM1", 
  "MAP1LC3B", "SIRT1", "PARP4", "ACTR2", "YWHAB", "KPNA3", "CASP8", 
  "HBP1", "ATP6V1A", "NAP1L1", "SRSF5", "SERBP1", "PTMA", "HMGB1", 
  "TET2", "RHOT1", "TNFSF10", "POLB", "NAMPT", "IFRD1", "LRRK2", 
  "SULT1B1", "SDCBP", "AQP9", "IQGAP2", "EVI5", "SEC24B", "VPS13B", 
  "TRPM7", "MAPK8", "PIK3R1", "CUL2", "YTHDF3", "MAP2K4", "ATF6", 
  "COMMD10", "MAT2B", "CASP3", "CHURC1", "KRAS", "RSL24D1", "METAP2", 
  "CALM2", "SUMO1", "SCP2", "ING3", "UBE2D2", "PTGES3", "TMEM230", 
  "HSP90AA1", "TARDBP", "PHF10")
gene_order_formatted <- tolower(gene_order)
gene_order_formatted <- paste0(toupper(substr(gene_order_formatted, 1, 1)), 
                               substr(gene_order_formatted, 2, nchar(gene_order_formatted)))
gene_order_final <- gene_order_formatted[gene_order_formatted %in% rownames(target)]
target_subset <- target[gene_order_final, ]
sample_order <- c("T73", "T74", "T75", "T76", "T77", "T78", "T79", "T80",
                  "D137", "D138", "D139", "D140", "D141", "D142", "D143", "D144")
target_subset <- target_subset[, sample_order]
target_subset=as.data.frame(target_subset)
dft <- target_subset %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  group_by(Gene) %>%  
  mutate(Zscore = scale(Expression)[,1]) %>%
  ungroup() %>%
  mutate(Gene = factor(Gene, levels = rev(gene_order_final)), 
         Sample = factor(Sample, levels = sample_order))
p <- ggplot(dft, aes(x = Sample, y = Gene, fill = Zscore)) +
  geom_tile() +
  labs(x = NULL, y = NULL)
mytheme <- theme(  panel.grid = element_blank(),  legend.position = "right",
  legend.text = element_text(size = 18),  legend.title = element_text(size = 18),
  axis.ticks.y = element_blank(),  axis.title = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
  axis.text.y = element_text(colour = 'black', size = 18))
p2 <- p + scale_fill_paletteer_c(
  "ggthemes::Classic Red-Blue",
  direction = -1,  name = "Expression level\n(Z-score)",
  limits = c(-10, 10)) + mytheme
p2 <- p2 + guides(  fill = guide_colorbar(    title.position = "top",
    title.hjust = 0.5,    barwidth = 1.5,    barheight = 15,  ticks = FALSE ))
Group <- data.frame(  Sample_id = sample_order,
  Group = c(rep("Ta", 8), rep("DU", 8)),   Yaxis = "Annotation")
mytheme2 <- theme(
  panel.grid = element_blank(),  legend.position = "right",
  legend.text = element_text(size = 18),  legend.title = element_text(size = 18),
  axis.ticks.y = element_blank(),  axis.title = element_blank(),
  axis.text.x = element_blank(),  axis.text.y = element_blank())
Block1 <- ggplot(data = Group, aes(Sample_id, Yaxis, fill = Group)) +
  geom_tile() +  scale_fill_manual(values = c("#3B9AB2", "#78B7C5")) + 
  theme_void() + mytheme2
p_final <- p2 %>% insert_top(Block1, height = 0.04)
print(p_final)
ggsave("GSE203624_80intersect.heatmap.pdf", p_final, width = 7, height = 7)
######
data.target=target_subset[rownames(target_subset)%in%c("Mapk8","Sirt1","Kras","Pik3r1","Mapk1"),]
library(tidyverse)
library(ggplot2)
library(ggpubr)
gene_expression_log <- data.target
##
target2=data.tpm[rownames(data.tpm)%in%c("Pam","Efemp1","Igfbp6","Igfbp2"),]
gene_expression_log <- target2
##
gene_expression_log <- as.data.frame(gene_expression_log)
zscore_data <- t(apply(gene_expression_log, 1, scale))
colnames(zscore_data) <- colnames(gene_expression_log)
rownames(zscore_data) <- rownames(gene_expression_log)
long_data_zscore <- as.data.frame(zscore_data) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(  cols = -Gene,    names_to = "Sample",
    values_to = "Zscore"  ) %>%
  mutate(    Group = ifelse(Sample %in% c("T73", "T74", "T75", "T76", "T77", "T78", "T79", "T80"), "Ta", "DU"),  Group = factor(Group, levels = c("Ta", "DU"))  )
geneorder <- rownames(gene_expression_log)
summary_data_zscore <- long_data_zscore %>%
  group_by(Gene, Group) %>%
  summarise(    mean_Zscore = mean(Zscore),    se = sd(Zscore) / sqrt(n())
  ) %>%
  ungroup() %>%
  mutate(Group = factor(Group, levels = c("Ta", "DU")))
stat_test_zscore <- long_data_zscore %>%
  group_by(Gene) %>%
  summarise(    p_value = t.test(      Zscore[Group == "Ta"],
      Zscore[Group == "DU"]    )$p.value
  ) %>%
  mutate(    signif = case_when(      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",      p_value < 0.05 ~ "*",  TRUE ~ "ns"    ),
    y.position = max(long_data_zscore$Zscore) * 1.1  )
summary_data_zscore$Gene <- factor(summary_data_zscore$Gene, levels = geneorder)
long_data_zscore$Gene <- factor(long_data_zscore$Gene, levels = geneorder)
stat_test_zscore$Gene <- factor(stat_test_zscore$Gene, levels = geneorder)
pdf("gene_expression_zscore_barplot.pdf", width = 7, height = 5)
ggplot() +
  geom_col(    data = summary_data_zscore,
    aes(x = Gene, y = mean_Zscore, fill = Group),
    position = position_dodge(width = 0.75),    width = 0.75,
    color = "black",    linewidth = 0.4,    alpha = 1  ) +
  geom_errorbar(    data = summary_data_zscore,
    aes(x = Gene, ymin = mean_Zscore - se, ymax = mean_Zscore + se, group = Group),
    position = position_dodge(width = 0.75),    width = 0.5,    color = "black",
    linewidth = 0.5  ) +
  geom_point(    data = long_data_zscore,
    aes(x = Gene, y = Zscore, group = Group, fill = Group),
    position = position_jitterdodge(      dodge.width = 0.75,  jitter.width = 0.2 ),
    size = 2,  shape = 21, color = "black", stroke = 0.4, alpha = 0.8  ) +
  geom_text( data = stat_test_zscore,  aes(x = Gene, y = y.position, label = signif),
    size = 6, position = position_dodge(width = 0.9)
  ) +  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0.2)),
    breaks = seq(-3, 3, by = 1)  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual( values = c("Ta" = "#3B9AB2", "DU" = "#78B7C5"),
    breaks = c("Ta", "DU")  ) +
  labs( x = "Gene Symbol", y = "Z-score (Normalized Expression)",
    title = "Gene Expression Z-scores (Ta group first)"  ) +
  theme_classic() +
  theme( axis.text.x = element_text( angle = 45,hjust = 1, vjust = 1, size = 14,face = "italic"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +  coord_cartesian(clip = "off")
dev.off()
