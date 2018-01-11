# pro_correlation.R by Kuan Huang @ WashU 201507
# VAF correlation plotting and analysis for WHIM manuscript figure 1

library(VennDiagram)
library(grDevices)
library("PerformanceAnalytics")

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/pro_correlation")
source("/Users/khuang/bin/LIB_exp.R")

ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
LFQ=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified.txt',header=TRUE, sep="\t", fill=T)
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)

ITRAQphoCol = read.table(row.names=1, header=TRUE, sep="\t", file="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt")

##### ITRAQ phosphoproteome #####
pdf(paste(pd,'ITRAQ_phosphoproteome_all_pairwise.pdf', sep="_"), width=20, height=20,useDingbats=FALSE)
chart.Correlation(ITRAQpho[,-c(1,2)])
dev.off()

pdf(paste(pd,'ITRAQ_phosphoproteome_WHIM13_dup.pdf', sep="_"), width=5, height=5,useDingbats=FALSE)
chart.Correlation(ITRAQpho[,c(18,22)])
dev.off()

##### ITRAQ proteome #####
## pre-process ##
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$datETcollapsed
# 12698 genes
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ_dup=ITRAQ.c[,c(16,20)]
pdf(paste(pd,'ITRAQ_proteome_WHIM13_dup.pdf', sep="_"), width=5, height=5,useDingbats=FALSE)
chart.Correlation(ITRAQ_dup)
dev.off()
ITRAQ_proteome = ITRAQ.c[,-c(17,18,20)]
pdf(paste(pd,'ITRAQ_proteome_all_pairwise.pdf', sep="_"), width=20, height=20,useDingbats=FALSE)
chart.Correlation(ITRAQ_proteome)
dev.off()
# ITRAQ_proteome_mad = smad(ITRAQ_proteome)
# ITRAQ_proteome_mad.m=melt(ITRAQ_proteome_mad)
ITRAQ_proteome_norm = normalize(ITRAQ_proteome)
ITRAQ_proteome_norm.m=melt(ITRAQ_proteome_norm)

### compare WHIM8 and WHIM14 proteome/phosphoproteome ###
ITRAQphoCol_2 = ITRAQphoCol[,c("WHIM8","WHIM14")]
ITRAQ.c_2 = ITRAQ.c[,c("WHIM8","WHIM14")]
colnames(ITRAQphoCol_2) = paste(colnames(ITRAQphoCol_2), "PHO", sep="_")
colnames(ITRAQ.c_2) = paste(colnames(ITRAQ.c_2), "PRO", sep="_")
ITRAQ_merge2 = merge(ITRAQphoCol_2,ITRAQ.c_2,by="row.names")
pdf(paste(pd,'ITRAQ_pro_pho_WHIM8_WHIM14.pdf', sep="_"), width=10, height=10,useDingbats=FALSE)
chart.Correlation(ITRAQ_merge2[,-1])
dev.off()

## compare all proteome'phospho
colnames(ITRAQphoCol) = paste(colnames(ITRAQphoCol), "PHO", sep="_")
colnames(ITRAQ.c) = paste(colnames(ITRAQ.c), "PRO", sep="_")
ITRAQ_merge = merge(ITRAQphoCol,ITRAQ.c,by="row.names")
ITRAQ_cor = cor(ITRAQ_merge[,-1],use="pairwise")
library(corrplot)
pdf(paste(pd,'ITRAQ_pro_pho_all.pdf', sep="_"), width=10, height=10,useDingbats=FALSE)
corrplot(ITRAQ_cor,order="AOE")
dev.off()
corrplot(ITRAQ_cor, order="AOE", method="circle", tl.pos="lt", type="upper",        
         tl.col="black", tl.cex=0.6, tl.srt=45, 
         addCoef.col="black", addCoefasPercent = TRUE,
         p.mat = 1-abs(ITRAQ_cor), sig.level=0.50, insig = "blank")

### compare WHIM8 and WHIM14 proteome/phosphoproteome ###
ITRAQphoCol_2 = ITRAQphoCol[,c("WHIM18_PHO","WHIM2_PHO")]
ITRAQ.c_2 = ITRAQ.c[,c("WHIM18_PRO","WHIM2_PRO")]
#colnames(ITRAQphoCol_2) = paste(colnames(ITRAQphoCol_2), "PHO", sep="_")
#colnames(ITRAQ.c_2) = paste(colnames(ITRAQ.c_2), "PRO", sep="_")
ITRAQ_merge2 = merge(ITRAQphoCol_2,ITRAQ.c_2,by="row.names")
pdf(paste(pd,'ITRAQ_pro_pho_WHIM18_WHIM2.pdf', sep="_"), width=10, height=10,useDingbats=FALSE)
chart.Correlation(ITRAQ_merge2[,-1])
dev.off()

##### LFQ proteome #####
row.names(LFQ) = LFQ$RefSeq_id
LFQ.c = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$datETcollapsed #note also got rid of rows seen in <=2 samples; 5802 genes
LFQ.NP = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$group2row
# get rid of WHIM16_2, WHIM2_2
LFQ_dup1=LFQ.c[,c(5,6)]
pdf(paste(pd,'LFQ_proteome_WHIM16_dup.pdf', sep="_"), width=5, height=5,useDingbats=FALSE)
chart.Correlation(LFQ_dup1)
dev.off()

LFQ_dup2=LFQ.c[,c(8,9)]
pdf(paste(pd,'LFQ_proteome_WHIM2_dup.pdf', sep="_"), width=5, height=5,useDingbats=FALSE)
chart.Correlation(LFQ_dup2)
dev.off()
# ITRAQ_dup=ITRAQ.c[,c(16,20)]
# pdf(paste(pd,'ITRAQ_proteome_WHIM13_dup.pdf', sep="_"), width=20, height=20)
# chart.Correlation(ITRAQ_dup)
# dev.off()
LFQ_proteome = LFQ.c[,-c(6,9)]
pdf(paste(pd,'LFQ_proteome_all_pairwise.pdf', sep="_"), width=20, height=20,useDingbats=FALSE)
chart.Correlation(LFQ_proteome)
dev.off()
LFQ_proteome.m = melt(LFQ_proteome)
#ggplot(LFQ_proteome.m, aes(x=value)) + geom_density() + xlim(-10,10) # check what's going on with the double peaks in LFQ
LFQ_proteome_mad = smad(LFQ_proteome)
LFQ_proteome_mad.m=melt(LFQ_proteome_mad)
LFQ_proteome_norm = normalize(LFQ_proteome)
LFQ_proteome_norm.m=melt(LFQ_proteome_norm)
#ggplot(LFQ_proteome_mad.m, aes(x=value)) + geom_density() + xlim(-10,10) # check what's going on with the triple peaks in LFQ
# ggplot(LFQ_proteome_norm.m, aes(x=value)) + geom_density() + xlim(-10,10) 
# ggplot(ITRAQ_proteome_norm.m, aes(x=value)) + geom_density() + xlim(-10,10)

##### Merge and compare ###
ITRAQ_LFQ_merge = merge(ITRAQ_proteome_norm.m,LFQ_proteome_norm.m, by=c("Var1","Var2"))
cor(ITRAQ_LFQ_merge$value.x,ITRAQ_LFQ_merge$value.y, use = 'pairwise.complete.obs', method="pearson")
#[1,] 0.614629
cor(ITRAQ_LFQ_merge$value.x,ITRAQ_LFQ_merge$value.y, use = 'pairwise.complete.obs', method="spearman")
#[1,] 0.6039114
colnames(ITRAQ_LFQ_merge)[1] = "gene"
colnames(ITRAQ_LFQ_merge)[2] = "sample"

fn=paste(pd,'ITRAQ_LFQ_norm_protein_correlation.pdf', sep="_")
p = ggplot(ITRAQ_LFQ_merge, aes(x=value.x, y=value.y, color=sample)) + 
  geom_point(alpha=0.1, size=1) + geom_smooth(method=lm, se=F, alpha=0.2) 
p = p + theme_bw() + #ggtitle("Normalized protein expression correlation between ITRAQ and LFQ proteome: r=0.62 (Pearson)") +
  xlab("ITRAQ protein expression") + ylab("LFQ protein expression") + xlim(-5,5) + ylim(-5,5) + coord_fixed() 
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=16), axis.text.y = element_text(colour="black", size=16))
p
ggsave(file=fn)

#chart.Correlation(ITRAQ_LFQ_merge[,c(3,4)])
pro_list = list(ITRAQ_refseq_proteins = ITRAQ$Description, LFQ_refseq_proteins = LFQ$RefSeq_id)

fn = paste(pd, '_ITRAQ_LFQ_refseq_overlap.pdf', sep = "")
venn.plot = venn.diagram(
  x = pro_list, 
  filename = NULL,
  lwd = 4,
  fill = c("royalblue1", "orange1"),
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  #cat.col = c("royalblue1", "orange1"),
  cat.cex = 1,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = c(0.16, 0.16),
  cat.pos = c(240, 120), 
  cat.default.pos = c("text"), 
  #height=2200, width=2200, 
  margin=0.2
)
fn = paste(figurePath, date, '_ITRAQ_LFQ_protein_overlap.pdf', sep = "")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

gene_list = list(ITRAQ_proteins = unique(ITRAQ$Gene), LFQ_proteins = LFQ$gene_name)
venn.plot = venn.diagram(
  x = gene_list, 
  filename = NULL,
  lwd = 4,
  fill = c("royalblue1", "orange1"),
  alpha = 0.2,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  #cat.col = c("royalblue1", "orange1"),
  cat.cex = 1,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
   cat.dist = c(0.16, 0.16),
   cat.pos = c(240, 120), 
  cat.default.pos = c("text"), 
  #height=2200, width=2200, 
  margin=0.2
)
fn = paste(figurePath, date, '_ITRAQ_LFQ_protein_overlap.pdf', sep = "")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

### add RNA and CNV

RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified',header=TRUE, sep="\t")

gene_list2 = list(ITRAQ_proteome_genes = unique(ITRAQ$Gene), LFQ_proteome_genes = unique(LFQ$gene_name), mRNA_genes=unique(RNA$gene))#, CNV_genes=unique(CNV$gene))
venn.plot = venn.diagram(
  x = gene_list2, 
  filename = NULL,
  lwd = 4,
  #fill = c("royalblue1", "orange1"),
  alpha = 0.2,
  scaled =T,
  label.col = "black",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  #cat.col = c("royalblue1", "orange1"),
  cat.cex = 1,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  #cat.dist = c(0.16, 0.16),
  #cat.pos = c(240, 120), 
  cat.default.pos = c("text"), 
  #height=2200, width=2200, 
  margin=0.2,
)

fn = paste(figurePath, date, '_3_level_overlap.pdf', sep = "")
pdf(file=fn)
grid.draw(venn.plot)
dev.off()

##### Phospho #####