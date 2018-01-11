# cluster.R by Kuan Huang @ WashU 201506
# cluster WHIMs based on (1) iTRAQ (2) LFQ, (1) proteome (2) phosphoproteome 
# cluster WHIMs based on RNA-Seq, CNV

# dependencies
library(ggplot2)
library(reshape)
library(RColorBrewer)
library("gplots")
library(matrixStats)

# others
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/clustering")
source("/Users/khuang/bin/LIB_exp.R")

## function

get.clinical.scale = function() {
  # Set1 colors
  #colors = c(NA, "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999") # positive is red
  # colors = c(NA, "#636363", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey
  # use Perou's intrinsic subtype colors instead
  colors = c(NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#662506", "#A65628", "#F781BF", "#999999") #positive is dark grey       
  color.names = c("negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB","Lymphoepithioma")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="status", values=colors)
  
  return(clinical.color.scale)
}

# plot_clin: plot clinical data based on order from the corresponding heatmap2 clustering plot
plot_clin = function(matrix, clus_order, clin, figure){
  order = as.data.frame(colnames(matrix[,clus_order]))
  colnames(order) = colnames(clin)[1]
  # add missing samples to clin through merge
  m = merge(order,clin, by = colnames(clin)[1], all.x = TRUE)
  row.names(m)=m[,colnames(clin)[1]]
  # reorder and melt for hmting
  m = m[colnames(matrix[,clus_order]),]
  m.m = melt(m, id=colnames(clin)[1])
  m.m$ExternalIdentifierName<-with(m.m,factor(ExternalIdentifierName,levels = colnames(matrix[,clus_order])))
  color.scale = get.clinical.scale()
  
  colourCount=length(unique(m.m$value))
  p = ggplot(data=m.m, aes(x=ExternalIdentifierName, y=variable)) + geom_tile(aes(fill = value)) + 
    scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(colourCount)) +
    xlab("Sample") + ylab("")  + theme_bw() + guides(fill=guide_legend(title="status",nrow=2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=14), axis.text.y = element_text(colour="black", size=14), legend.position="bottom")
  p = p + color.scale
  p
  ggsave(file=figure, useDingbats=FALSE)
} 

# proteome and phosphoproteome files to cluster
ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
LFQ=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
RSEM=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt',header=TRUE, sep="\t")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
#row.names(clin)=clin$ExternalIdentifierName
clin=clin[,-9]
#clin = t(clin[,1:5])
#write.table(clin, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified_cleaned.txt', quote=F, sep = '\t')

Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])

### process and cluster ITRAQ proteome
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ_m = as.matrix(ITRAQ[,-c(1,2,3)])
ITRAQ_m = ITRAQ_m[,-c(17,18,20)]

q99=quantile(ITRAQ_m, probs=0.99, na.rm=T)
q1=quantile(ITRAQ_m, probs=0.01, na.rm=T)
ITRAQ_m2 = matrix(,nrow=dim(ITRAQ_m)[1],ncol=dim(ITRAQ_m)[2])
colnames(ITRAQ_m2)=colnames(ITRAQ_m)
row.names(ITRAQ_m2)=row.names(ITRAQ_m)
for (i in 1:nrow(ITRAQ_m)){
   if ( (sum(ITRAQ_m[i,][!is.na(ITRAQ_m[i,])] > q99) + sum(ITRAQ_m[i,][!is.na(ITRAQ_m[i,])] < q1)) < 1){
     ITRAQ_m2[i,]=ITRAQ_m[i,]
   }
}
ITRAQ_m2 = ITRAQ_m2[rowSums(is.na(ITRAQ_m2)) <= 10,] # 7574 observations
SD=rowSds(ITRAQ_m2, na.rm=TRUE)
ITRAQ_m2 = ITRAQ_m2[SD>2.15,] # 344 observations

# from directly plotting the subtypes
# samples=colnames(ITRAQ_m2)
# samples[samples %in% Basal]="forestgreen"
# samples[samples %in% LumB]="orange"
# samples[samples %in% Her2]="purple"
# samples[samples %in% CLDN_low]="blue"
# 
# pdf(paste(date,'ITRAQ_proteome-ratio-norm_naMax10_SD1.5_mid98.pdf', sep="_"), height=15, width=10)
# par(oma=c(1,2,1,3))
# ITRAQ_m2_hm = heatmap.2(ITRAQ_m2, trace="none",na.color="white", notecol="black",
#                                cexRow=0.8,cexCol=1.5, ColSideColors = samples, 
#                                labRow=NA,col=getPalette, margins=c(5,5)) #
# par(lend = 1)  
# legend("topright",    # location of the legend on the heatmap plot
#        legend = c("Basal", "Luminal B", "Her2-E", "CLDN_low"), # category labels
#        col = c("forestgreen", "orange", "purple", "blue"),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )
# dev.off()

pdf(paste(pd,'ITRAQ_proteome-ratio-norm_naMax10_SD2.15.pdf', sep="_"))
par(oma=c(1,2,1,3))
ITRAQ_m2_hm = heatmap.2(ITRAQ_m2, trace="none", na.color="white", notecol="black",
                      cexRow=0.8,cexCol=1.5, scale="none",dendrogram='column', 
                      labRow=NA,col=getPalette, margins=c(5,5))
clus_order1 = ITRAQ_m2_hm$colInd
dev.off()

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_clin_proteome-ratio-norm_naMax10_SD2.15.pdf', sep="_")
plot_clin(ITRAQ_m2, clus_order1, clin, figure)

##### no basal #####
ITRAQ_m3 = ITRAQ_m2[,!(colnames(ITRAQ_m2) %in% Basal)]
pdf(paste(pd,'ITRAQ_proteome-ratio-norm_naMax10_SD2.15_noBasal.pdf', sep="_"))
par(oma=c(1,2,1,3))
ITRAQ_m3_hm = heatmap.2(ITRAQ_m3, trace="none", na.color="white", notecol="black",
                        cexRow=0.8,cexCol=1.5, scale="none",dendrogram='column', 
                        labRow=NA,col=getPalette, margins=c(5,5))
clus_order = ITRAQ_m3_hm$colInd
dev.off()

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_clin_proteome-ratio-norm_naMax10_SD2.15_noBasal.pdf', sep="_")
plot_clin(ITRAQ_m3, clus_order, clin, figure)

### process and cluster LFQ proteome
row.names(LFQ) = LFQ$RefSeq_id
LFQ = LFQ[,-c(21:26)]
LFQ= LFQ[,-c(6,9)]
LFQ_m = as.matrix(LFQ)
q99=quantile(LFQ_m, probs=0.99, na.rm=T)
q1=quantile(LFQ_m, probs=0.01, na.rm=T)
LFQ_m2 = matrix(,nrow=dim(LFQ_m)[1],ncol=dim(LFQ_m)[2])
colnames(LFQ_m2)=colnames(LFQ_m)
row.names(LFQ_m2)=row.names(LFQ_m)
for (i in 1:nrow(LFQ_m)){
  if ( (sum(LFQ_m[i,][!is.na(LFQ_m[i,])] > q99) + sum(LFQ_m[i,][!is.na(LFQ_m[i,])] < q1)) < 1){
    LFQ_m2[i,]=LFQ_m[i,]
  }
}

LFQ_m2 = LFQ_m2[rowSums(is.na(LFQ_m2)) <= 10,]
SD=rowSds(LFQ_m2, na.rm=TRUE)
LFQ_m2 = LFQ_m2[SD>1.5,] #312

pdf(paste(pd,'LFQ_proteome_naMax10_SD1.5.pdf.pdf', sep="_"))
par(oma=c(1,2,1,3))
LFQ_m2_hm = heatmap.2(LFQ_m2, trace="none",na.color="white", notecol="black",
                                              cexRow=0.8,cexCol=1.5, scale="none",dendrogram='column', 
                                              labRow=NA,col=getPalette, margins=c(5,5)) #
clus_order2 = LFQ_m2_hm$colInd
dev.off()

## plot corresponding clinical panel 
figure = paste(pd,'LFQ_clin_proteome_naMax10_SD1.5.pdf', sep="_")
plot_clin(LFQ_m2, clus_order2 , clin, figure)

### compare the two marker sets


# samples=colnames(LFQ_m)
# samples[samples %in% Basal]="forestgreen"
# samples[samples %in% LumB]="orange"
# samples[samples %in% Her2]="purple"
# samples[samples %in% CLDN_low]="blue"
# 
# pdf(paste(date,'LFQ_proteome_naMax10_SD2.pdf', sep="_"), height=15, width=10)
# par(oma=c(1,2,1,3))
# LFQ_m_hm = heatmap.2(LFQ_m, trace="none",na.color="white", notecol="black",
#                        cexRow=0.8,cexCol=1.5, ColSideColors = samples, 
#                        labRow=NA,col=getPalette, margins=c(5,5)) #
# par(lend = 1)  
# legend("topright",    # location of the legend on the heatmap plot
#        legend = c("Basal", "Luminal B", "Her2-E", "CLDN_low"), # category labels
#        col = c("forestgreen", "orange", "purple", "blue"),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )
# dev.off()


# ### process and cluster ITRAQ phosphoproteome
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho) = sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho = ITRAQpho[,-c(19,20,22)]
ITRAQpho_m = data.matrix(ITRAQpho[,-c(1,2)]) # 56651 sites
SD=rowSds(ITRAQpho_m, na.rm=TRUE)
ITRAQpho_m2 = ITRAQpho_m[SD>2.15,] # 8483 sites
ITRAQpho_m3 = ITRAQpho_m2[rowSums(is.na(ITRAQpho_m2)) <= 10,] #3279 sites

q99=quantile(ITRAQpho_m3, probs=0.99, na.rm=T)
q1=quantile(ITRAQpho_m3, probs=0.01, na.rm=T)
ITRAQpho_m4 = matrix(,nrow=dim(ITRAQpho_m3)[1],ncol=dim(ITRAQpho_m3)[2])
colnames(ITRAQpho_m4)=colnames(ITRAQpho_m3)
row.names(ITRAQpho_m4)=row.names(ITRAQpho_m3)
for (i in 1:nrow(ITRAQpho_m3)){
  if ( (sum(ITRAQpho_m3[i,][!is.na(ITRAQpho_m3[i,])] > q99) + sum(ITRAQpho_m3[i,][!is.na(ITRAQpho_m3[i,])] < q1)) < 1){
    ITRAQpho_m4[i,]=ITRAQpho_m3[i,]
  }
}

ITRAQpho_m4 = ITRAQpho_m4[rowSums(is.na(ITRAQpho_m4)) <= 10,] #2822 sites
Genes = sub("-NP.*","",row.names(ITRAQpho_m4))
length(unique(Genes)) # 1244 genes

pdf(paste(pd,'ITRAQ_phosphoproteome-ratio-norm_naMax10_SD2.15.pdf', sep="_"))
par(oma=c(1,2,1,3))
ITRAQpho_m4_hm = heatmap.2(ITRAQpho_m4, trace="none",na.color="white", notecol="black",
                           cexRow=0.8,cexCol=1.5, scale="none",dendrogram='column', 
                           labRow=NA,col=getPalette, margins=c(5,5))
clus_order3 = ITRAQpho_m4_hm$colInd
dev.off()

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_clin_phosphoproteome-ratio-norm_naMax10_SD2.15.pdf', sep="_")
plot_clin(ITRAQpho_m4, clus_order3, clin, figure)
# 
# ### process and cluster LFQ phosphoproteome
# row.names(LFQpho) = make.names(LFQpho$phospho_site, unique=T)
# LFQpho = LFQpho[,-c(19:24)]
# colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
# colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
# LFQpho_m = data.matrix(LFQpho)
# SD=rowSds(LFQpho_m, na.rm=TRUE)
# LFQpho_m2 = LFQpho_m[SD>2,]
# LFQpho_m3 = LFQpho_m2[rowSums(is.na(LFQpho_m2)+is.nan(LFQpho_m2)) <= 10,]
# 
# pdf(paste(date,'LFQ_phosphoproteome-ratio-norm_naMax10_SD2.pdf', sep="_"))
# par(oma=c(1,2,1,3))
# LFQpho_m3_hm = heatmap.2(LFQpho_m3, trace="none",na.color="white", notecol="black",cexRow=0.8,cexCol=0.8, col=col_palette2, margins=c(5,5)) #
# clus_order3 = LFQpho_m3_hm$colInd
# dev.off()
# 
# ## plot corresponding clinical panel 
# figure = paste(date,'LFQ_clin_phosphoproteome-ratio-norm_naMax10_SD2.pdf', sep="_")
# plot_clin(LFQpho_m3, clus_order3, cl in, figure)
# 
# ### process and cluster RSEM data
# row.names(RSEM) = RSEM$sample
# RSEM2 = RSEM[,colnames(RSEM) %in% colnames(ITRAQ)] 
# RSEM_m = data.matrix(RSEM2)
# #RSEM_m = log2(RSEM_m)
# SD=rowSds(RSEM_m, na.rm=TRUE)
# RSEM_m2 = RSEM_m[SD>2,]
# RSEM_m3 = RSEM_m2[rowSums(is.na(RSEM_m2)+is.nan(RSEM_m2)) <= 10,]
# 
# pdf(paste(date,'RSEM_naMax3_SD2.pdf', sep="_"))
# par(oma=c(1,2,1,3))
# RSEM_m3_hm = heatmap.2(RSEM_m3, trace="none",na.color="white", notecol="black",cexRow=0.8,cexCol=0.8, col=col_palette2, margins=c(5,5)) #
# clus_order4 = RSEM_m3_hm$colInd
# dev.off()
# 
# ## plot corresponding clinical panel 
# figure = paste(date,'RSEM_clin_naMax3_SD2.pdf', sep="_")
# plot_clin(RSEM_m3, clus_order4, clin, figure)
# 
# 
