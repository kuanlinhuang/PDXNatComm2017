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

## function\
# clinical color scale inherited from Matt
get.clinical.scale = function() {
  # Set1 colors
  #colors = c("#f0f0f0", "#636363", "#f0f0f0", "#636363", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  #colors = c(NA, "#101010", NA, "#636363", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
#   # use Perou's intrinsic subtype colors instead
#   colors = c(NA, "#101010", NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
#   
#   color.names = c("wt","mut","negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB")
  colors = c(NA, "#101010", NA, "#636363", "#CE2427","#EF5591","#FFFF33","#8FBCE5","#423996","#50A547") #positive is dark grey       
  
  color.names = c("wt","mut","negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB","LumA","Lymphoma")
  
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="mutation", values=colors)
  
  return(clinical.color.scale)
}
# plot_clin: plot clinical data based on order from the corresponding heatmap2 clustering plot
plot_clin = function(matrix, clus_order, clin, figure){
  color.scale = get.clinical.scale()
  
  order = as.data.frame(colnames(matrix[,clus_order]))
  colnames(order) = colnames(clin)[1]
  # add missing samples to clin through merge
  m = merge(order,clin, by = colnames(clin)[1], all.x = TRUE)
  row.names(m)=m[,colnames(clin)[1]]
  # reorder and melt for hmting
  m = m[colnames(matrix[,clus_order]),]
  m.m = melt(m, id=colnames(clin)[1])
  m.m$ExternalIdentifierName<-with(m.m,factor(ExternalIdentifierName,levels = colnames(matrix[,clus_order])))
  
  #geom_point(aes(x=sample, y=measurement, color=value), size=4, shape=16)
  colourCount=length(unique(m.m$value))
  #p = ggplot(data=m.m, aes(x=ExternalIdentifierName, y=variable)) + geom_tile(aes(fill = value)) 
    #scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(colourCount)) +
  p = ggplot(data=m.m, aes(x=ExternalIdentifierName, y=variable)) + geom_point(aes(color = value), size=6, shape=16)
  p = p + scale_y_discrete()
  p = p + color.scale
  p = p + theme_bw() 
  p = p + xlab("Sample") + ylab("")  
  #p = p + guides(fill=guide_legend(title="Status",nrow=2))
  p = p + guides(fill=FALSE,color=F)
  p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=14), axis.text.y = element_text(colour="black", size=14), legend.position="bottom") 
  p = p + theme(axis.title.x = element_text(colour="black", size=16, face="bold"))
  p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
  p = p + theme(panel.border = element_blank(), panel.background=element_blank())
  p = p + theme(axis.ticks = element_blank())
  p
  ggsave(file=figure, height=3.5, width = 9, useDingbats=FALSE)
} 


clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
#row.names(clin)=clin$ExternalIdentifierName
#clin=clin[,-9]
#clin$STAGE=as.factor(clin$STAGE)
#clin=clin[,c(1:8)]

Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])

##### find common LFQ and iTRAQ markers #####
ITRAQ = read.table(row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt',header=TRUE, sep="\t")
LFQ=read.table(row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt',header=TRUE, sep="\t", fill=T)

ITRAQ_m = as.matrix(ITRAQ)
LFQ_m = as.matrix(LFQ)
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
ITRAQ_m2 = ITRAQ_m2[rowSums(is.na(ITRAQ_m2)) <= 10,] 
LFQ_m2 = LFQ_m[rowSums(is.na(LFQ_m)) <= 10,] 

LFQ_m3 = LFQ_m2[rowSds(LFQ_m2, na.rm=TRUE)>0.75,] 
ITRAQ_m3 = ITRAQ_m2[rowSds(ITRAQ_m2, na.rm=TRUE)>2,] 

ITRAQ_mp = ITRAQ_m3#[row.names(ITRAQ_m3) %in% row.names(LFQ_m3),] #494 markers
LFQ_mp = LFQ_m3[row.names(LFQ_m3) %in% row.names(ITRAQ_m3),] #176 markers

# output for gene annotation
write.table(ITRAQ_mp, col.names=NA, quote=F, sep = '\t', file=paste(pd,"ITRAQ_pro_mp.txt",sep=""))

# read in gene annotation for row labeling
# gene_anno = read.table(row.names=1, sep = '\t', file=paste(pd,"ITRAQ_pro_mp_gene2path_top1.txt",sep=""))
# gene_anno_ord = gene_anno[row.names(ITRAQ_mp),,drop=F]
# gene_anno_index = order(gene_anno_ord$V2)
# gene_anno_ord = gene_anno_ord[gene_anno_index,]
# ITRAQ_mp = ITRAQ_mp[gene_anno_index,]


# num_colors = nlevels(gene_anno_ord)
# gene_anno_ord_colors = colorRampPalette(YlGnBu)(num_colors)
# gene_anno_color = gene_anno_ord_colors[gene_anno_ord]
# 
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# gene_anno_ord_colors=sample(color, num_colors)
# gene_anno_color = gene_anno_ord_colors[gene_anno_ord]

##### process and cluster LFQ proteome
#breaks = c(min(LFQ_mp, na.rm=T), seq(-5,5,length=1023), max(LFQ_mp, na.rm=T))
breaks = c(-5.1, seq(-5,5,length=1023), max(LFQ_mp, na.rm=T))
pdf(paste(pd,'LFQ_proteome_naMax10_SDi2.pdf.pdf', sep="_"))
par(oma=c(1,2,1,3))
LFQ_mp_hm = heatmap.2(LFQ_mp, trace="none",na.color="white", notecol="black", breaks = breaks, density.info = "density",
                      cexRow=0.8,cexCol=1.5, scale="none",dendrogram='column', 
                      labRow=NA,col=getPalette, margins=c(5,5)) #
clus_order2 = LFQ_mp_hm$colInd
dev.off()

LFQ_clus = LFQ_mp[rev(LFQ_mp_hm$rowInd), LFQ_mp_hm$colInd]
write.table(LFQ_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"LFQ_pro_clus.txt",sep=""))

## plot corresponding clinical panel 
figure = paste(pd,'LFQ_clin_proteome_naMax10_SDi2.pdf', sep="_")
plot_clin(LFQ_mp, clus_order2 , clin, figure)

##### process and cluster ITRAQ proteome
row.names(ITRAQ_mp)

breaks = c(min(ITRAQ_mp, na.rm=T), seq(-5,5,length=1023), 6)#max(ITRAQ_mp, na.rm=T))
pdf(paste(pd,'ITRAQ_proteome-ratio-norm_naMax10_SD2.pdf', sep="_"),h=10,w=10)
par(oma=c(1,2,1,3))
ITRAQ_mp_hm = heatmap.2(ITRAQ_mp, trace="none", na.color="white", notecol="black", breaks = breaks, density.info = "density",
                      cexRow=0.2,cexCol=1.5, scale="none",dendrogram='column',
                      Rowv = TRUE, Colv = TRUE,
                      labCol=NA, labRow=NA, #RowSideColors=gene_anno_color,#labRow=NA,
                      col=getPalette, margins=c(5,5))

#par(lend = 1)  
# legend("bottomleft",    # location of the legend on the heatmap plot
#        #legend = unique(gene_anno_ord), # category labels
#        #col = unique(gene_anno_color),  # color key
#        lty= 1,             # line style
#        lwd = 1,         # line width
#        cex=.2
# )

clus_order1 = ITRAQ_mp_hm$colInd
dev.off()

ITRAQ_clus = ITRAQ_mp[rev(ITRAQ_mp_hm$rowInd), ITRAQ_mp_hm$colInd]
write.table(ITRAQ_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"ITRAQ_pro_clus.txt",sep=""))

# add LFQ result
LFQ_proteome = matrix(c("A","A","A","B","B","B","A","A","A","B","N/A","B","B","A","B","A","N/A","B","N/A","N/A","N/A","B","N/A","B"), ncol=1)
clin2 = cbind(clin, LFQ_proteome)  

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_clin_proteome-ratio-norm_naMax10_SD2.pdf', sep="_")
plot_clin(ITRAQ_mp, clus_order1, clin, figure)


##### process and cluster ITRAQ phosphoproteome
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho) = sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho = ITRAQpho[,-c(19,20,22)]
ITRAQpho_m = data.matrix(ITRAQpho[,-c(1,2)]) # 56651 sites
SD=rowSds(ITRAQpho_m, na.rm=TRUE)
ITRAQpho_m2 = ITRAQpho_m[SD>2.5,] # 8483 sites
ITRAQpho_m3 = ITRAQpho_m2[rowSums(is.na(ITRAQpho_m2)) <= 10,] #1733 sites

Genes = sub("-NP.*","",row.names(ITRAQpho_m3))
length(unique(Genes)) # 877 genes

breaks = c(min(ITRAQpho_m3, na.rm=T), seq(-5,5,length=1023), max(ITRAQpho_m3, na.rm=T))
pdf(paste(pd,'ITRAQ_phosphoproteome-ratio-norm_naMax10_SD2.5.pdf', sep="_"))
par(oma=c(1,2,1,3))
ITRAQpho_m4_hm = heatmap.2(ITRAQpho_m3, trace="none",na.color="white", notecol="black",
                           breaks = breaks, density.info = "density",
                           cexRow=0.8,cexCol=1.5, scale="none",dendrogram='column',
                           Rowv = TRUE, Colv = TRUE,
                           labCol=NA,
                           labRow=NA,col=getPalette, margins=c(5,5))
clus_order3 = ITRAQpho_m4_hm$colInd
dev.off()

ITRAQpho_clus = ITRAQpho_m3[rev(ITRAQpho_m4_hm$rowInd), ITRAQpho_m4_hm$colInd]
write.table(ITRAQpho_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"ITRAQ_pho_clus.txt",sep=""))

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_clin_phosphoproteome-ratio-norm_naMax10_SD2.5.pdf', sep="_")
plot_clin(ITRAQpho_m3, clus_order3, clin, figure)

##### process and cluster RSEM
RSEM=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
pam50 = read.table(file="/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/pam50.list", header =F)
pam50 = t(pam50)
i1900 = read.table(file="/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/1900_intrinsic_gene.list", header =F)
i1900 = t(i1900)

row.names(RSEM) = make.names(RSEM[,1], unique=T)
RSEM=RSEM[,-1]
RNA = as.matrix(RSEM[,colnames(RSEM) %in% colnames(ITRAQ)])
# RNA2 = RNA[rowSds(RNA, na.rm=TRUE)>2.5,] # 8483 sites
# RNA3 = RNA2[rowSums(is.na(RNA2)) <= 10,]
RNA2 = RNA[row.names(RNA) %in% pam50,]

breaks = c(min(RNA2, na.rm=T), seq(-5,5,length=1023), max(RNA2, na.rm=T))
pdf(paste(pd,'pam50_RNA.pdf', sep="_"))
par(oma=c(1,2,1,3))
RNA2_hm = heatmap.2(RNA2, trace="none",na.color="white", notecol="black", breaks = breaks, density.info = "density",
                           cexRow=0.7,cexCol=1.5, scale="none",dendrogram='column', offsetCol=0, labCol=NA, labRow=NA,
                           col=getPalette, margins=c(5,5)) 
clus_order4 = RNA2_hm$colInd
dev.off()

RNA_clus = RNA2[rev(RNA2_hm$rowInd), RNA2_hm$colInd]
write.table(RNA_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"RNA_clus.txt",sep=""))

## plot corresponding clinical panel 
figure = paste(pd,'pam50_RNA_clin.pdf', sep="_")
plot_clin(RNA2, clus_order4, clin, figure)

#### RNA using proteome gene #####
RNA2 = RNA[row.names(RNA) %in% row.names(ITRAQ_m3),]

breaks = c(min(RNA2, na.rm=T), seq(-5,5,length=1023), max(RNA2, na.rm=T))
pdf(paste(pd,'RNA_iTRAQgenes.pdf', sep="_"))
par(oma=c(1,2,1,3))
RNA2_hm = heatmap.2(RNA2, trace="none",na.color="white", notecol="black", breaks = breaks, density.info = "density",
                    cexRow=0.7,cexCol=1.5, scale="none",dendrogram='column', offsetCol=0, labCol=NA, labRow=NA,
                    col=getPalette, margins=c(5,5)) 
clus_order4 = RNA2_hm$colInd
dev.off()

RNA_clus = RNA2[rev(RNA2_hm$rowInd), RNA2_hm$colInd]
write.table(RNA_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"RNA_iTRAQgenes_clus.txt",sep=""))

## plot corresponding clinical panel 
figure = paste(pd,'RNA_iTRAQgenes_clin.pdf', sep="_")
plot_clin(RNA2, clus_order4, clin, figure)

##### proteome using PAM50 #####
ITRAQ_pam50 = ITRAQ_m2[row.names(ITRAQ_m2) %in% pam50,]

breaks = c(min(ITRAQ_pam50, na.rm=T), seq(-5,5,length=1023), 6)
pdf(paste(pd,'pam50_iTRAQ_PRO.pdf', sep="_"))
par(oma=c(1,2,1,3))
ITRAQ_pam50_hm = heatmap.2(ITRAQ_pam50, trace="none",na.color="white", notecol="black", breaks = breaks, density.info = "density",
                    cexRow=0.7,cexCol=1.5, scale="none",dendrogram='column', offsetCol=0, labCol=NA, labRow=NA,
                    col=getPalette, margins=c(5,5)) 
clus_order_ITRAQ_pam50 = ITRAQ_pam50_hm$colInd
dev.off()

ITRAQ_pam50_clus = ITRAQ_pam50[rev(ITRAQ_pam50_hm$rowInd), ITRAQ_pam50_hm$colInd]
write.table(ITRAQ_pam50_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"ITRAQ_pam50_clus.txt",sep=""))

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_pam50_clin.pdf', sep="_")
plot_clin(ITRAQ_pam50, clus_order_ITRAQ_pam50, clin, figure)

##### process and cluster mouse proteins #####
mouseITRAQ = read.table(header=TRUE, quote="", row.names=1, file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/mouse-only/proteome-ratio-norm_collapsed.txt')
mouseITRAQ_m = as.matrix(mouseITRAQ) #8438 proteins
mouseITRAQ_m2 = mouseITRAQ_m[rowSums(is.na(mouseITRAQ_m)) <= 10,] #5622
SD=rowSds(mouseITRAQ_m2, na.rm=TRUE)
mouseITRAQ_m3 = mouseITRAQ_m2[SD>2,] # 8483 sites
#dim(mouseITRAQ_m3)

breaks = c(min(mouseITRAQ_m3, na.rm=T), seq(-5,5,length=1023), 6)#max(ITRAQ_mp, na.rm=T))
pdf(paste(pd,'ITRAQ_mouse-only-proteome-ratio-norm_naMax10_SD2.pdf', sep="_"),h=10,w=10)
par(oma=c(1,2,1,3))
mouseITRAQ_mp_hm = heatmap.2(mouseITRAQ_m3, trace="none", na.color="white", notecol="black", breaks = breaks, density.info = "density",
                        cexRow=0.2,cexCol=1.5, scale="none",dendrogram='column',
                        Rowv = TRUE, Colv = TRUE,
                        labCol=NA, labRow=NA, #RowSideColors=gene_anno_color,#labRow=NA,
                        col=getPalette, margins=c(5,5))

clus_order1 = mouseITRAQ_mp_hm$colInd
dev.off()

mouseITRAQ_clus = mouseITRAQ_m3[rev(mouseITRAQ_mp_hm$rowInd), mouseITRAQ_mp_hm$colInd]
write.table(mouseITRAQ_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"mouseITRAQ_pro_clus.txt",sep=""))
 

## plot corresponding clinical panel 
figure = paste(pd,'ITRAQ_clin_mouse-only-proteome-ratio-norm_naMax10_SD2.pdf', sep="_")
plot_clin(mouseITRAQ_m3, clus_order1, clin, figure)




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
