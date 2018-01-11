# diffExp.R by Kuan Huang @ WashU 201507
# find differentially expressed proteins

# dependencies
library(ggplot2)
library(reshape)
library(RColorBrewer)
library("gplots")
library(matrixStats)
library(samr)

# mis
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression")
# system("mkdir figures")
# date=Sys.time()
# date = sub(" .*","",date)
# date = paste(date, "KH", sep="_")
# date = paste("figures/",date,sep="")

source("/Users/khuang/bin/LIB_exp.R")

col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))

RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(rev(RdBu))

# kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
# drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/unionOfGeneVariantAADrug.tsv_hugoified_gene.list', header=FALSE, stringsAsFactors = F)

# functions
## to-do: try using samr and compare to t-test

## diff_exp: implement t-test, wilcoxon rank sum test, to define differential expressed genes in two vectors
diff_exp = function(m, g1, g2){
  g1.n = deparse(substitute(g1))
  g2.n = deparse(substitute(g2))
  m.n = deparse(substitute(m))
  x = m[,colnames(m) %in% g1] 
  y = m[,colnames(m) %in% g2]
  if (nrow(x) != nrow(y)){stop("The number of rows in x is different from that in y!")}
  # print header
  stats = matrix(,nrow=nrow(x),ncol=6)
  row.names(stats)=row.names(x)
  colnames(stats) = c(paste(m.n,g1.n,"mean"), paste(m.n,g2.n,"mean"), "t_test_Tstat", "t_test_p", "w_test_Wstat", "w_test_p")
  for (i in 1:nrow(x)){
    if (sum(!is.na(x[i,]))<2 | sum(!is.na(y[i,]))<2){stats[i,]=rep("NA",6)} else{
      # t-test
      t = t.test(x[i,],y[i,])
      t.p = t$p.value
      t.tstat = t$statistic
      t.conf1 = t$conf.int[1]
      t.conf2 = t$conf.int[2]
      t.meanx = t$estimate[1]
      t.meany = t$estimate[2]
      # Wilcoxon Rank Sum Test
      w = wilcox.test(x[i,],y[i,])
      w.p = w$p.value
      w.Wstat = w$statistic
    
      # return the results
      stats[i,] = c(t.meanx, t.meany, t.tstat, t.p, w.Wstat, w.p)
    }
  }
  t_test_fdr=p.adjust(stats[,"t_test_p"], method="BH")
  stats=cbind(stats, t_test_fdr)
  w_test_fdr=p.adjust(stats[,"w_test_p"], method="BH")
  stats=cbind(stats, w_test_fdr)
  stats=stats[order(stats[,"t_test_fdr"], decreasing=FALSE),]
  
  #sig.genes=row.names(stats[stats[,"t_test_fdr"]<0.01 & !is.na(stats[,"t_test_fdr"]),])
  sig.genes=row.names(stats[1:30,])
  # significant genes in g1 vs. g2 g1+g2 samples
  gs = c(g1,g2)
  m.gs=m[row.names(m) %in% sig.genes,colnames(m) %in% gs]
  
  pdf(paste(date,m.n,g1.n,"vs",g2.n,"top30_fdr.pdf", sep="_"))
  par(oma=c(1,2,5,3))
  mc=colnames(m.gs) %in% g1
  mc[mc]="forestgreen"
  mc[mc=="FALSE"]="orange"
  
  mhm = heatmap.2(m.gs, trace="none",na.color="white", notecol="black",
                      ColSideColors = mc,
                      cexRow=1.1,cexCol=1.4, col=getPalette) #
  par(lend = 1)  
  legend("topright",    # location of the legend on the heatmap plot
         legend = c(g1.n, g2.n), # category labels
         col = c("forestgreen", "orange"),  # color key
         lty= 1,             # line style
         lwd = 10            # line width
  )
  #clus_order1 = mhm$colInd
  dev.off()
  
  return(list("stats"=stats, "sig.genes"=sig.genes))
} 

# proteome and phosphoproteome files 
ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
ITRAQ_human = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
ITRAQpho_human = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
pam50 = read.table(file="/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/pam50.list", header =F)
pam50 = t(pam50)

# define pam50 based on clinical files
clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
clin_human = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/proteome-v2.22041_data/tumor-metadata.csv', header=T, sep=',')

Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])

Basal_h=as.vector(clin_human[clin_human$PAM50=="Basal",1])
LumB_h=as.vector(clin_human[clin_human$PAM50=="LumB",1])
LumA_all=as.vector(clin_human[clin_human$PAM50=="LumA",1])
Her2_h=as.vector(clin_human[clin_human$PAM50=="Her2",1])
Normal=as.vector(clin_human[clin_human$PAM50=="Normal",1])

Basal_all=c(Basal, Basal_h)
LumB_all=c(LumB, LumB_h)
Her2_all=c(Her2, Her2_h)

##### ITRAQ proteome #####
## pre-process ##
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$datETcollapsed
# 12698 genes
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ_proteome = ITRAQ.c[,-c(17,18,20)]
rm(ITRAQ.c)

ITRAQ_proteome = read.table(row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt',header=TRUE, sep="\t")

WHIM=colnames(ITRAQ_proteome)

##### human ITRAQ proteome #####
## pre-process ##
row.names(ITRAQ_human) = ITRAQ_human$Description
ITRAQ_human = ITRAQ_human[ITRAQ_human$Gene != "",]
ITRAQ_human.c = collapseRows(ITRAQ_human[,-c(1,2,3)], rowGroup=ITRAQ_human$Gene, rowID=ITRAQ_human$Description)$datETcollapsed
# 11349 genes
colnames(ITRAQ_human.c) = sub("\\.", "-", colnames(ITRAQ_human.c))
colnames(ITRAQ_human.c) = sub("\\..*", "", colnames(ITRAQ_human.c))

human=colnames(ITRAQ_human.c)

##### merged proteome #####
merged_proteome=merge(ITRAQ_proteome, ITRAQ_human.c, by = "row.names") #10502 genes
rownames(merged_proteome)=merged_proteome$Row.names
merged_proteome=merged_proteome[,-1]
merged_proteome_m=as.matrix(merged_proteome)

merged_proteome_m_human_vs_WHIM = diff_exp(merged_proteome_m, human, WHIM)$stats
tn = paste(pd, "ITRAQ_proteome_human_vs_xeno_diff_exp.txt")
write.table(merged_proteome_m_human_vs_WHIM, file=tn, quote=F, sep = '\t', row.names=T)

### de novo proteome approach ###
whim_genes = row.names(as.matrix(ITRAQ_proteome)[rowSds(as.matrix(ITRAQ_proteome), na.rm=TRUE)>2,])

#merged_proteome_m2 = merged_proteome_m
q99=quantile(merged_proteome_m, probs=0.998, na.rm=T)
q1=quantile(merged_proteome_m, probs=0.002, na.rm=T)
merged_proteome_m2 = matrix(,nrow=dim(merged_proteome_m)[1],ncol=dim(merged_proteome_m)[2])
colnames(merged_proteome_m2)=colnames(merged_proteome_m)
row.names(merged_proteome_m2)=row.names(merged_proteome_m)
for (i in 1:nrow(merged_proteome_m)){
  if ( (sum(merged_proteome_m[i,][!is.na(merged_proteome_m[i,])] > q99) + sum(merged_proteome_m[i,][!is.na(merged_proteome_m[i,])] < q1)) < 1){
    merged_proteome_m2[i,]=merged_proteome_m[i,]
  }
}
merged_proteome_m2 = merged_proteome_m2[row.names(merged_proteome_m2) %in% whim_genes,]
merged_proteome_m2 = merged_proteome_m2[rowSums(is.na(merged_proteome_m2)) <= 10,] # 5068 observations
SD=rowSds(merged_proteome_m2, na.rm=TRUE)
merged_proteome_m2 = merged_proteome_m2[SD>2,] # 478 proteins



# SD=rowSds(merged_proteome, na.rm=TRUE)
# merged_proteome = merged_proteome[SD>2,] #2531 left
# merged_proteome = merged_proteome[rowSums(is.na(merged_proteome)) <= 10,] #988


### conventional color theme for pam50
# Lum A dark blue, Lum B light Blue, HER2 Pink and Basal Red and Claudin low yellow?
samples=colnames(merged_proteome_m2)
samples=sub("\\.1", "", samples)
samples=sub("X263d3f", "263d3f", samples)
samples[samples %in% Basal_all]="red"
samples[samples %in% LumB_all]="#deebf7"
samples[samples %in% LumA_all]="#3182bd"
samples[samples %in% Her2_all]="pink"
samples[samples %in% Normal]="grey"
samples[samples %in% CLDN_low]="#ffeda0"

merged_proteome_m3 = merged_proteome_m2
colnames(merged_proteome_m3)[25:107] = " "
pdf(paste(pd,'merged_ITRAQ_proteome_naMax10_SD2.pdf', sep="_"), height=7, width=15)
par(oma=c(1,2,1,3))
merged_proteome_m2_hm = heatmap.2(merged_proteome_m3, trace="none",na.color="white", notecol="black",
                                  cexRow=0.8,cexCol=0.8, ColSideColors = samples, scale="none",dendrogram='column',
                                  labRow=NA,col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B", "Luminal A", "Her2-E", "Normal", "CLDN_low"), # category labels
       col = c("red", "#deebf7", "#3182bd", "pink", "grey", "#ffeda0"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()

##### Do differential expression between mouse and human ##### 

merged_proteome_m2_human_vs_WHIM = diff_exp(merged_proteome_m2, human, WHIM)$stats
tn = paste(pd, "ITRAQ_proteome_human_vs_xeno_diff_exp.txt")
write.table(merged_proteome_m2_human_vs_WHIM, file=tn, quote=F, sep = '\t', row.names=T)

merged_proteome_m2_human_vs_WHIM_nonsig = merged_proteome_m2_human_vs_WHIM[merged_proteome_m2_human_vs_WHIM[,"t_test_fdr"]>0.3,]
non_diff_genes = row.names(merged_proteome_m2_human_vs_WHIM_nonsig) # 211 proteins
merged_proteome_m2_nondiff = merged_proteome_m2[row.names(merged_proteome_m2) %in% non_diff_genes,]

merged_proteome_m3_nondiff = merged_proteome_m3[row.names(merged_proteome_m3) %in% non_diff_genes,]

# attach proteomic subtype from the CPTAC BRCA paper (Mertins et al., Nature 2016)
cptac_bc_cat = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_105BRCA/final_submitted/nature18003-s2/CPTAC_BC_SupplementaryTable01.txt"
cptac_human = read.table(file=cptac_bc_cat, header=T, sep='\t')

CPTAC_PRO_1 = as.vector(cptac_human[cptac_human$Proteome.Cluster..see.Fig..3b. %in% 1,1])
CPTAC_PRO_2 = as.vector(cptac_human[cptac_human$Proteome.Cluster..see.Fig..3b. %in% 2,1])
CPTAC_PRO_3 = as.vector(cptac_human[cptac_human$Proteome.Cluster..see.Fig..3b. %in% 3,1])

CPTAC_clus=colnames(merged_proteome_m2)
CPTAC_clus=sub("\\.1", "", CPTAC_clus)
CPTAC_clus=sub("X263d3f", "263d3f", CPTAC_clus)
CPTAC_clus[CPTAC_clus %in% CPTAC_PRO_1]="#b2182b"
CPTAC_clus[CPTAC_clus %in% CPTAC_PRO_2]="#2166ac"
CPTAC_clus[CPTAC_clus %in% CPTAC_PRO_3]="#1b7837"
CPTAC_clus[CPTAC_clus %in% c("263d3f-I","blcdb9-I","c4155b-C")]="#1b7837"
CPTAC_clus[!(CPTAC_clus %in% c("#b2182b","#2166ac","#1b7837"))]=NA


# ways to have two col side colors according to BioStar (https://www.biostars.org/p/18211/)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
clab=cbind(samples,CPTAC_clus)
colnames(clab)=c("Subtype","CPTAC proteomic cluster")

pdf(paste(pd,'merged_ITRAQ_proteome_nondiff_naMax10_SD2.pdf', sep="_"), height=7, width=15)
par(oma=c(1,2,1,3))
merged_proteome_m3_nondiff_hm = heatmap.2(merged_proteome_m3_nondiff, trace="none",na.color="white", notecol="black",
                                          cexRow=0.8,cexCol=0.8, ColSideColors = samples,labRow=NA,dendrogram='column', scale="none",
                                          col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B", "Luminal A", "Her2-E", "Normal", "CLDN_low"), # category labels
       col = c("red", "#deebf7", "#3182bd", "pink", "grey", "#ffeda0"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

BRCA77_WHIM24_clus = merged_proteome_m3_nondiff[rev(merged_proteome_m3_nondiff_hm$rowInd), merged_proteome_m3_nondiff_hm$colInd]
write.table(BRCA77_WHIM24_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"BRCA77_WHIM24_pro_clus.txt",sep=""))

dev.off()


### PAM50 based approach ###
merged_proteome_m2 = merged_proteome_m[row.names(merged_proteome_m) %in% pam50,]
# q99=quantile(merged_proteome_m, probs=0.998, na.rm=T)
# q1=quantile(merged_proteome_m, probs=0.002, na.rm=T)
# merged_proteome_m2 = matrix(,nrow=dim(merged_proteome_m)[1],ncol=dim(merged_proteome_m)[2])
# colnames(merged_proteome_m2)=colnames(merged_proteome_m)
# row.names(merged_proteome_m2)=row.names(merged_proteome_m)
# for (i in 1:nrow(merged_proteome_m)){
#   if ( (sum(merged_proteome_m[i,][!is.na(merged_proteome_m[i,])] > q99) + sum(merged_proteome_m[i,][!is.na(merged_proteome_m[i,])] < q1)) < 1){
#     merged_proteome_m2[i,]=merged_proteome_m[i,]
#   }
# }
merged_proteome_m2 = merged_proteome_m2[rowSums(is.na(merged_proteome_m2)) <= 10,] # 5068 observations
# SD=rowSds(merged_proteome_m2, na.rm=TRUE)
# merged_proteome_m2 = merged_proteome_m2[SD>2,] # 478 proteins



# SD=rowSds(merged_proteome, na.rm=TRUE)
# merged_proteome = merged_proteome[SD>2,] #2531 left
# merged_proteome = merged_proteome[rowSums(is.na(merged_proteome)) <= 10,] #988


### conventional color theme for pam50
# Lum A dark blue, Lum B light Blue, HER2 Pink and Basal Red and Claudin low yellow?
samples=colnames(merged_proteome_m2)
samples=sub("\\.1", "", samples)
samples=sub("X263d3f", "263d3f", samples)
samples[samples %in% Basal_all]="red"
samples[samples %in% LumB_all]="#deebf7"
samples[samples %in% LumA_all]="#3182bd"
samples[samples %in% Her2_all]="pink"
samples[samples %in% Normal]="grey"
samples[samples %in% CLDN_low]="#ffeda0"

merged_proteome_m3 = merged_proteome_m2
colnames(merged_proteome_m3)[25:107] = " "
pdf(paste(pd,'merged_ITRAQ_proteome_naMax10_PAM50.pdf', sep="_"), height=7, width=15)
par(oma=c(1,2,1,3))
merged_proteome_m2_hm = heatmap.2(merged_proteome_m3, trace="none",na.color="white", notecol="black",
                               cexRow=0.8,cexCol=0.8, ColSideColors = samples, scale="none",dendrogram='column',
                               labRow=NA,col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B", "Luminal A", "Her2-E", "Normal", "CLDN_low"), # category labels
       col = c("red", "#deebf7", "#3182bd", "pink", "grey", "#ffeda0"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()

##### Do differential expression between mouse and human ##### 

merged_proteome_m2_human_vs_WHIM = diff_exp(merged_proteome_m2, human, WHIM)$stats
tn = paste(pd, "ITRAQ_proteome_human_vs_xeno_diff_exp.txt")
write.table(merged_proteome_m2_human_vs_WHIM, file=tn, quote=F, sep = '\t', row.names=T)

merged_proteome_m2_human_vs_WHIM_nonsig = merged_proteome_m2_human_vs_WHIM[merged_proteome_m2_human_vs_WHIM[,"t_test_fdr"]>0.1,]
non_diff_genes = row.names(merged_proteome_m2_human_vs_WHIM_nonsig) # 211 proteins
merged_proteome_m2_nondiff = merged_proteome_m2[row.names(merged_proteome_m2) %in% non_diff_genes,]

merged_proteome_m3_nondiff = merged_proteome_m3[row.names(merged_proteome_m3) %in% non_diff_genes,]

# attach proteomic subtype from the CPTAC BRCA paper (Mertins et al., Nature 2016)
cptac_bc_cat = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_105BRCA/final_submitted/nature18003-s2/CPTAC_BC_SupplementaryTable01.txt"
cptac_human = read.table(file=cptac_bc_cat, header=T, sep='\t')

CPTAC_PRO_1 = as.vector(cptac_human[cptac_human$Proteome.Cluster..see.Fig..3b. %in% 1,1])
CPTAC_PRO_2 = as.vector(cptac_human[cptac_human$Proteome.Cluster..see.Fig..3b. %in% 2,1])
CPTAC_PRO_3 = as.vector(cptac_human[cptac_human$Proteome.Cluster..see.Fig..3b. %in% 3,1])

CPTAC_clus=colnames(merged_proteome_m2)
CPTAC_clus=sub("\\.1", "", CPTAC_clus)
CPTAC_clus=sub("X263d3f", "263d3f", CPTAC_clus)
CPTAC_clus[CPTAC_clus %in% CPTAC_PRO_1]="#b2182b"
CPTAC_clus[CPTAC_clus %in% CPTAC_PRO_2]="#2166ac"
CPTAC_clus[CPTAC_clus %in% CPTAC_PRO_3]="#1b7837"
CPTAC_clus[CPTAC_clus %in% c("263d3f-I","blcdb9-I","c4155b-C")]="#1b7837"
CPTAC_clus[!(CPTAC_clus %in% c("#b2182b","#2166ac","#1b7837"))]=NA


# ways to have two col side colors according to BioStar (https://www.biostars.org/p/18211/)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
clab=cbind(samples,CPTAC_clus)
colnames(clab)=c("Subtype","CPTAC proteomic cluster")

pdf(paste(pd,'merged_ITRAQ_proteome_nondiff_naMax10_SD2_PAM50.pdf', sep="_"), height=7, width=15)
par(oma=c(1,2,1,3))
merged_proteome_m3_nondiff_hm = heatmap.2(merged_proteome_m3_nondiff, trace="none",na.color="white", notecol="black",
                               cexRow=0.8,cexCol=0.8, ColSideColors = samples,labRow=NA,dendrogram='column', scale="none",
                               col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B", "Luminal A", "Her2-E", "Normal", "CLDN_low"), # category labels
       col = c("red", "#deebf7", "#3182bd", "pink", "grey", "#ffeda0"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

BRCA77_WHIM24_clus = merged_proteome_m3_nondiff[rev(merged_proteome_m3_nondiff_hm$rowInd), merged_proteome_m3_nondiff_hm$colInd]
write.table(BRCA77_WHIM24_clus, col.names=NA, quote=F, sep = '\t', file=paste(pd,"BRCA77_WHIM24_pro_clus_PAM50.txt",sep=""))

dev.off()

##### find signatures of WHIM17, WHIM46, and AR.A1AW #####
### look for the most up-regulated proteins in these samples ###
if (FALSE){
resist_proteome = merged_proteome_m2[non_diff_genes,c("WHIM17","WHIM46","AR-A1AW")]
resist_genes = row.names(resist_proteome[rowSums(resist_proteome>1.5, na.rm=T) == 3,])

# CNTRL
# GLRX
# IFIT3
# OAS2
# GBP4
# UBE2L6
# USP18
# MZB1
# HLA-E
# IKZF2
# IFI44
# HLA-DRB4
# GBP2
# GBP7
# TAP2
# TAPBP
# BST2
# HLA-F

##### process diff. expression and cluster ITRAQ phosphoproteome
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho = ITRAQpho[,-c(17,18,20)]
row.names(ITRAQpho) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho)),unique=T)
row.names(ITRAQpho) = make.names(sub(" _.*","",row.names(ITRAQpho)), unique=T)
row.names(ITRAQpho) = make.names(sub("_.*","",row.names(ITRAQpho)), unique=T)


### process, diff. expression and cluster LFQ phosphoproteome ###
LFQpho$sites=sub(".*\\(","",LFQpho$phospho_site)
LFQpho$sites=sub("\\)","",LFQpho$site)
LFQpho$sites=paste(LFQpho$gene_name,LFQpho$sites, sep=".")
row.names(LFQpho) = make.names(LFQpho$sites, unique=T)
colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
# 18229 phosphosites
LFQpho = LFQpho[,-c(19:25)]

LFQ_phosphoproteome=as.matrix(LFQpho)


### process, diff. expression and cluster RSEM data ###
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes 

mRNA=as.matrix(RNA)


### process, diff. expression and cluster CNV data ###
CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified',header=TRUE, sep="\t")
row.names(CNV)=CNV$gene
colnames(CNV) = sub("WHIM0","WHIM",colnames(CNV))
CNV=CNV[,-1] #16209 genes

CNV=as.matrix(CNV)
}
