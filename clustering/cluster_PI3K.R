# cluster_PI3K.R by Kuan Huang @ WashU 201506
# cluster WHIMs based on expressions in the PI3K pathway

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/clustering")
source("/Users/khuang/bin/LIB_exp.R")

cytokine = KEGG[["hsa04060\tCytokine-cytokine receptor interaction"]]
pi3k = as.vector(t(read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/PI3K_AKT_short_list")))
c_pi3k = as.vector(t(read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cynthia_pi3k.list")))
tri_neg=c("WHIM2","WHIM4","WHIM6","WHIM12","WHIM13","WHIM14","WHIM17","WHIM21","WHIM25","WHIM30","WHIM46")
lumb_er = c("WHIM9","WHIM18","WHIM16","WHIM20","WHIM24","WHIM26","WHIM27","WHIM37","WHIM43")
PI3mut = c("WHIM9","WHIM18","WHIM16","WHIM20","WHIM24")
exp = c("WHIM16","WHIM18","WHIM20")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])

## Proteome
ITRAQ = read.table(header=TRUE, sep="\t", row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
ITRAQ_pi3k = ITRAQ[row.names(ITRAQ) %in% pi3k,]
ITRAQ_pi3k = data.frame(ITRAQ_pi3k)
ITRAQ_pi3k$Dataset = "ITRAQ Protein"
ITRAQ_pi3k$Gene = row.names(ITRAQ_pi3k)

## Phosphoproteome

ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQpho_pi3k = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% pi3k,])

ITRAQpho_pi3k=ITRAQpho_pi3k[rowSums(is.na(ITRAQpho_pi3k))<3,]
ITRAQpho_pi3k=ITRAQpho_pi3k[order(row.names(ITRAQpho_pi3k)),]
ITRAQpho_pi3k = data.frame(ITRAQpho_pi3k)
ITRAQpho_pi3k$Dataset = "ITRAQ Phosphoprotein"
ITRAQpho_pi3k$Gene = row.names(ITRAQpho_pi3k)
#row.names(ITRAQpho_pi3k) = sub("(.*)","\\1_phospho",row.names(ITRAQpho_pi3k))

## RPPA
RPPA = read.table(file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RPPA/20151020_RPPA_curated.txt',header=TRUE, sep="\t")
RPPA = RPPA[,c(2,3,4,5,6)]
RPPA = RPPA[RPPA$Gene %in% pi3k,]
RPPA = RPPA[order(RPPA$Gene),]

ITRAQ_pi3k_exp = ITRAQ_pi3k[,colnames(ITRAQ_pi3k) %in% c(exp,"Gene","Dataset")]
ITRAQpho_pi3k_exp = ITRAQpho_pi3k[,colnames(ITRAQpho_pi3k) %in% c(exp,"Gene","Dataset")]

ITRAQ_pi3k_m = melt(ITRAQ_pi3k_exp, id.var=c("Dataset","Gene"))
ITRAQpho_pi3k_m = melt(ITRAQpho_pi3k_exp, id.var=c("Dataset","Gene"))
RPPA_m = melt(RPPA, id.var=c("Dataset","Gene"))

pro_lvls=rbind(ITRAQ_pi3k_m, ITRAQpho_pi3k_m, RPPA_m)
#pro_lvls = pro_lvls[order(pro_lvls$Sample),]
colnames(pro_lvls) = c("Dataset","Gene","Sample","Expression")

pro_lvls$Sample = factor(pro_lvls$Sample, levels = c("WHIM16", "WHIM18", "WHIM20"))
pro_lvls$Dataset = factor(pro_lvls$Dataset, levels = c("RPPA Protein", "ITRAQ Protein", "RPPA Phosphoprotein","ITRAQ Phosphoprotein"))

fn = paste(pd, "3WHIM_pi3k-pathway.pdf", sep="_")
# summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
# summary[summary$truncated_outlier_percentage>=22,]$truncated_outlier_percentage = 22
min_d = min(as.numeric(as.character(pro_lvls$Expression)), na.rm=T)
max_d = max(as.numeric(as.character(pro_lvls$Expression)), na.rm=T)
bound = max(c(max_d, -min_d))

p = ggplot(data=pro_lvls)
p = p + facet_grid(.~Sample)
p = p + coord_equal()
p = p + geom_tile(aes(x=Dataset, y=Gene, fill=Expression), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA, limits = c(-bound,bound))
p = p  + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14),
        axis.ticks = element_blank(), strip.text = element_text(size = 14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height=7.5, width=10, useDingbats=FALSE)


##### Both #####

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])

ITRAQ = read.table(header=TRUE, sep="\t", row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")

ITRAQ_pi3k = ITRAQ[row.names(ITRAQ) %in% c_pi3k,]
ITRAQ_pi3k = as.matrix(ITRAQ_pi3k)
ITRAQ_pi3k=ITRAQ_pi3k[order(row.names(ITRAQ_pi3k)),]
# pdf(paste(pd,'lumb_ER_ITRAQ_proteome_Cpi3k.pdf', sep="_"))
# par(oma=c(1,2,1,3))

samples=colnames(ITRAQ_pi3k)
samples[samples %in% Basal]="forestgreen"
samples[samples %in% LumB]="orange"
samples[samples %in% Her2]="purple"
samples[samples %in% CLDN_low]="blue"
samples[samples %in% c("WHIM17","WHIM46")] = "forestgreen"

##
ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQpho_pi3k = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% c_pi3k,])

ITRAQpho_pi3k=ITRAQpho_pi3k[rowSums(is.na(ITRAQpho_pi3k))<3,]
ITRAQpho_pi3k=ITRAQpho_pi3k[order(row.names(ITRAQpho_pi3k)),]
ITRAQpho_pi3k = ITRAQpho_pi3k[,colnames(ITRAQ_pi3k)]
row.names(ITRAQpho_pi3k) = sub("(.*)","\\1_PHO",row.names(ITRAQpho_pi3k))
row.names(ITRAQ_pi3k) = sub("(.*)","\\1_PRO",row.names(ITRAQ_pi3k))

ITRAQ_both=rbind(ITRAQ_pi3k, ITRAQpho_pi3k)
ITRAQ_both = ITRAQ_both[order(row.names(ITRAQ_both)),]

pdf(paste(pd,'whim_ITRAQ_pro_phospho_gene_pi3k.pdf', sep="_"), height=12)
par(oma=c(1,2,1,6))
ITRAQ_both_hm = heatmap.2(ITRAQ_both, trace="none",na.color="white", notecol="black", Rowv = "none",
                          cexRow=1.2,cexCol=1.4, ColSideColors = samples, scale="none",dendrogram='column',
                          col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B", "Her2-E", "CLDN_low"), # category labels
       col = c("forestgreen", "orange", "purple", "blue"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()

RNA=read.table(file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA) = make.names(RNA$gene, unique=T)
RNA_pi3k = RNA[row.names(RNA) %in% c_pi3k,colnames(RNA) %in% colnames(ITRAQ_both)]
RNA_pi3k = data.frame(RNA_pi3k)
#RNA_pi3k$Dataset = "RNA"
row.names(RNA_pi3k) = sub("(.*)","\\1_RNA",row.names(RNA_pi3k))
RNA_pi3k$WHIM46=NA
ITRAQ_3 = rbind(ITRAQ_both,RNA_pi3k)
ITRAQ_3 = ITRAQ_3[order(row.names(ITRAQ_3)),]
ITRAQ_3m = data.matrix(ITRAQ_3)

pdf(paste(pd,'whim_rna_pro_phospho_gene_pi3k.pdf', sep="_"), height=18)
par(oma=c(1,2,1,6))
ITRAQ_3_hm = heatmap.2(ITRAQ_3m, trace="none",na.color="white", notecol="black", Rowv = "none",
                          cexRow=1.2,cexCol=1.4, ColSideColors = samples, scale="none",dendrogram='column',
                          col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B", "Her2-E", "CLDN_low"), # category labels
       col = c("forestgreen", "orange", "purple", "blue"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


### only ER+ lumB ###

lumb_er = c("WHIM9","WHIM18","WHIM16","WHIM20","WHIM24","WHIM26","WHIM27","WHIM37","WHIM43")
PI3mut = c("WHIM9","WHIM18","WHIM16","WHIM20","WHIM24")
PI3wt = c("WHIM26","WHIM27","WHIM37","WHIM43")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')

ITRAQ = read.table(header=TRUE, sep="\t", row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')

ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")

ITRAQ_pi3k = ITRAQ[row.names(ITRAQ) %in% c_pi3k,]
ITRAQ_pi3k = as.matrix(ITRAQ_pi3k[,colnames(ITRAQ_pi3k) %in% lumb_er])
ITRAQ_pi3k=ITRAQ_pi3k[order(row.names(ITRAQ_pi3k)),]
# pdf(paste(pd,'lumb_ER_ITRAQ_proteome_Cpi3k.pdf', sep="_"))
# par(oma=c(1,2,1,3))

samples=colnames(ITRAQ_pi3k)
samples[samples %in% PI3mut]="forestgreen"
samples[samples %in% PI3wt]="orange"


# ITRAQ_pi3k_hm = heatmap.2(ITRAQ_pi3k, trace="none",na.color="white", notecol="black",Rowv = "none",
#                           cexRow=1.2,cexCol=1.2, ColSideColors = samples, scale="none",dendrogram='column',
#                           col=getPalette, margins=c(5,5)) #
# par(lend = 1)  
# legend("topright",    # location of the legend on the heatmap plot
#        legend = c("PIK3CA mutant", "PIK3CA wt"), # category labels
#        col = c("forestgreen", "orange"),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )
# 
# dev.off()

# ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
# ITRAQpho = ITRAQpho[,-c(19,20,22)]
# ITRAQpho_pi3k = ITRAQpho[ITRAQpho$gene %in% c_pi3k,]
# 
# row.names(ITRAQpho_pi3k) = ITRAQpho_pi3k$gene.site
# colnames(ITRAQpho_pi3k)<-sub("\\..*", "", colnames(ITRAQpho_pi3k))
# row.names(ITRAQpho_pi3k) = sub("-NP_\\d+_"," ",row.names(ITRAQpho_pi3k))
# row.names(ITRAQpho_pi3k) = make.names(sub(" _.*","",row.names(ITRAQpho_pi3k)), unique=T)
# row.names(ITRAQpho_pi3k) = make.names(sub("_.*","",row.names(ITRAQpho_pi3k)), unique=T)
# ITRAQpho_pi3k = as.matrix(ITRAQpho_pi3k[,-c(1,2)])
# 
# ITRAQpho_pi3k=ITRAQpho_pi3k[,colnames(ITRAQpho_pi3k) %in% lumb_er]
# ITRAQpho_pi3k=ITRAQpho_pi3k[rowSums(is.na(ITRAQpho_pi3k))<2,]
# ITRAQpho_pi3k=ITRAQpho_pi3k[order(row.names(ITRAQpho_pi3k)),]
# 
# pdf(paste(pd,'lumb_ER_ITRAQ_pro_pho_Cpi3k.pdf', sep="_"), height=10, width=5)
# par(oma=c(1,2,1,3))
# ITRAQpho_pi3k_hm = heatmap.2(ITRAQ_both, trace="none",na.color="white", notecol="black", Rowv = "none",
#                              cexRow=0.8,cexCol=1.0, ColSideColors = samples, scale="none",dendrogram='column',
#                              col=getPalette, margins=c(5,5)) #
# par(lend = 1)  
# legend("topright",    # location of the legend on the heatmap plot
#        legend = c("PIK3CA mutant", "PIK3CA wt"), # category labels
#        col = c("forestgreen", "orange"),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )
# dev.off()

##
ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQpho_pi3k = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% c_pi3k,])

ITRAQpho_pi3k=ITRAQpho_pi3k[,colnames(ITRAQpho_pi3k) %in% lumb_er]
ITRAQpho_pi3k=ITRAQpho_pi3k[rowSums(is.na(ITRAQpho_pi3k))<5,]
ITRAQpho_pi3k=ITRAQpho_pi3k[order(row.names(ITRAQpho_pi3k)),]
ITRAQpho_pi3k = ITRAQpho_pi3k[,colnames(ITRAQ_pi3k)]
row.names(ITRAQpho_pi3k) = sub("(.*)","\\1_phospho",row.names(ITRAQpho_pi3k))

ITRAQ_both=rbind(ITRAQ_pi3k, ITRAQpho_pi3k)
ITRAQ_both = ITRAQ_both[order(row.names(ITRAQ_both)),]

pdf(paste(pd,'lumB_ER_ITRAQ_pro_phospho_gene_Cpi3k.pdf', sep="_"), height=10)
par(oma=c(1,2,1,6))
ITRAQ_both_hm = heatmap.2(ITRAQ_both, trace="none",na.color="white", notecol="black", Rowv = "none",
                             cexRow=1.2,cexCol=1.4, ColSideColors = samples, scale="none",dendrogram='column',
                             col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("PIK3CA mutant", "PIK3CA wt"), # category labels
       col = c("forestgreen", "orange"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


### only tri_neg ###

tri_neg=c("WHIM2","WHIM4","WHIM6","WHIM12","WHIM13","WHIM14","WHIM17","WHIM21","WHIM25","WHIM30","WHIM46")
TP53mut = c("WHIM2","WHIM4","WHIM12","WHIM13","WHIM14","WHIM21","WHIM25","WHIM30")
TP53wt = c("WHIM6","WHIM17","WHIM46")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')

ITRAQ = read.table(header=TRUE, sep="\t", row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')

ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")

ITRAQ_tri_neg = ITRAQ[row.names(ITRAQ) %in% c_pi3k,]
ITRAQ_tri_neg = as.matrix(ITRAQ_tri_neg[,colnames(ITRAQ_tri_neg) %in% tri_neg])
ITRAQ_tri_neg=ITRAQ_tri_neg[order(row.names(ITRAQ_tri_neg)),]
# pdf(paste(pd,'lumb_ER_ITRAQ_proteome_Ctri_neg.pdf', sep="_"))
# par(oma=c(1,2,1,3))

samples=colnames(ITRAQ_tri_neg)
samples[samples %in% TP53mut]="forestgreen"
samples[samples %in% TP53wt]="orange"

##
ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQpho_pi3k = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% c_pi3k,])

ITRAQpho_tri_neg=ITRAQpho_pi3k[,colnames(ITRAQpho_pi3k) %in% tri_neg]
ITRAQpho_tri_neg=ITRAQpho_tri_neg[rowSums(is.na(ITRAQpho_tri_neg))<3,]
ITRAQpho_tri_neg=ITRAQpho_tri_neg[order(row.names(ITRAQpho_tri_neg)),]
ITRAQpho_tri_neg = ITRAQpho_tri_neg[,colnames(ITRAQ_tri_neg)]
row.names(ITRAQpho_tri_neg) = sub("(.*)","\\1_phospho",row.names(ITRAQpho_tri_neg))

ITRAQ_both= rbind(ITRAQ_tri_neg, ITRAQpho_tri_neg)
ITRAQ_both = ITRAQ_both[order(row.names(ITRAQ_both)),]

pdf(paste(pd,'tri_neg_ITRAQ_pro_phospho_gene_Cpi3k.pdf', sep="_"), height=10)
par(oma=c(1,2,1,6))
ITRAQ_both_hm = heatmap.2(ITRAQ_both, trace="none",na.color="white", notecol="black", Rowv = "none",
                          cexRow=1.2,cexCol=1.4, ColSideColors = samples, scale="none",dendrogram='column',
                          col=getPalette, margins=c(5,5)) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("TP53 mutant", "TP53 wt"), # category labels
       col = c("forestgreen", "orange"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()
