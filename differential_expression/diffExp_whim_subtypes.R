# diffExp.R by Kuan Huang @ WashU 201507
# find differentially expressed proteins

# dependencies
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression")
source("/Users/khuang/bin/LIB_exp.R")

# proteome and phosphoproteome files 
ITRAQ = read.table(row.names=1, header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
LFQ=read.table(row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt',header=TRUE, sep="\t", fill=T)
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
#LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
#RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')

Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])

##### ITRAQ proteome #####
ITRAQ_proteome_Basal_vs_LumB = diff_exp(m = ITRAQ, g1 = Basal, g2 = LumB, g1.n = "Basal", g2.n = "LuminalB", plot=T, pathwayC=F, pathwayT=TRUE)

##### LFQ proteome #####
Basal1 = Basal[Basal %in% colnames(LFQ)]
LumB1 = LumB[LumB %in% colnames(LFQ)]
LFQ_proteome_Basal_vs_LumB = diff_exp(m = LFQ, g1 = Basal1, g2 = LumB1, g1.n = "Basal", g2.n = "Luminal_B", plot=T, pathwayC=F, pathwayT=TRUE)

### merge and compare ### 
LFQ_ITRAQ_merge=merge(ITRAQ_proteome_Basal_vs_LumB$stats,LFQ_proteome_Basal_vs_LumB$stats, by = "row.names")

tn = paste(pd, "ITRAQ_LFQ_Basal_vs_LumB_diff_exp.txt")
write.table(LFQ_ITRAQ_merge, file=tn, quote=F, sep = '\t', row.names=F)

colnames(LFQ_ITRAQ_merge)[1] = "gene"
fn = paste(pd,'ITRAQ_vs_LFQ_basal_lumB_diff_logFDR.pdf', sep="_")
p = ggplot(data = LFQ_ITRAQ_merge, aes(x=-log10(as.numeric(as.character(t_test_fdr.x))) , y=-log10(as.numeric(as.character(t_test_fdr.y))), label=gene)) 
p = p + geom_point(alpha=0.5, size=1) + xlab("ITRAQ basal vs. luminal -log(FDR)") + ylab("LFQ basal vs. luminal -log(FDR)") + theme_bw()
p = p + xlim(0,5) + ylim(0,5) + coord_fixed() 
p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(t_test_fdr.x)))>1.8 & -log10(as.numeric(as.character(t_test_fdr.y)))>1.8,gene,"")),hjust=-0.05,vjust=0,size=3) # FDR 0.05
p
ggsave(file=fn)

# LFQ_ITRAQ_mergeP=merge(ITRAQ_proteome_Basal_vs_LumB$Pstats,LFQ_proteome_Basal_vs_LumB$Pstats, by = "hsa")
# colnames(LFQ_ITRAQ_mergeP)[2] = "PathwayName"
# fn = paste(pd,'ITRAQ_vs_LFQ_basal_lumB_diff_pathway_logFDR.pdf', sep="_")
# p = ggplot(data = LFQ_ITRAQ_mergeP, aes(x=-log10(as.numeric(as.character(Wilcox_fdr.x))) , y=-log10(as.numeric(as.character(Wilcox_fdr.y))), label=PathwayName)) 
# p = p + geom_point(alpha=0.5, size=1) + xlab("ITRAQ basal vs. luminal -log(FDR)") + ylab("LFQ basal vs. luminal -log(FDR)") + theme_bw()
# p = p + xlim(0,5) + ylim(0,5) + coord_fixed() 
# # for unknown reason this labeling is not working
# p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(Wilcox_fdr.x)))>1.30103 & -log10(as.numeric(as.character(Wilcox_fdr.y)))>1.30103, PathwayName, "")),hjust=-0.05,just=0,size=3) # FDR 0.05
# p
# ggsave(file=fn)

# # this result is hard to explain
# LFQ_ITRAQ_mergeP_sig = LFQ_ITRAQ_mergeP[-log10(as.numeric(as.character(LFQ_ITRAQ_mergeP$Wilcox_fdr.x)))>1.30103 & -log10(as.numeric(as.character(LFQ_ITRAQ_mergeP$Wilcox_fdr.y)))>1.30103 & LFQ_ITRAQ_mergeP$Direction.x == LFQ_ITRAQ_mergeP$Direction.y,]
# 
LFQ_ITRAQ_mergeR=merge(ITRAQ_proteome_Basal_vs_LumB$Rstats,LFQ_proteome_Basal_vs_LumB$Rstats, by = "ReactomeID")
LFQ_ITRAQ_mergeR_sig = LFQ_ITRAQ_mergeR[-log10(as.numeric(as.character(LFQ_ITRAQ_mergeR$Wilcox_fdr.x)))>2 & -log10(as.numeric(as.character(LFQ_ITRAQ_mergeR$Wilcox_fdr.y)))>2 & LFQ_ITRAQ_mergeR$Direction.x == LFQ_ITRAQ_mergeR$Direction.y,]
tn = paste(pd, "ITRAQ_vs_LFQ_REACTOMEpathway_sig_diff_exp.txt", sep="_")
write.table(LFQ_ITRAQ_mergeR_sig, file=tn, quote=F, sep = '\t', row.names=F)
# 
# p = ggplot(data = LFQ_ITRAQ_mergeGO, aes(x=-log10(as.numeric(as.character(FDR.Down.x))) , y=-log10(as.numeric(as.character(FDR.Down.y))), label=OntologyName)) 
# p = p + geom_point(alpha=0.5, size=1) + xlab("ITRAQ basal vs. luminal -log(FDR)") + ylab("LFQ basal vs. luminal -log(FDR)") + theme_bw()
# p = p + xlim(0,5) + ylim(0,5) + coord_fixed() 
# # for unknown reason this labeling is not working
# p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(FDR.Down.x)))>1.30103 & -log10(as.numeric(as.character(FDR.Down.y)))>1.30103, OntologyName, "")),hjust=-0.05,just=0,size=3) # FDR 0.05
# p


# # look at the count result
# # even harder to explain...
# LFQ_ITRAQ_mergeGO=merge(ITRAQ_proteome_Basal_vs_LumB$GOc,LFQ_proteome_Basal_vs_LumB$GOc, by = "row.names")
# colnames(LFQ_ITRAQ_mergeGO)[2] = "OntologyName"
# LFQ_ITRAQ_mergeGO_down_sig = LFQ_ITRAQ_mergeGO[-log10(as.numeric(as.character(LFQ_ITRAQ_mergeGO$FDR.Down.x)))>1.30103 & -log10(as.numeric(as.character(LFQ_ITRAQ_mergeGO$FDR.Down.y)))>1.30103,]
# LFQ_ITRAQ_mergeGO_up_sig = LFQ_ITRAQ_mergeGO[-log10(as.numeric(as.character(LFQ_ITRAQ_mergeGO$FDR.Up.x)))>1.30103 & -log10(as.numeric(as.character(LFQ_ITRAQ_mergeGO$FDR.Up.y)))>1.30103,]
# 
# #fn = paste(pd,'ITRAQ_vs_LFQ_basal_lumB_diff_pathway_logFDR.pdf', sep="_")
# p = ggplot(data = LFQ_ITRAQ_mergeGO, aes(x=-log10(as.numeric(as.character(FDR.Down.x))) , y=-log10(as.numeric(as.character(FDR.Down.y))), label=OntologyName)) 
# p = p + geom_point(alpha=0.5, size=1) + xlab("ITRAQ basal vs. luminal -log(FDR)") + ylab("LFQ basal vs. luminal -log(FDR)") + theme_bw()
# p = p + xlim(0,5) + ylim(0,5) + coord_fixed() 
# # for unknown reason this labeling is not working
# p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(FDR.Down.x)))>1.30103 & -log10(as.numeric(as.character(FDR.Down.y)))>1.30103, OntologyName, "")),hjust=-0.05,just=0,size=3) # FDR 0.05
# p
# #ggsave(file=fn)
# 
# p = ggplot(data = LFQ_ITRAQ_mergeGO, aes(x=-log10(as.numeric(as.character(FDR.Up.x))) , y=-log10(as.numeric(as.character(FDR.Up.y))), label=OntologyName)) 
# p = p + geom_point(alpha=0.5, size=1) + xlab("ITRAQ basal vs. luminal -log(FDR)") + ylab("LFQ basal vs. luminal -log(FDR)") + theme_bw()
# p = p + xlim(0,5) + ylim(0,5) + coord_fixed() 
# # for unknown reason this labeling is not working
# p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(FDR.Up.x)))>1.30103 & -log10(as.numeric(as.character(FDR.Up.y)))>1.30103, OntologyName, "")),hjust=-0.05,just=0,size=3) # FDR 0.05
# p
# 
# ##### ITRAQ phosphoproteome (gene-level) #####
# ITRAQpho_gene = read.table(row.names=1, header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt')
# ITRAQ_pho_Basal_vs_LumB = diff_exp(m = ITRAQpho_gene, g1 = Basal, g2 = LumB, g1.n = "Basal", g2.n = "Luminal B", plot=T, pathwayC=T, pathwayT=TRUE)

# ### process diff. expression and cluster ITRAQ phosphoproteome
# row.names(ITRAQpho) = ITRAQpho$gene.site
# colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# # 56651 phosphosites
# ITRAQpho=ITRAQpho[,-c(1,2)]
# # get rid of TaxIR, HumIR, WHIM13.1
# ITRAQpho = ITRAQpho[,-c(17,18,20)]
# row.names(ITRAQpho) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho)),unique=T)
# row.names(ITRAQpho) = make.names(sub(" _.*","",row.names(ITRAQpho)), unique=T)
# row.names(ITRAQpho) = make.names(sub("_.*","",row.names(ITRAQpho)), unique=T)
# 
# ITRAQ_phosphoproteome=as.matrix(ITRAQpho)
# ITRAQpho_Basal_vs_LumB = diff_exp(ITRAQ_phosphoproteome, Basal, LumB)$stats
# 
# ### process, diff. expression and cluster LFQ phosphoproteome ###
# LFQpho$sites=sub(".*\\(","",LFQpho$phospho_site)
# LFQpho$sites=sub("\\)","",LFQpho$site)
# LFQpho$sites=paste(LFQpho$gene_name,LFQpho$sites, sep=".")
# row.names(LFQpho) = make.names(LFQpho$sites, unique=T)
# colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
# colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
# # 18229 phosphosites
# LFQpho = LFQpho[,-c(19:25)]
# 
# LFQ_phosphoproteome=as.matrix(LFQpho)
# LFQpho_Basal_vs_LumB = diff_exp(LFQ_phosphoproteome, Basal, LumB)$stats
# 
# ### process, diff. expression and cluster RSEM data ###
# row.names(RNA)=make.names(RNA$gene, unique=T)
# RNA=RNA[,-1] #16209 genes 
# 
# mRNA=as.matrix(RNA)
# mRNA_Basal_vs_LumB = diff_exp(mRNA, Basal, LumB)$stats
# 
# ### process, diff. expression and cluster CNV data ###
# CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified',header=TRUE, sep="\t")
# row.names(CNV)=CNV$gene
# colnames(CNV) = sub("WHIM0","WHIM",colnames(CNV))
# CNV=CNV[,-1] #16209 genes
# 
# CNV=as.matrix(CNV)
# CNV_Basal_vs_LumB = diff_exp(CNV, Basal, LumB)$stats
