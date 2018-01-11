# heatmap_immuno.R by Kuan Huang @ WashU 201506
# cluster WHIMs based on expressions in the PI3K pathway

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/clustering")
source("/Users/khuang/bin/LIB_exp.R")

cytokine = KEGG[["hsa04060\tCytokine-cytokine receptor interaction"]]
pi3k = as.vector(t(read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/PI3K_AKT_short_list")))
#c_pi3k = as.vector(t(read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/cynthia_pi3k.list")))
tCell_f = read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_lists/Hugo_tCell_receptors.txt",quote="",header=T,sep="\t")
tCell = as.vector(t(tCell_f[,2]))
globulin_f = read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_lists/Hugo_Immunoglobulins.txt",quote="",header=T,sep="\t")
globulin = as.vector(t(globulin_f[,2]))
ck_f = read.table("/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_lists/immuno_genelist.txt_hugoified",quote="",header=T,sep="\t")
ck = as.vector(t(ck_f))

f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/reference_files/gene_lists/20160209_IHC_marker_KH.list"
r_f = read.table(f,quote="",header=F,sep="\t")
glist = as.vector(t(r_f))

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
RNA=read.table(file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA) = make.names(RNA$gene, unique=T)
ITRAQ_glist = ITRAQ[row.names(ITRAQ) %in% glist,]
mITRAQ_glist = melt(as.matrix(ITRAQ_glist))

fn = paste(pd, "WHIM_marker_ProteinExp_heatmap.pdf", sep="_")
# summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
# summary[summary$truncated_outlier_percentage>=22,]$truncated_outlier_percentage = 22
min_d = min(mITRAQ_glist$value, na.rm=T)
max_d = max(mITRAQ_glist$value, na.rm=T)
bound = max(c(max_d, -min_d))

p = ggplot(data=mITRAQ_glist)
#p = p + facet_grid(.~Sample)
p = p + coord_equal()
p = p + geom_tile(aes(x=Var2, y=Var1, fill=value), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA, limit=c(-bound,bound))
p = p  + theme_bw() + labs(y = "gene", x="sample") + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14),
        axis.ticks = element_blank(), strip.text = element_text(size = 14))#element_text(colour="black", ize=14))
p
ggsave(file=fn, height=3, width=10, useDingbats=FALSE)


RNA_glist = RNA[row.names(RNA) %in% glist,colnames(RNA) %in% colnames(ITRAQ)]
mRNA_glist = melt(as.matrix(RNA_glist))

fn = paste(pd, "WHIM_marker_RNA_Exp_heatmap.pdf", sep="_")
# summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
# summary[summary$truncated_outlier_percentage>=22,]$truncated_outlier_percentage = 22
min_d = min(mRNA_glist$value, na.rm=T)
max_d = max(mRNA_glist$value, na.rm=T)
bound = max(c(max_d, -min_d))

p = ggplot(data=mRNA_glist)
#p = p + facet_grid(.~Sample)
p = p + coord_equal()
p = p + geom_tile(aes(x=Var2, y=Var1, fill=value), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA, limit=c(-bound,bound))
p = p  + theme_bw() + labs(y = "gene", x="sample") + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14),
        axis.ticks = element_blank(), strip.text = element_text(size = 14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height=4, width=10, useDingbats=FALSE)


mITRAQ_glist$dataset = "Protein"
mRNA_glist$dataset = "mRNA"
combined = rbind(mITRAQ_glist, mRNA_glist)

fn = paste(pd, "WHIM_marker_RNA_Pro_Exp_heatmap.pdf", sep="_")
# summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
# summary[summary$truncated_outlier_percentage>=22,]$truncated_outlier_percentage = 22
min_d = min(combined$value, na.rm=T)
max_d = max(combined$value, na.rm=T)
bound = max(c(max_d, -min_d))

p = ggplot(data=combined)
p = p + facet_grid(dataset~.)
p = p + coord_equal()
p = p + geom_tile(aes(x=Var2, y=Var1, fill=value), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA, limit=c(-bound,bound))
p = p  + theme_bw() + labs(y = "Gene", x="Sample") + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90, hjust = 1, vjust =0.5), axis.text.y = element_text(colour="black", size=14),
        axis.ticks = element_blank(), strip.text = element_text(size = 14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height=8, width=10, useDingbats=FALSE)

###


ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQpho_glist = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% glist,])

ITRAQpho_glist$dataset = "Phosphoproteome"
combined = rbind(mITRAQ_glist, mRNA_glist)

fn = paste(pd, "WHIM_marker_RNA_Pro_Exp_heatmap.pdf", sep="_")
# summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
# summary[summary$truncated_outlier_percentage>=22,]$truncated_outlier_percentage = 22
min_d = min(combined$value, na.rm=T)
max_d = max(combined$value, na.rm=T)
bound = max(c(max_d, -min_d))

p = ggplot(data=combined)
p = p + facet_grid(dataset~.)
p = p + coord_equal()
p = p + geom_tile(aes(x=Var2, y=Var1, fill=value), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA, limit=c(-bound,bound))
p = p  + theme_bw() + labs(y = "Gene", x="Sample") + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90, hjust = 1, vjust =0.5), axis.text.y = element_text(colour="black", size=14),
        axis.ticks = element_blank(), strip.text = element_text(size = 14))#element_text(colour="black", size=14))
p
ggsave(file=fn, height=8, width=10, useDingbats=FALSE)




#####
RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes

RNA_tCell = RNA[row.names(RNA) %in% tCell,]
RNA_globulin = RNA[row.names(RNA) %in% globulin,]
RNA_CD20 = RNA[row.names(RNA) %in% "MS4A1",]


## Phosphoproteome

ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQpho_pi3k = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% pi3k,])

ITRAQpho_pi3k=ITRAQpho_pi3k[rowSums(is.na(ITRAQpho_pi3k))<5,]
ITRAQpho_pi3k=ITRAQpho_pi3k[order(row.names(ITRAQpho_pi3k)),]
ITRAQpho_pi3k = data.frame(ITRAQpho_pi3k)
ITRAQpho_pi3k$Dataset = "ITRAQ Phosphoprotein"
ITRAQpho_pi3k$Gene = row.names(ITRAQpho_pi3k)
#row.names(ITRAQpho_pi3k) = sub("(.*)","\\1_phospho",row.names(ITRAQpho_pi3k))

## RPPA
RPPA = read.table(file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RPPA/20151020_RPPA_curated.txt',header=TRUE, sep="\t")
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
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="blaglist", size=14, angle=90), axis.text.y = element_text(colour="blaglist", size=14),
        axis.tiglists = element_blank(), strip.text = element_text(size = 14))#element_text(colour="blaglist", size=14))
p
ggsave(file=fn, height=7.5, width=10, useDingbats=FALSE)