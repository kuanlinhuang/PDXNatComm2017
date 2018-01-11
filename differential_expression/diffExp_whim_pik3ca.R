# diffExp.R by Kuan Huang @ WashU 201507
# find differentially expressed proteins

# mis
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression")
source("/Users/khuang/bin/LIB_exp.R")

# proteome and phosphoproteome files 
ITRAQ_proteome = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt',header=TRUE, sep="\t")
LFQ_proteome =read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified',header=TRUE, sep="\t", fill=T)
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')

normal=as.vector(clin[clin$PIK3CA.mutation=="wt",1])
pik3ca=as.vector(clin[clin$PIK3CA.mutation=="mut",1])

##### ITRAQ phosphoproteome #####
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho = ITRAQpho[,-c(17,18,20)]
ITRAQpho.na = ITRAQpho[rowSums(is.na(ITRAQpho)) <= 14,] #35838 phosphosites
rm(ITRAQpho)

row.names(ITRAQpho.na) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub(" _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("_.*","",row.names(ITRAQpho.na)), unique=T)
ITRAQ_phosphoproteome = as.matrix(ITRAQpho.na)
ITRAQ_proteome_pik3ca = diff_exp(ITRAQ_phosphoproteome, pik3ca, normal)

##### ITRAQ proteome #####
## pre-process ##
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$datETcollapsed
# 12698 genes
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ_proteome = ITRAQ.c[,-c(17,18,20)]
rm(ITRAQ.c)

## differential expression and cluster ITRAQ proteome ##
ITRAQ_proteome_Basal_vs_LumB = diff_exp(ITRAQ_proteome, Basal, LumB)$stats

##### LFQ proteome #####
row.names(LFQ) = LFQ$RefSeq_id
LFQ.c = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$datETcollapsed #5802 genes
#LFQ.NP = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$group2row
rm(LFQ)
# get rid of WHIM16_2, WHIM2_2
LFQ_proteome = LFQ.c[,-c(6,9)]
rm(LFQ.c)

LFQ_proteome_Basal_vs_LumB = diff_exp(LFQ_proteome, Basal, LumB)$stats

### merge and compare ### 
LFQ_ITRAQ_merge=merge(basal_vs_lum,LFQ.basal_vs_lum, by = "row.names")
colnames(LFQ_ITRAQ_merge)[1] = "gene"
fn = paste(date,'ITRAQ_vs_LFQ_basal_lumB_diff_logP.pdf', sep="_")
p = ggplot(data = LFQ_ITRAQ_merge, aes(x=-log10(as.numeric(as.character(t_test_p.x))) , y=-log10(as.numeric(as.character(t_test_p.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("ITRAQ basal vs. luminal -logP") + ylab("LFQ basal vs. luminal -logP") + theme_bw()
p = p + xlim(0,9) + ylim(0,9)
p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(t_test_p.x)))>5 | -log10(as.numeric(as.character(t_test_p.y)))>5,gene,"")),hjust=-0.05,just=0,size=8)
p
ggsave(file=fn, width=12, height=8)

### process diff. expression and cluster ITRAQ phosphoproteome
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho = ITRAQpho[,-c(17,18,20)]
row.names(ITRAQpho) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho)),unique=T)
row.names(ITRAQpho) = make.names(sub(" _.*","",row.names(ITRAQpho)), unique=T)
row.names(ITRAQpho) = make.names(sub("_.*","",row.names(ITRAQpho)), unique=T)

ITRAQ_phosphoproteome=as.matrix(ITRAQpho)
ITRAQpho_Basal_vs_LumB = diff_exp(ITRAQ_phosphoproteome, Basal, LumB)$stats

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
LFQpho_Basal_vs_LumB = diff_exp(LFQ_phosphoproteome, Basal, LumB)$stats

### process, diff. expression and cluster RSEM data ###
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes 

mRNA=as.matrix(RNA)
mRNA_Basal_vs_LumB = diff_exp(mRNA, Basal, LumB)$stats

### process, diff. expression and cluster CNV data ###
CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified',header=TRUE, sep="\t")
row.names(CNV)=CNV$gene
colnames(CNV) = sub("WHIM0","WHIM",colnames(CNV))
CNV=CNV[,-1] #16209 genes

CNV=as.matrix(CNV)
CNV_Basal_vs_LumB = diff_exp(CNV, Basal, LumB)$stats
