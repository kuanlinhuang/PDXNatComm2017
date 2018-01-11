## prepare_exp_file.R
# 201508 by KH at TGI

### dependencies ###
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/")
source("/Users/khuang/bin/LIB_exp.R")
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/drug_list.txt', header=FALSE, stringsAsFactors = F)
druggable = as.vector(t(drugList))

plot_barplot = function(matrix, x_string, fill_string=NULL, fileName="data.pdf"){
  
  fn = paste(pd, fileName,sep ="_")
  
  if (is.null(fill_string)){
    p = ggplot(matrix,aes_string(x = x_string))
  } else{
    p = ggplot(matrix,aes_string(x = x_string, fill = fill_string))
  }
  
  p = p + geom_bar() + theme_bw() + theme_nogrid()
  p = p + labs(x = x_string, y="counts")
  p = p + geom_vline(xintercept = 9.5, alpha = 0.5)
  p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
  p
  ggsave(file=fn, useDingbats=FALSE)
}

##### add number of NAs to each protein #####
ITRAQ_pro = read.table(header=T, sep="\t", quote="", file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/whim-proteome-v3.01/20160818_WHIM_BC_SupplementaryTable_Proteome.SGS.1.v3.01011.txt')
ITRAQ_pho = read.table(header=T, sep="\t", quote="", file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/whim-phosphoproteome-v3.01011/20160729_WHIM_BC_SupplementaryTable_Phosphosite.txt')
LFQ_pro = read.table(header=T, sep="\t", quote="", file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/label_free_all_WHIM_LFQ_Global/20161201_PDX_supp_table_LFQ_proteome.txt')
ITRAQ_pro$numSamplesObserved = rowSums(!is.na(ITRAQ_pro[,c(14:40)]))#[,c(14:29,32,34:40)]))
ITRAQ_pho$numSamplesObserved = rowSums(!is.na(ITRAQ_pho[,c(32:58)]))#[,c(32:47,50,52:58)]))
LFQ_pro$numSamplesObserved = rowSums(!is.na(LFQ_pro[,c(7:26)]))#[,c(7:11,13,14,16:26)]))
#write.table(ITRAQ_pro, quote=F, col.names=NA,sep = '\t', file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/whim-proteome-v3.01/20160818_WHIM_BC_SupplementaryTable_Proteome.SGS.1.v3.01011_NAcount.txt')
#write.table(ITRAQ_pho, quote=F, col.names=NA,sep = '\t', file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/whim-phosphoproteome-v3.01011/20160729_WHIM_BC_SupplementaryTable_Phosphosite_NAcount.txt')
#write.table(LFQ_pro, quote=F, col.names=NA,sep = '\t', file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/label_free_all_WHIM_LFQ_Global/20161201_PDX_supp_table_LFQ_proteome_NAcount.txt')

### some stats on druggable ###
ITRAQ_pro_druggable = ITRAQ_pro[ITRAQ_pro$geneSymbol %in% druggable,]
ITRAQ_pro_druggable_24 = ITRAQ_pro_druggable[,c(14:40)]
table(rowSums(!is.na(ITRAQ_pro_druggable_24)))
colnames(LFQ_pro)[1] = "geneSymbol"
ITRAQ_pro$data = "ITRAQ_proteome"
ITRAQ_pho$data = "ITRAQ_phophosite"
LFQ_pro$data = "LFQ_proteome"
numSampleMatrix = rbind(ITRAQ_pro[,c("geneSymbol","numSamplesObserved","data")],ITRAQ_pho[,c("geneSymbol","numSamplesObserved","data")],LFQ_pro[,c("geneSymbol","numSamplesObserved","data")])
numSampleMatrix$druggable = "Other genes"
numSampleMatrix$druggable[numSampleMatrix$geneSymbol %in% druggable] = "Druggable genes"

p = ggplot(numSampleMatrix,aes(x = numSamplesObserved))
p = p + facet_grid(druggable~data, drop=T, scales = "free")
p = p + geom_bar(binwidth=1) + theme_bw() + theme_nogrid()
p = p + labs(y="Marker counts", x = "Number of samples observed")
p = p + geom_vline(xintercept = 9.5, alpha = 0.5)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12))#element_text(colour="black", size=14))
p
fn=paste(pd,"SampleNoNAcount_bar.pdf",sep="_")
ggsave(file=fn, useDingbats=FALSE)

# LFQ_pro_druggable = LFQ_pro[LFQ_pro$gene_name %in% druggable,]
# ITRAQ_pho_druggable = ITRAQ_pho[ITRAQ_pho$geneSymbol %in% druggable,]
# 
# plot_barplot(ITRAQ_pro,x_string = "numSamplesObserved",fileName="iTRAQpro_SampleNoNAcount_bar.pdf")
# plot_barplot(ITRAQ_pho,x_string = "numSamplesObserved",fileName="iTRAQpho_SampleNoNAcount_bar.pdf")
# plot_barplot(LFQ_pro,x_string = "numSamplesObserved",fileName="LFQ_SampleNoNAcount_bar.pdf")
# plot_barplot(ITRAQ_pro,x_string = "numSamplesObserved",fileName="iTRAQpro_SampleNoNAcount_bar.pdf")
# plot_barplot(ITRAQ_pho,x_string = "numSamplesObserved",fileName="iTRAQpho_SampleNoNAcount_bar.pdf")
# plot_barplot(LFQ_pro_druggable,x_string = "numSamplesObserved",fileName="LFQ_druggable_SampleNoNAcount_bar.pdf")

sum(ITRAQ_pro$numSamplesProteinObserved>=10)
sum(ITRAQ_pho$numSamplesPhosphositeObserved>=10)
sum(LFQ_pro$numSamplesProteinObserved>=10)

ITRAQ_pro$category = NA
ITRAQ_pro$category[ITRAQ_pro$geneSymbol %in% druggable] = "Druggable genes"

ITRAQ_pro$rep_diff = ITRAQ_pro$WHIM13.PDXBC06 - ITRAQ_pro$WHIM13.PDXBC07
ITRAQ_pro$rep_avg = (ITRAQ_pro$WHIM13.PDXBC06 + ITRAQ_pro$WHIM13.PDXBC07)/2
p = ggplot(data=ITRAQ_pro)
p = p + geom_point(aes(x=rep_avg,y=rep_diff, colour = category),alpha=0.3,pch=16) + theme_bw() + theme_nogrid()
p = p + labs(y = "Difference between WHIM13.PDXBC06 and WHIM13.PDXBC07", x="Average Protein Expression")
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust=0.5), 
              axis.text.y = element_text(colour="black", size=12), legend.position = "bottom")#element_text(colour="black", size=14))
p
fn = paste(pd, "technical_replicate_diff_and_avg.pdf",sep ="_")
ggsave(file=fn, useDingbats=FALSE)

##### collapse protein isoforms #####
### ITRAQ proteome-mouse only ###
# mouseITRAQ = read.table(header=F, sep="\t", fill=T,quote="",file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/mouse-only/proteome-ratio-norm.gct')
mouseITRAQ = read.table(header=F, sep="\t", fill=T,quote="",file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/whim-proteome-v3.01/whim-proteome-mouseonly-v3.01011-ratio-norm.gct')
colnames(mouseITRAQ) = as.character(unlist(mouseITRAQ[2,,drop=T]))
mouseITRAQ = mouseITRAQ[-c(1:2),]
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 13719 isoforms
mouseITRAQ[,-c(1:2)] = as.data.frame(lapply(mouseITRAQ[,-c(1:2)],function (x) as.numeric(as.character(x))))
row.names(mouseITRAQ) = mouseITRAQ$Description
colnames(mouseITRAQ) = sub("\\..*", "", colnames(mouseITRAQ))
mouseITRAQ$Name = as.character(mouseITRAQ$Name)
mouseITRAQ$Name = gsub("^ ","", mouseITRAQ$Name)
mouseITRAQ$Name = gsub(" .*","", mouseITRAQ$Name) 
mouseITRAQ.c = collapseRows(mouseITRAQ[,-c(1,2)], rowGroup=mouseITRAQ$Name, rowID=mouseITRAQ$Description, connectivityBasedCollapsing = TRUE)$datETcollapsed # 12698 genes
# get rid of TaxIR, HumIR, WHIM13.1
mouseITRAQ.c_1 = mouseITRAQ.c[,-c(17,18,20)]
#write.table(mouseITRAQ.c_1 , col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/mouse-only/proteome-ratio-norm_collapsed.txt')

### ITRAQ proteome ###
# ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
ITRAQ = read.table(file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/whim-proteome-v3.01/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 13719 isoforms
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description, connectivityBasedCollapsing = TRUE)$datETcollapsed # 12698 genes
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ.c_1 = ITRAQ.c[,-c(17,18,20)]
#write.table(ITRAQ.c_1 , col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
# ### LFQ proteome ###
# LFQ=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified',header=TRUE, sep="\t", fill=T)
# # collapse protein isoforms: use the maximum value if two; if more choose most representative
# row.names(LFQ) = LFQ$RefSeq_id #8648 protein isoforms
# LFQ.c = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id, connectivityBasedCollapsing = TRUE)$datETcollapsed #5802 genes
# # get rid of WHIM16_2, WHIM2_2
# LFQ.c_1 = LFQ.c[,-c(6,9)]
# #write.table(LFQ.c_1, col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt')

### OUTLIER IN ITRAQ PHOSPHO###
# average the values of all phophosites in that gene
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho = ITRAQpho[,-c(19,20,22)]
ITRAQpho.c = collapseRows(ITRAQpho[,-c(1,2)], rowGroup=ITRAQpho$gene, rowID=row.names(ITRAQpho), method= "Average")$datETcollapsed
#write.table(ITRAQpho.c, col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt')

# one representative phosphosite per gene
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho = ITRAQpho[,-c(19,20,22)]
ITRAQpho.c = collapseRows(ITRAQpho[,-c(1,2)], rowGroup=ITRAQpho$gene, rowID=row.names(ITRAQpho), method= "maxRowVariance", connectivityBasedCollapsing = TRUE)
ITRAQpho.collapsed = ITRAQpho.c$datETcollapsed
ITRAQpho.NP = ITRAQpho.c$group2row
# merge
dim(ITRAQpho.NP)
#[1] 8766    2
dim(ITRAQpho.collapsed)
#[1] 8766   27
ITRAQpho.c.m = merge(ITRAQpho.NP, ITRAQpho.collapsed, by="row.names")
dim(ITRAQpho.c.m)
#[1] 8766   30
row.names(ITRAQpho.c.m) = ITRAQpho.c.m$selectedRowID
ITRAQpho.c.m = ITRAQpho.c.m[,-c(1,2,3)]
row.names(ITRAQpho.c.m) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.c.m)), unique=T)
row.names(ITRAQpho.c.m) = make.names(sub("\\. _.*","",row.names(ITRAQpho.c.m)), unique=T)
row.names(ITRAQpho.c.m) = make.names(sub("\\._.*","",row.names(ITRAQpho.c.m)), unique=T)
row.names(ITRAQpho.c.m) = make.names(sub(" _.*","",row.names(ITRAQpho.c.m)), unique=T)
row.names(ITRAQpho.c.m) = make.names(sub("_.*","",row.names(ITRAQpho.c.m)), unique=T)
#write.table(ITRAQpho.c.m, col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_selected.txt')

# ### CNV ###
# #CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified',header=TRUE, sep="\t")
# CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM_WUPCC_shared/WHIM_CNV/201603/2016-03-24_WHIM_CNV.tsv',header=TRUE, sep="\t")
# row.names(CNV)=CNV$gene
# colnames(CNV) = sub("WHIM0","WHIM",colnames(CNV))
# CNV=CNV[,-1] #16209 genes
# 
# CNV.m=as.matrix(CNV) 
# # log transform: log(cn / <cn>), where <cn> is the mean
# CNV.n.m = CNV.m
# for (i in 1:nrow(CNV.n.m)){
#   CNV.n.m[i,]=log(CNV.n.m[i,]/mean(CNV.n.m[i,]), base=10)
# } 
# 
# #write.table(CNV.n.m, quote = F, col.names=NA, sep ="\t", file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM_WUPCC_shared/WHIM_CNV/201603/2016-03-24_WHIM_CNV_log10_normalized.tsv")
# 
# CNV.m2 = log(CNV.m, base = 2) - 1
# #write.table(CNV.m2, quote = F, col.names=NA, sep ="\t", file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM_WUPCC_shared/WHIM_CNV/201603/2016-03-24_WHIM_CNV_log2_transformed.tsv")
# 
# ### human proteome ###
# ITRAQ_human = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
# row.names(ITRAQ_human) = ITRAQ_human$Description
# ITRAQ_human.c = collapseRows(ITRAQ_human[,-c(1,2,3)], rowGroup=ITRAQ_human$Gene, rowID=ITRAQ_human$Description, connectivityBasedCollapsing = TRUE)$datETcollapsed
# colnames(ITRAQ_human.c) = sub("\\.", "-", colnames(ITRAQ_human.c))
# colnames(ITRAQ_human.c) = sub("\\..*", "", colnames(ITRAQ_human.c))
# #write.table(ITRAQ_human.c, col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt')
# 
# ### human phosphoproteome ###
# ITRAQpho_human = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# # average the values of all phophosites in that gene
# row.names(ITRAQpho_human) = ITRAQpho_human$Gene.site
# ITRAQpho_human.c = collapseRows(ITRAQpho_human[,-c(1,2)], rowGroup=ITRAQpho_human$Gene, rowID=row.names(ITRAQpho_human), method= "Average")$datETcollapsed
# colnames(ITRAQpho_human.c) = sub("\\.", "-", colnames(ITRAQpho_human.c))
# colnames(ITRAQpho_human.c) = sub("\\..*", "", colnames(ITRAQpho_human.c))
# #write.table(ITRAQpho_human.c, col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm_collapsed.txt')
# 
# # one representative phosphosite per gene
# ITRAQpho_human = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# row.names(ITRAQpho_human) = ITRAQpho_human$Gene.site
# ITRAQpho_human.c = collapseRows(ITRAQpho_human[,-c(1,2)], rowGroup=ITRAQpho_human$Gene, rowID=row.names(ITRAQpho_human), method= "maxRowVariance", connectivityBasedCollapsing = TRUE)
# ITRAQpho_human.collapsed = ITRAQpho_human.c$datETcollapsed
# colnames(ITRAQpho_human.collapsed) = sub("\\.", "-", colnames(ITRAQpho_human.collapsed))
# colnames(ITRAQpho_human.collapsed) = sub("\\..*", "", colnames(ITRAQpho_human.collapsed))
# ITRAQpho_human.NP = ITRAQpho_human.c$group2row
# # merge
# dim(ITRAQpho_human.NP)
# #[1] 6653    2
# dim(ITRAQpho_human.collapsed)
# ITRAQpho_human.c.m = merge(ITRAQpho_human.NP, ITRAQpho_human.collapsed, by="row.names")
# dim(ITRAQpho_human.c.m)
# row.names(ITRAQpho_human.c.m) = ITRAQpho_human.c.m$selectedRowID
# ITRAQpho_human.c.m = ITRAQpho_human.c.m[,-c(1,2,3)]
# row.names(ITRAQpho_human.c.m) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho_human.c.m)), unique=T)
# row.names(ITRAQpho_human.c.m) = make.names(sub("\\. _.*","",row.names(ITRAQpho_human.c.m)), unique=T)
# row.names(ITRAQpho_human.c.m) = make.names(sub("\\._.*","",row.names(ITRAQpho_human.c.m)), unique=T)
# row.names(ITRAQpho_human.c.m) = make.names(sub(" _.*","",row.names(ITRAQpho_human.c.m)), unique=T)
# row.names(ITRAQpho_human.c.m) = make.names(sub("_.*","",row.names(ITRAQpho_human.c.m)), unique=T)
# #write.table(ITRAQpho_human.c.m, col.names=NA, quote=F, sep = '\t', file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm_selected.txt')
# 
