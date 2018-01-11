# overlap
# find top from expression data --> applicable to any protein list of interest
# (1) Find proteins with NA in more than x % of samples & average expression lower than y # get rid of low expressors or ones that are not measured well
# Filter for proteins whose expression is 1SD away from the distribution # find outliers across samples (or use “inner” and “outer” fences, Let qr(i) be the rth percentile of the x i j values for gene i, and IQR(i) = q75(i)−q25(i), the interquartile range. Finally, note that values greater than the limit q75(i) + IQR(i) are defined to be outliers in the usual statistical sense.)
# (2) Filter for proteins whose expression is in the top 20% quantile # find outlier within sample, or whose rank changed??
# (3)	Rank by normalized expression values
# (4)	Downstream targets are phosphorylated? (differential expression FDR cut-off)
# nice note for outlier analysis: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.html

### dependencies ###
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/outlier")
source("/Users/khuang/bin/LIB_exp.R")

### outlier related functions ###
## to-do: find the probability density at the threshold point and then define p value?
# http://stats.stackexchange.com/questions/154015/r-kernel-discriminant-analysis-and-how-to-get-the-kernel-cumulative-density-est

### OUTLIER IN ITRAQ PROTEOME ###
ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 13719 isoforms
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$datETcollapsed
# 12698 genes
# ITRAQ.NP = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$group2row
#rm(ITRAQ)
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ.c_1 = ITRAQ.c[,-c(17,18,20)]
rm(ITRAQ.c)
ITRAQ.c.na = ITRAQ.c_1[rowSums(is.na(ITRAQ.c_1)) <= 14,] #10588 genes
rm(ITRAQ.c_1)
#write.table(ITRAQ.c_1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt', quote=F, sep = '\t')

ITRAQ.c.na.d = ITRAQ.c.na[row.names(ITRAQ.c.na) %in% t(drugList), ] # 116 genes
ITRAQ.c.na.d.m = as.matrix(ITRAQ.c.na.d)

# find outliers
outliers=find_outlier(ITRAQ.c.na.d.m, "ITRAQ druggable proteome")
outlier_zscore=outliers$outlier_zscore
outlier=outliers$outlier

# find modified z-score distribution
outlier_zscore_m = melt(outlier_zscore)


# # plot
# pdf(paste(pd,'WHIM_proteome-ratio-norm_ITRAQ_mzscore_distribution.pdf', sep="_"))
# ggplot(outlier_zscore_m, aes(x=value)) + geom_density(aes(colour=X2, fill=X2), alpha=0.1) + 
#   scale_colour_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(27)) +
#   xlab("modified Z-score in iTRAQ proteome") +
#   theme(text = element_text(colour="black", size=16),axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
# dev.off()
# 
# pdf(paste(pd,'all_WHIM_proteome-ratio-norm_ITRAQ_mzscore_distribution.pdf', sep = "_"))
# p = ggplot(outlier_zscore_m, aes(x=as.numeric(outlier_zscore_m$value))) + geom_density()#aes(colour=X2, fill=X2), alpha=0.1) + 
#   #scale_colour_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(27)) +
# p = p + geom_vline(xintercept=maxs, colour="red", alpha=0.2) + geom_vline(xintercept=mins, colour="blue", alpha=0.2) + geom_vline(xintercept=adj_peak, colour="red")
# #p = p + scale_fill_manual(name="Inflection point",values=c("red", "blue"), labels = c("Local maximum","Local minimum"))
# p = p + xlab("modified Z-score in iTRAQ proteome") + 
#   scale_x_continuous(breaks=maxs) + theme_bw() + 
#   theme(text = element_text(colour="black", size=16), axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
# p
# dev.off()

# find top 5% of ITRAQ outliers
top_outliers = find_top_outlier(outliers, "ITRAQ proteome", matrix=ITRAQ.c.na.d.m)
top_outlier_zscore = top_outliers$top_outlier_zscore
top_outlier = top_outliers$top_outlier
top_outlier_boolean = top_outliers$top_outlier_boolean

## MIS ##
# grep a gene list of choice and find enrichment of outlier in that list
# at the same time maybe get rid of all the members not on the list
# pro1_in_Klist = pro_1na14_2[pro_1na14$Gene %in% t(kinaseList), ] # 505 isoforms
# outlier_in_Klist = outlier[pro_1na14$Gene %in% t(kinaseList), ]
# 
# pro1_in_Dlist = pro_1na14_2[pro_1na14$Gene %in% t(drugList), ] # 119 isoforms
# outlier_in_Dlist = outlier[pro_1na14$Gene %in% t(drugList), ] 
# 
# # how many genes are outliers in each list
# 
# all = dim(outlier)[1]
# all_outlier = dim(outlier[rowSums(outlier, na.rm=T) >= 1,])[1] # 3197 isoforms that have at least 1 outliers, 4534 in total
# kinase_all = dim(pro1_in_Klist)[1]
# kinase_outlier = dim(outlier_in_Klist[rowSums(outlier_in_Klist, na.rm=T) >= 1,])[1] # 155 isoforms that have at least 1 outliers, 220 in total
# drug_all = dim(pro1_in_Dlist)[1]
# drug_outlier =dim(outlier_in_Dlist[rowSums(outlier_in_Dlist, na.rm=T) >= 1,])[1] # 38 isoforms that have at least 1 outliers, 64 in total
# # outlier enrichment: gene level
# t.test(matrix(c(all, all_outlier, kinase_all, kinase_outlier), nrow=2)) # gene level enrichment: p=0.2366
# t.test(matrix(c(all, all_outlier, drug_all, drug_outlier), nrow=2)) # gene level enrichment: p=0.2576
# 
# # outlier enrichment: gene & sample count level
# all_outlier_c = sum(rowSums(outlier, na.rm=T))
# kinase_outlier_c = sum(rowSums(outlier_in_Klist, na.rm=T))
# drug_outlier_c = sum(rowSums(outlier_in_Dlist, na.rm=T))
# 
# t.test(matrix(c(all, all_outlier_c, kinase_all, kinase_outlier_c), nrow=2)) # count level enrichment: p=0.2036
# t.test(matrix(c(all, all_outlier_c, drug_all, drug_outlier_c), nrow=2)) # count level enrichment: p=0.2245
# 
# rownames(outlier[rowSums(outlier, na.rm=T) >= 3,]) # this list enriched for GO protein binding P=1.5*10E-11

### OUTLIER IN LFQ###
LFQ=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified',header=TRUE, sep="\t", fill=T)

# collapse protein isoforms: use the maximum value if two; if more choose most representative
row.names(LFQ) = LFQ$RefSeq_id #8648 protein isoforms
LFQ.c = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$datETcollapsed #5802 genes
#LFQ.NP = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$group2row
rm(LFQ)
# get rid of WHIM16_2, WHIM2_2
LFQ.c_1 = LFQ.c[,-c(6,9)]
rm(LFQ.c)
#write.table(LFQ.c_1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt', quote=F, sep = '\t')
LFQ.c.na = LFQ.c_1[rowSums(is.na(LFQ.c_1)) <= 8,] #4439 genes
rm(LFQ.c_1)

LFQ.c.na.d = LFQ.c.na[row.names(LFQ.c.na) %in% t(drugList), ] #61 genes
LFQ.c.na.d.m = as.matrix(LFQ.c.na.d)

# find outliers
LFQ_outliers=find_outlier(LFQ.c.na.d.m, "LFQ proteome")
LFQ_outlier_zscore=LFQ_outliers$outlier_zscore
LFQ_outlier=LFQ_outliers$outlier

# find top 5% of LFQ proteome outliers
LFQ_top_outliers = find_top_outlier(LFQ_outliers, "LFQ proteome", matrix=LFQ.c.na.d.m)
LFQ_top_outlier_zscore = LFQ_top_outliers$top_outlier_zscore
LFQ_top_outlier = LFQ_top_outliers$top_outlier
LFQ_top_outlier_boolean = LFQ_top_outliers$top_outlier_boolean

### overlap between ITRAQ and LFQ
LFQ_top_outlier2=LFQ_top_outlier[,1:5]
LFQ_top_outlier_zscore2=LFQ_top_outlier_zscore[,1:5]
LFQ_top_outlier_boolean2=LFQ_top_outlier_boolean[,1:5]
top_outlier2=top_outlier[,1:5]
top_outlier_zscore2=top_outlier_zscore[,1:5]
top_outlier_boolean2=top_outlier_boolean[,1:5]
top_overlap = merge(LFQ_top_outlier2, top_outlier2, by = "row.names", all.x=T)
colnames(top_overlap)[1] = "WHIM"
top_zscore_overlap = merge(LFQ_top_outlier_zscore2, top_outlier_zscore2, by = "row.names", all.x=T)
colnames(top_zscore_overlap)[1] = "WHIM"
top_boolean_overlap = merge(LFQ_top_outlier_boolean2, top_outlier_boolean2, by = "row.names", all.x=T)
colnames(top_boolean_overlap)[1] = "WHIM"


# plot the overlap 
top_overlap.m <- melt(top_overlap, id.var = "WHIM")
top_zscore_overlap.m <- melt(top_zscore_overlap, id.var = "WHIM")
top_boolean_overlap.m <- melt(top_boolean_overlap, id.var = "WHIM")
colnames(top_boolean_overlap.m)[3]="outlier"
colnames(top_zscore_overlap.m)[3]="expression_score"

pdf(paste(pd,'ITRAQ_LFQ_top5percent_mzscore_proteome.pdf',sep ="_"), height=10, width=8)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)
outlier.colors=c("NA", "#000000")

p = ggplot()
p = p + geom_tile(data=top_zscore_overlap.m, aes(x=as.factor(WHIM), y=variable, fill=expression_score), linetype="blank") + scale_fill_gradientn(colours=getPalette(100))
p = p + geom_tile(data=top_boolean_overlap.m, aes(x=as.factor(WHIM), y=variable, color=outlier), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors)
p = p + geom_text(data=top_overlap.m,aes(x=as.factor(WHIM), y=variable, label = value), color="red", size=3, angle=90)
p = p + xlab("Sample") + ylab("Top druggable protein") + theme_bw() + 
  theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(,colour="black", size=12))
p
dev.off()

### mis: find correlation between ITRAQ and LFQ proteomes ###

# ITRAQ.c.na.m = melt(ITRAQ.c.na)
# LFQ.c.na.m = melt(LFQ.c.na)
# ITRAQ_LFQ_merge = merge(ITRAQ.c.na.m,LFQ.c.na.m, by=c("X1","X2"))
# cor(ITRAQ_LFQ_merge$value.x,ITRAQ_LFQ_merge$value.y, use = 'pairwise.complete.obs', method="pearson")
# # [1,] 0.2040447
# png('ITRAQ_LFQ_correlation.png')
# p = ggplot(ITRAQ_LFQ_merge, aes(x=value.x, y=value.y, color=X2)) + 
#   geom_point(alpha=0.2) + geom_smooth(method=lm, se=F, alpha=0.2) + ggtitle("Correlation between ITRAQ and LFQ proteome: 0.204") +
#   xlab("ITRAQ proteome level") + ylab("LFQ proteome level") + xlim(-15,15) + ylim(-15,15)
# p
# dev.off()
# 
# # correlation in druggable proteins
# ITRAQ.c.na.m = melt(ITRAQ.c.na.d)
# LFQ.c.na.m = melt(LFQ.c.na.d)
# ITRAQ_LFQ_d.merge = merge(ITRAQ.c.na.m,LFQ.c.na.m, by=c("X1","X2"), all=F)
# cor(ITRAQ_LFQ_d.merge$value.x,ITRAQ_LFQ_d.merge$value.y, use = 'pairwise.complete.obs', method="pearson")
# # [1,] 0.1628424
# png('ITRAQ_LFQ_druggable_correlation.png')
# p = ggplot(ITRAQ_LFQ_d.merge, aes(x=value.x, y=value.y, color=X2)) + 
#   geom_point(alpha=0.2) + geom_smooth(method=lm, se=F, alpha=0.2) + ggtitle("Correlation between ITRAQ and LFQ druggable proteome: 0.163") +
#   xlab("ITRAQ proteome level") + ylab("LFQ proteome level") + xlim(-15,15) + ylim(-15,15)
# p
# dev.off()

# correlation of z-scores would make sense since iTRAQ ratios are based to reference
# z_score_merge=merge(outlier_zscore_m,LFQ_outlier_zscore_m, by=c("X1","X2"), all=F)
# cor(z_score_merge$value.x,z_score_merge$value.y, use = 'pairwise.complete.obs', method="pearson")
# # [1,] 0.6374896
# cor(z_score_merge$value.x,z_score_merge$value.y, use = 'pairwise.complete.obs', method="spearman")
# # [1,] 0.6496502
# pdf(paste(pd,'ITRAQ_LFQ_druggable_zscore_correlation.pdf', sep="_"))
# p = ggplot(z_score_merge, aes(x=value.x, y=value.y, color=X2)) + 
#   geom_point(alpha=0.2) + geom_smooth(method=lm, se=F, alpha=0.2) 
# p = p + theme_bw() + ggtitle("Z-score correlation between ITRAQ and LFQ druggable proteome: 0.65") +
#   xlab("ITRAQ proteome z-score") + ylab("LFQ proteome z-score") + xlim(-5,5) + ylim(-5,5)
# p
# dev.off()

### OUTLIER IN ITRAQ PHOSPHO###
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))

# output a collapsed table for figure 1

#ITRAQpho.c = collapseRows(ITRAQpho[,-1], rowGroup=ITRAQpho$Gene, rowID=row.names(ITRAQpho))$datETcollapsed
# 12698 genes
#ITRAQpho.NP = collapseRows(ITRAQpho[,-1], rowGroup=ITRAQpho$Gene, rowID=row.names(ITRAQpho))$group2row
#ITRAQpho.c1 = ITRAQpho.c[,-c(17,18,20)]
#write.table(ITRAQpho.c1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt', quote=F, sep = '\t')


# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]

# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho = ITRAQpho[,-c(17,18,20)]

ITRAQpho.na = ITRAQpho[rowSums(is.na(ITRAQpho)) <= 14,] #35838 phosphosites
rm(ITRAQpho)

genes = sub("-NP.*", "", row.names(ITRAQpho.na))
ITRAQpho.na.d = ITRAQpho.na[genes %in% t(drugList), ] # 523 phosphosites

row.names(ITRAQpho.na.d) = sub("-NP_\\d+_"," ",row.names(ITRAQpho.na.d))
row.names(ITRAQpho.na.d) = make.names(sub(" _.*","",row.names(ITRAQpho.na.d)), unique=T)
row.names(ITRAQpho.na.d) = make.names(sub("_.*","",row.names(ITRAQpho.na.d)), unique=T)
ITRAQpho.na.d.m = as.matrix(ITRAQpho.na.d)

# use the adjusted z-score to find outliers
ITRAQpho_outliers=find_outlier(ITRAQpho.na.d.m, "ITRAQ phosphoproteome")
ITRAQpho_outlier_zscore=ITRAQpho_outliers$outlier_zscore
ITRAQpho_outlier=ITRAQpho_outliers$outlier

# # find modified z-score distribution
# ITRAQpho_outlier_zscore_m = melt(ITRAQpho_outlier_zscore)
# pdf(paste(pd,'WHIM_phosphoproteome-ratio-norm_ITRAQ_mzscore_distribution.pdf', sep="_"))
# ggplot(ITRAQpho_outlier_zscore_m, aes(x=value)) + geom_density(aes(colour=X2, fill=X2), alpha=0.1) + 
#   scale_colour_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(27)) +
#   xlab("modified Z-score in iTRAQ druggable phoshpoproteome") +
#   theme(text = element_text(colour="black", size=16),axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
# dev.off()
# 
# pdf(paste(pd,'all_WHIM_phosphoproteome-ratio-norm_ITRAQ_mzscore_distribution.pdf', sep="_"))
# p = ggplot(ITRAQpho_outlier_zscore_m, aes(x=as.numeric(ITRAQpho_outlier_zscore_m$value))) + geom_density()#aes(colour=X2, fill=X2), alpha=0.1) + 
# #scale_colour_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(27)) +
# p = p + xlab("modified Z-score in iTRAQ proteome") + 
#   scale_x_continuous(breaks=c(-4,-2,0,1.5,2,2.25,2.5,5)) + theme_bw() + 
#   theme(text = element_text(colour="black", size=16), axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
# p
# dev.off()

# find top 5% of ITRAQ phosphosite outliers
ITRAQpho_top_outliers = find_top_outlier(ITRAQpho_outliers, "ITRAQ phosphoproteome", matrix=ITRAQpho.na.d.m)
ITRAQpho_top_outlier_zscore = ITRAQpho_top_outliers$top_outlier_zscore
ITRAQpho_top_outlier = ITRAQpho_top_outliers$top_outlier
ITRAQpho_top_outlier_boolean = ITRAQpho_top_outliers$top_outlier_boolean

### OUTLIER IN LFQ phospho###
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)

LFQpho$sites=sub(".*\\(","",LFQpho$phospho_site)
LFQpho$sites=sub("\\)","",LFQpho$site)
LFQpho$sites=paste(LFQpho$gene_name,LFQpho$sites, sep=".")

row.names(LFQpho) = make.names(LFQpho$sites, unique=T)
colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
# system(paste("sed 's/NP//g' <<<", row.names(ITRAQpho))) doesn't work yet
# 18229 phosphosites
LFQpho = LFQpho[,-c(19:25)]

LFQpho.na = LFQpho[rowSums(is.na(LFQpho)) <= 8,] #3763 phosphosites! This is a great loss.
rm(LFQpho)

genes = sub("\\..*", "", row.names(LFQpho.na))
LFQpho.na.d = LFQpho.na[genes %in% t(drugList), ] #42 phosphosites
LFQpho.na.d.m = as.matrix(LFQpho.na.d)

# find outliers
LFQpho_outliers=find_outlier(LFQpho.na.d.m, "LFQ phosphoproteome")
LFQpho_outlier_zscore=LFQpho_outliers$outlier_zscore
LFQpho_outlier=LFQpho_outliers$outlier

# LFQpho_outlier_zscore = as.matrix(LFQpho_outlier_zscore)
# # find modified z-score distribution
# LFQpho_outlier_zscore_m = melt(LFQpho_outlier_zscore)
# pdf(paste(pd,'all_WHIM_LFQpho_Global_minimum1_norm_nameadded_human_cleaned_mzscore_distribution.pdf', sep = "_"))
# ggplot(LFQpho_outlier_zscore_m, aes(x=value)) + geom_density(aes(colour=X2, fill=X2), alpha=0.1) + 
#   scale_colour_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(27)) +
#   xlab("modified Z-score in LFQ phosphoproteome") +
#   theme(text = element_text(colour="black", size=14),axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))
# dev.off()

# find top 5% z-score and genes within each sample
LFQpho_top_outliers = find_top_outlier(LFQpho_outliers, "LFQ phosphoproteome", matrix=LFQpho.na.d.m)
LFQpho_top_outlier_zscore = LFQpho_top_outliers$top_outlier_zscore
LFQpho_top_outlier = LFQpho_top_outliers$top_outlier
LFQpho_top_outlier_boolean = LFQpho_top_outliers$top_outlier_boolean

### overlap between ITRAQ and LFQ ###

ITRAQpho_top_outlier=ITRAQpho_top_outlier[,1:5]
top_pho_overlap = merge(LFQpho_top_outlier, ITRAQpho_top_outlier, by = "row.names", all.x=T)
colnames(top_pho_overlap)[1] = "WHIM"

ITRAQpho_top_outlier_zscore=ITRAQpho_top_outlier_zscore[,1:5]
top_pho_zscore_overlap = merge(LFQpho_top_outlier_zscore, ITRAQpho_top_outlier_zscore, by = "row.names", all.x=T)
colnames(top_pho_zscore_overlap)[1] = "WHIM"

ITRAQpho_top_outlier_boolean=ITRAQpho_top_outlier_boolean[,1:5]
top_pho_boolean_overlap = merge(LFQpho_top_outlier_boolean, ITRAQpho_top_outlier_boolean, by = "row.names", all.x=T)
colnames(top_pho_boolean_overlap)[1] = "WHIM"

## plot the overlap 
top_pho_overlap.m <- melt(top_pho_overlap, id.var = "WHIM")
top_pho_zscore_overlap.m <- melt(top_pho_zscore_overlap, id.var = "WHIM")
top_pho_boolean_overlap.m <- melt(top_pho_boolean_overlap, id.var = "WHIM")

pdf(paste(pd,'ITRAQ_LFQ_top5percent_mzscore_phosphoproteome.pdf', sep="_"), height=6, width=12)
OrRd = brewer.pal(9, "OrRd") 
getPalette = colorRampPalette(OrRd)

p = ggplot()
p = p + geom_tile(data=top_pho_zscore_overlap.m, aes(x=as.factor(WHIM), y=variable, fill=value), linetype="blank") + scale_fill_gradientn(colours=getPalette(100))
p = p + geom_tile(data=top_pho_boolean_overlap.m, aes(x=as.factor(WHIM), y=variable, color=value), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors)
p = p + geom_text(data=top_pho_overlap.m,aes(x=as.factor(WHIM), y=variable, label = value), size=1.5, angle=45)
p = p + xlab("WHIM sample") + ylab("Top 5% druggable phosphosites") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(colour="black", size=12))
p
dev.off()

### find outlier only within ITRAQ proteome and phosphoproteome so we don't miss out on some samples### 
ITRAQpho_top_outlier=ITRAQpho_top_outlier[,1:5]
top_ITRAQ_overlap = merge(top_outlier, ITRAQpho_top_outlier, by = "row.names", all.x=T)
colnames(top_ITRAQ_overlap)[1] = "WHIM"

ITRAQpho_top_outlier_zscore=ITRAQpho_top_outlier_zscore[,1:5]
top_ITRAQ_zscore_overlap = merge(top_outlier_zscore, ITRAQpho_top_outlier_zscore, by = "row.names", all.x=T)
colnames(top_ITRAQ_zscore_overlap)[1] = "WHIM"

ITRAQpho_top_outlier_boolean=ITRAQpho_top_outlier_boolean[,1:5]
top_ITRAQ_boolean_overlap = merge(top_outlier_boolean, ITRAQpho_top_outlier_boolean, by = "row.names", all.x=T)
colnames(top_ITRAQ_boolean_overlap)[1] = "WHIM"

fn = paste(pd,'ITRAQ_top5percent_pro_phospho.txt', sep="_")
write.table(top_ITRAQ_overlap, file=fn, quote=F, row.names=F, sep="\t")

fn = paste(pd,'ITRAQ_top5percent_pro_phospho_zscore.txt', sep="_")
write.table(top_ITRAQ_zscore_overlap, file=fn, quote=F, row.names=F, sep="\t")

## plot the overlap 
top_ITRAQ_overlap.m <- melt(top_ITRAQ_overlap, id.var = "WHIM")
top_ITRAQ_zscore_overlap.m <- melt(top_ITRAQ_zscore_overlap, id.var = "WHIM")
top_ITRAQ_boolean_overlap.m <- melt(top_ITRAQ_boolean_overlap, id.var = "WHIM")

pdf(paste(pd,'ITRAQ_top5percent_mzscore_pro_phosphoproteome.pdf', sep="_"), height=6, width=12)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)

p = ggplot()
p = p + geom_tile(data=top_ITRAQ_zscore_overlap.m, aes(x=as.factor(WHIM), y=variable, fill=value), linetype="blank") + scale_fill_gradientn(colours=getPalette(100))
p = p + geom_tile(data=top_ITRAQ_boolean_overlap.m, aes(x=as.factor(WHIM), y=variable, color=value), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors)
p = p + geom_text(data=top_ITRAQ_overlap.m,aes(x=as.factor(WHIM), y=variable, label = value), size=1.5, angle=45, colour="red")
p = p + xlab("WHIM sample") + ylab("Top 5% druggable ITRAQsITRAQsites") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(colour="black", size=12))
p
dev.off()

### find extreme phosphos that are not indicated by pro 
# phospho (ITRAQpho_outlier_zscore) and pro (outlier_zscore) in iTRAQ of the same sample, see how they correlate in z-score
ITRAQpho_outlier_zscore.dt = as.data.frame(ITRAQpho_outlier_zscore)
ITRAQpho_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(ITRAQpho_outlier_zscore.dt))
ITRAQpho_outlier_zscore.dt$Site = row.names(ITRAQpho_outlier_zscore.dt)
outlier_zscore.dt = as.data.frame(outlier_zscore)
outlier_zscore.dt$Gene = sub("\\..*", "", row.names(outlier_zscore.dt))

ITRAQpho_outlier_zscore.m = melt(ITRAQpho_outlier_zscore.dt, id.var = c("Gene","Site"))
outlier_zscore.dt.m = melt(outlier_zscore.dt, id.var = "Gene")
pro_pho_overlap = merge(ITRAQpho_outlier_zscore.m, outlier_zscore.dt.m, by = c("Gene","variable"), all.x=T)

cor(pro_pho_overlap$value.x, pro_pho_overlap$value.y, use = 'pairwise.complete.obs')

pro_pho_overlap$gene = rep("Other", nrow(pro_pho_overlap))
pro_pho_overlap[pro_pho_overlap$Gene == "ERBB2",]$gene="ERBB2"
pro_pho_overlap[pro_pho_overlap$Gene == "PAK1",]$gene="PAK1"

pdf(paste(pd,'ITRAQ_phospho_vs_pro_mzscore.pdf', sep="_"))

p = ggplot(data = pro_pho_overlap, aes(x=value.y, y=value.x, colour=variable, shape=gene)) 
p = p + geom_point(alpha=0.3, size=1.5) + xlab("Protein expression modified z-score") + ylab("Phosphosite expression modified z-score") 
p = p + scale_shape_manual(values=c(0,16,2)) + guides(colour=FALSE)
p = p + theme_bw() + coord_fixed() #+ theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
p = p + xlim(-5,7.5) + ylim(-5,7.5) + geom_vline(xintercept=2.158469739,alpha=0.2) + geom_hline(yintercept=3.137479213,alpha=0.2)
#p = p + geom_text(aes(label=ifelse(value.x > 3 & value.x > value.y + 2,Site,'')),hjust=-0.05,just=0,size=3,alpha=0.7)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
dev.off()
# p2 = ggplot(data = pro_pho_overlap,aes(x=value.x)) + 
#   geom_density(alpha=0.5) + 
#   scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
#   scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
#   theme_bw() +
#   theme0(plot.margin = unit(c(1,0,0,2.2),"lines"))
# 
# p3 = ggplot(data = pro_pho_overlap,aes(x=value.y)) + 
#   geom_density(alpha=0.5) + 
#   scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
#   scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
#   theme_bw() + 
#   coord_flip()  +
#   theme0(plot.margin = unit(c(1,0,0,2.2),"lines"))
# 
# grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
#              arrangeGrob(p,p3,ncol=2,widths=c(3,1)),
#              heights=c(1,3))
#http://stackoverflow.com/questions/17370460/scatterplot-with-alpha-transparent-histograms-in-r



### RNA-Seq ###
RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes

RNA.d = RNA[row.names(RNA) %in% t(drugList), ] #142 genes
RNA.d.m=as.matrix(RNA.d) 

# find outliers
RNA_outliers=find_outlier(RNA.d.m, "mRNA expression")
RNA_outlier_zscore=RNA_outliers$outlier_zscore
RNA_outlier=RNA_outliers$outlier

# find top 5% of mRNA outliers
RNA_top_outliers = find_top_outlier(RNA_outliers, "mRNA", matrix=RNA.d.m)
RNA_top_outlier_zscore = RNA_top_outliers$top_outlier_zscore
RNA_top_outlier = RNA_top_outliers$top_outlier
RNA_top_outlier_boolean = RNA_top_outliers$top_outlier_boolean

### overlap between mRNA and ITRAQ proteome
top_overlap = merge(RNA_top_outlier, top_outlier, by = "row.names")
colnames(top_overlap)[1] = "WHIM"
top_zscore_overlap = merge(RNA_top_outlier_zscore, top_outlier_zscore, by = "row.names")
colnames(top_zscore_overlap)[1] = "WHIM"
top_boolean_overlap = merge(RNA_top_outlier_boolean, top_outlier_boolean, by = "row.names")
colnames(top_boolean_overlap)[1] = "WHIM"

top_overlap.m <- melt(top_overlap, id.var = "WHIM")
top_zscore_overlap.m <- melt(top_zscore_overlap, id.var = "WHIM")
top_boolean_overlap.m <- melt(top_boolean_overlap, id.var = "WHIM")

pdf(paste(pd,'mRNA_ITRAQ_top5percent_mzscore_proteome.pdf',sep ="_"), height=6, width=12)
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)
outlier.colors=c("NA", "#000000")

p = ggplot()
p = p + geom_tile(data=top_zscore_overlap.m, aes(x=as.factor(WHIM), y=variable, fill=value), linetype="blank") + scale_fill_gradientn(colours=getPalette(100))
p = p + geom_tile(data=top_boolean_overlap.m, aes(x=as.factor(WHIM), y=variable, color=value), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors)
p = p + geom_text(data=top_overlap.m,aes(x=as.factor(WHIM), y=variable, label = value), color="red", size=2.5, angle=45)
p = p + xlab("WHIM sample") + ylab("Top 5% druggable gene") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(colour="black", size=12))
p
dev.off()

### find outlier proteome expression not identified by transcriptome ###
outlier_zscore.dt = as.data.frame(unlist(outlier_zscore))
outlier_zscore.dt$Gene = row.names(outlier_zscore.dt)
RNA_outlier_zscore.dt = as.data.frame(unlist(RNA_outlier_zscore))
RNA_outlier_zscore.dt$Gene = row.names(RNA_outlier_zscore.dt)

outlier_zscore.dt.m = melt(outlier_zscore.dt, id.var = c("Gene"))
RNA_outlier_zscore.dt.m = melt(RNA_outlier_zscore.dt, id.var = "Gene")
pro_rna_overlap = merge(RNA_outlier_zscore.dt.m, outlier_zscore.dt.m, by = c("Gene","variable"))

cor(pro_rna_overlap$value.x, pro_rna_overlap$value.y, use = 'pairwise.complete.obs', method="spearman")

pro_rna_overlap$gene = rep("Other", nrow(pro_rna_overlap))
pro_rna_overlap[pro_rna_overlap$Gene == "ERBB2",]$gene="ERBB2"
pro_rna_overlap[pro_rna_overlap$Gene == "AURKA",]$gene="AURKA" #found to be amplified in basal breast tumor

pdf(paste(pd,'rna_vs_pro_mzscore.pdf', sep="_"))
p = ggplot(data = pro_rna_overlap, aes(x=value.x, y=value.y, colour=variable, shape=gene)) 
p = p + geom_point(alpha=0.3, size=1.5) + xlab("mRNA expression modified z-score") + ylab("Protein expression modified z-score") 
p = p + scale_shape_manual(values=c(0,2,16)) + guides(colour=FALSE)
p = p + theme_bw() + coord_fixed() #+ theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))
p = p + xlim(-5,7.5) + ylim(-5,7.5) + geom_vline(xintercept=2.436497542,alpha=0.2) + geom_hline(yintercept=2.158469739,alpha=0.2)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
dev.off()

### CNV ### 
CNV = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified',header=TRUE, sep="\t")
row.names(CNV)=CNV$gene
colnames(CNV) = sub("WHIM0","WHIM",colnames(CNV))
CNV=CNV[,-1] #16209 genes

CNV.d = CNV[row.names(CNV) %in% t(drugList), ] #144 genes
CNV.d.m=as.matrix(CNV.d) 
# log transform: log(cn / <cn>), where <cn> is the mean
CNV.d.n.m = CNV.d.m
for (i in 1:nrow(CNV.d.n.m)){
  CNV.d.n.m[i,]=log(CNV.d.n.m[i,]/mean(CNV.d.n.m[i,]), base=10)
} 

# find outliers
CNV_outliers=find_outlier(CNV.d.n.m, "CNV")
CNV_outlier_zscore=CNV_outliers$outlier_zscore
CNV_outlier=CNV_outliers$outlier

# find top 5% of mCNV outliers
CNV_top_outliers = find_top_outlier(CNV_outliers, "CNV", matrix=CNV.d.m)
CNV_top_outlier_zscore = CNV_top_outliers$top_outlier_zscore
CNV_top_outlier = CNV_top_outliers$top_outlier
CNV_top_outlier_boolean = CNV_top_outliers$top_outlier_boolean



### all levels: CNV, RNA, proteome, phosphoproteome ### 
# keep the 24 WHIMs in ITRAQ proteome
all_levels=list("top_outlier"=top_outlier,"LFQ_top_outlier"=LFQ_top_outlier,"RNA_top_outlier"=RNA_top_outlier,"CNV_top_outlier"=CNV_top_outlier)
top_overlap_all = merge(ITRAQpho_top_outlier, LFQpho_top_outlier, by = "row.names", all.x=T)
for (i in all_levels){
  rownames(top_overlap_all)=top_overlap_all[,1]
  top_overlap_all=top_overlap_all[,-1]
  top_overlap_all = merge(top_overlap_all, i, by = "row.names", all.x=T)
}
colnames(top_overlap_all)[1] = "sample"

all_levels2=list("top_outlier_zscore"=top_outlier_zscore,"LFQ_top_outlier_zscore"=LFQ_top_outlier_zscore,"RNA_top_outlier_zscore"=RNA_top_outlier_zscore,"CNV_top_outlier_zscore"=CNV_top_outlier_zscore)
ITRAQpho_top_outlier_zscore.n = ITRAQpho_top_outlier_zscore/mean(ITRAQpho_top_outlier_zscore)
LFQpho_top_outlier_zscore.n = LFQpho_top_outlier_zscore/mean(LFQpho_top_outlier_zscore)
top_overlap_all_zscore = merge(ITRAQpho_top_outlier_zscore.n, LFQpho_top_outlier_zscore.n, by = "row.names", all.x=T)
for (i in all_levels2){
  rownames(top_overlap_all_zscore)=top_overlap_all_zscore[,1]
  top_overlap_all_zscore=top_overlap_all_zscore[,-1]
  i.n=i/mean(i)
  top_overlap_all_zscore = merge(top_overlap_all_zscore, i.n, by = "row.names", all.x=T)
}
colnames(top_overlap_all_zscore)[1] = "sample"

all_levels3=list("top_outlier_boolean"=top_outlier_boolean,"LFQ_top_outlier_boolean"=LFQ_top_outlier_boolean,"RNA_top_outlier_boolean"=RNA_top_outlier_boolean,"CNV_top_outlier_boolean"=CNV_top_outlier_boolean)
top_overlap_all_boolean = merge(ITRAQpho_top_outlier_boolean, LFQpho_top_outlier_boolean, by = "row.names", all.x=T)
for (i in all_levels3){
  rownames(top_overlap_all_boolean)=top_overlap_all_boolean[,1]
  top_overlap_all_boolean=top_overlap_all_boolean[,-1]
  top_overlap_all_boolean = merge(top_overlap_all_boolean, i, by = "row.names", all.x=T)
}
colnames(top_overlap_all_boolean)[1] = "sample"

top_overlap_all.m <- melt(top_overlap_all, id.var = "sample")
top_overlap_all_zscore.m <- melt(top_overlap_all_zscore, id.var = "sample")
top_overlap_all_boolean.m <- melt(top_overlap_all_boolean, id.var = "sample")

fn=paste(pd,'all_level_top5_mzscore.pdf',sep ="_")
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)
outlier.colors=c("NA", "#000000")

p = ggplot()
p = p + geom_tile(data=top_overlap_all_zscore.m, aes(x=as.factor(sample), y=variable, fill=value), linetype="blank") + scale_fill_gradientn(colours=getPalette(100), na.value="white", name="relative modified z-score")# values=rescale(seq(0,6,by=12/99)))
p = p + geom_tile(data=top_overlap_all_boolean.m, aes(x=as.factor(sample), y=variable, color=value), fill=NA, size=0.5) + scale_colour_manual(values = outlier.colors, name="KDE outlier")
p = p + geom_text(data=top_overlap_all.m ,aes(x=as.factor(sample), y=variable, label = value), color="red", size=1.5, angle=90)
p = p + xlab("sample") + ylab("top druggable targets") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=12), axis.text.y = element_text(colour="black", size=12))
p
ggsave(file=fn, width = 210, height = 400, units = "mm")


##### BRCA77 human data ##### 
###proteome###
BRCA77 = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 12553 isoforms
row.names(BRCA77) = BRCA77$Description
BRCA77.c = collapseRows(BRCA77[,-c(1,2,3)], rowGroup=BRCA77$Gene, rowID=BRCA77$Description)$datETcollapsed
# 11349 genes
BRCA77.c.na = BRCA77.c[rowSums(is.na(BRCA77.c)) <= 67,] #11349 genes
rm(BRCA77)
rm(BRCA77.c)
BRCA77.c.na.d = BRCA77.c.na[row.names(BRCA77.c.na) %in% t(drugList), ] #136 genes
BRCA77.m = as.matrix(BRCA77.c.na.d)

# find outliers
BRCA77_outliers=find_outlier(BRCA77.m, "BRCA ITRAQ proteomes")
BRCA77_outlier_zscore=BRCA77_outliers$outlier_zscore
BRCA77_outlier=BRCA77_outliers$outlier

# find top 5% of mCNV outliers
BRCA77_top_outliers = find_top_outlier(BRCA77_outliers, "BRCA ITRAQ proteomes", h=6, w=25, matrix=BRCA77.m)
BRCA77_top_outlier_zscore = BRCA77_top_outliers$top_outlier_zscore
BRCA77_top_outlier = BRCA77_top_outliers$top_outlier
BRCA77_top_outlier_boolean = BRCA77_top_outliers$top_outlier_boolean


###phosphoproteome###
BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 33239 phosphosites
row.names(BRCA77pho) = BRCA77pho$Gene.site
BRCA77pho = BRCA77pho[,-c(1,2)]
BRCA77pho.na = BRCA77pho[rowSums(is.na(BRCA77pho)) <= 67,] #33239 phosphosites
rm(BRCA77pho)

genes = sub("-NP.*", "", row.names(BRCA77pho.na))
BRCA77pho.na.d = BRCA77pho.na[genes %in% t(drugList), ] #482 phosphosites

row.names(BRCA77pho.na.d) = sub("-NP_\\d+_"," ",row.names(BRCA77pho.na.d))
row.names(BRCA77pho.na.d) = make.names(sub(" _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("_.*","",row.names(BRCA77pho.na.d)), unique=T)

BRCA77pho.m = as.matrix(BRCA77pho.na.d)

# find outliers
BRCA77pho_outliers=find_outlier(BRCA77pho.m, "BRCA ITRAQ phosphoproteomes")
BRCA77pho_outlier_zscore=BRCA77pho_outliers$outlier_zscore
BRCA77pho_outlier=BRCA77pho_outliers$outlier

# find top 5% outliers
BRCA77pho_top_outliers = find_top_outlier(BRCA77pho_outliers, "BRCA ITRAQ phosphoproteomes", h=12, w=25, matrix=BRCA77pho.m)
BRCA77pho_top_outlier_zscore = BRCA77pho_top_outliers$top_outlier_zscore
BRCA77pho_top_outlier = BRCA77pho_top_outliers$top_outlier
BRCA77pho_top_outlier_boolean = BRCA77pho_top_outliers$top_outlier_boolean

### find extreme phosphos that are not indicated by pro 
BRCA77pho_outlier_zscore.dt = as.data.frame(BRCA77pho_outlier_zscore)
BRCA77pho_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(BRCA77pho_outlier_zscore.dt))
BRCA77pho_outlier_zscore.dt$Site = row.names(BRCA77pho_outlier_zscore.dt)
BRCA77_outlier_zscore.dt = as.data.frame(BRCA77_outlier_zscore)
BRCA77_outlier_zscore.dt$Gene = sub("\\..*", "", row.names(BRCA77_outlier_zscore.dt))

BRCA77pho_outlier_zscore.m = melt(BRCA77pho_outlier_zscore.dt, id.var = c("Gene","Site"))
BRCA77_outlier_zscore.dt.m = melt(BRCA77_outlier_zscore.dt, id.var = "Gene")
BRCA77_pro_pho_overlap = merge(BRCA77pho_outlier_zscore.m, BRCA77_outlier_zscore.dt.m, by = c("Gene","variable"), all.x=T)

cor(BRCA77_pro_pho_overlap$value.x, BRCA77_pro_pho_overlap$value.y, use = 'pairwise.complete.obs')

fn = paste(pd,'BRCA77_phospho_vs_pro_mzscore.pdf', sep="_")
p = ggplot(data = BRCA77_pro_pho_overlap, aes(x=value.y, y=value.x, colour=variable, label=Site)) 
p = p + geom_point(alpha=0.3) + xlab("druggable proteome z-score") + ylab("druggable phosphosites z-score") + theme_bw()
p = p + xlim(-5,7.5) + ylim(-5,7.5) + geom_vline(xintercept=2,alpha=0.2) + geom_hline(yintercept=2,alpha=0.2)
p = p + geom_text(aes(label=ifelse(value.x > 2.5 & value.x > value.y + 2 & Gene=="ERBB2",paste(variable,Site),'')),hjust=-0.05,just=0,size=3,alpha=0.7)
p
ggsave(file=fn, height=10, width=10)
