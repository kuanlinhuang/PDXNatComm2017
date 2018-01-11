# outlier analysis pipeline
# nice note for outlier analysis: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.html

### dependencies ###
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/outlier")
source("/Users/khuang/bin/LIB_exp.R")

# rewriting to look at sample-wise z-score distribution to define outliers
# function to find outliers using the modified zscores. Take a matrix as input
# requires function find_kde_inflection
# requires library reshape
# identify the adjacent local minimum as the threshold to define outliers
find_outlier = function(m, name, genes=NULL, minNum = 10, filter=TRUE, plot=FALSE, h=6, w=10){ 
  print("##### OUTLIER ANALYSIS #####")
  num = nrow(m)
  if (filter){
    m2 = m[rowSums(!is.na(m)) >= minNum, ]
    m2_genes = m2[row.names(m2) %in% genes, ]
  } else{ m2 = m; m2_genes=m}
  num_NA = nrow(m2)
  num_gene = nrow(m2_genes)
  m2_genes = as.matrix(m2_genes)
  print(paste("Looking for outliers in", deparse(substitute(genes)), "of", deparse(substitute(m)), sep=" "))
  print(paste("Original number of genes:", num, "; NA filtered:", num_NA, "; Gene list filtered:", num_gene, sep=" "))
  
  kde_dir_cmd = paste("mkdir figures/", date, "/KDE/", sep="")
  system(kde_dir_cmd)
  outlier = matrix(,nrow=dim(m2_genes)[1],ncol=dim(m2_genes)[2])
  row.names(outlier) = row.names(m2_genes)
  colnames(outlier) = colnames(m2_genes)
  outlier_zscore = outlier
  
  # gene-wise zscore
  for (i in 1:nrow(m2_genes)){
    m = median(m2_genes[i,], na.rm=TRUE)
    rowMAD=mad(m2_genes[i,], na.rm=TRUE)
    # modified z-score for outlier: Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    outlier_zscore[i,]  = 0.6745*(m2_genes[i,]-m)/rowMAD
  }
  
  # sample-wise outlier
  for (i in 1:ncol(outlier_zscore)) { 
    sample = colnames(outlier_zscore)[i]
    outlier_zscore_s = outlier_zscore[,i,drop=F]
    
    #     library(mixtools)
    #     mixMdl = normalmixEM(outlier_zscore_s, k=4)
    #     plot(mixMdl,which=2)
    #     library(fpc)
    #     pamk(outlier_zscore_s, krange=2:15)
    kde=find_kde_criticals(outlier_zscore_s)
    maxs=kde$maxs
    mins=kde$mins
    adj_max=kde$adj_max
    adj_min=kde$adj_min
    
    
    # make sure the adj_min threshold is greater than 95% quantile
    per95 = quantile(outlier_zscore_s, na.rm=T, probs=0.95)
    if (!is.na(adj_min)){ if (adj_min < per95){adj_min=mins[mins>per95][1]}}
    
    # make a high threshold if adjacent min is missing
    if (is.na(adj_min)){adj_min = kde$dmode + 3*mad(outlier_zscore_s, na.rm=T)
                        mins = c(mins,adj_min)
    }
    
    ## plot ## 
    if (plot){
      df = as.data.frame(outlier_zscore_s)
      fn = paste("figures/", date, "/KDE/",date, "_", sample, "_", name, '_KDE_mzscore_distribution.pdf', sep = "")
      p = ggplot(data=df, aes_string(x=sample)) + geom_density()
      p = p + geom_vline(xintercept=maxs, colour="red", alpha=0.2) + geom_vline(xintercept=mins, colour="blue", alpha=0.2) + geom_vline(xintercept=adj_min, colour="blue")
      p = p + xlab(paste("modified Z-score in", name, sep=" ")) + 
        scale_x_continuous(limits=c(-5,5),breaks=c(mins,maxs)) + theme_bw() + 
        theme(text = element_text(colour="black", size=16), axis.text.x = element_text(angle = 90,colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
      p
      ggsave(file=fn)
    }
    outlier[,i] = outlier_zscore[,i] >= adj_min
    print(paste(sample, name, "KDE_outlier_threshold:", adj_min, "Number_of_outliers", sum(outlier[,i], na.rm=T), sep = "  "))
  }
  #return(list("outlier_zscore"=outlier_zscore, "outlier"=outlier))
  ### get the returned list from function find_outlier() and find top 5% z-score and genes within each sample
  
  # set up return matrixes
  zscore=outlier_zscore
  outlier=outlier
  num_5percent = dim(zscore)[1]
  top_outlier_zscore = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
  top_outlier = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
  top_outlier_boolean = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
  top_outlier_raw = matrix(,nrow=dim(zscore)[2],ncol=num_5percent)
  row.names(top_outlier_zscore)=colnames(zscore)
  row.names(top_outlier)=colnames(zscore)
  row.names(top_outlier_boolean)=colnames(zscore)
  row.names(top_outlier_raw)=colnames(zscore)
  colnames(top_outlier_zscore)=c(1:num_5percent)
  for (i in 1:num_5percent){colnames(top_outlier_zscore)[i] = paste(name, "top", colnames(top_outlier_zscore)[i], sep=" ")}
  colnames(top_outlier)=colnames(top_outlier_zscore)
  colnames(top_outlier_boolean)=colnames(top_outlier_zscore)
  colnames(top_outlier_raw)= colnames(top_outlier_zscore)
  # find top 5
  for (i in colnames(zscore)){
    whim=zscore[,i]
    a = whim[order(whim, decreasing=TRUE)][1:num_5percent]
    top_outlier_zscore[i,] = a
    whim2 = outlier[,i]
    top_outlier_boolean[i,] = whim2[order(whim, decreasing=TRUE)][1:num_5percent]
    whim3 = m2_genes[,i]
    top_outlier_raw[i,] = whim3[order(whim, decreasing=TRUE)][1:num_5percent]
    top_outlier[i,] = names(a)
  }
  
  a=rbind(top_outlier, top_outlier_zscore)
  a = a[order(row.names(a)),]
  fn = paste(pd,name,'outlier_m-zscore.txt', sep="_")
  write.table(a, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  
  b=rbind(top_outlier, top_outlier_raw)
  b = b[order(row.names(b)),]
  fn = paste(pd,name,'outlier_raw_exp.txt', sep="_")
  write.table(b, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  if (plot){
    # plot
    top_outlier.m <- melt(top_outlier[,c(1:5)])
    top_outlier_zscore.m <- melt(top_outlier_zscore[,c(1:5)])
    top_outlier_boolean.m <- melt(top_outlier_boolean[,c(1:5)])
    
    fn = paste(pd, name, 'top5_mzscore.pdf',sep ="_")
    YlGnBu = brewer.pal(9, "YlGnBu") 
    getPalette = colorRampPalette(YlGnBu)
    outlier.colors=c("NA", "#000000")
    
    p = ggplot()
    p = p + geom_tile(data=top_outlier_zscore.m, aes(x=as.factor(X1), y=X2, fill=value), linetype="blank") + scale_fill_gradientn(name= "modified z-score", colours=getPalette(100))
    p = p + geom_tile(data=top_outlier_boolean.m, aes(x=as.factor(X1), y=X2, color=value), fill=NA, size=0.5) + scale_colour_manual(name="KDE outlier",values = outlier.colors)
    p = p + geom_text(data=top_outlier.m,aes(x=as.factor(X1), y=X2, label = value), color="red", size=4, angle=90)
    p = p + xlab("sample") + ylab("top outlier targets") + theme_bw() + 
      theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=16), axis.text.y = element_text(colour="black", size=16))
    p
    ggsave(file=fn, height=h, width=w)
  }
  # return the top outliers
  return(list("outlier_zscore"=outlier_zscore, "outlier"=outlier,
              "top_outlier_zscore"=top_outlier_zscore, "top_outlier"=top_outlier, "top_outlier_boolean"=top_outlier_boolean))
}

##### OUTLIER IN ITRAQ PROTEOME #####
ITRAQ = read.table(row.names = 1, file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt",header=TRUE, sep="\t")

ITRAQ_druggable = find_outlier(ITRAQ, "ITRAQ druggable proteome", genes = druggable)

##### OUTLIER IN LFQ#####
LFQ=read.table(row.names = 1 , file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified_collapsed',header=TRUE, sep="\t", fill=T)

LFQ_druggable = find_outlier(LFQ, "LFQ druggable proteome", genes = druggable)

##### OUTLIER IN ITRAQ PHOSPHO ##### 
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
# 56651 phosphosites
ITRAQpho=ITRAQpho[,-c(1,2)]
print("Original number of ITRAQ phosphosites: 56651")
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho = ITRAQpho[,-c(17,18,20)]
ITRAQpho.na = ITRAQpho[rowSums(is.na(ITRAQpho)) <= 14,] #35838 phosphosites
print("NA filtered ITRAQ phosphosites: 35838")
genes = sub("-NP.*", "", row.names(ITRAQpho.na))
row.names(ITRAQpho.na) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\. _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\._.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub(" _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("_.*","",row.names(ITRAQpho.na)), unique=T)
ITRAQpho.na.d = ITRAQpho.na[genes %in% t(drugList), ] # 523 phosphosites
print("Druggable list filtered ITRAQ phosphosites: 523")

# outlier analysis
ITRAQpho_druggable = find_outlier(ITRAQpho.na.d, "ITRAQ druggable phosphoproteome", filter=FALSE, h=12)


##### OUTLIER IN LFQ phospho #####
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
LFQpho$sites=sub(".*\\(","",LFQpho$phospho_site)
LFQpho$sites=sub("\\)","",LFQpho$site)
LFQpho$sites=paste(LFQpho$gene_name,LFQpho$sites, sep=".")
row.names(LFQpho) = make.names(LFQpho$sites, unique=T)
colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
print("Original number of LFQ phosphosites: 18229")
# 18229 phosphosites
LFQpho = LFQpho[,-c(19:25)]
LFQpho.na = LFQpho[rowSums(is.na(LFQpho)) <= 8,] #3763 phosphosites! This is a great loss.
print("NA filtered LFQ phosphosites: 3763")
genes = sub("\\..*", "", row.names(LFQpho.na))
LFQpho.na.d = LFQpho.na[genes %in% t(drugList), ] #42 phosphosites
print("Druggable list filtered LFQ phosphosites: 42")
LFQpho_druggable = find_outlier(LFQpho.na.d, "LFQ druggable phosphoproteome", filter=FALSE, h=12)

##### RNA-Seq #####
RNA = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified',header=TRUE, sep="\t")
row.names(RNA)=make.names(RNA$gene, unique=T)
RNA=RNA[,-1] #16209 genes

RNA_druggable = find_outlier(RNA, "RNA druggable proteome", genes = druggable, w=12)
RNA_kinome = find_outlier(RNA, "RNA kinome", genes = kinome, w=12)

##### CNV #####
CNV = read.table(row.names = 1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified_normalized.tsv',header=TRUE, sep="\t")
# find outliers
CNV_druggable = find_outlier(CNV.n.m, "CNV druggable proteome", genes = druggable, w=12)
CNV_kinome = find_outlier(CNV.n.m, "CNV kinome", genes = kinome, w=12)


##### BRCA77 human data ##### 
###proteome###
BRCA77 = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 12553 isoforms
row.names(BRCA77) = BRCA77$Description
BRCA77.c = collapseRows(BRCA77[,-c(1,2,3)], rowGroup=BRCA77$Gene, rowID=BRCA77$Description)$datETcollapsed
# 11349 genes
#write.table(BRCA77.c,file = "/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt", quote = F, col.names=NA, sep ="\t")

BRCA77_druggable = find_outlier(BRCA77.c, "BRCA77 druggable proteome", genes = druggable, h=6, w=24)

###phosphoproteome###
BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 33239 phosphosites
row.names(BRCA77pho) = BRCA77pho$Gene.site
BRCA77pho = BRCA77pho[,-c(1,2)]
BRCA77pho.na = BRCA77pho[rowSums(is.na(BRCA77pho)) <= 67,] #33239 phosphosites
print("Original number of BRCA77 ITRAQ phosphosites: 33239")
print("NA filtered BRCA77 ITRAQ phosphosites: 33239")
genes = sub("-NP.*", "", row.names(BRCA77pho.na))
BRCA77pho.na.d = BRCA77pho.na[genes %in% t(drugList), ] #482 phosphosites
row.names(BRCA77pho.na.d) = sub("-NP_\\d+_"," ",row.names(BRCA77pho.na.d))
row.names(BRCA77pho.na.d) = make.names(sub(" _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("_.*","",row.names(BRCA77pho.na.d)), unique=T)
print("Druggable list filtered BRCA77 ITRAQ phosphosites: 482")

# find outliers
BRCA77pho_druggable = find_outlier(BRCA77pho.na.d, "BRCA77 druggable phosphoproteome", filter=F, h=10, w=24)

##### output the outlier table ##### 
