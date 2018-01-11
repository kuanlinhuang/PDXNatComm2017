# pathway_activation.R
# find activated KEGG pathway through z-score of proteome and phosphoproteomes
# plot the protein and phosphoprotein data to the KEGG pathway

##### pathway_activation.R #####
# Kuan-lin Huang @ WashU 2015 Oct
# plot the protein and phosphoprotein data to the pi3k pathway

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/pi3k_pathway")
source("/Users/khuang/bin/LIB_exp.R")
#library("reactome.db")
library(KEGGprofile)
library(biomaRt)

# Get list KEGG and REACT for analysis
load("/Users/khuang/bin/2015-08-01_Gene_Set.RData")

KEGG_signaling = KEGG[grep("signaling", names(KEGG))]

##### Can also get other pathways from Broad's GSEA website #####
# test for gene set activation using z-score
KEGGpathway_activation = function (m){
  m.n = deparse(substitute(m))
  # print header
  cat("##### KEGG pathway differential regulation analysis #####\n")
  cat("Detecting differentially expressed KEGG pathways in", m.n, "\n")
  
  for (s in colnames(m)){
    sample = m[,s,drop=F]
    sample = sample[!is.na(sample),,drop=F]
    stats = matrix(,,ncol=10)
    colnames(stats) = c("hsa","Pathway name","Num_genes","Num_pathway_genes","All_gene_mean", "Pathway_mean", "T_P", "GeneSet_P", "K-S_P", "Wilcox_P")
    for (pathway in names(KEGG_signaling)){
      hsa = strsplit(pathway, split = "\t")[[1]][1]
      pathway_name=strsplit(pathway, split = "\t")[[1]][2]
      pathway_genes = KEGG[[pathway]]
      inSet = rownames(sample) %in% pathway_genes
      numAllGene = length(inSet)
      numPathwayGene = table(inSet)[2]
      if (is.na(numPathwayGene) || numPathwayGene < 3){next}
      geneSetP = geneSetTest(inSet, sample)[1]
      setExp = sample[rownames(sample) %in% pathway_genes,]
      ksP = ks.test(x=setExp, y=sample)$p
      WilcoxP = wilcox.test(setExp,sample)$p.value
      T_p = t.test(x=setExp, y=sample)$p.value
      allM = mean(sample)
      setM = mean(setExp)
      a = c(hsa, pathway_name, numAllGene,numPathwayGene,allM, setM, T_p, geneSetP, ksP, WilcoxP)
      stats=rbind(stats,a)
      row.names(stats)=NULL
    }
    Wilcox_fdr=p.adjust(stats[,"Wilcox_P"], method="BH")
    stats=as.data.frame(cbind(stats, Wilcox_fdr))
    stats=stats[order(as.numeric(stats$Wilcox_fdr), stats$Wilcox_P, decreasing=FALSE),]
    tn = paste(pd,m.n,s,"KEGGSig_pathway_activation.txt", sep="_")
    write.table(stats, file=tn, quote=F, sep = '\t', row.names=F)
    cat("Results for", s, "printed to", tn, "\n")
  }
}

# less entries (roughly 200, whereas the web API = 287)
# KEGGpathway_activation2 = function (m){
#   m.n = deparse(substitute(m))
#   # print header
#   cat("##### KEGG pathway differential regulation analysis #####\n")
#   cat("Detecting differentially expressed KEGG pathways in", m.n, "\n")
#   
#   for (s in colnames(m)){
#     sample = m[,s,drop=F]
#     sample = sample[!is.na(sample),,drop=F]
#     stats = matrix(,,ncol=10)
#     colnames(stats) = c("hsa","Pathway name","Num_genes","Num_pathway_genes","All_gene_mean", "Pathway_mean", "T_P", "GeneSet_P", "K-S_P", "Wilcox_P")
#     for (pathway in names(KEGG_path)){
#       pathname = KEGG_path[[pathway]]
#       h_path = paste("hsa",pathway,sep="")
#       pathgeneID = KEGG_path2gene[[h_path]]
#       if (length(pathgeneID)==0){next}
#       mapTab = getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "entrezgene", values = pathgeneID, 
#                      mart = ensembl, uniqueRows=FALSE)
#       dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
#       pathway_genes = mapTab[-dupRows, 2]
#       inSet = rownames(sample) %in% pathway_genes
#       numAllGene = length(inSet)
#       numPathwayGene = table(inSet)[2]
#       if (is.na(numPathwayGene) || numPathwayGene < 3){next}
#       geneSetP = geneSetTest(inSet, sample)[1]
#       setExp = sample[rownames(sample) %in% pathway_genes,]
#       ksP = ks.test(x=setExp, y=sample)$p
#       WilcoxP = wilcox.test(setExp,sample)$p.value
#       T_p = t.test(x=setExp, y=sample)$p.value
#       allM = mean(sample)
#       setM = mean(setExp)
#       a = c(h_path, pathname, numAllGene,numPathwayGene,allM, setM, T_p, geneSetP, ksP, WilcoxP)
#       stats=rbind(stats,a)
#       row.names(stats)=NULL
#       
#       sample2 = as.data.frame(cbind(sample,inSet))
#       colnames(sample2)[1]="S"
#       #       ### plot dens
#       #       p = ggplot(data=sample2, aes(x=factor(inSet), y =S))
#       #       p = p + geom_violin(aes(fill=factor(inSet)), alpha=0.5) + geom_jitter(height=0, alpha=0.2)
#       #       p = p + theme_bw() +ylab("normalized protein expression ratio")
#       #       p
#     }
#     Wilcox_fdr=p.adjust(stats[,"Wilcox_P"], method="BH")
#     stats=as.data.frame(cbind(stats, Wilcox_fdr))
#     stats = stats[order(as.numeric(stats$Wilcox_fdr),stats$Wilcox_P, decreasing=FALSE),]
#     tn = paste(pd,m.n,s,"KEGGpathway_activation2.txt", sep="_")
#     write.table(stats, file=tn, quote=F, sep = '\t', row.names=F)
#     cat("Results for", s, "printed to", tn, "\n")
#   }
# }

Reactome_pathway_activation = function (m){
  m.n = deparse(substitute(m))
  # print header
  cat("##### REACTOME pathway differential regulation analysis #####\n")
  cat("Detecting differentially expressed Reactome pathways in", m.n, "\n")
  
  for (s in colnames(m)){
    sample = m[,s,drop=F]
    sample = sample[!is.na(sample),,drop=F]
    stats = matrix(,,ncol=10)
    colnames(stats) = c("reactID","Pathway name","Num_genes","Num_pathway_genes","All_gene_mean", "Pathway_mean", "T_P", "GeneSet_P", "K-S_P", "Wilcox_P")
    for (pathway in names(REACT)){
      reactID = strsplit(pathway, split = "\t")[[1]][1]
      pathway_name=strsplit(pathway, split = "\t")[[1]][2]
      pathway_genes = REACT[[pathway]]
      inSet = rownames(sample) %in% pathway_genes
      numAllGene = length(inSet)
      numPathwayGene = table(inSet)[2]
      if (is.na(numPathwayGene) || numPathwayGene < 3){next}
      geneSetP = geneSetTest(inSet, sample)[1]
      setExp = sample[rownames(sample) %in% pathway_genes,]
      ksP = ks.test(x=setExp, y=sample)$p
      WilcoxP = wilcox.test(setExp,sample)$p.value
      T_p = t.test(x=setExp, y=sample)$p.value
      allM = mean(sample)
      setM = mean(setExp)
      a = c(reactID, pathway_name, numAllGene,numPathwayGene,allM, setM, T_p, geneSetP, ksP, WilcoxP)
      stats=rbind(stats,a)
      row.names(stats)=NULL
    }
    Wilcox_fdr=p.adjust(stats[,"Wilcox_P"], method="BH")
    stats=as.data.frame(cbind(stats, Wilcox_fdr))
    stats=stats[order(as.numeric(stats$Wilcox_fdr), stats$Wilcox_P, decreasing=FALSE),]
    tn = paste(pd,m.n,s,"REACTpathway_activation.txt", sep="_")
    write.table(stats, file=tn, quote=F, sep = '\t', row.names=F)
    cat("Results for", s, "printed to", tn, "\n")
  }
}


getKeggPathwayGenes<-function(hsaID){
  kegg_pathway=readLines(paste("http://rest.kegg.jp/get/",hsaID,sep=""))
  # KGML: showing the relations: http://rest.kegg.jp/get/hsa04012/kgml 
  kegg_pathway_genes=kegg_pathway[grep(";",kegg_pathway)]
  kegg_pathway_genes=sub(";.*", "",kegg_pathway_genes )
  kegg_pathway_genes=kegg_pathway_genes[-grep("CLASS",kegg_pathway_genes)]
  kegg_pathway_genes[grep("GENE", kegg_pathway_genes)]=gsub("GENE","",kegg_pathway_genes[grep("GENE", kegg_pathway_genes)])
  
#   genes =gsub("\\s","",gsub("\\s[0-9]+", "", perl = T, kegg_pathway_genes))
#   kegg_pathway_genes=paste("hsa",gsub("^\\s+","", perl=T, kegg_pathway_genes), sep="")
#   hsa = sub("\\s.*","",kegg_pathway_genes)
#   kegg_hsa_gene = cbind(hsa,genes)
#   return(kegg_hsa_gene)
  
  kegg_pathway_genes=gsub("\\s","",gsub("\\s[0-9]+", "", perl = T, kegg_pathway_genes))
  return(kegg_pathway_genes)
  
}

plotKEGG_pathway = function(exp,pathwayNum) {
  exp=exp[rowSums(is.na(exp))==0,,drop=FALSE]
  exp2=rep(0,nrow(exp))
  exp = cbind(exp,exp2)
  # convert hugo gene name rownames to hsaID; bug: this step requires two columns...
  exp = convertId(exp,filters="hgnc_symbol")  
  exp=exp[,1,drop=F]
  limit = max(-min(exp),max(exp))
  col = col_by_value(exp, col = RdBu1024, breaks=seq(-limit,limit,length.out=1025),showColorBar = T)
  
  temp = plot_pathway(exp, type = "bg", bg_col = col, text_col = "black",
                     magnify = 1.2, species = "hsa", database_dir = system.file("extdata", package = "KEGGprofile"),
                     pathway_id = pathwayNum)
}

plotKEGG_pathwayDev = function(exp,pathwayNum) {
  exp=exp[rowSums(is.na(exp))==0,,drop=FALSE]
  exp2=rep(0,nrow(exp))
  exp = cbind(exp,exp2)
  # convert hugo gene name rownames to hsaID; bug: this step requires two columns. That's why the exp2 column is added
  exp = convertId(exp,filters="hgnc_symbol")  
  exp=exp[,-ncol(exp),drop=F]
  limit = max(-min(exp),max(exp))
  col = col_by_value(exp, col = RdBu1024, breaks=seq(-limit,limit,length.out=1025),showColorBar = T)
  
  temp = plot_pathway(exp, type = "bg", bg_col = col, text_col = "black",
                      magnify = 1.2, species = "hsa", database_dir = system.file("extdata", package = "KEGGprofile"),
                      pathway_id = pathwayNum)
}

##### PLOTTING FUNCTION #####
### TODO: add proteomic (scale color?) 
# plot with color scale on all genes
plotWrap = function (sample, pathway){
  whim = ITRAQpho_outlier_zscore.c2[,sample,drop=F]
  plotKEGG_pathway(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}
# plot human data with color scale on genes only in the specific pathway
plotWrapH = function (sample, pathway){
  hsaID = paste("hsa",pathway,sep="")
  genes = getKeggPathwayGenes(hsaID)
  whim = BRCA77pho_outlier_zscore.c[,sample,drop=F]
  whim = whim[row.names(whim) %in% genes,,drop=F]
  plotKEGG_pathway(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}
# plot WHIM data with color scale on genes only in the specific pathway
plotWrapW = function (sample, pathway){
  hsaID = paste("hsa",pathway,sep="")
  genes = getKeggPathwayGenes(hsaID)
  whim = ITRAQpho_outlier_zscore.c[,sample,drop=F]
  whim = whim[row.names(whim) %in% genes,,drop=F]
  plotKEGG_pathway(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}
# plot WHIM data with color scale on genes only in the specific pathway; both proteome and phospho
# the resulting plot looks wrong, check or not use!
plotWrapWDev = function (sample, pathway){
  hsaID = paste("hsa",pathway,sep="")
  genes = getKeggPathwayGenes(hsaID)
  whim_pro = ITRAQ_druggable_pro_zscore[,sample,drop=F]
  whim_pho = ITRAQpho_outlier_zscore.c[,sample,drop=F]
  whim = merge(whim_pro, whim_pho, by = "row.names", all=T)
  row.names(whim) = whim[,1]
  whim = whim[,-1]
  whim = whim[row.names(whim) %in% genes,,drop=F]
  whim[is.na(whim)] = 0 # limitation: set NA to 0 otherwise they don't show up
  plotKEGG_pathwayDev(whim, pathway)
  command= paste("mv hsa",pathway,"_profile_bg.png ",pd,"_",sample,"_",pathway,"_bg.png", sep="")
  system(command)
}

##### PROTEOME DATA #####
ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 13719 isoforms
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$datETcollapsed # 12698 genes
# ITRAQ.NP = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$group2row
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ.c_1 = ITRAQ.c[,-c(17,18,20)]

# ITRAQ_outlier = find_outlier(ITRAQ.c_1, "ITRAQ proteome")
# ITRAQ_outlier_zscore = ITRAQ_outlier$outlier_zscore
#KEGGpathway_activation(ITRAQ_outlier_zscore)
#Reactome_pathway_activation(ITRAQ_outlier_zscore)

##### PHOHSPHO DATA ##### 
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
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
ITRAQpho.na = as.matrix(ITRAQpho.na)
# use the adjusted z-score to find outliers
ITRAQpho_outliers=find_outlier(ITRAQpho.na, "ITRAQ phosphoproteome")#, filter = F, plot=F)
ITRAQpho_outlier_zscore=ITRAQpho_outliers$outlier_zscore
genes = sub("\\..*", "", row.names(ITRAQpho_outlier_zscore))
ITRAQpho_outlier_zscore.c = collapseRows(ITRAQpho_outlier_zscore, rowGroup=genes, rowID=row.names(ITRAQpho_outlier_zscore))$datETcollapsed
KEGGpathway_activation(ITRAQpho_outlier_zscore.c)
#Reactome_pathway_activation(ITRAQpho_outlier_zscore.c)

##### plotting for WHIM #####
# for (whim in colnames(ITRAQpho_outlier_zscore.c)){
# #   plotWrapW(whim, "03460") # Fanconi anemia
# #   plotWrapW(whim, "04014") # Ras signaling
#   plotWrapW(whim, "04915") # ER signaling
# #   plotWrapW(whim, "04151") # PI3K-AKT signaling
# #   plotWrapW(whim, "04012") # ERBB2 signaling
# #   plotWrapW(whim, "04064") # NFkB signaling
# #   plotWrapW(whim, "04010") # MAPK signaling
# #   plotWrapW(whim, "04110") # cell cycle
# #   plotWrapW(whim, "04115") # p53 signaling
# #   plotWrapW(whim, "05200") # pathways in cancer
# }

##### BRCA77 phosphoproteome #####
BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
# collapse protein isoforms: use the maximum value if two; if more choose most representative (the representative row according to the least number of missing data, the highest sample mean, the highest sample variance, the highest connectivity)
# 33239 phosphosites
row.names(BRCA77pho) = BRCA77pho$Gene.site
BRCA77pho = BRCA77pho[,-c(1,2)]
BRCA77pho.na = BRCA77pho[rowSums(is.na(BRCA77pho)) <= 67,] #33239 phosphosites
rm(BRCA77pho)
row.names(BRCA77pho.na) = make.names(sub("-NP_\\d+_"," ",row.names(BRCA77pho.na)), unique=T)
row.names(BRCA77pho.na) = make.names(sub(" _.*","",row.names(BRCA77pho.na)), unique=T)
row.names(BRCA77pho.na) = make.names(sub("_.*","",row.names(BRCA77pho.na)), unique=T)
BRCA77pho.m = as.matrix(BRCA77pho.na)
BRCA77pho_outliers=find_outlier(BRCA77pho.m, "BRCA ITRAQ phosphoproteomes")
BRCA77pho_outlier_zscore=BRCA77pho_outliers$outlier_zscore
BRCA77pho_outlier=BRCA77pho_outliers$outlier
genes = sub("\\..*", "", row.names(BRCA77pho_outlier_zscore))
BRCA77pho_outlier_zscore.c = collapseRows(BRCA77pho_outlier_zscore, rowGroup=genes, rowID=row.names(BRCA77pho_outlier_zscore))$datETcollapsed
KEGGpathway_activation(BRCA77pho_outlier_zscore.c)
#Reactome_pathway_activation(BRCA77pho_outlier_zscore.c)

# for (h in colnames(BRCA77pho_outlier_zscore.c)){
#   plotWrapH(h, "04151") # PI3K-AKT signaling
#   plotWrapH(h, "04012") # ERBB2 signaling
#   plotWrapH(h, "04064") # NFkB signaling
#   plotWrapH(h, "04010") # MAPK signaling
#   plotWrapH(h, "04110") # cell cycle
#   plotWrapH(h, "04115") # p53 signaling
#   plotWrapH(h, "05200") # pathways in cancer
# }

# plot violin for top pathways?
#       ### plot dens
#       p = ggplot(data=sample2, aes(x=factor(inSet), y =S))
#       p = p + geom_violin(aes(fill=factor(inSet)), alpha=0.5) + geom_jitter(height=0, alpha=0.2)
#       p = p + theme_bw() +ylab("normalized protein expression ratio")
#       p





