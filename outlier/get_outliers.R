# get_outliers.R
# Takes a gene (row) by sample (column) expression table; and a gene list
# Perform outlier anlaysis to find the top outliers in the data
# usage: Rscript get_outliers.R input_table input_geneList

### dependencies ###
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/outlier")
source("/Users/khuang/bin/LIB_exp.R")

cat("##### Running get_outliers.R by KH... ...\n")
cat(paste("Date: ", date, "\n", sep=""))
args=commandArgs(TRUE)
inputTable=args[1]
cat(paste("Input expression table: ", args[1], "\n"))
if(length(args)>=2){
  inputGenes=args[2]
  cat(paste("Input gene list: ", args[2], "\n"))
}
if(length(args)>=3){
  inputName=args[3]
  cat(paste("Input analysis name: ", args[3], "\n"))
}

# rewriting to look at sample-wise z-score distribution to define outliers
# function to find outliers using the modified zscores. Take a matrix as input
# requires function find_kde_inflection
# requires library reshape
# identify the adjacent local minimum as the threshold to define outliers
find_outlier = function(m, name, genes=NULL, minNum = 10, filter=TRUE, plotKDE=TRUE, plot=TRUE){ 
  # set-ups
  cat("##### OUTLIER ANALYSIS #####\n")
  num = nrow(m)
  if (filter){
    m2 = m[rowSums(!is.na(m)) >= minNum, ]
    m2_genes = m2[row.names(m2) %in% genes, ]
  } else{ m2 = m; m2_genes=m}
  num_NA = nrow(m2)
  num_gene = nrow(m2_genes)
  m2_genes = as.matrix(m2_genes)
  cat(paste("Looking for outliers in", deparse(substitute(genes)), "of", deparse(substitute(m)), "\n", sep=" "))
  cat(paste("Original number of genes:", num, "; NA filtered:", num_NA, "; Gene list filtered:", num_gene, "\n", sep=" "))
  kde_dir_cmd = paste("mkdir figures/", date, "/KDE/", sep="")
  system(kde_dir_cmd)
  
  # set up outlier matrices
  outlier = matrix(,nrow=dim(m2_genes)[1],ncol=dim(m2_genes)[2])
  row.names(outlier) = row.names(m2_genes)
  colnames(outlier) = colnames(m2_genes)
  outlier_zscore = outlier
  NumOutliers = c()
  
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
    
    kde=find_kde_criticals(outlier_zscore_s)
    maxs=kde$maxs
    mins=kde$mins
    adj_max=kde$adj_max
    adj_min=kde$adj_min
    
    # make sure the adj_min threshold is greater than 95% quantile
    per95 = quantile(outlier_zscore_s, na.rm=T, probs=0.95)
    if (!is.na(adj_min)){ if (adj_min < per95){adj_min=mins[mins>per95][1]}}
    
    if (!is.na(adj_min)){ if (adj_min < 1.5){adj_min=mins[mins>1.5][1]}}
    
    # make a high threshold if adjacent min is missing
    if (is.na(adj_min)){adj_min = kde$dmode + 3*mad(outlier_zscore_s, na.rm=T)
                        mins = c(mins,adj_min)
    }
    
    ## plot ## 
    if (plotKDE){
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
    cat(paste(sample, name, "KDE_outlier_threshold:", adj_min, "Number_of_outliers", sum(outlier[,i], na.rm=T), "\n", sep = "  "))
    NumOutliers = c(NumOutliers, sum(outlier[,i], na.rm=T))
  }
  
  # plot the count of outliers
  fn = paste(pd, name, 'outlier_count.pdf',sep ="_")
  pdf(fn)
  hist(x=NumOutliers, breaks=seq(min(NumOutliers)-0.5, max(NumOutliers)+0.5, by=1))
  dev.off()
  
  # set up return matrixes that are ranked by outlier value
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
  # find top outliers
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
  
  # print out the tables
  a=rbind(top_outlier, top_outlier_zscore)
  a = a[order(row.names(a)),]
  fn = paste(pd,name,'outlier_m-zscore.txt', sep="_")
  write.table(a, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  b=rbind(top_outlier, top_outlier_raw)
  b = b[order(row.names(b)),]
  fn = paste(pd,name,'outlier_raw_exp.txt', sep="_")
  write.table(b, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
  c=rbind(top_outlier, top_outlier_boolean)
  c = c[order(row.names(c)),]
  fn = paste(pd,name,'outlier_boolean.txt', sep="_")
  write.table(c, file=fn, quote=F, row.names=T, sep="\t", col.names=NA)
  
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
    ggsave(file=fn, height=6, width=12)
  }
  # return the top outliers
  return(list("outlier_zscore"=outlier_zscore, "outlier"=outlier, "NumOutliers" = NumOutliers,
              "top_outlier_zscore"=top_outlier_zscore, "top_outlier"=top_outlier, "top_outlier_boolean"=top_outlier_boolean))
}

###### main code #####

exp = read.table(file=inputTable, header=TRUE, sep="\t", row.names=1)
#row.names(exp) = make.names(exp[,1], unique =T)
#exp = exp[,-1]
geneList = read.table(file=inputGenes, header=FALSE, stringsAsFactors = F)
geneL = as.vector(t(geneList))

outliers = find_outlier(m=exp, name=inputName, genes = geneL, plotKDE=T)

# summarize outliers
top = outliers$top_outlier[,c(1:3)]
top.m = melt(top)
fn = paste(pd, inputName, "top3_outlier_summary.pdf", sep="_")
p = ggplot(data=top.m, aes(x=factor(value), fill=X2))
p = p + geom_bar() + scale_fill_discrete(guide_legend(title="outlier rank"))
p = p + labs( x = "gene", y = "count of samples", title = paste(inputName, "top outlier counts")) + theme_bw()
p = p + theme(text = element_text(colour="black", size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=14), axis.text.y = element_text(colour="black", size=14))
p
ggsave(file=fn, width = 15, height=6)
