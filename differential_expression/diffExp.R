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
system("mkdir figures")
date=Sys.time()
date = sub(" .*","",date)
date = paste(date, "KH", sep="_")
date = paste("figures/",date,sep="")

col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))
options(digits = 3)

RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(RdBu)

kinaseList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/2014-12-05_CPTAC_Kinase.MATRIX.v3b5_sheet1_genes.list', header=FALSE, stringsAsFactors = F)
drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/unionOfGeneVariantAADrug.tsv_hugoified_gene.list', header=FALSE, stringsAsFactors = F)

# functions
## to-do: try using samr and compare to t-test

## diff_exp: implement t-test, wilcoxon rank sum test, to define differential expressed genes in two vectors
diff_exp = function(m, g1, g2, plot=FALSE, pathway=FALSE){
  g1.n = deparse(substitute(g1))
  g2.n = deparse(substitute(g2))
  m.n = deparse(substitute(m))
  
  # for limma
  g1g2 = c(g1,g2)
  m_g1g2 = m[,colnames(m) %in% g1g2]
  g2.b = as.numeric(colnames(m_g1g2) %in% g2)
  m_g1g2 = m_g1g2[rowSums(!is.na(m_g1g2[,g1]))>0 & rowSums(!is.na(m_g1g2[,g2]))>0,]
  n_genes = dim(m_g1g2)[1]
  fit = lmFit(m_g1g2, design=g2.b)
  fit = eBayes(fit)
  limma_t = topTable(fit, adjust.method="BH", number = n_genes, p.value=1.1)
  
  if (pathway){
    # pathway enrichment analysis: fisher's test on GO and KEGG
    library(biomaRt)
    ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    mapTab = getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = rownames(m), mart = ensembl, uniqueRows=FALSE)
    dupRows = union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
    if (length(dupRows)==0){entrezIds = mapTab} else {entrezIds = mapTab[-dupRows, ]}
    row.names(entrezIds) = entrezIds$hgnc_symbol
    mer = merge(m, entrezIds, by="row.names")
    mer = mer[!is.na(mer$entrezgene),]
    row.names(mer) = mer$entrezgene
    mer = mer[,-1]
    mer = mer[,-ncol(mer)]
    mer = mer[,-ncol(mer)]
    g1g2 = c(g1,g2)
    mer_g1g2 = mer[,colnames(mer) %in% g1g2]
    g2.b = as.numeric(colnames(mer_g1g2) %in% g2)
    mer_g1g2 = mer_g1g2[rowSums(!is.na(mer_g1g2[,g1]))>0 & rowSums(!is.na(mer_g1g2[,g2]))>0,]
    n_genes = dim(mer_g1g2)[1]
    fit2 = lmFit(mer_g1g2, design=g2.b)
    fit2 = eBayes(fit2)
    go.fisher = goana(fit2)
    go.fisher$FDR.Up=p.adjust(go.fisher$P.Up, method="BH")
    go.fisher$FDR.Down=p.adjust(go.fisher$P.Down, method="BH")
    go.fisher = go.fisher[order(go.fisher$FDR.Down),]
    tn = paste(pd,m.n,g1.n,"vs",g2.n,"diff_exp_GO.txt", sep="_")
    write.table(go.fisher, file=tn, quote=F, sep = '\t', col.names=NA)
    KEGG.fisher = kegga(fit2)
    KEGG.fisher$FDR.Up=p.adjust(KEGG.fisher$P.Up, method="BH")
    KEGG.fisher$FDR.Down=p.adjust(KEGG.fisher$P.Down, method="BH")
    KEGG.fisher = KEGG.fisher[order(KEGG.fisher$FDR.Down),]
    tn = paste(pd,m.n,g1.n,"vs",g2.n,"diff_exp_KEGG.txt", sep="_")
    write.table(KEGG.fisher, file=tn, quote=F, sep = '\t', col.names=NA)
  }
  
  # for t-test and Wilcoxon test
  # print header
  x = m_g1g2[,colnames(m_g1g2) %in% g1] 
  y = m_g1g2[,colnames(m_g1g2) %in% g2]
  
  stats = matrix(,nrow=nrow(x),ncol=6)
  row.names(stats)=row.names(x)
  colnames(stats) = c(paste(m.n,g1.n,"mean",sep="_"), paste(m.n,g2.n,"mean",sep="_"), "t_test_Tstat", "t_test_p", "w_test_Wstat", "w_test_p")
  for (i in 1:nrow(x)){
    if (sum(!is.na(x[i,]))<2 | sum(!is.na(y[i,]))<2){
      stats[i,]=c(mean(x[i,], na.rm=T),mean(y[i,], na.rm=T),rep("NA",4))} 
    else{
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
  stats = merge(limma_t, stats, by = "row.names")
  row.names(stats) = stats$Row.names
  stats=stats[,-1]
  stats=stats[order(stats$t_test_fdr, stats$adj.P.Val, decreasing=FALSE),]
  #colnames(stats)[1] = "gene"
  tn = paste(pd,m.n,g1.n,"vs",g2.n,"diff_exp.txt", sep="_")
  write.table(stats, file=tn, quote=F, sep = '\t', col.names=NA)
  
  #sig.genes=row.names(stats[stats[,"t_test_fdr"]<0.01 & !is.na(stats[,"t_test_fdr"]),])
  sig.genes=row.names(stats[1:30,])
  # significant genes in g1 vs. g2 g1+g2 samples
  gs = c(g1,g2)
  m.gs=m[row.names(m) %in% sig.genes,colnames(m) %in% gs]
  
  if (plot){
    pdf(paste(pd,m.n,g1.n,"vs",g2.n,"top30_fdr.pdf", sep="_"))
    par(oma=c(1,2,5,3))
    mc=colnames(m.gs) %in% g1
    mc[mc]="forestgreen"
    mc[mc=="FALSE"]="orange"
  
    mhm = heatmap.2(m.gs, trace="none",na.color="white", notecol="black",
                  ColSideColors = mc,
                  cexRow=1.1,cexCol=1.4, col=getPalette1) #
    par(lend = 1)  
    legend("topright",    # location of the legend on the heatmap plot
         legend = c(g1.n, g2.n), # category labels
         col = c("forestgreen", "orange"),  # color key
         lty= 1,             # line style
         lwd = 10            # line width
    )
    #clus_order1 = mhm$colInd
    dev.off()
  }
  
  return(list("stats"=stats, "sig.genes"=sig.genes, "GO"=go.fisher, "KEGG"=KEGG.fisher))
} 

# proteome and phosphoproteome files 
ITRAQ = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2.txt',header=TRUE, sep="\t")
LFQ=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified',header=TRUE, sep="\t", fill=T)
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
RSEM=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct.txt',header=TRUE, sep="\t")

clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')

Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
bl = c(Basal,LumB)

### process and cluster ITRAQ proteome
row.names(ITRAQ) = ITRAQ$Description
colnames(ITRAQ) = sub("\\..*", "", colnames(ITRAQ))
ITRAQ.c = collapseRows(ITRAQ[,-c(1,2,3)], rowGroup=ITRAQ$Gene, rowID=ITRAQ$Description)$datETcollapsed
# 12698 genes
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQ.c_1 = ITRAQ.c[,-c(17,18,20)]
rm(ITRAQ.c)
ITRAQ_proteome=as.matrix(ITRAQ.c_1)

ITRAQpho_Basal_vs_LumB = diff_exp(ITRAQ_proteome, Basal, LumB, plot=F)

ITRAQ.basal = ITRAQ.c_1[,colnames(ITRAQ.c_1) %in% basal] 
ITRAQ.LumB = ITRAQ.c_1[,colnames(ITRAQ.c_1) %in% LumB]
basal_vs_lum = diff_exp(ITRAQ.basal,ITRAQ.LumB)
t.test.sig.genes=row.names(basal_vs_lum[basal_vs_lum[,"t_test_fdr"]<0.005 & !is.na(basal_vs_lum[,"t_test_fdr"]),])
# significant genes in basal and luminal b in basal/luminbal B samples
ITRAQ.bl=ITRAQ.c_1[row.names(ITRAQ.c_1) %in% t.test.sig.genes,colnames(ITRAQ.c_1) %in% bl]

pdf(paste(date,'ITRAQ_proteome_basal_vs_luminal_fdr0.005.pdf', sep="_"))
par(oma=c(1,2,5,3))
basalc=colnames(ITRAQ.bl) %in% basal
basalc[basalc]="forestgreen"
basalc[basalc=="FALSE"]="orange"

ITRAQhm = heatmap.2(ITRAQ.bl, trace="none",na.color="white", notecol="black",
                    ColSideColors = basalc,
                    cexRow=1.1,cexCol=1.4, col=getPalette) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B"), # category labels
       col = c("forestgreen", "orange"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
#clus_order1 = ITRAQhm$colInd
dev.off()

### process and cluster LFQ proteome
row.names(LFQ) = LFQ$RefSeq_id
LFQ.c = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$datETcollapsed #5802 genes
#LFQ.NP = collapseRows(LFQ[,-c(21:26)], rowGroup=LFQ$gene_name, rowID=LFQ$RefSeq_id)$group2row
rm(LFQ)
# get rid of WHIM16_2, WHIM2_2
LFQ.c_1 = LFQ.c[,-c(6,9)]
rm(LFQ.c)

LFQ.basal = LFQ.c_1[,colnames(LFQ.c_1) %in% basal] 
LFQ.LumB = LFQ.c_1[,colnames(LFQ.c_1) %in% LumB]
LFQ.basal_vs_lum = diff_exp(LFQ.basal,LFQ.LumB)
t.test.sig.genes=row.names(LFQ.basal_vs_lum[LFQ.basal_vs_lum[,"t_test_fdr"]<0.03 & !is.na(LFQ.basal_vs_lum[,"t_test_fdr"]),])
# significant genes in basal and luminal b in basal/luminbal B samples
LFQ.bl=LFQ.c_1[row.names(LFQ.c_1) %in% t.test.sig.genes, colnames(LFQ.c_1) %in% bl]

pdf(paste(date,'LFQ_proteome_basal_vs_luminal_fdr0.03.pdf', sep="_"))
par(oma=c(1,2,5,3))
basalc=colnames(LFQ.bl) %in% basal
basalc[basalc]="forestgreen"
basalc[basalc=="FALSE"]="orange"

LFQhm = heatmap.2(LFQ.bl, trace="none",na.color="white", notecol="black",
                    ColSideColors = basalc,
                    cexRow=1.1,cexCol=1.4, col=getPalette) #
par(lend = 1)  
legend("topright",    # location of the legend on the heatmap plot
       legend = c("Basal", "Luminal B"), # category labels
       col = c("forestgreen", "orange"),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
#clus_order1 = LFQhm$colInd
dev.off()

### merge and compare ### 
LFQ_ITRAQ_merge=merge(basal_vs_lum,LFQ.basal_vs_lum, by = "row.names")
colnames(LFQ_ITRAQ_merge)[1] = "gene"
fn = paste(date,'ITRAQ_vs_LFQ_basal_lumB_diff_logP.pdf', sep="_")
p = ggplot(data = LFQ_ITRAQ_merge, aes(x=-log10(as.numeric(as.character(t_test_p.x))) , y=-log10(as.numeric(as.character(t_test_p.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("ITRAQ basal vs. luminal -logP") + ylab("LFQ basal vs. luminal -logP") + theme_bw()
p = p + xlim(0,9) + ylim(0,9)
p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(t_test_p.x)))>5 | -log10(as.numeric(as.character(t_test_p.y)))>5,gene,"")),hjust=-0.05,just=0,size=3,alpha=0.8)
p
ggsave(file=fn)

### TO-DO: compare to human, and then cluster with human###

### process and cluster ITRAQ phosphoproteome
row.names(ITRAQpho) = ITRAQpho$Gene
colnames(ITRAQpho) = sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho_m = data.matrix(ITRAQpho[,-c(1,2)])
SD=rowSds(ITRAQpho_m, na.rm=TRUE)
ITRAQpho_m2 = ITRAQpho_m[SD>2,]
ITRAQpho_m3 = ITRAQpho_m2[rowSums(is.na(ITRAQpho_m2)) <= 10,]



### process and cluster LFQ phosphoproteome
row.names(LFQpho) = make.names(LFQpho$phospho_site, unique=T)
LFQpho = LFQpho[,-c(19:24)]
colnames(LFQpho) = sub("_P", "", colnames(LFQpho))
colnames(LFQpho) = sub("Intensity.W", "WHIM", colnames(LFQpho))
LFQpho_m = data.matrix(LFQpho)
SD=rowSds(LFQpho_m, na.rm=TRUE)
LFQpho_m2 = LFQpho_m[SD>2,]
LFQpho_m3 = LFQpho_m2[rowSums(is.na(LFQpho_m2)+is.nan(LFQpho_m2)) <= 10,]



### process and cluster RSEM data
row.names(RSEM) = RSEM$sample
RSEM2 = RSEM[,colnames(RSEM) %in% colnames(ITRAQ)] 
RSEM_m = data.matrix(RSEM2)
#RSEM_m = log2(RSEM_m)
SD=rowSds(RSEM_m, na.rm=TRUE)
RSEM_m2 = RSEM_m[SD>2,]
RSEM_m3 = RSEM_m2[rowSums(is.na(RSEM_m2)+is.nan(RSEM_m2)) <= 10,]



