# diffExp.R by Kuan Huang @ WashU 201507
# find differentially expressed proteins

# mis
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression")
source("/Users/khuang/bin/LIB_exp.R")
load("/Users/khuang/bin/2015-08-01_Gene_Set.RData")

## diff_exp: implement t-test, wilcoxon rank sum test, to define differential expressed genes in two vectors
diff_exp = function(m, g1, g2, m.n="matrix", g1.n="group1", g2.n="group2", plot=FALSE, pathwayC=FALSE, pathwayT=FALSE){
  
  if (g1.n=="group1"){ g1.n = deparse(substitute(g1)) }
  if (g2.n=="group2"){ g2.n = deparse(substitute(g2)) }
  if (m.n=="matrix"){  m.n = deparse(substitute(m)) }
  cat("##### DIFFERENTIAL EXPRESSION ANALYSIS #####\n")
  cat(paste("Date: ", date, " By KH @ MGI\n", sep=""))
  cat(paste("Looking for differentially-expressed genes in", m.n, "; between", g1.n, "and", g2.n, "\n", sep=" "))
  
  m = as.matrix(m)
  # RUN limma
  g1g2 = c(g1,g2)
  m_g1g2 = m[,colnames(m) %in% g1g2]
  g2.b = as.numeric(colnames(m_g1g2) %in% g2)
  m_g1g2 = m_g1g2[rowSums(!is.na(m_g1g2[,g1]))>0 & rowSums(!is.na(m_g1g2[,g2]))>0,]
  n_genes = dim(m_g1g2)[1]
  fit = lmFit(m_g1g2, design=g2.b)
  fit = eBayes(fit)
  limma_t = topTable(fit, adjust.method="BH", number = n_genes, p.value=1.1)
  
  # pathway enrichment analysis based on limma result (count of sig genes)
  if (pathwayC){ 
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
    cat("GO differential expression results based on counts for ", m.n, g1.n,"vs",g2.n," printed to", tn, "\n")
    write.table(go.fisher, file=tn, quote=F, sep = '\t', col.names=NA)
    KEGG.fisher = kegga(fit2)
    KEGG.fisher$FDR.Up=p.adjust(KEGG.fisher$P.Up, method="BH")
    KEGG.fisher$FDR.Down=p.adjust(KEGG.fisher$P.Down, method="BH")
    KEGG.fisher = KEGG.fisher[order(KEGG.fisher$FDR.Down),]
    tn = paste(pd,m.n,g1.n,"vs",g2.n,"diff_exp_KEGG.txt", sep="_")
    write.table(KEGG.fisher, file=tn, quote=F, sep = '\t', col.names=NA)
    cat("KEGG differential expression results based on counts for ", m.n, g1.n,"vs",g2.n," printed to", tn, "\n")
  }
  
  ### RUN t-test and Wilcoxon test
  # set-up
  x = as.matrix(m_g1g2[,colnames(m_g1g2) %in% g1]) 
  y = as.matrix(m_g1g2[,colnames(m_g1g2) %in% g2])
  stats = matrix(,nrow=nrow(x),ncol=6)
  row.names(stats)=row.names(x)
  colnames(stats) = c(paste(m.n,g1.n,"mean",sep="_"), paste(m.n,g2.n,"mean",sep="_"), "t_test_Tstat", "t_test_p", "w_test_Wstat", "w_test_p")
  # TEST
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
  #stats=stats[order(stats$adj.P.Val, stats$t_test_fdr, decreasing=FALSE),] # limma top hits are segregated clustering worse than t-test... 
  #colnames(stats)[1] = "gene"
  tn = paste(pd,m.n,g1.n,"vs",g2.n,"diff_exp.txt", sep="_")
  write.table(stats, file=tn, quote=F, sep = '\t', col.names=NA)
  cat("KEGG differential expression results at gene level for ", m.n, g1.n,"vs",g2.n," printed to", tn, "\n")
  
  
  #sig.genes=row.names(stats[stats[,"t_test_fdr"]<0.01 & !is.na(stats[,"t_test_fdr"]),])
  sig.genes=row.names(stats[1:30,])
  m_g1g2.s=as.matrix(m_g1g2[row.names(m_g1g2) %in% sig.genes,])
  m_g1g2.s = m_g1g2.s[rowSums(is.na(m_g1g2.s))<=10,]
    
  if (plot){
    pdf(paste(pd,m.n,g1.n,"vs",g2.n,"top30_fdr.pdf", sep="_"))
    par(oma=c(1,2,5,3))
    mc=colnames(m_g1g2) %in% g1
    mc[mc]="forestgreen"
    mc[mc=="FALSE"]="orange"
    
    mhm = heatmap.2(m_g1g2.s, trace="none",na.color="white", notecol="black",
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
  
  # pathway analysis based on T statistics
  if (pathwayT){
    sample = as.data.frame(stats[,"t_test_Tstat",drop=F])
    sample$t_test_Tstat = as.numeric(levels(sample$t_test_Tstat))[sample$t_test_Tstat]
    sample = sample[!is.na(sample),,drop=F]
    
    stats = matrix(,,ncol=10)
    # Tstat positive means greater expression in mutant; negative means vise versa
    colnames(stats) = c("hsa","Pathway name","Num_genes","Num_pathway_genes","All_gene_mean", "Pathway_mean", "T_P", "GeneSet_P", "K-S_P", "Wilcox_P")
    for (pathway in names(KEGG)){
      hsa = strsplit(pathway, split = "\t")[[1]][1]
      pathway_name=strsplit(pathway, split = "\t")[[1]][2]
      pathway_genes = KEGG[[pathway]]
      inSet = rownames(sample) %in% pathway_genes
      numAllGene = length(inSet)
      numPathwayGene = table(inSet)[2]
      if (is.na(numPathwayGene) || numPathwayGene < 3){next}
      geneSetP = geneSetTest(inSet, sample$t_test_Tstat)[1]
      setExp = sample[rownames(sample) %in% pathway_genes,]
      ksP = ks.test(x=setExp, y=sample$t_test_Tstat)$p
      WilcoxP = wilcox.test(setExp,sample$t_test_Tstat)$p.value
      T_p = t.test(x=setExp, y=sample$t_test_Tstat)$p.value
      allM = mean(sample$t_test_Tstat)
      setM = mean(setExp)
      a = c(hsa, pathway_name, numAllGene,numPathwayGene,allM, setM, T_p, geneSetP, ksP, WilcoxP)
      stats=rbind(stats,a)
      row.names(stats)=NULL
    }
    Wilcox_fdr=p.adjust(stats[,"Wilcox_P"], method="BH")
    stats=as.data.frame(cbind(stats, Wilcox_fdr))
    stats=stats[order(as.numeric(stats$Wilcox_fdr), stats$Wilcox_P, decreasing=FALSE),]
    tn = paste(pd,m.n,g1.n,"vs",g2.n,"KEGGpathway_diff_exp.txt", sep="_")
    write.table(stats, file=tn, quote=F, sep = '\t', row.names=F)
    cat("KEGG differential expression results for ", m.n, g1.n,"vs",g2.n," printed to", tn, "\n")
  }
  
  return(list("stats"=stats, "sig.genes"=sig.genes)) #, "GO"=go.fisher, "KEGG"=KEGG.fisher))
} 
# proteome and phosphoproteome files 
ITRAQ_proteome = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt',header=TRUE, sep="\t", row.names=1)
LFQ_proteome =read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned.txt_hugoified_collapsed',header=TRUE, sep="\t", row.names=1)
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
LFQpho=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/label_free_all_WHIM_phospho_LFQ/all_WHIM_phospho_LFQ_minimum1_norm_nameadded_human_cleaned.txt',header=TRUE, sep="\t", fill=T)
ITRAQ_human_proteome = read.table(file="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt", header =T, sep="\t")

mut_matrix = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression/20150806_somatic_sample_gene_matrix.txt', header=T, sep='\t', row.names=1)

geneList = c("TP53", "PIK3CA", "ESR1")

##### ITRAQ proteome #####
mut = mut_matrix[rownames(mut_matrix) %in% colnames(ITRAQ_proteome),]
#for (gene in colnames(mut)){
for (gene in geneList){
  if (sum(is.na(mut[,gene])) + sum(mut[,gene]=="silent", na.rm=T) + sum(mut[,gene]=="rna", na.rm=T) < 19){
    wt = as.vector(c(rownames(mut[is.na(mut[,gene]),]),rownames(subset(mut, mut[,gene] == "silent")),rownames(subset(mut, mut[,gene] == "rna"))))
    mutant= as.vector(rownames(subset(mut, ! rownames(mut) %in% wt)))
    diff_x = diff_exp(m = ITRAQ_proteome, g1 = mutant, g2 = wt, g1.n = gene, g2.n = "wt",plot=T,pathwayC=T, pathwayT=TRUE)
  }
}

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
for (gene in geneList){
  if (sum(is.na(mut[,gene])) + sum(mut[,gene]=="silent", na.rm=T) + sum(mut[,gene]=="rna", na.rm=T) < 19){
    wt = as.vector(c(rownames(mut[is.na(mut[,gene]),]),rownames(subset(mut, mut[,gene] == "silent")),rownames(subset(mut, mut[,gene] == "rna"))))
    mutant= as.vector(rownames(subset(mut, ! rownames(mut) %in% wt)))
    diff = diff_exp(m = ITRAQ_phosphoproteome, g1 = mutant, g2 = wt, g1.n = gene, g2.n = "wt",plot=T,pathwayC=T, pathwayT=TRUE)
  }
}

### collapsed, gene level ###
ITRAQphosphoproteome = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
ITRAQ_phosphoproteome_gene = as.matrix(ITRAQphosphoproteome)
for (gene in geneList){
  if (sum(is.na(mut[,gene])) + sum(mut[,gene]=="silent", na.rm=T) + sum(mut[,gene]=="rna", na.rm=T) < 19){
    wt = as.vector(c(rownames(mut[is.na(mut[,gene]),]),rownames(subset(mut, mut[,gene] == "silent")),rownames(subset(mut, mut[,gene] == "rna"))))
    mutant= as.vector(rownames(subset(mut, ! rownames(mut) %in% wt)))
    diff = diff_exp(m = ITRAQ_phosphoproteome_gene, g1 = mutant, g2 = wt, g1.n = gene, g2.n = "wt",plot=T,pathwayC=T, pathwayT=TRUE)
  }
}

##### LFQ proteome #####
mut = mut_matrix[rownames(mut_matrix) %in% colnames(LFQ_proteome),]
#for (gene in colnames(mut)){
for (gene in geneList){
  if (sum(is.na(mut[,gene])) + sum(mut[,gene]=="silent", na.rm=T) + sum(mut[,gene]=="rna", na.rm=T) < 19){
    wt = as.vector(c(rownames(mut[is.na(mut[,gene]),]),rownames(subset(mut, mut[,gene] == "silent")),rownames(subset(mut, mut[,gene] == "rna"))))
    mutant= as.vector(rownames(subset(mut, ! rownames(mut) %in% wt)))
    diff = diff_exp(m = LFQ_proteome, g1 = mutant, g2 = wt, g1.n = gene, g2.n = "wt",plot=T,pathway=T)
  }
}

### merge and compare ### 
# LFQ_ITRAQ_merge=merge(basal_vs_lum,LFQ.basal_vs_lum, by = "row.names")
# colnames(LFQ_ITRAQ_merge)[1] = "gene"
# fn = paste(date,'ITRAQ_vs_LFQ_basal_lumB_diff_logP.pdf', sep="_")
# p = ggplot(data = LFQ_ITRAQ_merge, aes(x=-log10(as.numeric(as.character(t_test_p.x))) , y=-log10(as.numeric(as.character(t_test_p.y))), label=gene)) 
# p = p + geom_point(alpha=0.7) + xlab("ITRAQ basal vs. luminal -logP") + ylab("LFQ basal vs. luminal -logP") + theme_bw()
# p = p + xlim(0,9) + ylim(0,9)
# p = p + geom_text(aes(label=ifelse(-log10(as.numeric(as.character(t_test_p.x)))>5 | -log10(as.numeric(as.character(t_test_p.y)))>5,gene,"")),hjust=-0.05,just=0,size=8)
# p
# ggsave(file=fn, width=12, height=8)



