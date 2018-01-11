### pdx_vs_human_pro.R
# compare outlier and protein expression status between PDX and human
# outlier analysis pipeline
# nice note for outlier analysis: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.html

### dependencies ###
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/outlier")
source("/Users/khuang/bin/LIB_exp.R")

drugList = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/reference_files/gene_drug_list/Premed_raw_databases/drugBank/all_target_ids_all.txt.human.tsv_hugoified.tsv_list.txt_list', header=FALSE, stringsAsFactors = F)
druggable = as.vector(t(drugList))

##### data #####
### proteome comparison
ITRAQ = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
ITRAQ.d = ITRAQ[row.names(ITRAQ) %in% druggable,]
ITRAQ.d.m = melt(as.matrix(ITRAQ.d))
colnames(ITRAQ.d.m) = c("Gene","Sample","Pro")
ITRAQ.d.m$Species = "PDX"

### phospho
ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
gene.site=ITRAQpho$gene.site
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho=ITRAQpho[,-c(1,2)]
# get rid of TaxIR, HumIR, WHIM13.1
ITRAQpho.na = ITRAQpho[,-c(17,18,20)]
genes = sub("-NP.*", "", row.names(ITRAQpho.na))
row.names(ITRAQpho.na) = make.names(sub("-NP_\\d+_"," ",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\. _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("\\._.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub(" _.*","",row.names(ITRAQpho.na)), unique=T)
row.names(ITRAQpho.na) = make.names(sub("_.*","",row.names(ITRAQpho.na)), unique=T)
ITRAQpho.na.d = ITRAQpho.na[genes %in% druggable, ] 
#look up table
genesite2rowname = data.frame(cbind(as.character(gene.site),row.names(ITRAQpho.na)))
colnames(genesite2rowname) = c("Full_phosphosite","Phosphosite")

ITRAQpho.na.d.m = melt(as.matrix(ITRAQpho.na.d))
colnames(ITRAQpho.na.d.m) = c("Site","Sample","Pho")
ITRAQpho.na.d.m$Species = "PDX"

##### BRCA77 human data ##### 
###proteome###
BRCA77 = read.table(row.names=1,header=TRUE, sep="\t", file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_proteome/201507/BRCA77_unimodal_proteome-ratio-norm_exp_collapsed.txt')
BRCA77.d = BRCA77[row.names(BRCA77) %in% druggable,]
# selected only the sample without duplicates and not normal
BRCA77.d.s =  BRCA77.d[,nchar(colnames(BRCA77.d)) == 7]

BRCA77.d.m = melt(as.matrix(BRCA77.d.s))
colnames(BRCA77.d.m) = c("Gene","Sample","Pro")
BRCA77.d.m$Species = "Human"
BRCA77.d.m$Sample = paste("TCGA",BRCA77.d.m$Sample,sep=".")
BRCA77.d.m$Sample = gsub("\\.","-",BRCA77.d.m$Sample)

###phosphoproteome###
BRCA77pho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_105BRCA/BRCA105_phosphoproteome/201507/BRCA77_unimodal_phosphoproteome-ratio-norm.txt',header=TRUE, sep="\t")
row.names(BRCA77pho) = make.names(BRCA77pho$Gene.site,unique=T)
BRCA77pho.na = BRCA77pho[,-c(1,2)]
genes = sub(".NP.*", "", row.names(BRCA77pho.na))
BRCA77pho.na.d = BRCA77pho.na[genes %in% druggable, ]
row.names(BRCA77pho.na.d) = make.names(sub(".NP_\\d+_"," ",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("\\. _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("\\._.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub(" _.*","",row.names(BRCA77pho.na.d)), unique=T)
row.names(BRCA77pho.na.d) = make.names(sub("_.*","",row.names(BRCA77pho.na.d)), unique=T)
colnames(BRCA77pho.na.d) = gsub("\\...TCGA","",colnames(BRCA77pho.na.d))
# remove duplicated ones
BRCA77pho.na.d = BRCA77pho.na.d[,!duplicated(colnames(BRCA77pho.na.d))]
# remove normal samples
BRCA77pho.na.d =  BRCA77pho.na.d[,nchar(colnames(BRCA77pho.na.d)) == 7]
colnames(BRCA77pho.na.d) = paste("TCGA",colnames(BRCA77pho.na.d),sep=".")
colnames(BRCA77pho.na.d) = gsub("\\.","-",colnames(BRCA77pho.na.d))

BRCA77pho.na.d.m = melt(as.matrix(BRCA77pho.na.d))
colnames(BRCA77pho.na.d.m) = c("Site","Sample","Pho")
BRCA77pho.na.d.m$Species = "Human"

# clinical data
PDX_clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
colnames(PDX_clin)[1] = "Sample"
PDX_clin$Intrinsic.subtype = factor(PDX_clin$Intrinsic.subtype, levels = c("Basal", "CLDN low", "HER2-E", "LumB"))
PDX_clin2 = PDX_clin[,c(1,2)]
ITRAQ.d.m = merge(ITRAQ.d.m,PDX_clin2,by="Sample")
ITRAQpho.na.d.m = merge(ITRAQpho.na.d.m,PDX_clin2,by="Sample")

h_pam50_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_105BRCA/proteome-v2.22041_data/tumor-metadata.csv"
h_pam50 = read.csv(file=h_pam50_f, header=T)
colnames(h_pam50)[3] = "Intrinsic.subtype"
h_pam50$Intrinsic.subtype = gsub("Her2","HER2-E",h_pam50$Intrinsic.subtype)
h_pam50$Intrinsic.subtype = factor(h_pam50$Intrinsic.subtype, levels = c("Basal", "HER2-E", "LumA","LumB"))
h_pam50.2 = h_pam50[,c(1,3)]
h_pam50.2$Sample = paste("TCGA",h_pam50.2$Sample,sep="-")
BRCA77.d.m = merge(BRCA77.d.m,h_pam50.2,by="Sample")
BRCA77pho.na.d.m = merge(BRCA77pho.na.d.m,h_pam50.2,by="Sample")

# make combined datasets
cat_pro = rbind(ITRAQ.d.m,BRCA77.d.m)
cat_pho = rbind(ITRAQpho.na.d.m,BRCA77pho.na.d.m)

cat_pro = cat_pro[order(cat_pro$Gene),]
cat_pho = cat_pho[order(cat_pho$Site),]
cat_pho$Site = factor(cat_pho$Site, levels= unique(cat_pho$Site)[order(unique(cat_pho$Site))])
cat_pho$Gene = gsub("\\..*","",cat_pho$Site)

# takes out TCGA label
cat_pro$Sample = gsub("TCGA-","",cat_pro$Sample)
cat_pho$Sample = gsub("TCGA-","",cat_pho$Sample)

# outlier genes
inputTable="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/druggable_events/results/2016-08-09_filtered_outliers_b.txt"
outlier=read.table(inputTable, header = T, sep = "\t")
#outlier_s = outlier[outlier$outlier_score >= 2.5,]
outlier_pro = outlier[outlier$level=="PRO_I",]
s_gene = c(as.character(unique(outlier_pro$gene)),"FGFR2","RAF1","EGFR")

outlier_pho = outlier[outlier$level=="PHO",]
outlier_phosites = unique(outlier_pho$site)

### functions
get.clinical.scale = function() {
  # Set1 colors
  colors = c(NA, "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  # use Perou's intrinsic subtype colors instead
  colors = c(NA, "#101010", NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#3182bd", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
  
  color.names = c("wt","mut","negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB", "LumA") # add lumA and check color
  #color.names = c("negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="Intrinsic Subtype", values=colors)
  
  return(clinical.color.scale)
}

plot_marker = function(geneList=geneList){
  cat_pro_g = cat_pro[cat_pro$Gene %in% geneList,]
  cat_pro_g$outlier = FALSE
  # quick outlier analysis for each gene
  for (gene in geneList){
    IQR = quantile(cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "PDX",]$Pro, probs=0.75, na.rm=T) - quantile(cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "PDX",]$Pro, probs=0.25, na.rm=T) 
    cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "PDX",]$outlier = (cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "PDX",]$Pro >= quantile(cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "PDX",]$Pro, probs=0.75, na.rm=T) + 1.5*IQR) #inner fences
    
    IQR = quantile(cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "Human",]$Pro, probs=0.75, na.rm=T) - quantile(cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "Human",]$Pro, probs=0.25, na.rm=T) 
    cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "Human",]$outlier = (cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "Human",]$Pro >= quantile(cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "Human",]$Pro, probs=0.75, na.rm=T) + 1.5*IQR) #inner fences 
    
    if (!((TRUE %in% cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "PDX",]$outlier) & (TRUE %in% cat_pro_g[cat_pro_g$Gene == gene & cat_pro_g$Species == "Human",]$outlier))){
      if (!(gene %in% c("EGFR","AKT2","AKT3"))){
        cat_pro_g = cat_pro_g[cat_pro_g$Gene != gene,]
      }  
    }
  }
  
  # make the functionally validated one apparent
  cat_pro_g[cat_pro_g$Sample=="WHIM14" & cat_pro_g$Gene == "ERBB2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM35" & cat_pro_g$Gene == "ERBB2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM8" & cat_pro_g$Gene == "ERBB2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM6" & cat_pro_g$Gene == "ERBB2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM16" & cat_pro_g$Gene == "AKT1",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM18" & cat_pro_g$Gene == "AKT1",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM20" & cat_pro_g$Gene == "AKT1",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM16" & cat_pro_g$Gene == "AKT2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM18" & cat_pro_g$Gene == "AKT2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM20" & cat_pro_g$Gene == "AKT2",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM16" & cat_pro_g$Gene == "AKT3",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM18" & cat_pro_g$Gene == "AKT3",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM20" & cat_pro_g$Gene == "AKT3",]$outlier=TRUE
  cat_pro_g[cat_pro_g$Sample=="WHIM14" & cat_pro_g$Gene == "EGFR",]$outlier=TRUE
  
  # plot violin plots faceted by marker genes
  
  p = ggplot(data=cat_pro_g)
  p = p + facet_grid(.~Gene)
  p = p + geom_boxplot(aes(x=Species, y=Pro, fill=NULL),alpha=0.1, outlier.shape = NA) 
  p = p + geom_jitter(aes(x=Species, y=Pro, color=Intrinsic.subtype), size=1) #+ geom_point(aes(x=Status, y=value)) 
  p = p + geom_text(aes(x=Species, y=Pro, label = ifelse(outlier,as.character(Sample),NA)),size=2)
  p = p + labs(x = "", y = "Protein expression")
  p = p + theme_nogrid() + guides(fill=FALSE) 
  p = p + get.clinical.scale()
  p = p + theme(text = element_text(colour="black", size=16), axis.ticks.x = element_blank(),
                axis.text.x = element_text(colour="black", size=14,angle=90,vjust=0.5),#element_blank(), axis.title.x = element_blank(),  
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p = p + theme(legend.position="bottom") + ylim(-8,8)
  p
  fn = paste(pd, "outlier_genes_pdx_human_box.pdf", sep="_")
  ggsave(file=fn, width=12, height=6,useDingbats=FALSE)
}

plot_marker(s_gene)

SiteList = c("MAP2K1.T386t","AKT2.T451t","RAF1.S295s")
plot_phospho = function(SiteList=SiteList){
  cat_pho_g = cat_pho[cat_pho$Site %in% SiteList,]
  cat_pho_g$outlier = FALSE
  # quick outlier analysis for each Site
  for (Site in SiteList){
    IQR = quantile(cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "PDX",]$Pho, probs=0.75, na.rm=T) - quantile(cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "PDX",]$Pho, probs=0.25, na.rm=T) 
    cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "PDX",]$outlier = (cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "PDX",]$Pho >= quantile(cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "PDX",]$Pho, probs=0.75, na.rm=T) + 1.5*IQR) #inner fences
    
    IQR = quantile(cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "Human",]$Pho, probs=0.75, na.rm=T) - quantile(cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "Human",]$Pho, probs=0.25, na.rm=T) 
    cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "Human",]$outlier = (cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "Human",]$Pho >= quantile(cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "Human",]$Pho, probs=0.75, na.rm=T) + 1.5*IQR) #inner fences 
    
    if (!((TRUE %in% cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "PDX",]$outlier) & (TRUE %in% cat_pho_g[cat_pho_g$Site == Site & cat_pho_g$Species == "Human",]$outlier))){
      #if (!(cat_pho_g$Gene %in% c("EGFR","AKT2","AKT3"))){
      cat_pho_g = cat_pho_g[cat_pho_g$Site != Site,]
      #}  
    }
  }
  
  #   # make the functionally validated one apparent
    cat_pho_g[cat_pho_g$Sample=="WHIM14" & cat_pho_g$Site == "ERBB2.T701t",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM35" & cat_pho_g$Site == "ERBB2.T701t",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM8" & cat_pho_g$Site == "ERBB2.T701t",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM6" & cat_pho_g$Site == "ERBB2.T701t",]$outlier=TRUE
  
  cat_pho_g[cat_pho_g$Sample=="WHIM14" & cat_pho_g$Site == "ERBB2.T1240t",]$outlier=TRUE
  cat_pho_g[cat_pho_g$Sample=="WHIM35" & cat_pho_g$Site == "ERBB2.T1240t",]$outlier=TRUE
  cat_pho_g[cat_pho_g$Sample=="WHIM8" & cat_pho_g$Site == "ERBB2.T1240t",]$outlier=TRUE
  cat_pho_g[cat_pho_g$Sample=="WHIM6" & cat_pho_g$Site == "ERBB2.T1240t",]$outlier=TRUE
  
  cat_pho_g[cat_pho_g$Sample=="WHIM14" & cat_pho_g$Site == "ERBB2.Y1248y",]$outlier=TRUE
  cat_pho_g[cat_pho_g$Sample=="WHIM35" & cat_pho_g$Site == "ERBB2.Y1248y",]$outlier=TRUE
  cat_pho_g[cat_pho_g$Sample=="WHIM8" & cat_pho_g$Site == "ERBB2.Y1248y",]$outlier=TRUE
  cat_pho_g[cat_pho_g$Sample=="WHIM6" & cat_pho_g$Site == "ERBB2.Y1248y",]$outlier=TRUE
  
    cat_pho_g[cat_pho_g$Sample=="WHIM16" & cat_pho_g$Site == "AKT1.S122s",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM18" & cat_pho_g$Site == "AKT1.S122s",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM20" & cat_pho_g$Site == "AKT1.S122s",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM16" & cat_pho_g$Site == "AKT1.S475s",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM18" & cat_pho_g$Site == "AKT1.S475s",]$outlier=TRUE
    cat_pho_g[cat_pho_g$Sample=="WHIM20" & cat_pho_g$Site == "AKT1.S475s",]$outlier=TRUE
  #   cat_pho_g[cat_pho_g$Sample=="WHIM16" & cat_pho_g$Site == "AKT3",]$outlier=TRUE
  #   cat_pho_g[cat_pho_g$Sample=="WHIM18" & cat_pho_g$Site == "AKT3",]$outlier=TRUE
  #   cat_pho_g[cat_pho_g$Sample=="WHIM20" & cat_pho_g$Site == "AKT3",]$outlier=TRUE
  #   cat_pho_g[cat_pho_g$Sample=="WHIM14" & cat_pho_g$Site == "EGFR",]$outlier=TRUE
  
  
  # get clear phospho annotations
  cat_pho_g$Site = gsub("[a-z]","",cat_pho_g$Site)
  cat_pho_g$Site = gsub("\\."," p.",cat_pho_g$Site)
  
  # plot violin plots faceted by marker Sites
  p = ggplot(data=cat_pho_g)
  #p = p + facet_wrap(Gene~Site, nrow=5)
  p = p + facet_grid(.~Site)
  p = p + geom_boxplot(aes(x=Species, y=Pho, fill=NULL),alpha=0.1, outlier.shape = NA) 
  p = p + geom_jitter(aes(x=Species, y=Pho, color=Intrinsic.subtype), size=1) #+ geom_point(aes(x=Status, y=value)) 
  p = p + geom_text(aes(x=Species, y=Pho, label = ifelse(outlier,as.character(Sample),NA)),size=2)
  p = p + theme_nogrid() + guides(fill=FALSE) 
  p = p + labs(x = "", y = "Phosphosite expression")
  p = p + get.clinical.scale()
  p = p + theme(text = element_text(colour="black", size=16), axis.ticks.x = element_blank(),
                axis.text.x = element_text(colour="black", size=14,angle=90,vjust=0.5),#element_blank(), axis.title.x = element_blank(),  
                axis.text.y = element_text(colour="black", size=14), strip.text = element_text(size = 8))
  p = p + theme(legend.position="bottom") + ylim(-7.2,7.2)
  p
  fn = paste(pd, "outlier_Sites_pdx_human_box.pdf", sep="_")
  ggsave(file=fn, width=12, height=6,useDingbats=FALSE)
}
#plot_phospho(outlier_phosites)
selected_sites = c("AKT1.S122s","AKT1.S475s","ARAF.S299s","BRAF.S447s","BRAF.S750s",
                   "ERBB2.T701t","ERBB2.T1240t","ERBB2.Y1248y","HSP90AB1.Y56y","HSP90B1.S169s",
                   "MTOR.S2481s","PLK1.T210t","TOP2A.S1213s")
plot_phospho(selected_sites)
