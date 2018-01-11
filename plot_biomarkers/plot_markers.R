# heatmap_immuno.R by Kuan Huang @ WashU 201506
# cluster WHIMs based on expressions in the PI3K pathway

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/plot_biomarkers")
source("/Users/khuang/bin/LIB_exp.R")

# dependencies that can potentially be helpful
tri_neg=c("WHIM2","WHIM4","WHIM6","WHIM12","WHIM13","WHIM14","WHIM21","WHIM25","WHIM30")
lym = c("WHIM17","WHIM46")
lumb_er = c("WHIM9","WHIM18","WHIM16","WHIM20","WHIM24","WHIM26","WHIM27","WHIM37","WHIM43")
PI3mut = c("WHIM9","WHIM18","WHIM16","WHIM20","WHIM24")
exp = c("WHIM16","WHIM18","WHIM20")
clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
Basal=as.vector(clin[clin$Intrinsic.subtype=="Basal",1])
LumB=as.vector(clin[clin$Intrinsic.subtype=="LumB",1])
Her2=as.vector(clin[clin$Intrinsic.subtype=="HER2-E",1])
CLDN_low=as.vector(clin[clin$Intrinsic.subtype=="CLDN low",1])

# Datasets
ITRAQ = read.table(header=TRUE, sep="\t", row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')

# function
plot_pro_marker = function(marker=NULL){
  PRO_marker = ITRAQ[row.names(ITRAQ) %in% marker,]
  PRO_marker_p = data.frame(t(PRO_marker))
  PRO_marker_p$Sample="WHIM"
  PRO_marker_p$WHIM = row.names(PRO_marker_p)
  colnames(PRO_marker_p)[1] = "Value"
  
  fn = paste(pd, marker,"WHIM_Protein_exp.pdf", sep="_")
  p = ggplot(data=PRO_marker_p)
  p = p + geom_violin(aes(x=Sample, y=Value, fill="Sample"),alpha=0.1, linetype="blank") + guides(fill=FALSE) 
  p = p + geom_text(aes(x=Sample, y=Value, label=WHIM),position = position_jitter(w = 0.2, h = 0), size=3) + guides(label=FALSE) 
  p = p + labs(title = marker, y = "Protein expression") + theme_bw()
  p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_blank(), 
                axis.text.y = element_text(colour="black", size=12), strip.text = element_text(size = 8, angle = 90))
  p
  ggsave(file=fn, height=3, width=5, useDingbats=FALSE)
}


ITRAQpho = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp.txt',header=TRUE, sep="\t")
row.names(ITRAQpho) = ITRAQpho$gene.site
colnames(ITRAQpho)<-sub("\\..*", "", colnames(ITRAQpho))
ITRAQpho.na=ITRAQpho[,-c(1,2,19,20,22)]

plot_pho_marker = function(marker=NULL){
  PHO_marker = ITRAQpho.na[row.names(ITRAQpho) %in% marker,]
  PHO_marker_p = data.frame(t(PHO_marker))
  PHO_marker_p$Sample="WHIM"
  PHO_marker_p$WHIM = row.names(PHO_marker_p)
  colnames(PHO_marker_p)[1] = "Value"
  marker_name = sub("_\\d_.*","",marker)
  
  fn = paste(pd, marker,"WHIM_phosphosite_exp.pdf", sep="_")
  p = ggplot(data=PHO_marker_p)
  #p = p + facet_grid(variable~.)
  p = p + geom_violin(aes(x=Sample, y=Value, fill="Sample"),alpha=0.1, linetype="blank") + guides(fill=FALSE) 
  p = p + geom_text(aes(x=Sample, y=Value, label=WHIM),position = position_jitter(w = 0.2, h = 0), size=3) + guides(label=FALSE) 
  p = p + labs(title = marker_name, y = "Phosphosite expression") + theme_bw()
  p = p + theme(text = element_text(colour="black", size=14), axis.text.x = element_blank(), 
                axis.text.y = element_text(colour="black", size=12), strip.text = element_text(size = 8, angle = 90))
  p
  ggsave(file=fn, height=5, width=5, useDingbats=FALSE)
}

plot_pro_marker("MS4A1")
#plot_pro_marker("TRAF")
plot_pro_marker("TBX21")
plot_pho_marker("DPYSL3-NP_001184223_Y146y _1_1_146_146")
plot_pho_marker("TRAF1-NP_005649_S170s _1_1_170_170")
plot_pho_marker("TBX21-NP_037483_S503s _1_1_503_503")
