# vaf_correlation.R by Kuan Huang @ WashU 201507
# VAF correlation plotting and analysis for WHIM manuscript figure 1

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/vaf_correlation")
source("/Users/khuang/bin/plotting_essentials.R")

vaf=read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_variants/2013_cell_report_vafs.txt', header=TRUE, sep = "\t")
#vaf$Gene_aa = paste(vaf$Gene,vaf$Amino_acid_change)
cor(vaf$VAF_Tumor, vaf$VAF_Xeno, method ="pearson")
vaf$gene = rep("Other", nrow(vaf))
vaf[vaf$Gene == "PIK3CA",]$gene="PIK3CA"
vaf[vaf$Gene == "TP53",]$gene="TP53"
  
fn = paste(pd,'2013_cellReport_tumor_vs_xeno_VAF.pdf', sep="_")
p = ggplot(data = vaf, aes(x=VAF_Tumor , y=VAF_Xeno, color=WHIM, shape=gene)) 
p = p + geom_point(alpha=0.7, size=3) + xlab("VAF in human tumor sample") + ylab("VAF in xenograft") + theme_bw()
p = p + scale_shape_manual(values=c(16,0,2), guide="legend") +guides(colour=FALSE)
p = p + xlim(0,1) + ylim(0,1) + coord_fixed() 
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(angle = 90, vjust = 0.5, colour="black", size=16), 
              axis.text.y = element_text(colour="black", size=16))
p
ggsave(file=fn)