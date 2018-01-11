# heatmap_immuno.R by Kuan Huang @ WashU 201506
# cluster WHIMs based on expressions in the PI3K pathway

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/RPPA")
source("/Users/khuang/bin/LIB_exp.R")

# exp = c("WHIM16","WHIM18","WHIM20")
# 
ITRAQ = read.table(header=TRUE, sep="\t", row.names=1, file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt')
samples_24 = colnames(ITRAQ)
# 
# ## Phosphoproteome
# 
# ITRAQpho = read.table(row.names=1,file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt',header=TRUE, sep="\t")
# ITRAQpho_pi3k = as.matrix(ITRAQpho[row.names(ITRAQpho) %in% pi3k,])
# 
# ITRAQpho_pi3k=ITRAQpho_pi3k[rowSums(is.na(ITRAQpho_pi3k))<5,]
# ITRAQpho_pi3k=ITRAQpho_pi3k[order(row.names(ITRAQpho_pi3k)),]
# ITRAQpho_pi3k = data.frame(ITRAQpho_pi3k)
# ITRAQpho_pi3k$Dataset = "ITRAQ Phosphoprotein"
# ITRAQpho_pi3k$Gene = row.names(ITRAQpho_pi3k)

## RPPA
RPPA_file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_RPPA/ana55_169slides/Ana55_169\ Slides_Tissue\ ID_RPPA\ DATA_Set43_curated.txt"
#RPPA = read.table(file=RPPA_file,header=TRUE, sep="\t") # sample name by passage name
RPPA = read.table(file=RPPA_file,header=TRUE, skip=1, sep="\t") # sample name with only the latest passage
colnames(RPPA)[1] = "Antibody"
colnames(RPPA)[2] = "Genes"
RPPA_24 = RPPA[,colnames(RPPA) %in% c(samples_24,"Antibody","Genes")]

RPPA_24_m = melt(RPPA_24,by=c("Antibody","Genes"))
RPPA_24_m_p = RPPA_24_m[RPPA_24_m$Genes %in% c("AKT1, AKT2, AKT3,","PIK3CA,","MTOR,","ERBB2, "),]

fn = paste(pd, "RPPA_test.pdf", sep="_")
p = ggplot(data=RPPA_24_m_p)
p = p + facet_grid(Genes~.)
p = p + coord_equal()
p = p + geom_tile(aes(x=variable, y=Antibody, fill=value), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA)#, limits = c(-bound,bound))
p = p  + theme_bw() + 
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=14, angle=90), axis.text.y = element_text(colour="black", size=14),
        strip.text = element_text(size = 14))  #element_text(colour="blaglist", size=14))
p
ggsave(file=fn, height=7.5, width=10, useDingbats=FALSE)


#RPPA = read.table(file='/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/WHIM_RPPA/20151020_RPPA_curated.txt',header=TRUE, sep="\t")
# RPPA = RPPA[,c(2,3,4,5,6)]
# RPPA = RPPA[RPPA$Gene %in% pi3k,]
# RPPA = RPPA[order(RPPA$Gene),]





# ### merge with proteomic and phospho dataset ###
# ITRAQ_pi3k_exp = ITRAQ_pi3k[,colnames(ITRAQ_pi3k) %in% c(exp,"Gene","Dataset")]
# ITRAQpho_pi3k_exp = ITRAQpho_pi3k[,colnames(ITRAQpho_pi3k) %in% c(exp,"Gene","Dataset")]
# 
# ITRAQ_pi3k_m = melt(ITRAQ_pi3k_exp, id.var=c("Dataset","Gene"))
# ITRAQpho_pi3k_m = melt(ITRAQpho_pi3k_exp, id.var=c("Dataset","Gene"))
# RPPA_m = melt(RPPA, id.var=c("Dataset","Gene"))
# 
# pro_lvls=rbind(ITRAQ_pi3k_m, ITRAQpho_pi3k_m, RPPA_m)
# #pro_lvls = pro_lvls[order(pro_lvls$Sample),]
# colnames(pro_lvls) = c("Dataset","Gene","Sample","Expression")
# 
# pro_lvls$Sample = factor(pro_lvls$Sample, levels = c("WHIM16", "WHIM18", "WHIM20"))
# pro_lvls$Dataset = factor(pro_lvls$Dataset, levels = c("RPPA Protein", "ITRAQ Protein", "RPPA Phosphoprotein","ITRAQ Phosphoprotein"))
# 
# fn = paste(pd, "3WHIM_pi3k-pathway.pdf", sep="_")
# # summary$truncated_outlier_percentage = summary$rounded_outlier_percentage
# # summary[summary$truncated_outlier_percentage>=22,]$truncated_outlier_percentage = 22
# min_d = min(as.numeric(as.character(pro_lvls$Expression)), na.rm=T)
# max_d = max(as.numeric(as.character(pro_lvls$Expression)), na.rm=T)
# bound = max(c(max_d, -min_d))
# 
# p = ggplot(data=pro_lvls)
# p = p + facet_grid(.~Sample)
# p = p + coord_equal()
# p = p + geom_tile(aes(x=Dataset, y=Gene, fill=Expression), linetype="blank") + scale_fill_gradientn(name= "Expression", colours=getPalette(100), na.value=NA, limits = c(-bound,bound))
# p = p  + theme_bw() + 
#   theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="blaglist", size=14, angle=90), axis.text.y = element_text(colour="blaglist", size=14),
#         axis.tiglists = element_blank(), strip.text = element_text(size = 14))#element_text(colour="blaglist", size=14))
# p
# ggsave(file=fn, height=7.5, width=10, useDingbats=FALSE)