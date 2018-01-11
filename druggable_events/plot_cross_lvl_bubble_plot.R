# plot_mut_exp.R by Kuan Huang @ WashU 201507
# find differentially expressed genes due to effect of input mutations

# mis
library(plyr)
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/druggable_events/")
source("/Users/khuang/bin/plotting_essentials.R") 

cat(paste("Date: ", date, "\n", sep=""))
#args=commandArgs(TRUE)
inputTable="/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/druggable_events/results/2016-08-09_filtered_outliers_b.txt"


file=read.table(inputTable, header = T, sep = "\t")

mut.data = c("WHIM12","MUT","PIK3CA","PIK3CA p.V105_E109delinsT",2,
             "WHIM16","MUT","PIK3CA","PIK3CA p.H1047R",2,
             "WHIM18","MUT","PIK3CA","PIK3CA p.E545K",2,
             "WHIM20","MUT","PIK3CA","PIK3CA p.H542K",2,
             "WHIM24","MUT","PIK3CA","PIK3CA p.H1047R",2,
             "WHIM26","MUT","SF3B1","SF3B1 p.K700E",2,
             "WHIM9","MUT","PIK3CA","PIK3CA p.H1047R",2,
             "WHIM9","MUT","KRAS","KRAS p.A146V",2)
mutation=data.frame(matrix(data=mut.data,ncol = 5, byrow=T))
colnames(mutation) = colnames(file)
mutation$outlier_score = as.numeric(as.character(mutation$outlier_score))
outlier = rbind(mutation,file)

### data level ###
outlier = outlier[outlier$level != "PRO_L",]

outlier$level = factor(outlier$level, levels = c("MUT","CNV", "RNA", "PRO_I", "PHO"))
outlier$level = revalue(outlier$level, c("MUT"="Mutation", "PRO_I"="Protein", "PHO"="Phosphosite"))

### sample order###
clin = read.table(file='/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/data/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified.txt', header=T, sep='\t')
colnames(clin)[1] = "sample"
clin$Intrinsic.subtype = factor(clin$Intrinsic.subtype, levels = c("Basal", "CLDN low", "HER2-E", "LumB"))
clin2 = clin[,c(1,2)]
top_overlap2 = merge(clin2, outlier, by="sample",all.x=T)

sample.order =  unique(top_overlap2[order(top_overlap2$Intrinsic.subtype),c("sample")])
top_overlap2$sample = factor(top_overlap2$sample, levels=sample.order)

top_overlap2_s = top_overlap2[top_overlap2$outlier_score >= 2.5,]
s_gene = unique(top_overlap2_s$gene)

top_overlap2_gs = top_overlap2[top_overlap2$gene %in% c(as.character(s_gene),"FGFR2","RAF1"),]

# add WHIM14 as a sample level as it has no hits at all...
#levels(top_overlap2_gs$sample) <- c(levels(top_overlap2_gs$sample), "WHIM14")
# reorder
top_overlap2_gs$sample = factor(top_overlap2_gs$sample, levels=sample.order)


### bubble plot ###
YlGnBu = brewer.pal(9, "YlGnBu") 
getPalette = colorRampPalette(YlGnBu)
#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
#+ guides(fill=FALSE)
top_overlap2$sample = factor(top_overlap2$sample, levels=sample.order)

fn = paste(pd, "_all_lvl_outlier.pdf", sep="_")
p = ggplot(data=top_overlap2, aes(x=gene, y=sample, size=outlier_score, color = level))
#p = ggplot(data=outlier, aes(x=gene, y=sample, size=outlier_score))#, color=level))
#p = p + facet_grid(.~level)
p = p + facet_grid(.~level,drop=F,scales = "free", space = "free")
p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
#p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
p = p + scale_colour_manual(values = brewer.pal(5,"Set1"), guide=F)
p = p + labs( x = "Druggable Target", y = "Sample") + theme_bw()
p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=8, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
ggsave(file=fn, height=5, width=10, useDingbats=FALSE)

fn = paste(pd, "_all_lvl_outlier_geneView.pdf", sep="_")
p = ggplot(data=top_overlap2, aes(y=gene, x=level, size=outlier_score, color = level))
#p = p + facet_grid(.~sample,drop=T,scales = "free", space = "free")
p = p + facet_grid(.~sample)
p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
#p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
p = p + scale_colour_manual(values = brewer.pal(5,"Set1"), guide=F)
p = p + labs( x = "Druggable Target", y = "Sample") + theme_bw()
p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=8, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
ggsave(file=fn, height=8, width=16, useDingbats=FALSE)

fn = paste(pd, "_all_lvl_outlier_geneView_byPAM50.pdf", sep="_")
p = ggplot(data=top_overlap2, aes(y=gene, x=level, size=outlier_score, color = level))
p = p + facet_grid(Intrinsic.subtype~sample, drop=T, scales="free_y", space="free_y")
p = p + geom_point()
p = p + scale_colour_manual(values = brewer.pal(5,"Set1"), guide=F)
p = p + labs( x = "Druggable Target", y = "Sample") + theme_bw()
p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=8, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
ggsave(file=fn, height=10, width=16, useDingbats=FALSE)

fn = paste(pd, "_all_lvl_outlier_geneView_v2.pdf", sep="_")
p = ggplot(data=top_overlap2, aes(y=level, x=gene, size=outlier_score, color = level))
#p = p + facet_grid(.~sample,drop=T,scales = "free", space = "free")
p = p + facet_grid(sample~.)
p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
#p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
p = p + scale_colour_manual(values = brewer.pal(5,"Set1"), guide=F)
p = p + labs( x = "Druggable Target", y = "Sample") + theme_bw()
p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=8, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
ggsave(file=fn, height=16, width=8, useDingbats=FALSE)

### selected gene ###
top_overlap2_gs$cat = "Other"
top_overlap2_gs[top_overlap2_gs$gene %in% c("AKT1","AKT2","AKT3","PIK3CA","PIK3CG","MTOR"),]$cat = "PI3K"
top_overlap2_gs[top_overlap2_gs$gene %in% c("PDGFRA","PDGFRB","ERBB2","FGFR4","AURKA","FGFR2"),]$cat = "RTK-Signaling"
top_overlap2_gs[top_overlap2_gs$gene %in% c("MAPK1","ARAF","BRAF","RAF1"),]$cat = "MAPK-Signaling"
top_overlap2_gs[top_overlap2_gs$gene %in% c("MDM2","PRKDC","TOP2A"),]$cat = "Genome Integrity"
top_overlap2_gs$cat = factor(top_overlap2_gs$cat, levels = c("RTK-Signaling", "PI3K", "MAPK-Signaling", "Genome Integrity","Other"))
top_overlap2_gs$gene = factor(top_overlap2_gs$gene, levels = unique(top_overlap2_gs$gene)[order(as.character((unique(top_overlap2_gs$gene))))])

#sample.order = factor(sample.order, levels = sample.order)
#top_overlap2_gs$sample = factor(top_overlap2_gs$sample, levels=sample.order)
#top_overlap2_gs$sample_s = factor(top_overlap2_gs$sample, levels=sample.order)
# top_overlap2_gs$sample = factor(top_overlap2_gs$sample, levels= c("WHIM13", "WHIM14", "WHIM17", "WHIM2","WHIM21" , "WHIM25" , "WHIM30" , "WHIM4"  , "WHIM46" , "WHIM6"  , "WHIM12" , "WHIM11" , "WHIM35" , "WHIM47"
#                                                                     , "WHIM8"  , "WHIM16" , "WHIM18" , "WHIM20" , "WHIM24" , "WHIM26" , "WHIM27" , "WHIM37" , "WHIM43" , "WHIM9" ))
#top_overlap2_gs = top_overlap2_gs[order(top_overlap2_gs$sample_s),]

top_overlap2_gs[top_overlap2_gs$sample=="WHIM14",] = c("WHIM14","Basal","Mutation","PIK3CA",NA,0,"PI3K")
top_overlap2_gs$outlier_score = as.numeric(top_overlap2_gs$outlier_score)

fn = paste(pd, "all_lvl_outlier_geneView_selected.pdf", sep="_")
p = ggplot(data=top_overlap2_gs, aes(y=gene, x=level, size=outlier_score, color = level))
#p = p + facet_grid(.~sample,drop=T,scales = "free", space = "free")

p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
#p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
p = p + scale_colour_manual(values = brewer.pal(5,"Set1"))
p = p + facet_grid(cat~sample,drop=F,scales = "free_y", space = "free_y")
p = p + labs( x = "Druggable Target", y = "Sample") + theme_nogrid()
p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_blank(),#element_text(colour="black", size=8, angle=90, vjust=0.5), 
              axis.text.y = element_text(colour="black", size=8), axis.ticks = element_blank(),
              legend.position="bottom")
p
ggsave(file=fn, height=6, width=16, useDingbats=FALSE)

### pi3k pathway for the PTRC grant ### 
pi3k = c("PIK3CA","AKT1","AKT2","AKT3","MTOR")
outlier_pi3k = outlier[outlier$gene %in% pi3k,]
fn = paste(pd, "_all_lvl_outlier_pi3k.pdf", sep="_")
p = ggplot(data=outlier_pi3k, aes(x=gene, y=sample, size=outlier_score, color = level))
#p = ggplot(data=outlier, aes(x=gene, y=sample, size=outlier_score))#, color=level))
#p = p + facet_grid(.~level)
p = p + facet_grid(.~level,drop=T,scales = "free", space = "free")
p = p + geom_point()#alpha=0.5) #+ scale_colour_brewer(palette = co)#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
#p = p + scale_colour_manual(values = col_vector[1:length(unique(outlier$gene))],guide = FALSE)
p = p + scale_colour_manual(values = brewer.pal(5,"Set1"), guide=F)
p = p + labs( x = "Druggable Target", y = "Sample") + theme_bw()
p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=8, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=8))
p
ggsave(file=fn, height=3, width=4, useDingbats=FALSE)

# fn = paste(pd, "_all_lvl_outlier.pdf", sep="_")
# p = ggplot(data=outlier, aes(x=gene, fill=level))
# p = p + geom_bar()  + scale_fill_brewer(palette = "Set1")#+ scale_color_gradientn(name= "Outlier Score", colours=getPalette(100), na.value=NA, limits=c(0,6))
# p = p + labs( x = "Druggable Target", y = "Sample") + theme_bw()
# p = p + theme(text = element_text(colour="black", size=8), axis.text.x = element_text(colour="black", size=8, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=14))
# p
# ggsave(file=fn, height=7, width=7)

