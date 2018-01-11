# plot_treatment.R by Kuan Huang @ WashU 201604

### dependencies ###

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/")
source("/Users/khuang/bin/plotting_essentials.R") 
library("agricolae")
library(ggplot2)
library(reshape2)

### functions ###

anova_tukey = function(data){
  t_model = lm(data=data, formula= Tumor_Volume ~ Days + Group)
  summary(t_model)
  anova(t_model)
  a1 = aov(data=data, formula= Tumor_Volume ~ Group + Days)
  posthoc = TukeyHSD(x=a1, 'Group', conf.level=0.95)
  posthoc$Group
}

plot_curve = function(melted_m, name){
  fn = paste(pd, name, sep="_")
  p = ggplot(melted_m, aes(as.numeric(Days),as.numeric(Tumor_Volume), color=as.factor(Group))) +
    stat_summary(fun.data=mean_se, geom="pointrange")+ 
    aes(colour = as.factor(Group)) + stat_summary(fun.y = mean, geom="line") +
    labs(y=as.character("Tumor Volume (mm3)"), x=as.character("Days Post Treatment")) +
    theme_nogrid() + theme(legend.title=element_blank(), legend.position = "bottom") + 
    theme(axis.title.y = element_text(size = 14, angle = 90)) + 
    theme(axis.title.x = element_text(size = 14, angle = 0))
  p = p + facet_grid(.~Sample,drop=T,scales = "free", space = "free")
  p = p + expand_limits(x = 0, y = 0)
  p
  ggsave(file=fn, useDingbats=FALSE)
}

### all ### 
WHIM14 = read.table(header=T, sep="\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/20160907_WHIM14_lapatinib_exps_FC.txt")
wilcox.test(WHIM14[WHIM14$Sample=="WHIM14" & WHIM14$Group =="Vehicle","Fold_change"],WHIM14[WHIM14$Sample=="WHIM14" & WHIM14$Group =="Lapatinib","Fold_change"])
wilcox.test(WHIM14[WHIM14$Sample=="WHIM14" & WHIM14$Group =="Vehicle","Fold_change_24"],WHIM14[WHIM14$Sample=="WHIM14" & WHIM14$Group =="Lapatinib","Fold_change_24"])

WHIM14_t_p = read.table(header=T, sep="\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/20160907_WHIM14_herceptin_pertuzumab_exps_FC.txt")
wilcox.test(WHIM14_t_p[WHIM14_t_p$Sample=="WHIM14" & WHIM14_t_p$Group =="Mock","Fold_change"],WHIM14_t_p[WHIM14_t_p$Sample=="WHIM14" & WHIM14_t_p$Group =="Trastuzumab","Fold_change"])
wilcox.test(WHIM14_t_p[WHIM14_t_p$Sample=="WHIM14" & WHIM14_t_p$Group =="Mock","Fold_change"],WHIM14_t_p[WHIM14_t_p$Sample=="WHIM14" & WHIM14_t_p$Group =="Pertuzumab","Fold_change"])


her2 = read.table(header=T, sep="\t","/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/20160524_HER2_exps_FC.txt")

wilcox.test(her2[her2$Sample=="WHIM14-2" & her2$Group =="Vehicle","Fold_change_end"],her2[her2$Sample=="WHIM14-2" & her2$Group =="Neratinib","Fold_change_end"])
wilcox.test(her2[her2$Sample=="WHIM8-2" & her2$Group =="Vehicle","Fold_change_end"],her2[her2$Sample=="WHIM8-2" & her2$Group =="Herceptin","Fold_change_end"])

# her2$FC_g = "Slow_grower"
# her2[her2$Sample %in% c("WHIM2","WHIM6","WHIM30"),]$FC_g = "Fast_grower"
# log2(her2$Fold_change)
# 
# wilcox.test(her2[her2$Sample=="WHIM14-2" & her2$Group =="Vehicle","Fold_change"],her2[her2$Sample=="WHIM14-2" & her2$Group =="Neratinib","Fold_change"])
# wilcox.test(her2[her2$Sample=="WHIM30" & her2$Group =="Vehicle","Fold_change"],her2[her2$Sample=="WHIM30" & her2$Group =="Lapatinib","Fold_change"])
# t.test(her2[her2$Sample=="WHIM30" & her2$Group =="Vehicle","Fold_change"],her2[her2$Sample=="WHIM30" & her2$Group =="Lapatinib","Fold_change"])
# 
# #t.test(her2[her2$Sample=="WHIM14-2" & her2$Group =="Vehicle","Fold_change"],her2[her2$Sample=="WHIM14-2" & her2$Group =="Neratinib","Fold_change"])
# 
# #p = ggplot(her2, aes(x = Group , y = log2(Fold_change)))
# p = ggplot(her2, aes(x = Group , y = Fold_change, fill=Center))
# p = p + facet_grid(FC_g~Sample,drop=T,scales = "free")#, space = "free")
# p = p + geom_boxplot()#her2,aes(x = Group , y = as.numeric(Fold_change)))
# p = p + geom_jitter(alpha=0.3)
# p = p + theme_bw() + theme_nogrid()
# p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
# p
# 
# fn = paste(pd, "HER2_inhibition_treatment_FC.pdf", sep="_")
# ggsave(file=fn, width=8,useDingbats=FALSE)
# 
# ### 
# p = ggplot(her2, aes(x = Group , y = Fold_change, fill=Center))
# p = p + facet_grid(FC_g~Sample,drop=T,scales = "free")#, space = "free")
# p = p + geom_boxplot()#her2,aes(x = Group , y = as.numeric(Fold_change)))
# p = p + geom_jitter(alpha=0.3)
# p = p + theme_bw() + theme_nogrid()
# p = p + theme(axis.title = element_text(size=12), axis.text.x = element_text(colour="black", size=10, angle=90, vjust=0.5), axis.text.y = element_text(colour="black", size=10))#element_text(colour="black", size=14))
# p
# 
# fn = paste(pd, "HER2_inhibition_treatment_FC.pdf", sep="_")
# ggsave(file=fn, width=8,useDingbats=FALSE)


her2$g = "HER2-E"
her2[her2$Sample %in% c("WHIM2","WHIM6","WHIM30"),]$g = "HER2-negative"
her2[her2$Sample %in% c("WHIM14-1","WHIM14-2"),]$g = "WHIM14"
her2_2 = her2[her2$Sample %in% c("WHIM6","WHIM8-1","WHIM35","WHIM14-1"),]
her2_2$Sample = as.character(her2_2$Sample)
her2_2[her2_2$Sample == "WHIM8-1",]$Sample = "WHIM8"
her2_2[her2_2$Sample == "WHIM14-1",]$Sample = "WHIM14"
her2_2$Sample = factor(her2_2$Sample, levels=c("WHIM6","WHIM8","WHIM35","WHIM14"))
her2_2$Group = factor(her2_2$Group, levels=c("Vehicle","Lapatinib"))

# statistical test for lapatinib exp
wilcox.test(her2_2[her2_2$Sample=="WHIM6" & her2_2$Group =="Vehicle","Fold_change"],her2_2[her2_2$Sample=="WHIM6" & her2_2$Group =="Lapatinib","Fold_change"])
wilcox.test(her2_2[her2_2$Sample=="WHIM8" & her2_2$Group =="Vehicle","Fold_change"],her2_2[her2_2$Sample=="WHIM8" & her2_2$Group =="Lapatinib","Fold_change"])
wilcox.test(her2_2[her2_2$Sample=="WHIM35" & her2_2$Group =="Vehicle","Fold_change"],her2_2[her2_2$Sample=="WHIM35" & her2_2$Group =="Lapatinib","Fold_change"])
wilcox.test(her2_2[her2_2$Sample=="WHIM14" & her2_2$Group =="Vehicle","Fold_change"],her2_2[her2_2$Sample=="WHIM14" & her2_2$Group =="Lapatinib","Fold_change"])

p = ggplot(her2_2, aes(x = Group , y = Fold_change, fill=Group))
p = p + facet_wrap(~Sample,drop=T,nrow=1, scales = "free_x")#, space = "free")
p = p + geom_violin(alpha = 0.3,color=NA)#her2,aes(x = Group , y = as.numeric(Fold_change)))
p = p + geom_jitter(alpha=0.3)
p = p + stat_summary(color="gray25", fun.ymin=median, fun.ymax=median, geom="errorbar", size=0.5)
p = p + theme_bw() + theme_nogrid()
p = p + labs(x="",y="Fold change (%)") #+ ylim(-200,1100)
p = p + theme(axis.title = element_text(size=12), axis.text.x = element_blank(),#_text(colour="black", size=10, angle=90, vjust=0.5), 
              axis.text.y = element_text(colour="black", size=10), axis.ticks.x=element_blank())#element_text(colour="black", size=14))
p

fn = paste(pd, "HER2_inhibition_treatment_FC_organized.pdf", sep="_")
ggsave(file=fn, height=3, useDingbats=FALSE)
