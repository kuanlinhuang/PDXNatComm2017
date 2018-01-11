# plot_treatment.R by Kuan Huang @ WashU 201604

### dependencies ###

setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/")
source("/Users/khuang/bin/plotting_essentials.R") 
library("agricolae")
library(ggplot2)
library(reshape2)

### functions ###

preprocess_data = function(data){
  data_m = melt(data, id.var=c("Sample","Group", "Replicate"))
  colnames(data_m)[4] = "Day"
  colnames(data_m)[5] = "Tumor_Volume"
  data_m$Days = as.numeric(gsub("X","",data_m$Day))
  data_m$Group = relevel(data_m$Group, ref="Vehicle") #base the reference off the Vehicle group
}

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank())

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

### WHIM8 ### 
whim8 = read.csv("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/HER2_WHIM8_experiments.txt", 
                 header=T, sep="\t")

whim8_m = melt(whim8, id.var=c("Sample","Group", "Replicate"))
colnames(whim8_m)[4] = "Day"
colnames(whim8_m)[5] = "Tumor_Volume"
#whim8_m$Tumor_Volume = whim8_m$Tumor_Volume_KG * 1000
whim8_m$Days = as.numeric(gsub("X","",whim8_m$Day))
whim8_m$Group = relevel(whim8_m$Group, ref="Vehicle")

t_model = lm(data=whim8_m, formula= Tumor_Volume ~ Days + Group)
summary(t_model)
anova(t_model)

### WHIM14 ### 
whim14 = read.csv("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/HER2_WHIM14_experiments.txt", 
                header=T, sep="\t")

whim14_m = melt(whim14, id.var=c("Sample","Group", "Replicate"))
colnames(whim14_m)[4] = "Day"
colnames(whim14_m)[5] = "Tumor_Volume_KG"
whim14_m$Tumor_Volume = whim14_m$Tumor_Volume_KG * 1000
whim14_m$Days = as.numeric(gsub("X","",whim14_m$Day))
whim14_m$Group = relevel(whim14_m$Group, ref="Vehicle")

t_model = lm(data=whim14_m, formula= Tumor_Volume ~ Days + Group)
summary(t_model)
anova(t_model)

### WHIM16,18,20 ###
akt_whims = read.csv("/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/targeted_treatment/AKT_WHIM16.18.20_experiments.txt", 
                  header=T, sep="\t")

akt_whims_m = melt(akt_whims, id.var=c("Sample","Group", "Replicate"))
colnames(akt_whims_m)[4] = "Day"
colnames(akt_whims_m)[5] = "Tumor_Volume"
#akt_whims_m$Tumor_Volume = akt_whims_m$Tumor_Volume_KG * 1000
akt_whims_m$Days = as.numeric(gsub("X","",akt_whims_m$Day))
akt_whims_m$Group = relevel(akt_whims_m$Group, ref="Vehicle")

plot_curve(akt_whims_m, "3whims_PI3Kinhibitors.pdf")
plot_curve(akt_whims_m[akt_whims_m$Sample == "WHIM16",], "whim16_PI3Kinhibitors.pdf")
plot_curve(akt_whims_m[akt_whims_m$Sample == "WHIM18",], "whim18_PI3Kinhibitors.pdf")
plot_curve(akt_whims_m[akt_whims_m$Sample == "WHIM20",], "whim20_PI3Kinhibitors.pdf")
plot_curve(whim14_m,"WHIM14_Neratinib.pdf")
plot_curve(whim8_m,"WHIM8_Herceptin.pdf")

akt_whims_m$Days = as.factor(akt_whims_m$Days) #Tukey's test require all variables to be in factors
anova_tukey(akt_whims_m[akt_whims_m$Sample == "WHIM16",])
anova_tukey(akt_whims_m[akt_whims_m$Sample == "WHIM18",])
anova_tukey(akt_whims_m[akt_whims_m$Sample == "WHIM20",])

# ### previous plotting (unused now)
# fn = paste(pd, "akt_whims_inhibitors.pdf", sep="_")
# p = ggplot(akt_whims_m, aes(as.numeric(Days),as.numeric(Tumor_Volume), color=as.factor(Group))) +
#   stat_summary(fun.data=mean_se, geom="pointrange")+ 
#   aes(colour = as.factor(Group)) + stat_summary(fun.y = mean, geom="line") +
#   labs(y=as.character("Tumor Volume (mm3)"), x=as.character("Days Post Treatment")) +
#   theme_bw() + theme_nogrid() + theme(legend.title=element_blank(), legend.position = "bottom") + 
#   theme(axis.title.y = element_text(size = 14, angle = 90)) + 
#   theme(axis.title.x = element_text(size = 14, angle = 0))
# p = p + facet_grid(.~Sample,drop=T,scales = "free", space = "free")
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# fn = paste(pd, "whim8_Herceptin.pdf", sep="_")
# p = ggplot(whim8_m, aes(as.numeric(Days),as.numeric(Tumor_Volume), color=as.factor(Group))) +
#   stat_summary(fun.data=mean_se, geom="pointrange")+ 
#   aes(colour = as.factor(Group)) + stat_summary(fun.y = mean, geom="line") +
#   labs(y=as.character("Tumor Volume (mm3)"), x=as.character("Days Post Treatment")) +
#   theme_bw() + theme_nogrid() + theme(legend.title=element_blank(), legend.position = "bottom") + 
#   theme(axis.title.y = element_text(size = 14, angle = 90)) + 
#   theme(axis.title.x = element_text(size = 14, angle = 0))
# p
# ggsave(file=fn, useDingbats=FALSE)
# 
# fn = paste(pd, "WHIM14_Herceptin.pdf", sep="_")
# p = ggplot(whim14_m, aes(as.numeric(Days),as.numeric(Tumor_Volume), color=as.factor(Group))) +
#   stat_summary(fun.data=mean_se, geom="pointrange")+ 
#   aes(colour = as.factor(Group)) + stat_summary(fun.y = mean, geom="line") +
#   labs(y=as.character("Tumor Volume (mm3)"), x=as.character("Days Post Treatment")) +
#   theme_bw() + theme_nogrid() + theme(legend.title=element_blank(), legend.position = "bottom") + 
#   theme(axis.title.y = element_text(size = 14, angle = 90)) + 
#   theme(axis.title.x = element_text(size = 14, angle = 0))
# p
# ggsave(file=fn, useDingbats=FALSE)