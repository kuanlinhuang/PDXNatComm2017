# compareDiffExp.R by Kuan Huang @ WashU 201508
# compare differentially expressed genes found in human and WHIM

# do it for TP53 and PIK3CA for now; there is not a ESR1 mutant in human
# limma seem to work better for human

# mis
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression")
source("/Users/khuang/bin/LIB_exp.R")

##### proteome #####
### Xeno ###
xeno_PIK3CA_pro = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression/figures/2015-08-11/2015-08-11_KH_ITRAQ_proteome_PIK3CA_vs_wt_diff_exp.txt")
xeno_TP53_pro = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression/figures/2015-08-11/2015-08-11_KH_ITRAQ_proteome_TP53_vs_wt_diff_exp.txt")
### Human ###
human_PIK3CA_pro = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/analysis/differential_expression/figures/2015-08-13/2015-08-13_KH_ITRAQ_human_proteome_PIK3CA_vs_wt_diff_exp.txt")
human_TP53_pro = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/analysis/differential_expression/figures/2015-08-13/2015-08-13_KH_ITRAQ_human_proteome_TP53_vs_wt_diff_exp.txt")
### merge ###
merge_PIK3CA_pro = merge(xeno_PIK3CA_pro, human_PIK3CA_pro, by = "row.names")
merge_TP53_pro = merge(xeno_TP53_pro, human_TP53_pro, by = "row.names")
colnames(merge_PIK3CA_pro)[1] = "gene"
colnames(merge_TP53_pro)[1] = "gene"

##### phosphoprotoeme #####
### Xeno ###
xeno_PIK3CA_pho = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression/figures/2015-08-08/2015-08-08_KH_ITRAQ_phosphoproteome_PIK3CA_vs_wt_diff_exp.txt")
xeno_TP53_pho = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/differential_expression/figures/2015-08-08/2015-08-08_KH_ITRAQ_phosphoproteome_TP53_vs_wt_diff_exp.txt")
### Human ###
human_PIK3CA_pho = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/analysis/differential_expression/figures/2015-08-13/2015-08-13_KH_ITRAQ_human_phosphoproteome_PIK3CA_vs_wt_diff_exp.txt")
human_TP53_pho = read.table(header = T, sep = "\t", row.names=1, file = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/analysis/differential_expression/figures/2015-08-13/2015-08-13_KH_ITRAQ_human_phosphoproteome_TP53_vs_wt_diff_exp.txt")
### merge ###
merge_PIK3CA_pho = merge(xeno_PIK3CA_pho, human_PIK3CA_pho, by = "row.names")
merge_TP53_pho = merge(xeno_TP53_pho, human_TP53_pho, by = "row.names")
colnames(merge_PIK3CA_pho)[1] = "gene"
colnames(merge_TP53_pho)[1] = "gene"

### make comparison plot based on limma FDR ###
# merges = list('merge_PIK3CA_pro' = merge_PIK3CA_pro, "merge_TP53_pro"=merge_TP53_pro, "merge_PIK3CA_pho" = merge_PIK3CA_pho, "merge_TP53_pho" = merge_TP53_pho)
# for (merged in merges){
#   mergedName = deparse(substitute(merged))
#   maxlim = max(-log10(as.numeric(as.character(merged[,"adj.P.Val.y"]))))
#   colnames(merged)[1] = "gene"
#   fn = paste(pd,mergedName,'limma_diff_exp_QQplot.pdf', sep="_")
#   p = ggplot(data = merged, aes(x=-log10(as.numeric(as.character(adj.P.Val.x))) , y=-log10(as.numeric(as.character(adj.P.Val.y))), label=gene)) 
#   p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
#   p = p + xlim(0,maxlim) + ylim(0,maxlim)
#   p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(adj.P.Val.x)))*log10(as.numeric(as.character(adj.P.Val.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
#   p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
#   p
#   ggsave(file=fn)
# }

maxlim = max(-log10(as.numeric(as.character(merge_PIK3CA_pro[,"adj.P.Val.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_PIK3CA_pro",'limma_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_PIK3CA_pro, aes(x=-log10(as.numeric(as.character(adj.P.Val.x))) , y=-log10(as.numeric(as.character(adj.P.Val.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(adj.P.Val.x)))*log10(as.numeric(as.character(adj.P.Val.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)

maxlim = max(-log10(as.numeric(as.character(merge_TP53_pro[,"adj.P.Val.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_TP53_pro",'limma_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_TP53_pro, aes(x=-log10(as.numeric(as.character(adj.P.Val.x))) , y=-log10(as.numeric(as.character(adj.P.Val.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(adj.P.Val.x)))*log10(as.numeric(as.character(adj.P.Val.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)

maxlim = max(-log10(as.numeric(as.character(merge_PIK3CA_pho[,"adj.P.Val.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_PIK3CA_pho",'limma_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_PIK3CA_pho, aes(x=-log10(as.numeric(as.character(adj.P.Val.x))) , y=-log10(as.numeric(as.character(adj.P.Val.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(adj.P.Val.x)))*log10(as.numeric(as.character(adj.P.Val.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)

maxlim = max(-log10(as.numeric(as.character(merge_TP53_pho[,"adj.P.Val.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_TP53_pho",'limma_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_TP53_pho, aes(x=-log10(as.numeric(as.character(adj.P.Val.x))) , y=-log10(as.numeric(as.character(adj.P.Val.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(adj.P.Val.x)))*log10(as.numeric(as.character(adj.P.Val.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)


### t-test ###
maxlim = max(-log10(as.numeric(as.character(merge_PIK3CA_pro[,"t_test_fdr.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_PIK3CA_pro",'ttest_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_PIK3CA_pro, aes(x=-log10(as.numeric(as.character(t_test_fdr.x))) , y=-log10(as.numeric(as.character(t_test_fdr.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(t_test_fdr.x)))*log10(as.numeric(as.character(t_test_fdr.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)

maxlim = max(-log10(as.numeric(as.character(merge_TP53_pro[,"t_test_fdr.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_TP53_pro",'ttest_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_TP53_pro, aes(x=-log10(as.numeric(as.character(t_test_fdr.x))) , y=-log10(as.numeric(as.character(t_test_fdr.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(t_test_fdr.x)))*log10(as.numeric(as.character(t_test_fdr.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)

maxlim = max(-log10(as.numeric(as.character(merge_PIK3CA_pho[,"t_test_fdr.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_PIK3CA_pho",'ttest_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_PIK3CA_pho, aes(x=-log10(as.numeric(as.character(t_test_fdr.x))) , y=-log10(as.numeric(as.character(t_test_fdr.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(t_test_fdr.x)))*log10(as.numeric(as.character(t_test_fdr.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)

maxlim = max(-log10(as.numeric(as.character(merge_TP53_pho[,"t_test_fdr.y"]))))
colnames(merged)[1] = "gene"
fn = paste(pd,"merge_TP53_pho",'ttest_diff_exp_QQplot.pdf', sep="_")
p = ggplot(data = merge_TP53_pho, aes(x=-log10(as.numeric(as.character(t_test_fdr.x))) , y=-log10(as.numeric(as.character(t_test_fdr.y))), label=gene)) 
p = p + geom_point(alpha=0.7) + xlab("xenograft -log(FDR)") + ylab("human -log(FDR)") + theme_bw()
p = p + xlim(0,maxlim) + ylim(0,maxlim)
p = p + geom_text(aes(label=ifelse(log10(as.numeric(as.character(t_test_fdr.x)))*log10(as.numeric(as.character(t_test_fdr.y)))>5,gene,"")),hjust=-0.05,just=0,size=5)
p = p + theme(axis.title = element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
p
ggsave(file=fn)
