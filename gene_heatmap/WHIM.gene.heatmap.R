# WHIM.gene.heatmap.R
# Kuan Huang and Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript WHIM.merger.R [-v] [-G genes.dat] clinical.fn rna.fn proteome.itraq.fn proteome.lfq.fn phospho.fn cnv.fn mutation.fn out.fn
#
# melt, combine, and normalize various proteogenomic data for use in heatmap figure
# -G: list of genes of interest.  Only these will be retained.

library("reshape2")
library("plyr")
#library("ggplot2")

# input files
setwd("/Users/khuang/Box Sync/PhD/proteogenomics/CPTAC_24WHIM/analysis/gene_heatmap")
source("/Users/khuang/bin/LIB_exp.R")
D="/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_24WHIM/data/"

clinical.fn=paste(sep="",D,"/WHIM_other/20141125_CPTAC_WHIM_samples-reviewed-by-Jeremy-Shun_simplified_cleaned.txt")
rna.fn=paste(sep="",D,"/WHIM_RNASeq/CPTAC_WHIMS_n40_rnaSeq_log_mdcnt_80pct_cleaned.txt_hugoified")
proteome.itraq.fn=paste(sep="",D,"/WHIM_proteome/proteome-ratio-norm/proteome-ratio-norm_exp_v2_collapsed.txt")
proteome.lfq.fn=paste(sep="",D,"/WHIM_proteome/label_free_all_WHIM_LFQ_Global/all_WHIM_LFQ_Global_minimum1_norm_nameadded_human_cleaned_collapsed.txt")
phospho.fn=paste(sep="",D,"/WHIM_phosphoproteome/phosphoproteome-ratio-norm/phosphoproteome-ratio-norm_exp_collapsed.txt")
cnv.fn=paste(sep="",D,"/WHIM_CNV/cnv_whims_7_8_2015.tsv_hugoified_normalized.tsv")
mutation.fn=paste(sep="",D,"/WHIM_variants/201509_variants/20150902_somatic_whims_combined_08_06_2015_matrix.txt")

options("width"=270) # useful for debugging
options(warn=2)

get_val_arg = function(args, flag, default) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = args[ix+1] } else { val = default }
    return(val)
}

get_bool_arg = function(args, flag) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = TRUE } else { val = FALSE }
    return(val)
}

parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    genes.fn = get_val_arg(args, "-G", NULL)

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.fn = args[length(args)];               args = args[-length(args)]
    mutation.fn = args[length(args)]; args = args[-length(args)]
    cnv.fn = args[length(args)]; args = args[-length(args)]
    phospho.fn = args[length(args)]; args = args[-length(args)]
    proteome.lfq.fn = args[length(args)];      args = args[-length(args)]
    proteome.itraq.fn = args[length(args)];      args = args[-length(args)]
    rna.fn = args[length(args)];               args = args[-length(args)]
    clinical.fn = args[length(args)];               args = args[-length(args)]

    # this is non-functional code.

    val = list( 'verbose'=verbose, 'out.fn' = out.fn, 'genes.fn'=genes.fn, 'phospho.fn' = phospho.fn, 
        'proteome.itraq.fn' = proteome.itraq.fn, 'proteome.lfq.fn' = proteome.lfq.fn, 
        'rna.fn' = rna.fn, 'clinical.fn' = clinical.fn, 'cnv.fn'=cnv.fn, 'mutation.fn'=mutation.fn)
    if (val$verbose) { print(val) }

    return (val)
}


# Normalization is done for each gene individually.
normalize.mad = function(x) {
    m = median(x, na.rm=TRUE)
    rowMAD=mad(x, na.rm=TRUE)
    return((x-m)/rowMAD)
}

normalize.logmean = function(x) {
    m = mean(x, na.rm=TRUE)
    print(m)
    print(x)
    return( log( x / m ) )
}

get.genes = function(genes.fn) {
    genes.dat = scan(genes.fn, what="character", sep="\n", quiet=TRUE)
    return(genes.dat)
}

parse.data = function(data.fn, value.name, genes) {
    data = read.table(data.fn, row.names=NULL, header=TRUE, sep = "\t", comment.char="#",stringsAsFactors=F)
    names(data)[1] = 'gene'
    data = melt(data, id.vars = c("gene"), variable.name="sample", value.name=value.name)

    genes.of.interest = genes
    data = data[data$gene %in% genes.of.interest,]
    return(data)
}

# clinical has per-sample information
# according to Kuan, get rid of PIK3CAdetails
parse.clinical = function(data.fn) {
    clinical = read.table(clinical.fn, row.names=NULL, header=TRUE, sep = "\t", comment.char="#",stringsAsFactors=F)
    # clinical is melted and cast again, so that it takes the form,
    #  sample       ER     HER2 Intrinsic.subtype PIK3CAdetails       PR
    #  1  WHIM2 negative negative             Basal            wt negative

    names(clinical)[1] = "type"
    clinical = melt(clinical, id.vars = c("type"), variable.name="sample")
    PIK = which(clinical$type == "PIK3CAdetails")
    clinical = clinical[-PIK,]
    clinical = dcast(clinical, sample ~ type)
    return(clinical)
}

parse.mutation = function(data.fn, genes) {
    data = parse.data(data.fn, 'mutation', genes)
    # Rename mutation 'wt' to NA according to Kuan email 9/3/15
    if (nrow(data)>0){
      is.wt = which(data$mutation == "wt")
      data[is.wt,'mutation'] = NA
    }
    return(data)
}

# We also normalize WHIM names by converting "WHIM02" to "WHIM2"
parse.cnv = function(data.fn, genes) {
    data = parse.data(data.fn, 'cnv', genes)
    data$sample = gsub("WHIM0", "WHIM", data$sample)
#    data = ddply(data, "gene", transform, cnv.mad = normalize.mad( cnv )) # normalize per gene
    return(data)
}

parse.rna = function(data.fn, genes) {
    data = parse.data(data.fn, 'rna', genes)
#    data = ddply(data, "gene", transform, rna.mad = normalize.mad( rna )) # normalize per gene
    return(data)
}

parse.itraq = function(data.fn, genes) {
    data = parse.data(data.fn, 'itraq', genes)
#    data = ddply(data, "gene", transform, itraq.mad = normalize.mad( itraq )) # normalize per gene
    return(data)
}

parse.lfq = function(data.fn, genes) {
    data = parse.data(data.fn, 'lfq', genes)
#    data = ddply(data, "gene", transform, lfq.mad = normalize.mad( lfq )) # normalize per gene
    return(data)
}

parse.phospho = function(data.fn, genes) {
    data = parse.data(data.fn, 'phospho', genes)
#    data = ddply(data, "gene", transform, phospho.mad = normalize.mad( phospho )) # normalize per gene
    return(data)
}

get.data = function(data) {
  #   sample   gene      itraq   itraq.mad        lfq    lfq.mad       rna    rna.mad      cnv cnv.logmean mutation       ER     HER2 Intrinsic.subtype       PR    phospho phospho.mad
  # 1 WHIM11   EGFR -3.7457486 -0.20715912 -2.7140330 -1.2388351 -4.886062 -1.2744750 1.789554  -0.7706783     <NA> positive negative            HER2-E negative -3.2510627 -0.24278237
  
  # ordering by intrinsic subtype
  data$Intrinsic.subtype = factor(data$Intrinsic.subtype, levels = c("Basal", "CLDN low", "HER2-E", "LumB", "LumA","Lymphoma"))
  #data$Intrinsic.subtype = factor(data$Intrinsic.subtype, levels = c("LumB", "CLDN low", "HER2-E", "Basal"))
  sample.order =  unique(data[order(data$Intrinsic.subtype),c("sample")])
  data$sample = factor(data$sample, levels=sample.order)
  
  #data.numeric = melt(data[,c("sample", "gene", "itraq.mad", "lfq.mad", "rna.mad", "phospho.mad", "cnv.logmean")], id.vars=c("sample", "gene"), variable.name="measurement")
  #data.numeric = melt(data[,c("sample", "gene", "itraq.mad", "lfq.mad", "rna.mad", "phospho.mad", "cnv.mad")], id.vars=c("sample", "gene"), variable.name="measurement")
  data.numeric = melt(data[,c("sample", "gene", "itraq", "rna", "phospho", "cnv")], id.vars=c("sample", "gene"), variable.name="measurement")
  data.subtype = melt(data[,c("sample", "gene", "Intrinsic.subtype")], id.vars=c("sample", "gene"), variable.name="measurement")
  
  # take out just the data we'll be using for data.status, and make the levels of ER, HER2, PR, and Intrinsic.subtype be the
  # same.  This is done to avoid a warning, "attributes are not identical across measure variables; they will be dropped" in the melt step.
  # see discussion here: http://stackoverflow.com/questions/25688897/reshape2-melt-warning-message
  data.s = data[,c("sample", "gene", "ER", "HER2", "Intrinsic.subtype", "PR")]
  all.levels = c(levels(data.s$ER), levels(data.s$Intrinsic.subtype))
  data.s[,c("ER", "HER2", "Intrinsic.subtype", "PR")] = lapply(data.s[,c("ER", "HER2", "Intrinsic.subtype", "PR")], factor, levels=all.levels)
  data.status = melt(data.s, id.vars=c("sample", "gene"), variable.name="measurement")
  # data.status = melt(data[,c("sample", "gene", "ER", "HER2", "Intrinsic.subtype", "PR")], id.vars=c("sample", "gene"), variable.name="measurement")  # this yields warning
  
  return(list('data.numeric'=data.numeric, 'data.subtype'=data.subtype, 'data.status'=data.status))
}
get.heatmap.scales = function() {
  RdBu.11 = brewer.pal(11, "RdBu") # red -> white -> blue
  RdBu.11[6]="gray90" # change midpoint color
  getPalette = colorRampPalette(RdBu.11)
  #fill.scale = scale_fill_gradientn(colours=rev(getPalette(100)), limits=c(-10,10), breaks=c(-10, 10), labels=c("Co-occurring", "Mutually Exclusive"), name="")
  fill.color.scale = scale_fill_gradientn(colours=rev(getPalette(100)), limits=c(-10,10), breaks=c(-10, 10),na.value="white"
                                          , labels=c("Low", "High"), name="Expression")
  
  colors = brewer.pal(5,"Dark2") #c("#1B9E77", "#D95F02", "#7570B3", "#E7298A") # based on Dark2
  color.names = c("missense", "frame_shift_ins", "nonsense", "splice_site", "complex_insertion")
  names(colors) = color.names
  mutation.color.scale = scale_color_manual(name="mutation", values=colors, na.value=NA)#, guide=TRUE)
  
  return(list('fill.color.scale'=fill.color.scale, 'mutation.color.scale'=mutation.color.scale))
}

plot.heatmap = function(data.numeric, data.subtype, plot.sample.name) {
  
  data.numeric$value_plot = data.numeric$value
  data.numeric$value_plot[data.numeric$value_plot>=10] = 10
  data.numeric$value_plot[data.numeric$value_plot<=-10] = -10
  
  p = ggplot(data=data.numeric)
  p = p + facet_grid(gene~.)
  
  p = p + geom_tile(aes(x=sample, y=measurement, fill=value_plot))
  p = p + geom_point(data=data.subtype, aes(x=sample, y=measurement, color=value), size=4, shape=16) #shape 15 is square, 16 is circle
  
  p = p + theme_bw() 
  p = p + coord_equal() # makes tiles square
  #     p = p + scale_y_discrete(breaks=rev(c("mutation", "cnv.logmean", "rna.mad", "itraq.mad", "lfq.mad", "phospho.mad")), 
  #                              labels=rev(c("Somatic Mutation", "CNV", "RNA-Seq", "MS ITRAQ proteome", "MS LFQ proteome", "MS ITRAQ Phosphoproteome")),
  #                              limits=rev(c("mutation", "cnv.logmean", "rna.mad", "itraq.mad", "lfq.mad", "phospho.mad")))
  
  p = p + scale_y_discrete(breaks=rev(c("Intrinsic.subtype", "cnv", "rna", "itraq", "phospho")), 
                           labels=rev(c("Intrinsic subtype", "CNV", "RNA-Seq", "MS ITRAQ proteome", "MS ITRAQ Phosphoproteome")),
                           limits=rev(c("Intrinsic.subtype", "cnv", "rna", "itraq", "phospho")))
  
  scales = get.heatmap.scales() 
  p = p + scales$fill.color.scale
  color.scale = get.clinical.scale()
  p = p + color.scale
  
  if (plot.sample.name) {
    p = p + theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=9))  
  } 
  p = p + theme( axis.text.y = element_text(hjust = 1, vjust=0.5, size=9))  
  p = p + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
  p = p + theme(panel.border = element_blank(), panel.background=element_blank())
  p = p + theme(axis.ticks = element_blank())
  p = p + theme(legend.position = "bottom")
  return(p)
}

get.clinical.scale = function() {
  # Set1 colors
  colors = c(NA, "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  # use Perou's intrinsic subtype colors instead
  colors = c(NA, "#101010", NA, "#636363", "red", "pink", "#ffeda0", "#deebf7", "#FFFF33", "#A65628", "#F781BF", "#999999") #positive is dark grey       
  # use better colors
  colors = c(NA, "#101010", NA, "#636363", "#CE2427","#EF5591","#FFFF33","#8FBCE5","#423996","#50A547") #positive is dark grey       
  
  color.names = c("wt","mut","negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB","LumA","Lymphoma")
  #color.names = c("negative", "positive", "Basal", "HER2-E", "CLDN low", "LumB")
  names(colors) = color.names
  clinical.color.scale = scale_color_manual(name="mutation", values=colors)
  
  return(clinical.color.scale)
}


plot.status = function(data.status) {
  color.scale = get.clinical.scale()
  sample.order = levels(data.status$sample)
  
  p = ggplot(data=data.status)
  
  p = p + geom_point(aes(x=sample, y=measurement, color=value), size=4, shape=16)
  
  p = p + theme_bw() 
  #    p = p + coord_equal() # makes tiles square
  
  p = p + scale_y_discrete(breaks=rev(c("Intrinsic.subtype", "ER", "HER2", "PR")),
                           labels=rev(c("Subtype", "ER", "HER2", "PR")),
                           limits=rev(c("Intrinsic.subtype", "ER", "HER2", "PR")))
  
  p = p + color.scale
  #    p = p + theme(legend.position="none")
  p = p + theme(axis.text.x = element_blank())
  #p = p + theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=6))    # this is for testing only
  p = p + theme(axis.text.y = element_text(hjust = 1, vjust=0.5, size=6, color="gray10"))    # this is for testing only
  p = p + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  p = p + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
  p = p + theme(panel.border = element_blank(), panel.background=element_blank())
  p = p + theme(axis.ticks = element_blank())
  p = p + theme(legend.position="top")
  
  return(p)
}

##### MAIN #####

genes = "MSLN"

mutation = parse.mutation(mutation.fn, genes)
clinical = parse.clinical(clinical.fn)
cnv = parse.cnv(cnv.fn, genes)
rna = parse.rna(rna.fn, genes)
itraq = parse.itraq(proteome.itraq.fn, genes)
lfq = parse.lfq(proteome.lfq.fn, genes)
phospho = parse.phospho(phospho.fn, genes)

# note, LQF has a subset of samples that ITRAQ does
#data = merge(itraq, lfq, by=c("sample", "gene"), all=TRUE)
#data = merge(data, rna, by=c("sample", "gene"), all=TRUE)
data = merge(itraq, rna, by=c("sample", "gene"), all=TRUE)
data = merge(data, cnv, by=c("sample", "gene"), all=TRUE)
data = merge(data, mutation, by=c("sample", "gene"), all=TRUE)
data = merge(data, clinical, by="sample", all=TRUE)
data = merge(data, phospho, by=c("sample", "gene"), all=FALSE)  # we want only those which have phospho data - according to Kuan.

data = get.data(data)
ph = plot.heatmap(data$data.numeric, data$data.subtype, T)
ph
fn = paste(pd, genes,"24_PDX_omics.pdf", sep="_")
ggsave(file=fn, width=6, useDingbats=FALSE)
# con = file(out.fn, open="wt")
# timestamp = paste0("# Created ", Sys.time())
# writeLines(paste("# Normalized WHIM data", timestamp, sep="\n"), con)
# write.table(data, con, sep="\t", quote=FALSE, row.names=FALSE)
# close(con)
# cat(sprintf("    Saved to %s\n", out.fn))
