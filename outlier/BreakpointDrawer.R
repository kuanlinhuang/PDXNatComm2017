# Matthew Wyczalkowski
# mwyczalk@genome.wustl.edu
# The Genome Institute
#
# Usage: Rscript BreakpointDrawer.R [-v] [-a chrom.annotation.ggp] [-A virus.annotation.ggp] [-H histogram.ggp] 
#               [-t title] [-S] [-h height] [-w width] [-L] [-b] discordant.ggp chrom.depth.ggp virus.depth.ggp out.pdf
#
# Read in various GGP objects and write to PDF
# in cases where there are no features to annotate (e.g., no genes in range of interest) policy is that 
# no annotation.ggp file is passed.  
#
# Input arguments:
# -v: verbose
# -S: make simplified structure plot, with only discordant reads and read depth
# -L: do not print axis labels
# -b: make big axis text
#

options("width"=180) # useful for debugging
suppressPackageStartupMessages(library("ggplot2"))
library('grid')
library('gridExtra', quietly=TRUE)
library('gridBase', quietly=TRUE)

#library('plyr', quietly=TRUE)  
#library('zoo', quietly=TRUE)
#library("reshape2")

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
    plot.structure = get_bool_arg(args, "-S")
    chrom.annotation.fn = get_val_arg(args, "-a", NULL)
    virus.annotation.fn = get_val_arg(args, "-A", NULL)
    histogram.fn = get_val_arg(args, "-H", NULL)
    title = get_val_arg(args, "-t", NULL)
    height = as.numeric(get_val_arg(args, "-h", 8))
    width = as.numeric(get_val_arg(args, "-w", 8))
    no.label = get_bool_arg(args, "-L")
    big.font = get_bool_arg(args, "-b")

    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    out.fn = args[length(args)];             args = args[-length(args)]
    virus.depth.fn = args[length(args)];             args = args[-length(args)]
    chrom.depth.fn = args[length(args)];             args = args[-length(args)]
    discordant.fn = args[length(args)];      args = args[-length(args)]

    val = list('verbose'=verbose, 'discordant.fn'= discordant.fn, 'histogram.fn'= histogram.fn, 'chrom.depth.fn'= chrom.depth.fn,
               'virus.depth.fn'= virus.depth.fn, 'out.fn'= out.fn, 'chrom.annotation.fn'=chrom.annotation.fn,
               'virus.annotation.fn'=virus.annotation.fn, 'title'=title, 'plot.structure'=plot.structure, 'height'=height, 'width'=width,
               'no.label'=no.label, 'big.font'=big.font)
    if (val$verbose) { print(val) }

    return (val)
}

# Precise alignment of panels is based on work here: 
    # /Users/mwyczalk/Data/Virus/Virus_2013.9a/RSEM-Exon/RPKM-Scatter/src/RPKM_scatter_plotter.R
# Read more about approach here: http://stackoverflow.com/questions/17370853/align-ggplot2-plots-vertically

# from http://zevross.com/blog/2014/11/20/under-the-hood-of-ggplot2-graphics-in-r/
    #The ggplot_build function outputs a list of data frames (one for each layer)
    #and a panel object with information about axis limits among other things. The
    #ggplot_gtable function, which takes the ggplot_build object as input, builds
    #all grid graphical objects (known as “grobs”) necessary for displaying the
    #plot. You can manipulate the output from ggplot_build.

# The general idea is to set the ranges of .gb objects first, then the widths/heights of .gt objects
# details of .gb and .gt objects can be obtained with print(str(foo))
# .ggp and .gt objects can be visualized with grid.draw(foo)

assemble_ggp = function(discordant.ggp, chrom.depth.ggp, virus.depth.ggp, histogram.ggp, chrom.annotation.ggp, virus.annotation.ggp, title.ggp) {

    # the following lines will demonstrate that depth, annotation, and discordant plots are aligned precisely 
    # the numbers below are for BA-4077, but can be generated for all figures by reading .gb$x.range, .gb$yrange values
        # discordant.gb: x.range: 1 7905  y.range: 68633064 68792035
#    virus.marks = c(100, 4000, 7000)
#    chrom.marks = c(68633064 + 1000, 68633064 + (68792035 - 68633064) / 2, 68792035 - 1000)
#    virus.annotation.ggp = virus.annotation.ggp + geom_vline(xintercept=virus.marks)
#    virus.depth.ggp = virus.depth.ggp + geom_vline(xintercept=virus.marks)
#    discordant.ggp = discordant.ggp + geom_hline(yintercept=virus.marks)
#    chrom.annotation.ggp = chrom.annotation.ggp + geom_vline(xintercept=chrom.marks)
#    chrom.depth.ggp = chrom.depth.ggp + geom_vline(xintercept=chrom.marks)
#    discordant.ggp = discordant.ggp + geom_vline(xintercept=chrom.marks)

    discordant.gb = ggplot_build(discordant.ggp)
    chrom.depth.gb = ggplot_build(chrom.depth.ggp)
    virus.depth.gb = ggplot_build(virus.depth.ggp)

    chrom.depth.gb$panel$ranges[[1]]$y.range = discordant.gb$panel$ranges[[1]]$y.range   
    virus.depth.gb$panel$ranges[[1]]$x.range = discordant.gb$panel$ranges[[1]]$x.range

    discordant.gt = ggplot_gtable(discordant.gb)
    chrom.depth.gt = ggplot_gtable(chrom.depth.gb)
    virus.depth.gt = ggplot_gtable(virus.depth.gb)

    chrom.depth.gt$heights = discordant.gt$heights
    virus.depth.gt$widths = discordant.gt$widths

    if (!is.null(chrom.annotation.ggp)) {
#    # this fails in cases where chrom.annotation.ggp is a blank grobRect, such as when no genes in range
#    # for this reason, for now, implement policy of no chrom.annotation.ggp in cases where no features to annotate.
        chrom.annotation.gb = ggplot_build(chrom.annotation.ggp)
        chrom.annotation.gb$panel$ranges[[1]]$y.range = chrom.depth.gb$panel$ranges[[1]]$y.range
        chrom.annotation.gt = ggplot_gtable(chrom.annotation.gb)     
        chrom.annotation.gt$heights= chrom.depth.gt$heights
        chrom.gt = arrangeGrob(chrom.depth.gt, chrom.annotation.gt, widths=c(0.5,0.5), ncol=2, nrow=1)
    } else 
        chrom.gt = chrom.depth.gt

    if (!is.null(virus.annotation.ggp)) {
        virus.annotation.gb = ggplot_build(virus.annotation.ggp)
        virus.annotation.gb$panel$ranges[[1]]$x.range = virus.depth.gb$panel$ranges[[1]]$x.range
        virus.annotation.gt = ggplot_gtable(virus.annotation.gb)
        virus.annotation.gt$widths = virus.depth.gt$widths
        virus.gt = arrangeGrob(virus.annotation.gt, virus.depth.gt, heights=c(0.5,0.5), nrow=2)
    } else
        virus.gt = virus.depth.gt

    main.grob = arrangeGrob(chrom.gt, discordant.gt, histogram.ggp, virus.gt, widths=c(0.3,0.7), heights=c(0.7,0.3), ncol=2, nrow=2)

    annotated.grob = grid.arrange(title.ggp, main.grob, heights = c(0.025, 0.975), ncol=1, nrow=2, newpage=FALSE)
}

# This figure includes only discordant plot and read depth.
assemble_structure_ggp = function(discordant.ggp, chrom.depth.ggp, virus.depth.ggp, histogram.ggp, title.ggp, no.label=FALSE) {

    # get rid of axis titles because we'll draw them ourselves.
    chrom.depth.ggp = chrom.depth.ggp + theme(axis.title = element_blank())
    virus.depth.ggp = virus.depth.ggp + theme(axis.title = element_blank())

    discordant.gb = ggplot_build(discordant.ggp)
    chrom.depth.gb = ggplot_build(chrom.depth.ggp)
    virus.depth.gb = ggplot_build(virus.depth.ggp)

    chrom.depth.gb$panel$ranges[[1]]$y.range = discordant.gb$panel$ranges[[1]]$y.range
    virus.depth.gb$panel$ranges[[1]]$x.range = discordant.gb$panel$ranges[[1]]$x.range

    discordant.gt = ggplot_gtable(discordant.gb)
    chrom.depth.gt = ggplot_gtable(chrom.depth.gb)
    virus.depth.gt = ggplot_gtable(virus.depth.gb)

    chrom.depth.gt$heights = discordant.gt$heights
    virus.depth.gt$widths = discordant.gt$widths

    # make sure histogram dimensions match
    if (is.null(histogram.ggp)) {
        histogram.ggp = grid.rect(gp=gpar(col="white"))
    } else {
        histogram.gt = ggplot_gtable(ggplot_build(histogram.ggp))
        histogram.gt$heights = virus.depth.gt$heights
        histogram.gt$widths = chrom.depth.gt$widths
    }


    main.grob = arrangeGrob(chrom.depth.gt, discordant.gt, histogram.gt, virus.depth.gt, widths=c(0.3,0.7), heights=c(0.7,0.3), ncol=2, nrow=2)

    if (!is.null(title.ggp)) {
        annotated.grob = grid.arrange(title.ggp, main.grob, heights = c(0.025, 0.975), ncol=1, nrow=2, newpage=FALSE)
    } else {
        grid.draw(main.grob)
    }

    if (!no.label) { 
        # hard positioning of labels.  Ick.
        grid.text("Chromosome Copy Number", x = unit(0.160, "npc"), y = unit(0.285, "npc"), gp=gpar(fontsize=10))
        grid.text("Virus Copy Number", x = unit(0.295, "npc"), y = unit(0.150, "npc"), gp=gpar(fontsize=10), rot=-90)
    }
}


args = parse_args()
discordant.ggp = readRDS(args$discordant.fn)   # http://www.fromthebottomoftheheap.net/2012/04/01/saving-and-loading-r-objects/

chrom.depth.ggp = readRDS(args$chrom.depth.fn)
virus.depth.ggp = readRDS(args$virus.depth.fn)

chrom.annotation.ggp = NULL
virus.annotation.ggp = NULL
histogram.ggp = NULL
title.ggp = NULL

if (!is.null(args$chrom.annotation.fn))
    chrom.annotation.ggp = readRDS(args$chrom.annotation.fn)
if (!is.null(args$virus.annotation.fn))
    virus.annotation.ggp = readRDS(args$virus.annotation.fn)
if (!is.null(args$histogram.fn)) 
    histogram.ggp = readRDS(args$histogram.fn)

title_text = args$title
if (!is.null(title_text))  {
    title.ggp = textGrob(title_text)
} else {
    title.ggp = NULL
}

chrom.depth.ggp = chrom.depth.ggp + coord_flip() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab(sprintf("Chr %s Position", args$chrname))
virus.depth.ggp = virus.depth.ggp+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y=element_text(angle=-90, hjust=0.5)) + xlab(sprintf("%s Position", args$virname))

if (args$no.label) {
    histogram.ggp = histogram.ggp + theme(axis.title = element_blank(), legend.position="none")
}

# we're expanding meaning of big.font to do many customizaitons for -special figures for manuscript
if (args$big.font) {
    target.size = 12
    no.grid = theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
    resize.text = theme(axis.text = element_text(size=target.size))
    histogram.ggp = histogram.ggp + no.grid + resize.text
    chrom.depth.ggp = chrom.depth.ggp + no.grid + resize.text
    virus.depth.ggp = virus.depth.ggp + no.grid + resize.text
    discordant.ggp = discordant.ggp + no.grid + theme(axis.text = element_text(size=6, color="gray50"))

    histogram.ggp = histogram.ggp + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

cat(sprintf("Saving to %s\n", args$out.fn))
pdf(file=args$out.fn, width=args$width, height=args$height, useDingbats=FALSE)


if (args$plot.structure) {
    assemble_structure_ggp(discordant.ggp, chrom.depth.ggp, virus.depth.ggp, histogram.ggp, title.ggp, args$no.label)
} else {
    assemble_ggp(discordant.ggp, chrom.depth.ggp, virus.depth.ggp, histogram.ggp, chrom.annotation.ggp, virus.annotation.ggp, title.ggp, args$no.label)
}
