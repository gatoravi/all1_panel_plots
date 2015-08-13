require(reshape2)
library(scales)
library(gridExtra)
source("~/src/Scripts/R/reorder_chr.R")
source("~/src/Scripts/R/multiPlot.R")

usage <- function() {
# Returns the usage for this script
    return("\nUsage:\n\tRscript plot_cnv_panel.R all_loh.infile")
}

parse_args <-function() {
# Parses the different inputs to this script
    args <- commandArgs(trailingOnly = T)
    if(length(args) < 1) {
        stop(usage())
    }
    p1_hq <- args[1]
    return(c(p1_hq))
}

read_file <- function(f1) {
# Read the file into an object
    return(read.table(f1, sep = "\t", head = T))
}

plot_panel <- function(all_hq, primary1_segments, primary2_segments,
                       relapse1_segments, relapse2_segments) {
# Plotting happens here
    print(head(all_hq))
    print(head(primary1_segments))
    p1 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Primary1, color = Primary1_outlier, size = 1)) +
          scale_size_continuous(range = c(1,1)) + scale_y_continuous(limits = c(0, 100)) +
          #geom_segment(data = primary1_segments, aes(x = chr_start, xend = chr_stop, y = seg_mean * 100, yend = seg_mean * 100, group = chr_start), color = "blue") +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          xlab("Primary1") + ylab("") + theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red")) +
          theme(strip.text.y = element_text("Primary1"))+ theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"))
    p2 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Primary2, color = Primary2_outlier, size = 1)) +
          scale_size_continuous(range = c(1,1)) + scale_y_continuous(limits = c(0, 100)) +
          #geom_segment(data = primary2_segments, aes(x = chr_start, xend = chr_stop, y = seg_mean * 100, yend = seg_mean * 100, group = chr_start), color = "blue") +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) +
          xlab("Primary2") + ylab("") +theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red")) + theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"))
    p3 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Relapse1, color = Relapse1_outlier, size = 1)) +
          scale_size_continuous(range = c(1,1)) + scale_y_continuous(limits = c(0, 100)) +
          #geom_segment(data = relapse1_segments, aes(x = chr_start, xend = chr_stop, y = seg_mean * 100, yend = seg_mean * 100, group = chr_start), color = "blue") +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) +
          xlab("Relapse1") + ylab("") +theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red"))+ theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"))
    p4 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Relapse2, color = Relapse2_outlier, size = 1)) +
          scale_size_continuous(range = c(1,1)) + scale_y_continuous(limits = c(0, 100)) +
          geom_segment(data = relapse2_segments, aes(x = chr_start, xend = chr_stop, y = seg_mean * 100, yend = seg_mean * 100, group = chr_start, size = 3), color = "black") +
          geom_segment(data = relapse2_segments, aes(x = chr_start, xend = chr_stop, y = 100 - seg_mean * 100, yend = 100 - seg_mean * 100, group = chr_start, size = 3), color = "black") +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) +
          xlab("Relapse2") + ylab("") +theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red")) +
          theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"))
    png("./all1_loh_panel.png", width = 16, height = 10, units = "in", res = 300)
    grid.arrange(p1, p2, p3, p4, left = textGrob("Variant allele frequency",  gp = gpar(fontsize=18, fontface="bold"), rot = 90), nrow = 4, ncol = 1)
    dev.off()
}


main <- function() {
# Everything starts here
    require(ggplot2)
    loh_robject = "all_loh.robject"
    if(file.exists(loh_robject)) {
        print("Loading R object")
        load(loh_robject)
    } else {
        all_hq <- parse_args()
        all <- read_file(all_hq)
        #zscore of 2 or more
        all$Primary1_outlier <- all$Primary1>mean(all$Primary1, na.rm = T)+3*sd(all$Primary1, na.rm = T)|all$Primary1<mean(all$Primary1, na.rm = T)-3*sd(all$Primary1, na.rm = T)
        all$Primary2_outlier <- all$Primary2>mean(all$Primary2, na.rm = T)+3*sd(all$Primary2, na.rm = T)|all$Primary2<mean(all$Primary2, na.rm = T)-3*sd(all$Primary2, na.rm = T)
        all$Relapse1_outlier <- all$Relapse1>mean(all$Relapse1, na.rm = T)+3*sd(all$Relapse1, na.rm = T)|all$Relapse1<mean(all$Relapse1, na.rm = T)-3*sd(all$Relapse1, na.rm = T)
        all$Relapse2_outlier <- all$Relapse2>mean(all$Relapse2, na.rm = T)+3*sd(all$Relapse2, na.rm = T)|all$Relapse2<mean(all$Relapse2, na.rm = T)-3*sd(all$Relapse2, na.rm = T)
        print(head(all))
        all <- reorder_chr(all, "chr")
        save(all, file = loh_robject)
    }
    primary1_segments <- read_segments("primary1_segments.cbs")
    primary2_segments <- read_segments("primary2_segments.cbs")
    relapse1_segments <- read_segments("relapse1_segments.cbs")
    relapse2_segments <- read_segments("relapse2_segments.cbs")
    plot_panel(all, primary1_segments, primary2_segments,
               relapse1_segments, relapse2_segments)
}

read_segments <- function(segments_file) {
    segments <- read.table(segments_file, head = T)
    segments <- segments[segments$seg_mean > 0.05, ]
    segments$chr <- segments$chrom
    segments
}

main()


