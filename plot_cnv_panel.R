require(reshape2)
library(scales)
library(gridExtra)
source("~/src/Scripts/R/reorder_chr.R")
source("~/src/Scripts/R/multiPlot.R")

usage <- function() {
# Returns the usage for this script
    return("\nUsage:\n\tRscript plot_cnv_panel.R merged_cnvs.hq")
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

plot_panel <- function(all_hq, gains, losses) {
# Plotting happens here
    print(head(all_hq))
    print(head(gains))
    p1 <- ggplot(all_hq, aes(x = pos, y = Primary1)) + stat_binhex(bins = 200) +
          scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          xlab("Primary1") + ylab("") +
          theme_bw() +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"), legend.position = "none",
                axis.title = element_text(size = 10, face = "bold"),
                plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"), strip.text = element_text(face = "bold", size = 10))
    p2 <- ggplot(all_hq, aes(x = pos, y = Primary2)) + stat_binhex(bins = 200) +
          scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          xlab("Primary2") + ylab("") +
          theme_bw() +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"), legend.position = "none",
                strip.text.x = element_blank(),
                axis.title = element_text(size = 10, face = "bold"),
                plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"), strip.background = element_blank())
    p3 <- ggplot(all_hq, aes(x = pos, y = Relapse1)) + stat_binhex(bins = 200) +
          scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          xlab("Relapse1") + ylab("") +
          theme_bw() +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"), legend.position = "none",
                strip.text.x = element_blank(),
                axis.title = element_text(size = 10, face = "bold"),
                plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"), strip.background = element_blank())
    p4 <- ggplot(all_hq, aes(x = pos, y = Relapse2)) + stat_binhex(bins = 200) +
          scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          xlab("Relapse2") + ylab("") +
          geom_segment(data = gains, aes(x = START, xend = END,
                                         y = Adjusted_CN1 - Adjusted_CN2,
                                         yend = Adjusted_CN1 - Adjusted_CN2,
                                         size = 1, group = START),
                       color = "forestgreen") +
          geom_segment(data = losses, aes(x = START, xend = END,
                                          y = Adjusted_CN1 - Adjusted_CN2,
                                          yend = Adjusted_CN1 - Adjusted_CN2,
                                          size = 1, group = START),
                       color = "forestgreen") +
          theme_bw() +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                axis.text.y = element_text(size = 10, face = "bold"), legend.position = "none",
                strip.text.x = element_blank(),
                axis.title = element_text(size = 10, face = "bold"),
                plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"), strip.background = element_blank())
    png("./all1_cnv_panel.png", width = 16, height = 10, units = "in", res = 300)
    grid.arrange(p1, p2, p3, p4, left = textGrob("Copy number difference",
                 gp = gpar(fontsize=18, fontface="bold"), rot = 90), nrow = 4, ncol = 1)
    dev.off()
}

melt_df <- function(df) {
# Melt the cnvs.hq values
    melt(df, id.vars=c("chr", "pos"),
        measure.vars = c("Primary1", "Primary2", "Relapse1", "Relapse2",
                         "Primary1_outlier"))
}

read_reviewed_segments <- function(gains, losses) {
    gains$chr <- gains[["CHR"]]
    gains$variable <- "Relapse2"
    gains <- reorder_chr(gains, "chr")
    #gains <- melt(gains, id.vars=c("chr", "nMarkers"), measure.vars = c("START", "END"))
    losses$chr <- losses[["CHR"]]
    losses$variable <- "Relapse2"
    losses <- reorder_chr(losses, "chr")
    return(list(gains, losses))
}

main <- function() {
# Everything starts here
    require(ggplot2)
    if(file.exists("all_cnvs.Robject")) {
        print("Loading R object")
        load("all_cnvs.Robject")
    } else {
        all_hq <- parse_args()
        all <- read_file(all_hq)
        #zscore of 3 or more
        all$Primary1_outlier <- factor(1, levels = c(1, 2, 3))
        all$Primary2_outlier <- factor(1, levels = c(1, 2, 3))
        all$Relapse1_outlier <- factor(1, levels = c(1, 2, 3))
        all$Relapse2_outlier <- factor(1, levels = c(1, 2, 3))
        all[which(all$Primary1 > (median(all$Primary1)+3*mad(all$Primary1))), "Primary1_outlier"] <- as.factor(2)
        all[which(all$Primary1 < (median(all$Primary1)-3*mad(all$Primary1))), "Primary1_outlier"] <- as.factor(3)
        all[which(all$Primary2 > (median(all$Primary2)+3*mad(all$Primary2))), "Primary2_outlier"] <- as.factor(2)
        all[which(all$Primary2 < (median(all$Primary2)-3*mad(all$Primary2))), "Primary2_outlier"] <- as.factor(3)
        all[which(all$Relapse1 > (median(all$Relapse1)+3*mad(all$Relapse1))), "Relapse1_outlier"] <- as.factor(2)
        all[which(all$Relapse1 < (median(all$Relapse1)-3*mad(all$Relapse1))), "Relapse1_outlier"] <- as.factor(3)
        all[which(all$Relapse2 > (median(all$Relapse2)+3*mad(all$Relapse2))), "Relapse2_outlier"] <- as.factor(2)
        all[which(all$Relapse2 < (median(all$Relapse2)-3*mad(all$Relapse2))), "Relapse2_outlier"] <- as.factor(3)
        all <- reorder_chr(all, "chr")
        #all <- melt_df(all)
        save(all, file = "all_cnvs.Robject")
    }
    #TODO - factor out
    gains <- read.table("cnaseq.cnvhmm.gains.merged.final.txt", head = T)
    losses <- read.table("cnaseq.cnvhmm.losses.merged.final.txt", head = T)
    gains_losses <- read_reviewed_segments(gains, losses)
    gains <- gains_losses[[1]]
    losses <- gains_losses[[2]]
    plot_panel(all, gains, losses)
}

main()
warnings()

