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
    p1 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Primary1, color = Primary1_outlier)) +
          scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2, 4)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          xlab("Primary1") + ylab("") + theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red")) +
          theme(strip.text.y = element_text("Primary1"))+ theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"))
    p2 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Primary2, color = Primary2_outlier)) +
          scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2, 4)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) +
          xlab("Primary2") + ylab("") +theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red")) + theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"))
    p3 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Relapse1, color = Relapse1_outlier)) +
          scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2, 4)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) +
          xlab("Relapse1") + ylab("") +theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red"))+ theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"))
    p4 <- ggplot(all_hq) +  geom_point(aes(x = pos, y = Relapse2, color = Relapse2_outlier)) +
          scale_y_continuous(limits = c(-4, 4), breaks = c(-2, 0, 2, 4)) +
          facet_grid(.~chr, scales = "free", space = "free") +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(strip.background = element_blank(), strip.text.x = element_blank()) +
          xlab("Relapse2") + ylab("") +theme(legend.position="none") +
          scale_color_manual(values = c("grey", "red")) +
          geom_segment(data = gains, aes(x = START, xend = END, y = Adjusted_CN1 - Adjusted_CN2, yend = Adjusted_CN1 - Adjusted_CN2, group = START), color = "blue") +
          geom_segment(data = losses, aes(x = START, xend = END, y = Adjusted_CN1 - Adjusted_CN2, yend = Adjusted_CN1 - Adjusted_CN2, group = START), color = "blue") +
          theme(plot.margin = unit(c(-0.5, 0.5, 0.5, 0.5), "cm"))
    #p1 <- ggplot(primary1) + geom_point(aes(x = pos, y = value, color = value)) +
    #      facet_grid(variable~chr, space = "free", scales = "free") + ylim(-2, 5) +
    #      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    #      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    #      xlab("") + ylab("Copy number difference") +
    #      scale_color_gradient2(limits = c(-2, 2), high = muted("red"), low = muted("blue"), oob = squish) +
    #      geom_segment(data = gains, aes(x = START, xend = END, y = Adjusted_CN1 - Adjusted_CN2, yend = Adjusted_CN1 - Adjusted_CN2, group = START), color = "green") +
    #      geom_segment(data = losses, aes(x = START, xend = END, y = Adjusted_CN1 - Adjusted_CN2, yend = Adjusted_CN1 - Adjusted_CN2, group = START), color = "green")
    pdf("~/all1_cnv_panel.pdf", width = 14)
    grid.arrange(p1, p2, p3, p4, left = textGrob("Copy number difference",  gp = gpar(fontsize=18, fontface="bold"), rot = 90), nrow = 4, ncol = 1)
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
        #zscore of 2 or more
        all$Primary1_outlier <- all$Primary1>mean(all$Primary1)+2*sd(all$Primary1)|all$Primary1<mean(all$Primary1)-2*sd(all$Primary1)
        all$Primary2_outlier <- all$Primary2>mean(all$Primary2)+2*sd(all$Primary2)|all$Primary2<mean(all$Primary2)-2*sd(all$Primary2)
        all$Relapse1_outlier <- all$Relapse1>mean(all$Relapse1)+2*sd(all$Relapse1)|all$Relapse1<mean(all$Relapse1)-2*sd(all$Relapse1)
        all$Relapse2_outlier <- all$Relapse2>mean(all$Relapse2)+2*sd(all$Relapse2)|all$Relapse2<mean(all$Relapse2)-2*sd(all$Relapse2)
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


