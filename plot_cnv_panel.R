require(reshape2)
source("~/src/Scripts/R/reorder_chr.R")

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

plot_panel <- function(primary1) {
# Plotting happens here
    print(summary(primary1))
    p1 <- ggplot(primary1) + geom_point(aes(x = POS, y = value)) +
          facet_grid(variable~CHR) + ylim(-2, 5) + theme(aspect.ratio = 1) +
          theme(axis.ticks = element_blank(), axis.text = element_blank()) +
          xlab("Position") + ylab("Copy number")
#    p1 <- p1 + geom_point(aes(x = POS, y = DIFF2)) +
#          facet_grid(variable~CHR) + ylim(-2, 5) + theme(aspect.ratio = 1) +
#          theme(axis.ticks = element_blank(), axis.text = element_blank())
    print(p1)
    ggsave("~/test_p2.pdf")
}

melt_df <- function(df) {
# Melt the cnvs.hq values
    melt(df, id.vars=c("CHR", "POS"))
}

main <- function() {
# Everything starts here
    require(ggplot2)
    p1_hq <- parse_args()
    primary1 <- read_file(p1_hq)
    primary1 <- reorder_chr(primary1, "CHR")
    primary1 <- melt_df(primary1)
    print(head(primary1))
    plot_panel(primary1)
}

main()


