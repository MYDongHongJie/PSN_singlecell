#!/usr/bin/env Rscript
library("optparse")
library("oeRtools")
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character", default = NULL, help = "input file name",
        metavar = "character"
    ),
    make_option(
        c("-m", "--mark"),
        type = "character", default = NULL, help = "select Up, Total or Down",
        metavar = "character"
    ),
    make_option(
        c("-o", "--outpath"),
        type = "character", default = NULL, help = "outfile directory",
        metavar = "character"
    )
)

opt_parser <- OptionParser(
    option_list = option_list,
    epilogue = paste0(
        "Rscript top20_KEGG.r -i enrichment-kegg-Group1-vs-Group2-Down.txt ",
        "-m Down -o outdir/"
    )
)
opt <- parse_args(opt_parser)
if (is.null(opt$input) | is.null(opt$outpath) | is.null(opt$mark)) {
    print_help(opt_parser)
    stop("--input --outpath --mark must be supplied", call. = FALSE)
}
if (!file.exists(opt$input)) {
    print("input is not exists")
    q()
}
if (!file.exists(opt$outpath)) {
    dir.create(opt$outpath, recursive = T)
}
opt$outpath <- gsub("/$", "", opt$outpath)

library(ggplot2)
library(RColorBrewer)

load_table_from_str <- function(string) {
    if (string == "") return(data.frame())
    df <- read.table(
        text = string, header = T, sep = "\t",
        quote = "", check.names = F
    )
    return(df)
}

enrich_str <- readChar(opt$input, file.info(opt$input)$size)
enrich <- load_table_from_str(enrich_str)

enrich <- enrich[which(enrich["ListHits"] > 2), ]

groupname <- gsub(
    "\\.(txt|xls)$", "", gsub("^enrichment-kegg-", "", basename(opt$input))
)
if (grepl(
    "-Down$", groupname,
    ignore.case = F, perl = F, fixed = F, useBytes = F
)) {
    groupname <- gsub("-Down$", "(Down)", groupname)
}
if (grepl(
    "-Up$", groupname,
    ignore.case = F, perl = F, fixed = F, useBytes = F
)) {
    groupname <- gsub("-Up$", "(Up)", groupname)
}
if (grepl(
    "-Total$", groupname,
    ignore.case = F, perl = F, fixed = F, useBytes = F
)) {
    groupname <- gsub("-Total$", "(Total)", groupname)
}

if ("p.value" %in% colnames(enrich)) {
    top20 <- head(enrich[order(enrich[, "p.value"]), ], 20)
} else {
    top20 <- enrich
}

if (nrow(top20) == 0) {
    print("top20 items = 0, program exit!")
    sink(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".pdf"))
    sink()
    sink(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png"))
    sink()
    sink(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".xls"))
    sink()
    print("No diff genes")
    q()
}
write.table(
    top20[, c(1, 2, 3, 7, 9, 10)],
    paste0(opt$outpath, "/KEGG.top.", opt$mark, ".xls"),
    sep = "\t",
    quote = FALSE, col.names = TRUE, row.names = FALSE
)

# top20[, "Term"] <- sub(
#     "^path:", "", paste(top20[, 1], top20[, 2], sep = ": ")
# )

modify_name <- function(x, thres = 40) {
    x <- as.character(x)
    if(nchar(x) <= thres) {
        xx <- paste0(paste0(rep(" ", thres - nchar(x)), collapse = ""), x)
    } else {
        xx <- paste0(substr(x, 1, thres - 3) , "...")
    }
    return(xx)
}

top20[, "Term"] <- vapply(top20[, "Term"], modify_name, "a", USE.NAMES = F)
top20[, "Term"] <- make.unique(top20$Term)

xlab <- "Enrichment Score"
title <- paste0(groupname, "\n", "KEGG Enrichment top 20")
size.lab <- "Number"

p <- qplot(Enrichment_score, Term,
    data = top20, size = ListHits,
    colour = p.value, xlab = xlab, ylab = ""
) +
    facet_grid(vars(Classification_level1), scales = "free", space = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(
            color = "black", size = 1, fill = NA
    )) +
    theme(axis.line = element_blank(), strip.background = element_blank()) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradientn(colours = brewer.pal(11, "RdYlBu")) +
    labs(size = size.lab) +
    labs(colour = "p.value") +
    guides(size = guide_legend(order = 1))

ggsave(
    paste0(opt$outpath, "/KEGG.top.", opt$mark, ".pdf"),
    height = 8, width = 10, plot = p
)
ggsave(
    paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png"),
    type = "cairo-png", height = 8, width = 10, plot = p
)

print(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png(pdf) is OK"))
