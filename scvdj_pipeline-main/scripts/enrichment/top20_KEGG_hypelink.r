#!/usr/bin/env Rscript
library("optparse")
library(scales)

# library("oeRtools")
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character", default = NULL,
        help = "input file name", metavar = "character"
    ),
    make_option(
        c("-m", "--mark"),
        type = "character", default = NULL,
        help = "select Up, Total or Down", metavar = "character"
    ),
    make_option(
        c("-o", "--outpath"),
        type = "character", default = NULL,
        help = "outfile directory", metavar = "character"
    )
)

opt_parser <- OptionParser(
    option_list = option_list,
    epilogue = paste0(
        "Rscript top20_KEGG.r ",
        "-i enrichment-kegg-Group1-vs-Group2-Down.txt -m Down -o outdir/"
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
if (
    grepl(
        "-Down$", groupname,
        ignore.case = F, perl = F, fixed = F,
        useBytes = F
    )) {
    groupname <- gsub("-Down$", "(Down)", groupname)
}
if (
    grepl(
        "-Up$", groupname,
        ignore.case = F, perl = F, fixed = F, useBytes = F
    )) {
    groupname <- gsub("-Up$", "(Up)", groupname)
}
if (
    grepl(
        "-Total$", groupname,
        ignore.case = F, perl = F, fixed = F, useBytes = F
    )) {
    groupname <- gsub("-Total$", "(Total)", groupname)
}

if ("p-value" %in% colnames(enrich)) {
    top20 <- head(enrich[order(enrich[, "p-value"]), ], 20)
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
    top20[, c(1, 2, 3, 4, 5, 9, 11, 12)],
    paste0(opt$outpath, "/KEGG.top.", opt$mark, ".xls"),
    sep = "\t",
    quote = FALSE, col.names = TRUE, row.names = FALSE
)

# top20[, "Term"] <- sub(
#     "^path:", "", paste(top20[, 1], top20[, 2], sep = ": ")
# )

modify_name <- function(x, thres = 3) {
    # par()
    x <- as.character(x)
    sw <- strwidth(x, family = "sans", units = "inches")
    bw <- strwidth(" ", family = "sans", units = "inches")
    dw <- strwidth(".", family = "sans", units = "inches")
    if (sw <= thres) {
        xx <- paste0(paste0(rep(" ", (thres - sw) / bw), collapse = ""), x)
    } else {
        xs <- substring(x, seq(1, nchar(x), 1), seq(1, nchar(x), 1))
        for (i in seq_len(length(xs))) {
            if (
                sum(
                    strwidth(xs[1:i], family = "sans", units = "inches")
                ) > thres - 3 * dw) {
                xx <- paste0(substr(x, 1, i - 1), "...")
                break
            }
        }
    }
    # dev.off()
    return(xx)
}

pdf(NULL)
top20[, "Term"] <- vapply(top20[, "Term"], modify_name, "a", USE.NAMES = F)
dev.off()

top20[, "Term"] <- make.unique(top20$Term)

top20$Classification_level1 <- abbreviate(
    top20$Classification_level1,
    use.classes = F, dot = T, minlength = 5
)

xlab <- "Enrichment Score"
title <- paste0(groupname, "\n", "KEGG Enrichment top 20")
size.lab <- "Number"
colnames(top20)[9] <- "pvalue"
top20 <- top20[order(top20$Enrichment_score, decreasing = F), ]
top20$Term <- factor(top20$Term, levels = top20$Term)

p <- qplot(Enrichment_score, Term,
    data = top20, size = ListHits,
    colour = pvalue, xlab = xlab, ylab = ""
) +
    facet_grid(vars(Classification_level1), scales = "free", space = "free_y") +
    theme_bw() +
    # theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(
        color = "black", size = 1, fill = NA
    )) +
    theme(
        text = element_text(family = "sans"),
        plot.title = element_text(face = "bold"),
        axis.line = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 7),
        strip.background = element_rect(color = "#e0e0e0", fill = "#e0e0e0")
    ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradientn(
        colours = brewer.pal(11, "RdYlBu"),
        labels = scientific
    ) +
    labs(size = size.lab) +
    labs(colour = "pvalue") +
    guides(size = guide_legend(order = 1))

ggsave(
    paste0(opt$outpath, "/KEGG.top.", opt$mark, ".pdf"),
    height = 7 - (20 - nrow(top20)) * 0.2, width = 8, plot = p
)
ggsave(
    paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png"),
    type = "cairo-png",
    height = 7 - (20 - nrow(top20)) * 0.2, width = 8, plot = p
)

print(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png(pdf) is OK"))
