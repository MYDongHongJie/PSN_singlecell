#!/usr/bin/env Rscript


# by The Coder, 20160527
library("optparse")
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
        "Rscript top10X3_GO.r ",
        "-i enrichment-go-Group1-vs-Group2-Down.txt -m Down -o outdir/"
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
library(stringr)
library(grid)
library(RColorBrewer)

groupname <- gsub(
    "\\.(txt|xls)$", "", gsub("^enrichment-go-", "", basename(opt$input))
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

limit70 <- function(s, df) {
    k <- as.character(s)
    if (str_length(s) > 55) {
        go <- subset(df, df["Term"] == s)
        go <- go$id[1]
        term <- sub("[^ ]+$", "...", substr(k, 1, 55))
        k <- paste0(go, ":", term)
    }
    return(k)
}

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

top10 <- function(i) {
    return(i[order(head(i, 10)["p-value"]), ])
}

load_table_from_str <- function(string) {
    if (string == "") return(data.frame())
    df <- read.table(
        text = string, header = T, sep = "\t",
        quote = "", check.names = F
    )
    return(df)
}

d_str <- readChar(opt$input, file.info(opt$input)$size)
d <- load_table_from_str(d_str)

d <- d[which(d[, "ListHits"] > 2), ]
if ("p-value" %in% colnames(d)) d <- d[order(d[, "p-value"]), ]

# stopifnot(nrow(d)>0)
if (nrow(d) == 0) {
    print("d items = 0, program exit!")
    sink(paste0(opt$outpath, "/GO.top.", opt$mark, ".pdf"))
    sink()
    sink(paste0(opt$outpath, "/GO.top.", opt$mark, ".png"))
    sink()
    sink(paste0(opt$outpath, "/GO.top.", opt$mark, ".xls"))
    sink()
    print("No diff genes")
    q()
}

dp <- rbind(
    top10(d[d["Category"] == "biological_process", ]),
    top10(d[d["Category"] == "cellular_component", ]),
    top10(d[d["Category"] == "molecular_function", ])
)
write.table(
    dp[, c(1, 2, 3, 4, 8, 10, 11)],
    paste0(opt$outpath, "/GO.top.", opt$mark, ".xls"),
    sep = "\t", quote = FALSE,
    col.names = TRUE, row.names = FALSE
)

d1 <- dp[which(dp[, 3] == "biological_process"), ]
d2 <- dp[which(dp[, 3] == "cellular_component"), ]
d3 <- dp[which(dp[, 3] == "molecular_function"), ]

d_l <- rbind(d1, d2, d3)
d_l[which(dp$Category == "biological_process"), "color"] <- "#6ea647"
d_l[which(dp$Category == "cellular_component"), "color"] <- "#e37933"
d_l[which(dp$Category == "molecular_function"), "color"] <- "#5893c9"

# d_l["Term"] <- apply(d_l["Term"], 1, limit70, d_l)

pdf(NULL)
d_l[, "Term"] <- vapply(d_l[, "Term"], modify_name, "a", USE.NAMES = F)
dev.off()

d_l[, "Term"] <- make.unique(d_l$Term)

d_l$Term <- factor(d_l$Term, levels = d_l$Term)
colnames(d_l)[8] <- "pvalue"

firstup <- function(x) {
    x <- as.character(x)
    words <- unlist(strsplit(x, "_"))
    for (i in seq_along(words)) {
        substr(words[i], 1, 1) <- toupper(substr(words[i], 1, 1))
    }
    paste(words, collapse = " ")
}

d_l$Category <- vapply(d_l$Category, firstup, "a", USE.NAMES = F)

if (max(d_l$pvalue) == "0") {
    min_y <- 0
    max_y <- 50
    p <- ggplot(
        data = d_l,
        aes(
            x = Term, y = -log(pvalue, 10), width = 0.6,
            fill = Category, space = 0.6, cex.main = 3, cex.lab = 2
        )
    ) +
        geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.5) +
        ylim(min_y, max_y) +
        labs(x = "", y = expression("-log"[10] * " p-value"), title = paste0(groupname, "\n", "Top 30 GO Term") ) +
        theme_bw() +
        scale_fill_manual(values = unique(d_l$color)) +
        theme( axis.text = element_text(# angle = 45, hjust = 1,
                size = 14,
                color = "black"
            ) ) +
        theme(legend.text = element_text(size = 14), legend.position = "none") +
        facet_grid(vars(Category), scales = "free", space = "free_y") +
        theme(text = element_text(size = 14), strip.text = element_text(size = 14) ) +
        theme( plot.title = element_text( hjust = 0.5, vjust = 4, size = 18, family = "ArialMT" ) ) +
        theme(plot.margin = unit(c(2, 2, 2, 2), "lines")) +
        theme(panel.grid = element_blank()) +
        theme(panel.grid = element_blank()) +
        theme(panel.border = element_rect( color = "black", size = 1, fill = NA ) ) +
        theme(
            text = element_text(family = "sans"),
            plot.title = element_text(face = "bold"),
            axis.text = element_text(size = 15, color = "black"),
            axis.line = element_blank(),
            strip.background = element_rect(color = "#e0e0e0", fill = "#e0e0e0")
        ) +
        coord_flip()
} else {
    p <- ggplot(
        data = d_l,
        aes(
            x = Term, y = -log(pvalue, 10), width = 0.6,
            fill = Category, space = 0.6, cex.main = 3, cex.lab = 2
        )
    ) +
        geom_bar(
            stat = "identity", position = position_dodge(0.7), width = 0.5
        ) +
        labs(
            x = "",
            y = expression("-log"[10] * " p-value"),
            title = paste0(groupname, "\n", "Top 30 GO Term")
        ) +
        theme_bw() +
        scale_fill_manual(values = unique(d_l$color)) +
        theme(
            axis.text = element_text(
                # angle = 45, hjust = 1,
                size = 14,
                color = "black"
            )
        ) +
        theme(
            legend.position = "none"
            # legend.title = element_blank(),
            # legend.position = "top",
            # legend.justification = 1.8,
            # legend.key.width = unit(1, "lines")
        ) +
        theme(legend.text = element_text(size = 14)) +
        facet_grid(vars(Category), scales = "free", space = "free_y") +
        theme(
            text = element_text(size = 14),
            strip.text = element_text(size = 14)
        ) +
        theme(
            plot.title = element_text(
                hjust = 0.5, vjust = 4, size = 16, family = "ArialMT"
            )
        ) +
        theme(plot.margin = unit(c(2, 2, 2, 2), "lines")) +
        theme(panel.grid = element_blank()) +
        theme(panel.grid = element_blank()) +
        theme(
            panel.border = element_rect(
                color = "black", size = 1, fill = NA
            )
        ) +
        theme(
            text = element_text(family = "sans"),
            plot.title = element_text(face = "bold"),
            axis.text = element_text(size = 15, color = "black"),
            axis.line = element_blank(),
            strip.background = element_rect(color = "#e0e0e0", fill = "#e0e0e0")
        ) +
        coord_flip(ylim = c(0, max(-log(d_l$pvalue, 10)) * 1.1), expand = FALSE)
}

pdf(NULL)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl("strip-r", g$layout$name))

fills <- c()
if (nrow(d1) > 0){
  fills <- append(x=fills,"#6ea647")
}
if (nrow(d2) > 0){
  fills <- append(x=fills,"#e37933")
}
if (nrow(d3) > 0){
  fills <- append(x=fills,"#5893c9")
}
# fills <- c("#6ea647", "#e37933", "#5893c9")

k <- 1
for (i in stripr) {
    j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
}
dev.off()

pdf(
    paste0(opt$outpath, "/GO.top.", opt$mark, ".pdf"),
    height = 11, width = 9
)
grid::grid.draw(g)
dev.off()

png(
    paste0(opt$outpath, "/GO.top.", opt$mark, ".png"),
    height = 11, width = 9, units = "in", type = "cairo-png", res = 300
)
grid::grid.draw(g)
dev.off()

print(paste0(opt$outpath, "/GO.top.", opt$mark, ".png(pdf) is OK"))
