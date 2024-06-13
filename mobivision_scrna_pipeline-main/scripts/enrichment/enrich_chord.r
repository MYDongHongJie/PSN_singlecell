#!/usr/bin/env Rscript

### 使用 GOplot R 包绘制圈图 ### 
library(GOplot,lib.loc = "/home/lipeng/miniconda3/lib/R/library")
library(stringr)
library(optparse)

################## Color ###################
#                                          #
#   +----------------------------------+   #
#   | up:   #febfca | down:    #86cdf9 |   #
#   +----------------------------------+   #
#                                          #
############################################

option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "input top terms file name",
        metavar = "character"
    ),
    make_option(c("-d", "--deg"),
        type = "character", default = NULL,
        help = "deg file name", metavar = "character"
    ),
    make_option(c("-s", "--sort"),
        type = "character", default = "p-value",
        help = "sort terms, default[p-value]", metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory",
        metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("--input must be supplied", call. = FALSE)
}


main <- function() {
    topFile <- opt$input
    degFile <- opt$deg
    outputDir <- opt$outdir
    sortBy <- opt$sort # topfile's column name

    topString <- readChar(topFile, file.info(topFile)$size)
    degString <- readChar(degFile, file.info(degFile)$size)

    fileName <- tools::file_path_sans_ext(basename(topFile))
    outputName <- sub("\\.", ".chord.", fileName)
    # get dbname from filename
    db <- unlist(strsplit(fileName, "[.]"))[1]

    topDf <- load_table_from_str(topString)
    degDf <- load_table_from_str(degString)

    if (nrow(topDf) < 1 || nrow(degDf) < 1) {
        sink(file.path(outputDir, paste(outputName, "pdf", sep = ".")))
        sink()
        return("No diff genes")
    }

    # sort by one column
    topDf <- topDf[order(topDf[[sortBy]], decreasing = T), ]

    chord <- make_chord_data(topDf, degDf, db = db)

    dir.create(file.path(outputDir), showWarnings = FALSE)

    # for (i in names(chordList)) {
    #     chordPlot <- enrich_chord(chordList[[i]], gene.order = "logFC")
    #     ggsave(file.path(outputDir, paste(fileName, i, "pdf", sep = ".")),
    #         chordPlot, "pdf",
    #         width = 7, height = 10
    #     )
    # }
    chordPlot <- enrich_chord(chord, gene.order = "logFC")
    write.table(
        chord,
        file.path(outputDir, paste(outputName, "xls", sep = ".")),
        sep = "\t", col.names = NA, quote = F
    )
    ggsave(file.path(outputDir, paste(outputName, "pdf", sep = ".")),
        chordPlot, "pdf",
        width = 7, height = 10
    )
    ggsave(file.path(outputDir, paste(outputName, "png", sep = ".")),
        chordPlot, "png",
        width = 7, height = 10, dpi=300
    )
}

### 修改 GOplot circle_dat 源码，不对基因名称转换成大写  
circle_dat <- function (terms, genes)
{
    colnames(terms) <- tolower(colnames(terms))
    #terms$genes <- toupper(terms$genes)
    #genes$ID <- toupper(genes$ID)
    tgenes <- strsplit(as.vector(terms$genes), ", ")
    if (length(tgenes[[1]]) == 1)
        tgenes <- strsplit(as.vector(terms$genes), ",")
    count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
    logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x,
        genes$ID)])
    if (class(logFC) == "factor") {
        logFC <- gsub(",", ".", gsub("\\.", "", logFC))
        logFC <- as.numeric(logFC)
    }
    s <- 1
    zsc <- c()
    for (c in 1:length(count)) {
        value <- 0
        e <- s + count[c] - 1
        value <- sapply(logFC[s:e], function(x) ifelse(x > 0,
            1, -1))
        zsc <- c(zsc, sum(value)/sqrt(count[c]))
        s <- e + 1
    }
    if (is.null(terms$id)) {
        df <- data.frame(category = rep(as.character(terms$category),
            count), term = rep(as.character(terms$term), count),
            count = rep(count, count), genes = as.character(unlist(tgenes)),
            logFC = logFC, adj_pval = rep(terms$adj_pval, count),
            zscore = rep(zsc, count), stringsAsFactors = FALSE)
    }
    else {
        df <- data.frame(category = rep(as.character(terms$category),
            count), ID = rep(as.character(terms$id), count),
            term = rep(as.character(terms$term), count), count = rep(count,
                count), genes = as.character(unlist(tgenes)),
            logFC = logFC, adj_pval = rep(terms$adj_pval, count),
            zscore = rep(zsc, count), stringsAsFactors = FALSE)
    }
    return(df)
}


# tab separated, first line as a header
load_table_from_str <- function(string) {
    if (string == "") return(data.frame())
    df <- read.table(
        text = string, header = T, sep = "\t",
        quote = "", check.names = F, comment.char = ""
    )
    return(df)
}

# parse input table, make chord data for GOplot
make_chord_data <- function(topDf, degDf, db = "GO") {
    if (db == "GO") {
        ctg <- "Category"
        n <- unique(topDf[, ctg, drop = T])
        term <- "Term"
    } else if (db == "KEGG") {
        ctg <- "Classification_level1"
        n <- gsub(" ", "_", unique(topDf[, ctg, drop = T]))
        term <- "Term"
    } else {
        ctg <- "Pathway"
        n <- unique(topDf[, ctg, drop = T])
        term <- "Pathway"
    }

    geneDf <- data.frame(ID = degDf$gene, logFC = degDf$log2FoldChange)
    termDf <- topDf[, c(ctg, "id", term, "p-value", "Gene")]
    colnames(termDf) <- c(
        "category", "ID", "term", "adj_pval", "genes"
    ) # adj_pval
    termDf$genes <- gsub(";", ", ", termDf$genes)
    termDf$category <- gsub(" ", "_", termDf$category) # KEGG

    # termDfList <- lapply(
    #     setNames(n, n),
    #     function(x) termDf[termDf$category == x, , drop = F]
    # )

    # p值排序，最多 10 term
    termDf <- head(termDf[order(termDf$adj_pval, decreasing=F), , drop = F], 10)

    # circDfList <- lapply(
    #     setNames(n, n),
    #     function(x) circle_dat(termDfList[[x]], geneDf)
    # )

    ### 修改 circle_dat 源码
    circDf <- circle_dat(termDf, geneDf)
    #circDf <- Goplot::circle_dat(termDf, geneDf)

    # fix Inf value
    geneDf$logFC[geneDf$logFC == -Inf] <-
        min(geneDf$logFC[!is.infinite(geneDf$logFC)])
    geneDf$logFC[geneDf$logFC == Inf] <-
        max(geneDf$logFC[!is.infinite(geneDf$logFC)])

    ### 修改 circle_dat 源码
    ## geneDf$ID <- toupper(geneDf$ID) # bug fix, same with circle_dat

    # most 10 genes in each term
    subgene <- function(geneDf, circDf) {
        # subg <- geneDf[which(geneDf$ID %in% circDf$genes), , drop = F]
        # subg <- subg[
        #     head(order(abs(subg$logFC), decreasing = T), 50), , drop = F
        # ]
        select_genes <- function(x) {
            subg <- geneDf[
                which(geneDf$ID %in% circDf$genes[circDf$term == x]), ,
                drop = F
            ]
            subg <- subg[
                head(order(abs(subg$logFC), decreasing = T), 10), , drop = F
            ]
        }

        subg <- do.call(rbind, lapply(unique(circDf$term), select_genes))
        subg <- subg[!duplicated(subg), ]

        subc <- circDf[which(circDf$genes %in% subg$ID), , drop = F]

        return(list(subg, subc))
    }

    # chordDfList <- lapply(
    #     setNames(n, n),
    #     function(x) {
    #         chord_dat(
    #             data = subgene(geneDf, circDfList[[x]])[[2]],
    #             genes = subgene(geneDf, circDfList[[x]])[[1]],
    #             process = unique(
    #                 subgene(geneDf, circDfList[[x]])[[2]]$term
    #             )
    #         )
    #     }
    # )

    chordDf <- chord_dat(
        data = subgene(geneDf, circDf)[[2]],
        genes = subgene(geneDf, circDf)[[1]],
        process = unique(
            subgene(geneDf, circDf)[[2]]$term
        )
    )

    return(chordDf)
}

# modified GOChord. colour, border, label default size ...
enrich_chord <- function(data, title, space, gene.order, gene.size, gene.space,
                         nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col,
                         border.size, process.label, limit) {
    y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
    Ncol <- dim(data)[2]
    if (missing(title)) {
        title <- ""
    }
    if (missing(space)) {
        space <- 0
    }
    if (missing(gene.order)) {
        gene.order <- "none"
    }
    if (missing(gene.size)) {
        gene.size <- 2
    }
    if (missing(gene.space)) {
        gene.space <- 0.2
    }
    if (missing(lfc.col)) {
        lfc.col <- c("#febfca", "azure", "#86cdf9")
    }
    if (missing(lfc.min)) {
        lfc.min <- -Inf
    }
    if (missing(lfc.max)) {
        lfc.max <- Inf
    }
    if (missing(border.size)) {
        border.size <- 0
    }
    if (missing(process.label)) {
        process.label <- 7
    }
    if (missing(limit)) {
        limit <- c(0, 0)
    }
    if (gene.order == "logFC") {
        data <- data[order(data[, Ncol], decreasing = T), , drop = F]
    }
    if (gene.order == "alphabetical") {
        data <- data[order(rownames(data)), , drop = F]
    }
    if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
        if (nlfc == 1) {
            # bug fix, do not drop
            cdata <- check_chord(data[, 1:(Ncol - 1), drop = F], limit)
            lfc <- sapply(rownames(cdata), function(x) {
                data[match(
                    x,
                    rownames(data)
                ), Ncol, drop = F]
            })
        } else {
            cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
            lfc <- sapply(rownames(cdata), function(x) {
                data[
                    ,
                    (Ncol - nlfc + 1), drop = F
                ]
            })
        }
    } else {
        cdata <- check_chord(data, limit)
        lfc <- 0
    }
    if (missing(ribbon.col)) {
        # 36 colors
        default.col <- c(
            "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
            "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
            "#71A99F", "#CCCC8F", "#9895AE", "#C9665B", "#668EA9", "#CA904E",
            "#8FB254", "#CAA4B7", "#AEAEAE", "#966697", "#A3BC9E", "#CCBE59",
            "#557F77", "#99996B", "#727083", "#974D44", "#4D6A7F", "#986C3B",
            "#6B853F", "#977B89", "#828282", "#714D71", "#7A8D76", "#998E43"
        )
        colRib <- default.col[1:dim(cdata)[2]]
    } else {
        colRib <- ribbon.col
    }
    nrib <- colSums(cdata)
    ngen <- rowSums(cdata)
    Ncol <- dim(cdata)[2]
    Nrow <- dim(cdata)[1]
    colRibb <- c()
    for (b in seq_len(length(nrib))) {
        colRibb <- c(colRibb, rep(
            colRib[b],
            202 * nrib[b]
        ))
    }
    r1 <- 1
    r2 <- r1 + 0.1
    xmax <- c()
    x <- 0
    for (r in seq_len(length(nrib))) {
        perc <- nrib[r] / sum(nrib)
        xmax <- c(xmax, (pi * perc) - space)
        if (length(x) <= Ncol - 1) {
            x <- c(x, x[r] + pi * perc)
        }
    }
    xp <- c()
    yp <- c()
    l <- 50
    for (s in 1:Ncol) {
        xh <- seq(x[s], x[s] + xmax[s], length = l)
        xp <- c(
            xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] +
                xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)),
            r2 * sin(x[s])
        )
        yp <- c(
            yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] +
                xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)),
            r2 * cos(x[s])
        )
    }
    df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol),
        each = 4 + 2 * l
    ))
    xp <- c()
    yp <- c()
    logs <- NULL
    x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
    xmax2 <- rep(-pi / Nrow + space, length = Nrow)
    for (s in 1:Nrow) {
        xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
        if (nlfc <= 1) {
            xp <- c(
                xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) *
                    sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]),
                r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)),
                r2 * sin(x2[s])
            )
            yp <- c(
                yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) *
                    cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]),
                r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)),
                r2 * cos(x2[s])
            )
        } else {
            tmp <- seq(r1, r2, length = nlfc + 1)
            for (t in 1:nlfc) {
                logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
                xp <- c(
                    xp, (tmp[t]) * sin(x2[s]), (tmp[t]) *
                        sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]),
                    tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t +
                        1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s])
                )
                yp <- c(
                    yp, (tmp[t]) * cos(x2[s]), (tmp[t]) *
                        cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]),
                    tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t +
                        1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s])
                )
            }
        }
    }
    if (lfc[1] != 0) {
        if (nlfc == 1) {
            df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow),
                each = 4 + 2 * l
            ), logFC = rep(lfc, each = 4 +
                2 * l))
        } else {
            df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc *
                Nrow)), each = 4 + 2 * l), logFC = rep(logs,
                each = 4 + 2 * l
            ))
        }
    } else {
        df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow),
            each = 4 + 2 * l
        ))
    }
    aseq <- seq(0, 180, length = length(x2))
    angle <- c()
    for (o in aseq) {
        if ((o + 270) <= 360) {
            angle <- c(angle, o + 270)
        } else {
            angle <- c(angle, o - 90)
        }
    }
    df_texg <- data.frame(
        xgen = (r1 + gene.space) * sin(x2 +
            xmax2 / 2), ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),
        labels = rownames(cdata), angle = angle
    )
    df_texp <- data.frame(
        xpro = (r1 + 0.15) * sin(x + xmax / 2),
        ypro = (r1 + 0.15) * cos(x + xmax / 2),
        labels = str_wrap(colnames(cdata), 70),
        stringsAsFactors = FALSE
    )
    cols <- rep(colRib, each = 4 + 2 * l)
    x.end <- c()
    y.end <- c()
    processID <- c()
    for (gs in seq_len(length(x2))) {
        val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] +
            1)
        pros <- which((cdata[gs, ] != 0) == T)
        for (v in 1:(length(val) - 1)) {
            x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
            y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
            processID <- c(processID, rep(pros[v], 2))
        }
    }
    df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
    df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end), ]
    x.start <- c()
    y.start <- c()
    for (rs in seq_len(length(x))) {
        val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] +
            1)
        for (v in 1:(length(val) - 1)) {
            x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
            y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
        }
    }
    df_bezier$x.start <- x.start
    df_bezier$y.start <- y.start
    df_path <- bezier(df_bezier, colRib)
    if (length(df_genes$logFC) != 0) {
        tmp <- sapply(df_genes$logFC, function(x) {
            ifelse(x >
                lfc.max, lfc.max, x)
        })
        logFC <- sapply(tmp, function(x) {
            ifelse(x < lfc.min,
                lfc.min, x
            )
        })
        df_genes$logFC <- logFC
    }
    g <- ggplot() +
        geom_polygon(data = df_process, aes(x, y,
            group = id
        ), fill = "gray70", inherit.aes = F, color = NA, size = border.size) +
        geom_polygon(
            data = df_process, aes(x, y, group = id),
            fill = cols, inherit.aes = F, alpha = 0.6, color = NA,
            size = border.size
        ) +
        ## 增加 term 图例
        geom_point(aes(x = xpro, y = ypro, size = factor(labels,
            levels = labels
        ), shape = NA), data = df_texp) +
        guides(size = guide_legend("",
            ncol = 2, byrow = T,
            override.aes = list(
                shape = 22, fill = unique(cols), colour = NA,
                size = 4
            )
        )) +
        theme(legend.text = element_text(
            size = process.label,
            margin = margin(r = 0, b = 0, l = 0, t = 0, unit = "pt")
        )) +
        geom_text(aes(xgen * .95, ygen * .95, label = labels, angle = angle),
            hjust = 1,
            data = df_texg, size = gene.size
        ) +
        geom_polygon(aes(
            x = lx,
            y = ly, group = ID
        ),
        data = df_path, fill = colRibb, alpha = 0.8,
        color = NA, size = border.size, inherit.aes = F
        ) +
        coord_fixed(clip = "off") +
        labs(title = title) +
        theme_blank
    if (nlfc >= 1) {
        bks <- unique(c(min(df_genes$logFC), max(df_genes$logFC)))
        lbs <- unique(c(round(min(df_genes$logFC)), round(max(df_genes$logFC))))
        if (length(bks) > length(lbs)) {
            lbs <- unique(
                c(round(min(df_genes$logFC), 3), round(max(df_genes$logFC), 3))
            )
        }
        g + geom_polygon(data = df_genes, aes(x, y,
            group = id,
            fill = logFC
        ), inherit.aes = F, color = NA, size = border.size) +
            scale_fill_gradient2("logFC",
                space = "Lab", low = lfc.col[3],
                mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(
                    title = expression(paste(log[2], "FC")), 
                    title.position = "top",
                    title.hjust = 0.5,
                    order = 1
                ), breaks = bks, labels = lbs
            ) + theme(
                legend.position = "bottom",
                # legend.background = element_blank(),
                legend.box = "vertical", legend.direction = "horizontal",
                legend.spacing.y = unit(0, "line"),
                legend.spacing.x = unit(0, "line"),
                legend.key = element_rect(fill = NA, color = NA),
                legend.key.height = unit( # symbol size
                    1.2,
                    "line"
                ), legend.key.width = unit(
                    .6,
                    "line"
                )
            ) + theme(
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(2, 2, 2, 2), 'cm'),
                legend.box.margin=margin(50, 0, 0, 0)
            )
    } else {
        g + geom_polygon(
            data = df_genes, aes(x, y, group = id),
            fill = "gray50", inherit.aes = F, color = NA, size = border.size
        ) +
            theme(
                legend.position = "bottom",
                # legend.background = element_blank(),
                legend.box = "vertical", legend.direction = "horizontal",
                legend.spacing.y = unit(0, "pt"),
                legend.spacing.x = unit(0, "line"),
                legend.key = element_rect(fill = NA, color = NA),
                legend.key.height = unit( # symbol size
                    1.2,
                    "line"
                ), legend.key.width = unit(
                    .6,
                    "line"
                )
            ) + theme(
                plot.background = element_rect(fill = "white"),
                plot.margin = unit(c(2, 2, 2, 2), 'cm'),
                legend.box.margin=margin(50, 0, 0, 0)
            )
    }
}

environment(enrich_chord) <- environment(GOChord)

##########################################################
argv <- commandArgs(T)
if (length(argv > 0)) main()
