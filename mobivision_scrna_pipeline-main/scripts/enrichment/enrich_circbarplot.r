#!/usr/bin/env Rscript
ofile <- parent.frame(2)$ofile

get_file_path <- function() {
  if (is.null(ofile)) {
    args <- commandArgs(FALSE)
    file <- args[grepl("^--file=", args)][1]
    ofile <- gsub("^--file=", "", file)
  }
  dirname(ofile)
}

basepath <- ifelse(
    length(getSrcDirectory(get_file_path)) > 0,
    getSrcDirectory(get_file_path),
    get_file_path()
)

source(file.path(basepath, "color.R"))

suppressPackageStartupMessages(library("optparse"))
# 加载 circlize 包
suppressPackageStartupMessages(library(circlize, lib.loc = "/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/R/library"))
## 使用 ComplexHeatmap 包额外添加图例
suppressPackageStartupMessages(library(ComplexHeatmap))


option_list <- list(
    make_option(c("-i", "--diff"),
        type = "character", default = NULL,
        help = "A-vs-B-diff-p-val-0.05-FC-2.gene.xls", metavar = "character"
    ),
    make_option(c("-e", "--enrich"),
        type = "character", default = NULL,
        help = "enrichment-[go|kegg]-A-vs-B-Total.xls", metavar = "character"
    ),
    make_option(c("-s", "--sort"),
        type = "character", default = "p-value",
        help = "columns for sorting [term/pathway]. default: p-value",
        metavar = "character"
    ),
    make_option(c("-o", "--outputdir"),
        type = "character", default = "./",
        help = "output directory for results", metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

##########################################################################
# https://github.com/lyao222lll/sheng-xin-xiao-bai-yu-2021/blob/main/2021/

mycol <- c(
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    "#CCEBC5", "#FFED6F", "#71A99F", "#CCCC8F", "#9895AE",
    "#C9665B", "#668EA9", "#CA904E", "#8FB254", "#CAA4B7"
)

load_table_from_str <- function(string) {
    if (string == "") return(data.frame())
    df <- read.table(
        text = string, header = T, sep = "\t",
        quote = "", check.names = F, comment.char = ""
    )
    return(df)
}

enrich_dat <- read.delim(opt$enrich, sep="\t",check.names=FALSE)
# enrich_str <- readChar(opt$enrich, file.info(opt$enrich)$size)
# enrich_dat <- load_table_from_str(enrich_str)

filename <- tools::file_path_sans_ext(basename(opt$enrich))


deg_dat <- read.delim(
    opt$diff,
    sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F,
    comment.char = ""
)

if (nrow(deg_dat) < 1 || nrow(enrich_dat) < 1) {
    sink(file.path(opt$outputdir, paste0(filename, ".circos.pdf")))
    sink()
    print("No diff genes")
    quit("no")
}

## 筛选对应差异基因数目大于 2 的 GO 条目
sub_data <- enrich_dat[which(enrich_dat[,"ListHits"]>2),]
enrich_dat <- head(sub_data[order(sub_data[[opt$sort]]), , drop = F], 20)

enrich_dat <- enrich_dat[
    order(enrich_dat$ListHits / enrich_dat$PopHits, decreasing = T), , drop = F
]

reg_count <- function(x, deg_dat, flag) {
    length(
        which(
            unlist(strsplit(x, ";")) %in%
                rownames(deg_dat)[deg_dat$Regulation == flag]
        )
    )
}

# Gene 转为字符型 
enrich_dat$Gene <- as.vector(enrich_dat$Gene)

dat <- data.frame(
    category = enrich_dat[
        ,
        grepl("Category|Classification_level1", colnames(enrich_dat)),
        drop = T
    ],
    gene_num.min = 0,
    gene_num.max = max(enrich_dat$PopHits),
    gene_num.rich = enrich_dat$PopHits,
    log.p = -log10(enrich_dat[["p-value"]]),
    rich.factor = enrich_dat$ListHits / enrich_dat$PopHits,
    up.regulated = vapply(
        enrich_dat$Gene,
        reg_count,
        deg_dat = deg_dat, flag = "Up",
        0,
        USE.NAMES = F
    ),
    down.regulated = vapply(
        enrich_dat$Gene,
        reg_count,
        deg_dat = deg_dat, flag = "Down",
        0,
        USE.NAMES = F
    )
)


# 默认按照原表格中的排列顺序
dat$id <- factor(enrich_dat[["id"]], levels = enrich_dat[["id"]])
rownames(dat) <- dat$id

#### 创建一个 pdf 画板
pdf(
    file.path(opt$outputdir, paste0(filename, ".circos.pdf")),
    width = 8, height = 6
)
# circle_size <- unit(1, "snpc")

## 整体布局
circos.par(
    gap.degree = 2, start.degree = 90, circle.margin = c(1e-5, .8, 1e-5, 1e-5)
)
# circos.par(
#     gap.degree = 2, start.degree = 90, circle.margin = c(1e-5, .8, 1e-5, 1e-5)
# )


## 第一圈
# 选择作图数据集，定义区块的基因总数量范围
plot_data <- dat[c("id", "gene_num.min", "gene_num.max")]
# 定义分组颜色
#cat_col <- mycol[seq_len(length(unique(dat$category)))]
#names(cat_col) <- unique(dat$category)
cat_col <- fix_color(unique(dat$category))
cat_color <- cat_col[dat$category]

 # 一个总布局
circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)
circos.track(
    # 圈图的高度、颜色等设置
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = cat_color,
    panel.fun = function(x, y) {
        # ylim、xlim 用于指定 id 文字标签添加的合适坐标
        ylim <- get.cell.meta.data("ycenter")
        xlim <- get.cell.meta.data("xcenter")
        # sector.name 用于提取 id 名称
        sector.name <- get.cell.meta.data("sector.index")
        # 绘制外周的刻度线
        circos.axis(
            h = "top", labels.cex = 0.4, major.tick.length = 0.4,
            labels.niceFacing = FALSE
        )
        # 将 id 文字标签添加在图中指定位置处
        circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
    }
)

## 第二圈，绘制富集的基因和富集 p 值
# 选择作图数据集，包括富集基因数量以及 p 值等信息
plot_data <- dat[c("id", "gene_num.min", "gene_num.rich", "log.p")]
# 标签数据集，仅便于作图时添加相应的文字标识用
label_data <- dat["gene_num.rich"]
# 定义一个 p 值的极值，以方便后续作图
p_max <- round(max(dat$log.p)) + 1
# 这两句用于定义 p 值的渐变颜色
RdYlBu <- c(
    "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
    "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4",
    "#313695"
)
colorsChoice <- colorRampPalette(rev(RdYlBu))
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
    plot_data,
    # 圈图的高度、颜色等设置
    track.height = 0.08, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
        # 区块的长度反映了富集基因的数量，颜色与 p 值有关
        circos.genomicRect(
            region, value,
            col = color_assign(value[[1]]), border = NA, ...
        )
        # 指定文字标签（富集基因数量）添加的合适坐标
        ylim <- get.cell.meta.data("ycenter")
        xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
        sector.name <- label_data[get.cell.meta.data("sector.index"), 1]
        # 将文字标签添（富集基因数量）加在图中指定位置处
        circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
    }
)

## 第三圈，绘制上下调基因
# 首先基于表格中上下调基因的数量，计算它们的占比
dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

# 随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c("id", "gene_num.min", "up")]
names(plot_data_up) <- c("id", "start", "end")
plot_data_up$type <- 1 # 分配 1 指代上调基因

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c("id", "up", "down")]
names(plot_data_down) <- c("id", "start", "end")
plot_data_down$type <- 2 # 分配 2 指代下调基因

# 选择作图数据集、标签数据集，并分别为上下调基因赋值不同颜色
plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c("up", "down", "up.regulated", "down.regulated")]
color_assign <- colorRamp2(breaks = c(1, 2), col = c("#febfca", "#86cdf9"))

# 继续绘制圈图
suppressMessages(circos.genomicTrackPlotRegion(
    plot_data,
    # 圈图的高度、颜色等设置
    track.height = 0.08, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(
            # 区块的长度反映了上下调基因的相对占比
            region, value,
            col = color_assign(value[[1]]), border = NA, ...
        )
        # 指定文字标签（上调基因数量）添加的合适坐标
        ylim <- get.cell.meta.data(
            "cell.bottom.radius"
        ) - 0.5
        xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
        sector.name <- label_data[get.cell.meta.data("sector.index"), 3]
        # 将文字标签（上调基因数量）添加在图中指定位置处
        circos.text(
            xlim, ylim, sector.name,
            cex = 0.4, niceFacing = FALSE
        )
        xlim <- (label_data[get.cell.meta.data("sector.index"), 2] +
            label_data[get.cell.meta.data("sector.index"), 1]) / 2
        sector.name <- label_data[get.cell.meta.data("sector.index"), 4]
        # 类似的操作，将下调基因数量的标签也添加在图中
        circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
    }
))

## 第四圈，绘制富集得分
# 选择作图数据集，标准化后的富集得分
plot_data <- dat[
    c("id", "gene_num.min", "gene_num.max", "rich.factor")
]
# 将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
label_data <- dat["category"]
# color_assign <- c(
#     "Metabolism" = "#F7CC13",
#     "Environment Information Processing " = "#954572",
#     "Organismal Systems" = "#0796E0"
# )

circos.genomicTrack(
    plot_data,
    # 圈图的高度、颜色等设置
    ylim = c(0, 1), track.height = 0.4, bg.col = "gray95", bg.border = NA,
    panel.fun = function(region, value, ...) {
        # sector.name 用于提取 id 名称，并添加在下一句中匹配的高级分类，以分配颜色
        sector.name <- get.cell.meta.data("sector.index")
        # 等位线
        # circos.lines(c(0, max(region)), c(0.1, 0.1), col = "gray", lwd = 0.3)
        circos.lines(c(0, max(region)), c(0.2, 0.2), col = "gray", lwd = 0.3)
        # circos.lines(c(0, max(region)), c(0.3, 0.3), col = "gray", lwd = 0.3)
        circos.lines(c(0, max(region)), c(0.4, 0.4), col = "gray", lwd = 0.3)
        # circos.lines(c(0, max(region)), c(0.5, 0.5), col = "gray", lwd = 0.3)
        circos.lines(c(0, max(region)), c(0.6, 0.6), col = "gray", lwd = 0.3)
        # circos.lines(c(0, max(region)), c(0.7, 0.7), col = "gray", lwd = 0.3)
        circos.lines(c(0, max(region)), c(0.8, 0.8), col = "gray", lwd = 0.3)
        # circos.lines(c(0, max(region)), c(0.9, 0.9), col = "gray", lwd = 0.3)
        # 绘制矩形区块，高度代表富集得分，颜色代表分类
        circos.genomicRect(
            region, value,
            col = cat_col[label_data[sector.name, 1]],
            border = NA, ytop.column = 1, ybottom = 0, ...
        )
    }
)

polygon(
    x = c(-0.2, -0.2, -0.05, -0.05),
    y = c(0.12, 0.18, 0.18, 0.12),
    col = "grey", border = NA
)
text(-0.075, 0.15, "Number of Genes", cex = .5, pos = 4)
text(-0.125, 0.15, "100", cex = .5)

polygon(
    x = c(-0.2, -0.2, -0.05, -0.05),
    y = c(0.02, 0.08, 0.08, 0.02),
    col = "#febfca", border = NA
)
text(-0.075, 0.05, "Up-regulated", cex = .5, pos = 4)

polygon(
    x = c(-0.2, -0.2, -0.05, -0.05),
    y = c(-0.08, -0.02, -0.02, -0.08),
    col = "#86cdf9", border = NA
)
text(-0.075, -0.05, "Down-regulated", cex = .5, pos = 4)

polygon(
    x = c(-0.2, -0.2, -0.05, -0.05),
    y = c(-0.125, -0.175, -0.19, -0.11),
    col = "grey", border = NA
)
text(-0.075, -0.15, "Rich Factor [0-1)", cex = .5, pos = 4)

## 绘图完毕后，不要忘了清除痕迹，以免影响下一次作图
circos.clear()

category_legend <- Legend(
    labels = names(cat_col),
    type = "points", pch = NA,
    background = cat_col,
    labels_gp = gpar(fontsize = 7),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
)

# updown_legend <- Legend(
#     labels = c("Up-regulated", "Down-regulated"),
#     type = "boxplot", pch = 32,
#     legend_gp = gpar(width = 8),
#     background = c("#febfca", "#86cdf9"),
#     labels_gp = gpar(fontsize = 7),
#     grid_height = unit(0.5, "cm"),
#     grid_width = unit(0.5, "cm")
# )

pvalue_legend <- Legend(
    col_fun = colorRamp2(
        round(seq(0, p_max, length.out = 6), 0),
        colorRampPalette(rev(RdYlBu))(6)
    ),
    legend_height = unit(3, "cm"),
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 8),
    title_position = "topleft",
    title = expression(paste("-", log[10], "p-value"))
)

lgd_list_vertical <- packLegend(category_legend, pvalue_legend)
pushViewport(viewport(x = 0.85, y = 0.5))
grid.draw(lgd_list_vertical)
upViewport()

# lgd_cent <- packLegend(updown_legend)
# pushViewport(viewport(x = 0.35, y = 0.5))
# grid.draw(lgd_cent)
# upViewport()

## 关闭 pdf 画板
dev.off()

# pdf convert png 
setwd(opt$outputdir)
system("for i in  `ls *circos.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")
