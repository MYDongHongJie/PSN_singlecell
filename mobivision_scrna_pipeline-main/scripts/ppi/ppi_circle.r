
library("optparse")
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character", default = NULL,
        help = "top300 ppi", metavar = "character"
    ),
    make_option(
        c("-d", "--diff"),
        type = "character", default = NULL,
        help = "diff file", metavar = "character"
    ),
    # make_option(
        # c("-n", "--number"),
        # type = "character", default = "25",
        # help = "show number", metavar = "character"
    # ),
    make_option(
        c("-o", "--output"),
        type = "character", default = NULL,
        help = "outfile directory", metavar = "character"
    )
)

opt_parser <- OptionParser(
    option_list = option_list,
    epilogue = paste0(
        "Rscript ppi_circle.r ",
        "-i group_case-vs-control-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.tsv -d /public/scRNA_works/works/ziqingzhen/test/report/result/Diffexp/group_case-vs-control-diff-pval-0.05-FC-1.5_anno.xls -o outdir/"
    )
)

opt <- parse_args(opt_parser)
if (is.null(opt$input) | is.null(opt$diff)) {
    print_help(opt_parser)
    stop("--input --diff must be supplied", call. = FALSE)
}
if (!file.exists(opt$input)) {
    print("input is not exists")
    q()
}
if ( is.null(opt$output)){
    print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
    output_dir = getwd()
}else{
    output_dir = opt$output
    if ( !file.exists(output_dir) ){
        dir.create(output_dir)
    }
}

library(tidyverse)
library(igraph)
library(dplyr)
library(ggraph)

df_ppi <- read.table(opt$input, sep = "\t", header = T)
df_fc <- read.delim(opt$diff, sep = "\t", header = T)

#diffname <- gsub(
#    "top_300_diff-", "",
#    gsub("_gene2gene_network.txt", "", basename(opt$input))
#)
diffname <- gsub("\\.(tsv)$", "", gsub(".ppi_network", "", basename(opt$input)))
#filename = strsplit(diffname, split='_' )[[1]][2]
# if (dim(df_ppi)[1] < 1) {
    # # "ppi fill is empty"
    # sink(
        # paste0(
            # opt$outpath, "/top_",25, "_diff-", diffname, "_network.pdf"
        # )
    # )
    # sink()
    # q()
# }


# get links
#df_ppi <-df_ppi
#df_ppi <- head(df_ppi[order(df_ppi[,"score"]),], 50)



# cal gene nodes and weight,add FC

list1 <- as.character(df_ppi$preferredName_A)
list2 <- as.character(df_ppi$preferredName_B)
list <- as.character(cbind(list1, list2))
df_num <- as.data.frame(table(list))
names(df_num) <- c("gene", "num")
gene2Num <- df_num[order(df_num[, 2], decreasing = T), ]
gene2FC <- df_fc[c("gene", "FoldChange")]
mergr_nodes <- merge(
    gene2Num, gene2FC,
    by.x = "gene", by.y = "gene", all.x = T
)
# get color
#min_fc = max(!is.infinite(mergr_nodes$FoldChange))

tmp <- mergr_nodes[which(mergr_nodes$FoldChange != "Inf" & mergr_nodes$FoldChange != "-Inf"), ]
max <- max(tmp$FoldChange)
min <- min(tmp$FoldChange)
mergr_nodes$FoldChange <- as.numeric(sub("-Inf", min, mergr_nodes$FoldChange))
mergr_nodes$FoldChange <- as.numeric(sub("Inf", max, mergr_nodes$FoldChange))

df_Up <- mergr_nodes[which(mergr_nodes[, 3] > 1), ]
df_Down <- mergr_nodes[which(mergr_nodes[, 3] < 1), ]
df_Up <- df_Up[order(df_Up[, 2], decreasing = F), ]
df_Down <- df_Down[order(df_Down[, 2], decreasing = F), ]
palette_Up <- colorRampPalette(
    c("#febfca", "#ff4d40")
)(n = length(rownames(df_Up)))
palette_Down <- colorRampPalette(
    c("steelblue", "#86cdf9")
)(n = length(rownames(df_Down)))
df_Up$color <- palette_Up
df_Down$color <- palette_Down
my_nodes <- rbind(df_Up, df_Down)

my_link <- df_ppi[c(3, 4, 6)]
# rm NA if less than top num
my_link <- my_link %>% drop_na()

colnames(my_link)[3] <- "weight"

color_num <- seq(length(rownames(df_num)), 1)
# get network
net <- graph.data.frame(my_link, my_nodes, directed = F)
# nodes size &color
deg <- igraph::degree(net, mode = "all")
V(net)$size <- log2(deg * 50 + 50)
V(net)$deg <- deg

list_sort <- as.data.frame(V(net)$name)
names(list_sort) <- "gene"
re_col <- left_join(list_sort, my_nodes[, c(1, 4)], by = "gene")
V(net)$color <- re_col$color

print(max(abs(log2(V(net)$FoldChange+0.001))))

p0 <- ggraph(net, layout = "circle")

p <- p0 + geom_edge_arc(
    color = "grey",
    strength = 0.1, alpha = 0.5,
    end_cap = circle(3, "mm"), start_cap = circle(3, "mm")
) +
    geom_node_point(
        pch = 19,
        aes(color = log2(V(net)$FoldChange+0.001), size = deg)
    )+
    geom_node_text(
        aes(label = name),
        size = 4,
        repel = TRUE,
        #max.overlaps = Inf
    ) +
    # scale_edge_width(range = c(0.5, 1.5), guide = "none") +
    scale_size(
        range = c(5, 8),
        'Degree',
    ) +
    scale_colour_gradientn(
        expression(paste(log[2], " FC")),
        limits = c(
            -max(abs(log2(V(net)$FoldChange+0.001))),
            max(abs(log2(V(net)$FoldChange+0.001)))
        ),
        colours = c("steelblue", "#86cdf9", "azure", "#febfca", "#ff4d40")
    ) +
    guides(size = guide_legend(order = 1)) +
    theme_void() + theme(plot.background = element_rect(color = "black",
                                        size = 0.1),
         plot.margin = margin(t = 15,  # 顶部边缘距离
                              r = 15,  # 右边边缘距离
                              b = 15,  # 底部边缘距离
                              l = 15))

ggsave(file.path(output_dir,paste0(diffname,".top_25_ppi_network.png",collapse = "")), plot = p, dpi=300,height = 7, width = 7.5,bg = "white")
ggsave(file.path(output_dir,paste0(diffname,".top_25_ppi_network.pdf",collapse = "")), plot = p, height = 7, width = 7.5)

