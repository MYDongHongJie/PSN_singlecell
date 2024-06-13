suppressWarnings({
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("argparse") )
    suppressPackageStartupMessages(library("ggplot2"))
    suppressPackageStartupMessages(library("RColorBrewer"))
})

parser = ArgumentParser(description = "VlnPlot.",
                        usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input", type = "character",
             help = "The input seurat object in several possible format.")
parser$add_argument("-f", "--informat", type = "character", default = "h5seurat",
             help = "The format of data object, the possible choices can be:h5seurat,(seurat)rds,(sce)rds, loom.[default: %(default)s]")
parser$add_argument("-x", "--extraGene", type = "character", default = NULL,
             help = "The extra gene list of interest to visualize specified by the user.")
parser$add_argument("-g", "--groupby", type = "character", default = NULL,
             help = "The grouppinig variable in the metadata for separate the cells to visulize marker genes.")
parser$add_argument("-y", "--splitby", type = "character", default = NULL,
             help = "The variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
parser$add_argument("-s", "--pointsize", type = "double", default = NULL,
             help = "the point size in the plot.")
parser$add_argument("-o", "--output", type = "character", default = "./",
             help = "the output directory of results."  )
opt = parser$parse_args()

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive = T)
    }
}

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

data_ob = OESingleCell::ReadX(input = opt$input, informat = opt$informat,verbose = F)
gene_name = read.delim(opt$extraGene, sep = "\t", header = T)
match = CaseMatch(search = as.vector(gene_name$gene),match = rownames(data_ob))

ymaxs = ceiling(max(data_ob@assays$RNA@data[match,]))

modify_vlnplot <- function(object, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),color.use = NULL,cols = color.use,...) {
       p <- VlnPlot(object, features = features, pt.size = pt.size, cols = color.use,... ) +
               xlab("") + ylab("Expression Levels") + ggtitle(features) +
               theme(legend.position = "right",
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_line(),
               axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5),
               plot.margin = plot.margin ) +
               guides(fill=guide_legend(title="Group"))
       return(p)
}

## main function
StackedVlnPlot <- function(object, features, split.by = NULL, angle.x = 45, vjust.x = NULL, hjust.x = NULL, 
    pt.size = 0, plot.margin = margin(0, 0, 0, 0, "cm"), ymaxs = NULL,color.use=NULL,ncols=ncols,cols = color.use,...) {
    if (is.null(vjust.x) | is.null(hjust.x)) {
        angle = c(0, 45, 90)
        hjust = c(0, 1, 1)
        vjust = c(0, 1, 0.5)
        vjust.x = vjust[angle == angle.x]
        hjust.x = hjust[angle == angle.x]
    }
    color.use <- CustomCol2(1:ncols)
    plot_list <- purrr::map(features, function(x) modify_vlnplot(object = object, 
        features = x, split.by = split.by, pt.size = pt.size, cols = color.use,...))
    plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] + 
        theme(axis.text.x = element_text(), axis.ticks.x = element_line()) + 
        theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
            vjust = vjust.x)) + theme(axis.text.x = element_text(size = 10))
    plot_list <- purrr::map2(plot_list, ymaxs, function(x, y) x + 
        scale_y_continuous(breaks = seq(0,ymaxs,1)) + expand_limits(y = y)) 
    p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
}

ncols <- length(unique(data_ob@meta.data$group))

for (i in match){
      sum_all <- StackedVlnPlot(object = data_ob, 
                                features = i, 
                                group.by = opt$groupby,
                                split.by= opt$splitby,
                                pt.size= opt$pointsize ,
                                ymaxs=ymaxs,
                                ncols=ncols)
      ggsave(file.path('./',paste0(i,"_vlnplot.pdf",collapse="")),width=5 + 0.5*(length(unique(data_ob@meta.data[,opt$groupby]))),height=4,plot=sum_all)
      ggsave(file.path('./',paste0(i,"_vlnplot.png",collapse="")),width=5 + 0.5*(length(unique(data_ob@meta.data[,opt$groupby]))),height=4,plot=sum_all)
}
