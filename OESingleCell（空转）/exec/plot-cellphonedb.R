sub_cellphonedb_pt <- subparsers$add_parser("cellphonedb_pt", help = "plot for cellphonedb result.")
sub_cellphonedb_pt$add_argument("-d", "--plot", type = "character", default = "dotplot",
                                help = "[REQUIRED]the comma seperated plot type list to visualize the cell communications,currently dotplot, network and circos are supported.[default: clusters]")
sub_cellphonedb_pt$add_argument("--pvalue", type = "double", default = 0.05,
                                help = "[REQUIRED]the comma seperated plot type list to visualize the cell communications,currently dotplot, network and circos are supported.[default: %(default)s]")
sub_cellphonedb_pt$add_argument("-n", "--topn", type = "integer", default = 5,
                                help = "[OPTIONAL]the top N interacting gene pairs for each cell pairs. This is only setted for circos plot.[default: %(default)s]")
sub_cellphonedb_pt$add_argument("--topby", type = "character", default = "expr",
                                help = "[OPTIONAL]the variable used to order the interacting pairs for circos plot.[default: %(default)s]")
sub_cellphonedb_pt$add_argument("-v", "--onlySigs", type = "character", default = "true",
                                help = "[OPTIONAL]only keep the significant gene pairs in all interacting cell pairs for dot plot.[default: %(default)s]")
sub_cellphonedb_pt$add_argument("-c", "--colorschema", type = "character", default = "customecol2",
                                help = "[OPTIONAL]the supported color schema currently:seurat,blindless,paried,colx22,jet,tableau20,tableau10medium,colorblind10,cyclic.[default: %(default)s]")
sub_cellphonedb_pt$add_argument("--heatmapcol", type = "character", default = NULL,
                                help = "[OPTIONAL]if assign the color of heatmap,can do like this:'60,40,0+#fb4f4f,#fbc93d,#6cc0e5'[default: %(default)s]")

args <- commandArgs(TRUE)
if ("cellphonedb_pt" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "cellphonedb_pt") {
    if (is.null(opt$output)) {
      print("NO output directory specified,the current directory will be used!")
      output_dir <- getwd()
    }else {
      if (file.exists(opt$output)) {
        output_dir <- opt$output
      }else {
        output_dir <- opt$output
        dir.create(output_dir)
      }
    }

    if (tolower(opt$onlySigs) == "true") {
      is_onlySig <- TRUE
    }else {
      is_onlySig <- FALSE
    }

    output_dir <- normalizePath(output_dir)


    if (is.null(opt$informat) || opt$informat == "rds") {
      cellp <- readRDS(opt$input)
      data <- cellp$ligrec
    }else if (opt$informat == "matrix") {
      data <- read.table(opt$input, header = T, sep = "\t", stringsAsFactors = F)
    }

    if (is.null(opt$colorschema)) {
      colorschema <- "customecol2"
    }else {
      colorschema <- opt$colorschema
    }

    if (is.null(opt$pvalue)) {
      pvalue <- 0.05
    }else {
      pvalue <- opt$pvalue
    }

    if (is.null(opt$topn)) {
      topn <- 5
    }else {
      topn <- opt$topn
    }

    if (is.null(opt$topby)) {
      topby <- "expr"
    }else {
      topby <- opt$topby
    }

    for (plotx in unlist(strsplit(opt$plot, ",", perl = T))) {

      if (plotx == "dotplot") {
        ggdot <- OESingleCell::LRDotplot(data = data, is_onlySig = is_onlySig, xsize = 10, xangle = 45)
        png_width <- length(unique(data$receptor_cell))
        png_height <- length(unique(data$receptor))
        ggplot2::ggsave(file.path(output_dir, "cell_comm_dotplot.pdf"), limitsize = F, plot = ggdot, width = png_width * 4, height = png_height / png_width * 5)
        ggplot2::ggsave(file.path(output_dir, "cell_comm_dotplot.png"), limitsize = F, plot = ggdot, width = png_width * 4, height = png_height / png_width * 5)
      }

      if (plotx == "network") {
        suppressPackageStartupMessages(library("network"))
        suppressPackageStartupMessages(library("igraph"))
        celltypes <- unique(c(data$receptor_cell, data$ligand_cell))
        colx <- OESingleCell::SelectColors(palette = colorschema, n = length(celltypes))
        pdf(file.path(output_dir, "cell_comm_network.pdf"))
        OESingleCell::LRNetwork(data %>% dplyr::filter(pval < pvalue),
                                col = colx,
                                edge.label.cex = 0.3,
                                edge.max.width = 5,
                                arrow.width = 0.5,
                                vertex.label.cex = 0.5,
                                vertex.size = 10)
        dev.off()
        # summary
        net <- data %>%
          dplyr::filter(pval < pvalue) %>%
          dplyr::group_by(ligand_cell, receptor_cell) %>%
          dplyr::summarize(n = dplyr::n()) %>%
          dplyr::rename(significant_pairs = n)
        write.table(net, file.path(output_dir, "cell_comm_network_summary.xls"), quote = F, sep = "\t", row.names = F)
      }

      if (plotx == "circos") {
        suppressPackageStartupMessages(library("circlize"))
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        celltypes <- unique(c(data$receptor_cell, data$ligand_cell))
        colx <- OESingleCell::SelectColors(palette = colorschema, n = length(celltypes))
        pdf(file.path(output_dir, "cell_comm_circos_plot.pdf"), width = 20, height = 10)
        OESingleCell::LRCircos(data, gap.degree = 0.05, cell_col = colx, screenvar = topby,
                               topn = opt$topn, labels.cex = 0.35, link.lwd = 1, arr.length = 0.1)
        dev.off()
        png(file.path(output_dir, "cell_comm_circos_plot.png"), width = 1600, height = 800)
        OESingleCell::LRCircos(data, gap.degree = 0.05, cell_col = colx, screenvar = topby,
                               topn = opt$topn, labels.cex = 0.35, link.lwd = 1, arr.length = 0.1)
        dev.off()
      }
      if (plotx == "chorddiagram") {
        suppressPackageStartupMessages(library("circlize"))
        pdf(file.path(output_dir, "cell_comm_chorddiagram_plot.pdf"))
        OESingleCell::LRChorddiagram(data, grid_orbit_h = 0.02, label_h = 0.04, diffHeight = 4, output_dir = output_dir)
        dev.off()
        png(file.path(output_dir, "cell_comm_chorddiagram_plot.png"), width = 7, height = 7, res = 96, units = "in")
        OESingleCell::LRChorddiagram(data, grid_orbit_h = 0.02, label_h = 0.04, diffHeight = 4, output_dir = output_dir)
        dev.off()
        circos.clear()
      }

      if (plotx == "heatmap") {
        suppressPackageStartupMessages(library("ComplexHeatmap"))
        suppressPackageStartupMessages(library("RColorBrewer"))
        pdf(file.path(output_dir, "cell_comm_heatmap_plot.pdf"))
        print(OESingleCell::LRHeatmap(data=data,col4heatmap=opt$heatmapcol))
        dev.off()
        png(file.path(output_dir, "cell_comm_heatmap_plot.png"), width = 7, height = 7, res = 96, units = "in")
        print(OESingleCell::LRHeatmap(data=data,col4heatmap=opt$heatmapcol))
        dev.off()
      }

      if (plotx == "barplot") {
        celltypes <- unique(c(data$receptor_cell, data$ligand_cell))
        colx <- OESingleCell::SelectColors(palette = opt$colorschema, n = length(celltypes))
        names(colx) <- celltypes
        out <- OESingleCell::LRBarplot(data = data, bar_width = 0.6)
        ggsave(file.path(output_dir, "cell_comm_histogram_plot.pdf"), plot = out, width = length(celltypes) + 2)
        ggsave(file.path(output_dir, "cell_comm_histogram_plot.png"), plot = out, width = length(celltypes) + 2, dpi = 1000)
      }
    }
  }
}


