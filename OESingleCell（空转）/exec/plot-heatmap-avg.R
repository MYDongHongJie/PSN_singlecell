### example
docstring <- "example1:\\n\\n\\
  scVis -i  seurat.rds  -f rds --dataslots counts --assay SCT -o results avg_heatmap -g gene.list  -s TRUE  -c clusters -f FALSE "
sub_avg_heatmap <- subparsers$add_parser("avg_heatmap",
                                         description = docstring,
                                         formatter_class = 'argparse.RawTextHelpFormatter',
                                         #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                         argument_default = "True",
                                         help = "visualize the bunch of features in a heatmap using ComplexHeatmap.")
sub_avg_heatmap$add_argument("-g", "--genelist", type = "character",
                             help = "The genelist used to generate the heatmap with header.")
sub_avg_heatmap$add_argument( "--filter", type = "logical", default = FALSE,
                             help = "filter gene list by log2FoldChange&topn or not: TRUE or FALSE.[default: FALSE]")
sub_avg_heatmap$add_argument("-c", "--averageby", type = "character",
                             help = "The variable level of heatmap. e.g. clusters or sampleid or celltype")
sub_avg_heatmap$add_argument("-r", "--rowcluster", type = "logical", default = FALSE,
                             help = "boolean values determining if rows should be clustered or hclust object.[default: FALSE]")
sub_avg_heatmap$add_argument("-l", "--colcluster", type = "logical", default = FALSE,
                             help = "boolean values determining if columns should be clustered or hclust object.[default: FALSE]")
sub_avg_heatmap$add_argument("-s", "--scale", type = "logical", default = TRUE,
                             help = "boolean values determining if row scale or not,if the dataslots is scale.data, no need.[default: TRUE]")
sub_avg_heatmap$add_argument("-n", "--showname", type = "logical", default =TRUE,
                             help = "Whether to display the row name,TRUE or FALSE.[default: TRUE]")
# sub_avg_heatmap$add_argument("-w", "--gaps_row",type = "character", default = NULL,
#            help="Whether to display gaps by row")
sub_avg_heatmap$add_argument("--topn", type = "integer", default = 25,
                             help = "[OPTIONAL] the number of top DEGs to visualizse.[default: %(default)s]")
sub_avg_heatmap$add_argument("--order", type = "character",default=NULL,
                             help = paste0("[OPTIONAL] specify the order of colnames to present from left to right.",
                                           "OR can also be used to show subset of levels for factors specifyed by --collapedby.[default:NULL]"))
sub_avg_heatmap$add_argument("--var2use", "-q", type = "character", default = "clusters",
                             help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with levels4var.[default: %(default)s]")
sub_avg_heatmap$add_argument("--levels4var", "-u", type = "character", default = NULL,
                             help = "[OPTIONAL] subset of factor levels for the specified factor by --var2use.[default:NULL]")

args <- commandArgs(TRUE)
if ("avg_heatmap" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "avg_heatmap") {
    #### load query object data=========================================================================================
    seurat_ob <- OESingleCell::ReadX(input = opt$input,
                                     informat = opt$informat,
                                     assays = assays[1],
                                     data.use = dataslots,
                                     verbose = F)
    Seurat::DefaultAssay(seurat_ob) <- assays
    ###get the subset of cells used for visualization if necessary =====================================================
    if (!is.null(opt$levels4var)) {
      seurat_ob <- base::subset(seurat_ob, subset = opt$levels4var == unlist(strsplit(opt$levels4var, ",", perl = T)))
      seurat_ob@meta.data[[ident2use]] <- factor(seurat_ob@meta.data[[opt$var2use]],
                                                 levels = sort(unique(seurat_ob@meta.data[[opt$var2use]])))
    }
    ###=================================================================================================================
    ## filter gene by some additional conditions
    genelist <- readr::read_tsv(opt$genelist,show_col_types=FALSE)
    if (opt$filter) {
      up <- dplyr::filter(genelist, FoldChange > 1) %>%
            dplyr::arrange(desc(log2FoldChange)) %>%
            dplyr::top_n(opt$topn, log2FoldChange) %>%
            dplyr::select(gene)
      down <- dplyr::filter(genelist, FoldChange < 1) %>%
              dplyr::arrange(log2FoldChange) %>%
              dplyr::top_n(as.numeric(paste0("-", opt$topn)), log2FoldChange) %>%
              dplyr::select(gene)
      genelist <- rbind(up, down)
      #genelist[,1] <- factor(as.character(genelist[, 1]))
    }
    ## check gene
    match <- Seurat::CaseMatch(search = as.vector(genelist$gene), match = rownames(seurat_ob))
    matched_gene <- genelist %>% dplyr::filter(gene %in% names(match))
    unmatched_gene <- genelist %>% dplyr::filter(!gene %in% names(match))
    if (dim(unmatched_gene)[1] > 0) {
      futile.logger::flog.warn("There are some mismatched gene symbol, Please check filtered_gene.xls for the gene name.")
      unmatched_gene %>%
        as.data.frame() %>%
        readr::write_tsv(file.path(output_dir, "unmatched_gene.xls"))
    }
    if (dim(matched_gene)[1] < 1) stop("The number of matching gene is smaller than 2, please check the input gene list.")
    ### get average counts expression data from seurat object===========================================================

    gene_data <- Seurat::AverageExpression(seurat_ob,
                                           assays = assays[1],
                                           features = matched_gene$gene,
                                           return.seurat = FALSE,
                                           group.by = opt$averageby,
                                           slot = dataslots,
                                           verbose = TRUE) %>%
                 as.data.frame   #%>%
    # dplyr::rename_at( vars(contains( assays[1])),
    #                         dplyr::funs(stringr::str_replace(., paste0(assays[1],"."), "")))
    colnames(gene_data) <- gsub(paste0(assays[1], "."), "", colnames(gene_data))
    if (!is.null(opt$order)) {
      gene_data <- gene_data %>%
                   tibble::rownames_to_column("gene") %>%
                   dplyr::select( c("gene",unlist(strsplit(opt$order, ",", perl = T)))) %>%
                   tibble::column_to_rownames("gene")
    }
    ###
    colnames(gene_data) <- paste0(opt$averageby, "_", colnames(gene_data))
    readr::write_tsv(gene_data %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column("gene"),
                     file.path(output_dir, paste0("heatmap_", opt$averageby, "_counts.xls")))
    ##scale data by row===============================================================================================
    if (opt$scale) {
      futile.logger::flog.appender("scaling data")
      if (assays[1] == "SCT") {
        gene_data <- gene_data %>%
          as.matrix() %>%
          log1p() %>%
          t() %>%
          scale(center = TRUE, scale = TRUE) %>%
          t()
      }else {
        gene_data <- gene_data %>%
          as.matrix() %>%
          +0.00001 %>%
          log(base = exp(1)) %>%
          t() %>%
          scale(center = TRUE, scale = TRUE) %>%
          t()
      }
      readr::write_tsv(gene_data %>%
                         as.data.frame() %>%
                         tibble::rownames_to_column("gene"),
                       file.path(output_dir, paste0("heatmap_", opt$averageby, "_log_counts_scaled.xls")))
    }

    ## ploting =========================================================================================================
    pdf(NULL) ##for removing  Rplot.pdf files
    showname <- ifelse(dim(gene_data)[1] > 100 && opt$showname != TRUE, FALSE, TRUE)
    plot <- ComplexHeatmap::Heatmap(as.matrix(gene_data), # 输入数据为矩阵
                                    name = 'heatmap for gene', # 热图图例名称
                                    col = grDevices::colorRampPalette(c("blue", "white", "red"))(n = 299), #配色方案
                                    border = FALSE, # 显示边框
                                    cluster_rows = opt$rowcluster, #按照行进行聚类
                                    cluster_columns = opt$colcluster, #按照列进行聚类
                                    #column_km = 2, # 划分列聚类
                                    #row_km = 2, # 划分行聚类
                                    show_column_names = T, #显示列名
                                    show_row_names = showname, #显示行名与否
                                    column_names_side = ifelse(opt$colcluster == "FALSE", "top", "bottom"), , #列名位置 top bottom
                                    row_names_rot = 0, ##行名旋转角度
                                    column_names_rot = 45, ##列名旋转角度
                                    row_names_side = ifelse(opt$rowcluster == "FALSE", "left", "right"), #行名位置 left right
                                    column_names_gp = grid::gpar(fontsize = 10), ##列字体大小
                                    row_names_gp = grid::gpar(fontsize = 10), ##行字体大小
                                    #alpha=0.5, #设置图片颜色透明度
                                    # width = ncol(gene_data)*unit(5, "mm"), # 格子的宽度
                                    # height = nrow(gene_data)*unit(2, "mm"),# 格子的高度
                                    top_annotation = NULL, # 添加左侧注释信息
                                    right_annotation = NULL, # 添加右侧注释信息
                                    left_annotation = NULL, # 添加左侧注释信息
                                    bottom_annotation = NULL, # 添加下方注释信息
                                    heatmap_legend_param = list(title = "",
                                                                #title_position = "topcenter", # 标题相对图例的位置 topcenter, topleft, leftcenter, lefttop.
                                                                #at=c(0,10), #图例范围
                                                                labels_gp = grid::gpar(fontsize = 10), ##图例字体
                                                                legend_direction = "vertical", ##图例放置
                                                                legend_height = unit(5, "cm") #图例长度
                                    ))
    ###save plot =======================================================================================================
    OESingleCell::save_ggplots(file.path(output_dir, paste0("heatmap_", opt$averageby, "_average")),
                               plot = ggplotify::as.ggplot(plot),
                               width = 36 * (dim(gene_data)[2] +
                                 1 +
                                 ceiling(max(nchar(rownames(gene_data))) / 5)) / 72,
                               height = ifelse(showname,
                                               12 * (dim(gene_data)[1] +
                                                 1 +
                                                 ceiling(max(nchar(colnames(gene_data))) / 5)) / 72,
                                               (8 * dim(gene_data)[1] + 108) / 72))
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)

    quit()
  }
}
