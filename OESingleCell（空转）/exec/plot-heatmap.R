sub_heatmap <- subparsers$add_parser("heatmap", help = "visualize the bunch of features in a heatmap using ComplexHeatmap.")
sub_heatmap$add_argument("-l", "--markers", type = "character",
                         help = "the file of marker genes table to be visulized.")
sub_heatmap$add_argument("-g", "--groupby", type = "character", default = "ident",
                         help = "the groupping variable name for columns in column metadata.")
sub_heatmap$add_argument("-n", "--topn", type = "integer", default = 5,
                         help = "the number of top markers for each cluster to visualizse.")
sub_heatmap$add_argument("-c", "--topby", type = "character", default = "avg_log2FC",
                         help = "the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.")
sub_heatmap$add_argument("-x", "--extraGene", type = "character", default = NULL,
                         help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.")
sub_heatmap$add_argument("--group-colors", type = "character", default = "blindless",
                         help = "one or more color schema names list for each groupping variables including the additional groups if supplied.[default: %(default)s]")
sub_heatmap$add_argument("--extra-groups", type = "character", default = NULL,
                         help = "the comma separated names list of additional groupping variables for columns in column metadata.")
sub_heatmap$add_argument("-r", "--group-orderby", type = "character", default = NULL,
                         help = "one or more groups(comma separated) names to order the columns of matrix.")
sub_heatmap$add_argument("-b", "--col-splitby", type = "character", default = "ident",
                         help = "[OPTIONAL]one of the groupping variables in '--group-by' and 'extra-groups' used to split the columns of matrix to display.")
sub_heatmap$add_argument("--row-annovar", type = "character", default = NULL,
                         help = "[OPTIONAL]the row/features annotation varibale name in the row/features annotation metadata.")
sub_heatmap$add_argument("--row-splitby", type = "character", default = NULL,
                         help = "[OPTIONAL]the row annotation variable used to split the row of the matrix to display.")
sub_heatmap$add_argument("--row-colors", type = "character", default = NULL,
                         help = "[OPTIONAL]the color schema name for row annotations.")
sub_heatmap$add_argument("--row-marks", type = "character", default = NULL,
                         help = "[OPTIONAL]the additional row/features names list to highlight with padding style.")
sub_heatmap$add_argument("--sample-ratio", type = "double", default = NULL,
                         help = "[OPTIONAL]the ratio of random subsample for each group when drawing heatmap.")
sub_heatmap$add_argument("--style", type = "character", default = "seurat",
                         help = "[OPTIOANL]the plot style for heatmap, seurat or complexheatmap.[default: %(default)s]")
sub_heatmap$add_argument("--row-txtsize", type = "character", default = 5.5,
                         help = "the size of row text labels.")
sub_heatmap$add_argument("--value-colors", type = "character", default = "#FF08BD,#000000,#FFF42F",
                         help = "the color scale used mapping the value of each cell in the matrix.[default: %(default)s] Format can be:low_color, middle_color, high_color.")

args <- commandArgs(TRUE)
if ("heatmap" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "heatmap" ){
    if (is.null(opt$markers) & is.null(opt$extraGene)) { stop("NO marker genes is AVAILABLE!") }
    # prepare markers to be visualized
    topn_markers <- data.frame()
    if (!is.null(opt$markers)) {
      markers2vis <- read.table(opt$markers, sep = "\t", header = T)
      topn_markers <- markers2vis %>%
                      dplyr::group_by(cluster) %>%
                      # dplyr::arrange(p_val,desc(avg_log2FC),desc(gene_diff)) %>%
                      dplyr::arrange(p_val, desc(avg_log2FC)) %>%
                      dplyr::top_n(n = as.integer(opt$topn), wt = !!rlang::sym(opt$topby)) %>%
                      dplyr::mutate(folder_suffix = paste0("cluster", cluster)) %>%
                      dplyr::ungroup() %>%
                      # dplyr::arrange(cluster) %>%
                      dplyr::select(cluster, gene, folder_suffix)
                      # dplyr::select(cluster,gene)
    }
    # add extra customized genes if existis
    if (!is.null(opt$extraGene)) {
      extra_gene <- read.table(opt$extraGene, sep = "\t", header = T)
      if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene") { colnames(extra_gene)[1] <- "extra" }
      formated_extra_gene <- as.data.frame(tidyr::gather(extra_gene, key = "cluster", value = "GENE"))
      match <- OESingleCell::CaseMatch(search = as.vector(formated_extra_gene$GENE), match = rownames(data_ob))
      formated_extra_gene <-  formated_extra_gene %>%
                              dplyr::filter(GENE %in% names(match)) %>%
                              dplyr::mutate(gene = match, folder_suffix = cluster) %>%
                              dplyr::select(cluster, gene, folder_suffix)
      topn_markers <- rbind(topn_markers, formated_extra_gene)
    }

    markers2vis4heatmap <- unique(as.vector(topn_markers$gene))
    # read the specified assay and data slot in data object into memory
    # for visualization usually only "data" slot is loaded,which can save a lot of time and memory
    data_ob <- OESingleCell::ReadX( input = opt$input,
                                    informat = opt$informat,
                                    assays = assays,
                                    data.use = dataslots, # "data","scale.data"
                                    reductions = opt$reduct,
                                    graphs = FALSE,
                                    images = FALSE, # no graph object needed here
                                    verbose = FALSE)
    #get subset of cells used col.name = col.name = or visualization if necessay
    if (!is.null(opt$predicate)) {
      df <- OESingleCell::colData(data_ob)
      desired_cells <- subset(df, eval(parse(text = opt$predicate)))
      print(desired_cells)
      data_ob <- subset(data_ob,cells=rownames(desired_cells))
    }
    # downsample to specified ratio in case of too many cells to plot
    if (!is.null(opt$sample_ratio)) {
      sampled_cellmeta <- cellmeta %>%
                          tibble::rownames_to_column() %>%
                          dplyr::group_by(.dots = opt$groupby) %>%
                          dplyr::sample_frac(size = opt$sample_ratio, replace = F) %>%
                          dplyr::column_to_rownames()
      data_ob <- data_ob[, rownames(sampled_cellmeta)]
    }

    cellmeta <- OESingleCell::colData(data_ob)

    if (!is.null(opt$group_orderby)) {
      opt$group_orderby <- opt$groupby
    }
    clusters <- cellmeta[[opt$groupby]]
    # data_ob[[opt$groupby]] = factor(clusters, levels = sort( unique(as.numeric(clusters))))
    # data_ob[[groupby]] = as.factor(data_ob[[groupby]])

    # color settings
    group_colors <- unlist(strsplit(opt$group_colors, ",", perl = T))
    if (!is.null(opt$extra_groups)) {
      extra_groups <- unlist(strsplit(opt$extra_groups, ",", perl = T))
    }else {
      extra_groups <- NULL
    }
    names(group_colors) <- union(opt$groupby, extra_groups)

    # if ( !is.null(opt$row_annovar) ){
    # }

    if (!is.null(opt$row_marks)) {
      row_markers <- read.table(opt$extraGene, sep = "\t", header = T)
      row_markers <- as.vector(extra_gene[, 1])
    }else {
      row_markers <- NULL
    }

    value_colors <- unlist(strsplit(opt$value_colors, ",", perl = T))
    if (style == "seurat") {
      subseted_seurat <- subseted_seurat %>% Seurat::ScaleData(assay = "SCT", features = markers2vis4heatmap)
      ggheat <- DoHeatmap(object = subseted_seurat,
                          features = markers2vis4heatmap,
                          group.colors = colors2use,
                          group.by = opt$groupby, group.bar = T, label = F) +
                theme(axis.text.y = element_text(size = 4))

    }else {
      ggcheat <- OESingleCell::DoComplexHeatmap(data_ob,
                                                features = markers2vis4heatmap,
                                                group.by = opt$groupby,
                                                cell.annotation = extra_groups,
                                                group.colors = as.list(group_colors),
                                                group.order = opt$group_orerby,
                                                col.split.by = opt$col_splitby,
                                                x.size = opt$row_txtsize,
                                                gene.marks = row_markers,
                                                value.colors = circlize::colorRamp2(c(-1.5, 0, 1.5), value_colors),
                                                column_title = NULL)
      ggheat <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(ggheat, padding = grid::unit(c(2, 3, 2, 20), "mm"))))
    }
    ###
    OESingleCell::save_ggplots(file.path(output_dir, "marker_gene_heatmap"),
                               plot = ggheat,
                               dpi = 1000,
                               limitsize = F)
    ##
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }}