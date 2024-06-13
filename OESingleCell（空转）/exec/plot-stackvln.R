sub_stackvln = subparsers$add_parser("stackvln", help = "visualize the marker using stacked violin plot.")
sub_stackvln$add_argument("-l", "--markers",type ="character",
        help="the file of marker genes table to be visulized.")
sub_stackvln$add_argument("-n", "--topn", type="integer", default = 1,
        help = "the number of top ranked markers for each cluster to visualizse.")
sub_stackvln$add_argument("-c", "--topby", type = "character", default = "avg_LogFC",
        help="the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.")
sub_stackvln$add_argument("-x", "--extraGene", type = "character", default = NULL,
        help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.")
sub_stackvln$add_argument("-g", "--groupby", type = "character", default = "clusters",
        help = "[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize marker genes.")
sub_stackvln$add_argument("-y", "--splitby", type = "character", default = NULL,
        help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
sub_stackvln$add_argument("-s", "--pointsize", type = "double", default = 1,
        help = "[OPTIONAL]the point size in the plot. If you prefer no points, please set it as 0.")
sub_stackvln$add_argument("-a", "--alpha2use", type = "double", default = 0,
        help = "[OPTIONAL]the opacity of the pooints on the violin plot.")
sub_stackvln$add_argument("--flip", type = "character", default = "TRUE",
        help = "[OPTIONAL]wether to flip the coordination of the plot.")
sub_stackvln$add_argument("--dodge", type = "character", default = "TRUE",
        help = "[OPTIONAL]visualize the feature between the contrast groups separately for each level in variable specified by --groupby.")
args <- commandArgs(TRUE)
if ( "stackvln"  %in% args ){
    opt<-intial_setting()
    # prepare markers to be visualized
    topn_markers = data.frame()
    if ( !is.null(opt$markers) ){
        markers2vis = read.table( opt$markers, sep="\t", header = T )
        topn_markers  = markers2vis %>%
            dplyr::group_by(cluster) %>%
            dplyr::arrange(p_val,desc(avg_log2FC),desc(gene_diff)) %>%
            dplyr::top_n(opt$topn,opt$topby) %>%
            dplyr::arrange(cluster) %>%
            dplyr::mutate(folder_suffix = paste0("cluster",cluster)) %>%
            dplyr::select(cluster,gene,folder_suffix)
    }

    # add extra customized genes if existis
    if ( !is.null(opt$extraGene) ){
      extra_gene = read.table(opt$extraGene, sep="\t", header = T)
      if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
      formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
      match = OESingleCell::CaseMatch(search = as.vector(formated_extra_gene$GENE),match = rownames(data_ob))
      formated_extra_gene = formated_extra_gene %>%
                            dplyr::filter(GENE %in% names(match)) %>%
                            dplyr::mutate(gene = match,folder_suffix = cluster) %>%
                            dplyr::select(cluster, gene,folder_suffix)
      topn_markers = rbind(topn_markers, formated_extra_gene)
    }
    cellmeta = OESingleCell::colData(data_ob)

    # set color for visualizarion
    colors2use = OESingleCell::SelectColors(1:length(unique(cellmeta[[opt$groupby]])), palette = opt$palette)
    gg_stackvln = Seurat::VlnPlot( data_ob, features = topn_markers2vis, group.by = opt$groupby,
                                  split.plot = as.logical(opt$dodge),
                                  cols = colors2use, split.by = opt$splitby, pt_size = opt$pointsize,
                                  flip = as.logical(opt$flip), stack = T, fill.by = groupby)
    OESingleCell::save_ggplots(file.path(output_dir,"marker_gene_stacked_violin_plot"),
                            plot = ggdot, dpi = 1000 ,
                            height = length(topn_markers2vis) * 0.6, width = length(unique(cellmeta[[groupby]]))*0.5+1,
                            limitsize = F)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}