#
sub_dotplot = subparsers$add_parser("dotplot", help = "visualize the marker using dotplot")
sub_dotplot$add_argument("-l", "--markers",type ="character",
        help="the file of marker genes table to be visulized.")
sub_dotplot$add_argument("-n", "--topn", type="integer", default = 1,
        help = "the number of top ranked markers for each cluster to visualizse.")
sub_dotplot$add_argument("-c", "--topby", type = "character", default = "avg_log2FC",
        help="the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.")
sub_dotplot$add_argument("-x", "--extraGene", type = "character", default = NULL,
        help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.")
sub_dotplot$add_argument("-g", "--groupby", type = "character", default = "clusters",
        help = "[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize marker genes.")
sub_dotplot$add_argument("-y", "--splitby", type = "character", default = NULL,
        help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
sub_dotplot$add_argument("-s", "--pointsize", type = "double", default = 1,
        help = "[OPTIONAL]the point size in the plot. If you prefer no points, please set it as 0.")
sub_dotplot$add_argument( "--row-splitby", type = "character", default = NULL,
          help = "the row annotation variable used to split the row of the matrix to display.")
args <- commandArgs(TRUE)
if ( "dotplot" %in% args ){
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
  # data_ob = Seurat::SetIdent( data_ob, value = groupby )
  Seurat::Idents(data_ob) = groupby
  markers2vis4dotplot = unique(as.vector(topn_markers$gene))
  ggdot = Seurat::DotPlot(object = data_ob, features = markers2vis4dotplot ) + RotatedAxis()
  OESingleCell::save_ggplots(file.path(output_dir,"marker_gene_dotplot.png"),
                            plot = ggdot, dpi = 1000 ,limitsize = F)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}
