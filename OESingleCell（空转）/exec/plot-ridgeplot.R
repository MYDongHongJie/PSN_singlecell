## usage help:
docstring<-"example1(Select the topn marker gene from findAllmarker results and Visiual):\\n\\n\\
  scVis -i seurat.rds -f rds -o results --assay Spatial featureplot  --pt.alpha TRUE  "
sub_ridgeplot<- subparsers$add_parser("ridgeplot",
                                       description = docstring,
                                       formatter_class= 'argparse.RawTextHelpFormatter' ,
                                       #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                       argument_default = "True" ,
                                       help = "visualize the bunch of features in a featureplot.")
sub_ridgeplot$add_argument("-l", "--markers", type = "character",
                          help = "the file of marker genes table to be visulized.")
sub_ridgeplot$add_argument("-n", "--topn", type = "integer", default = 1,
                          help = "the number of top markers for each cluster to visualizse.")
sub_ridgeplot$add_argument("-c", "--topby", type = "character", default = "avg_log2FC",
                          help = "the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.")
sub_ridgeplot$add_argument("-x", "--extraGene", type = "character", default = NULL,
                          help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.")
sub_ridgeplot$add_argument("-g", "--groupby", type = "character", default = "clusters",
                          help = "[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize marker genes.")
sub_ridgeplot$add_argument("-y", "--splitby", type = "character", default = NULL,
                          help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
sub_ridgeplot$add_argument("--colorscale", type = "character", default = "grey,red",
                          help = "[OPTIONAL]the color scale used to map the gene expression value,format:low_color,high_color.[default: %(default)s]")
sub_ridgeplot$add_argument("-s", "--pointsize", type = "double", default = 1,
                          help = "[OPTIONAL]the point size in the plot.")
# sub_ridgeplot$add_argument("--reduct", type = "character", default = "umap",
#                           help = "[OPTIONAL]the previous calculated reduction result used in the featureplot,.")
args <- commandArgs(TRUE)
if ("ridgeplot" %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "ridgeplot" ){
  # step1.prepare markers to be visualized==============================================================================
    features <- data.frame()
    if (!is.null(opt$features)) {
      ## for feature list，eg: diff_markers
      if(file.exists(opt$features & ! is.null(opt$topn) & is.null((opt$topby) ) )){
        features2vis <- read.table(opt$features, sep = "\t", header = T)
        topn_markers <- markers2vis %>%
                        dplyr::group_by(cluster) %>%
                        dplyr::arrange(p_val, desc(avg_log2FC), desc(gene_diff)) %>%
                        dplyr::top_n(opt$topn, opt$topby) %>%
                        dplyr::arrange(cluster) %>%
                        dplyr::mutate(folder_suffix = paste0("cluster", cluster)) %>%
                        dplyr::select(cluster, gene, folder_suffix)
      } else{

      }}

    }
    # step2.add extra customized genes if existis=======================================================================
    if (!is.null(opt$extraGene)) {
      extra_gene <- read.table(opt$extraGene, sep = "\t", header = T)
      if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene") { colnames(extra_gene)[1] <- "extra" }
      formated_extra_gene <- as.data.frame(tidyr::gather(extra_gene, key = "cluster", value = "GENE"))
      match <- OESingleCell::CaseMatch(search = as.vector(formated_extra_gene$GENE), match = rownames(data_ob))
      formated_extra_gene <- formated_extra_gene %>%
                             dplyr::filter(GENE %in% names(match)) %>%
                             dplyr::mutate(gene = match, folder_suffix = cluster) %>%
                             dplyr::select(cluster, gene, folder_suffix)
      topn_markers <- rbind(topn_markers, formated_extra_gene)
    }
    cellmeta <- OESingleCell::colData(data_ob)
    # step3.set color for visualizarion=================================================================================
    colors2use <- OESingleCell::SelectColors(1:length(unique(cellmeta[[groupby]])), palette = opt$palette)
    for (clusterx in unique(topn_markers$folder_suffix)) {
      topn <- topn_markers %>%
              filter(folder_suffix == clusterx) %>%
              select(cluster, gene, folder_suffix)
      topn_markers2vis <- as.vector(topn$gene)

      path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
      if (file.exists(path4vis)) {
        output_dir <- path4vis
      }else {
        output_dir <- path4vis
        dir.create(output_dir)
      }
      # data_ob = Seurat::SetIdent( data_ob, value = groupby )===========================================================
      Seurat::Idents(data_ob) <- groupby
      ggridge <- Seurat::RidgePlot(object = data_ob,
                                   features = topn_markers2vis,
                                   cols = colors2use,
                                   nCol = 3)
      OESingleCell::save_ggplots(file.path(output_dir, "marker_gene_ridgeplot"),
                                 plot = ggridge,
                                 dpi = 1000,
                                 limitsize = F)
    }
  ###
  write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
  quit()
}