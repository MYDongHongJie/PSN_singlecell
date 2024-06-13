sub_vlnplot = subparsers$add_parser("vlnplot", help = "visualize the marker using violin plot.")
sub_vlnplot$add_argument("-l", "--markers",type ="character",
        help="the file of marker genes table to be visulized.")
sub_vlnplot$add_argument("-n", "--topn", type="integer", default = 1,
        help = "the number of top ranked markers for each cluster to visualizse.")
sub_vlnplot$add_argument("-c", "--topby", type = "character", default = "avg_LogFC",
        help="the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.")
sub_vlnplot$add_argument("-x", "--extraGene", type = "character", default = NULL,
        help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.")
sub_vlnplot$add_argument("-g", "--groupby", type = "character", default = "clusters",
        help = "[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize marker genes.")
sub_vlnplot$add_argument("-y", "--splitby", type = "character", default = NULL,
        help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
sub_vlnplot$add_argument("-s", "--pointsize", type = "double", default = 1,
        help = "[OPTIONAL]the point size in the plot. If you prefer no points, please set it as 0.")
sub_vlnplot$add_argument("-a", "--alpha2use", type = "double", default = 0,
        help = "[OPTIONAL]the opacity of the pooints on the violin plot.")
sub_vlnplot$add_argument("--dodge", type = "character", default = "TRUE",
        help = "[OPTIONAL]visualize the feature between the contrast groups separately for each level in variable specified by --groupby.")
args <- commandArgs(TRUE)
if ( "vlonplot"  %in% args ){
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
  # Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)
  for ( clusterx in unique(topn_markers$folder_suffix) ){
    topn = topn_markers %>%
      dplyr::filter( folder_suffix == clusterx) %>%
      dplyr::select(cluster,gene,folder_suffix)
    topn_markers2vis = as.vector(topn$gene)

    path4vis = file.path(root_dir,paste0("markers_vis4",clusterx,collapse = ""))
    if ( file.exists( path4vis ) ){
        output_dir = path4vis
    }else{
        output_dir = path4vis
        dir.create(output_dir)
    }

    colors2use = OESingleCell::SelectColors(1:length(unique(OESingleCell::colData(data_ob)[[groupby]])))
    gs = lapply(topn_markers2vis,
                function(x) Seurat::VlnPlot(data_ob, features = x, cols = colors2use,
                                    pt.size = opt$pointsize, alpha = opt$alpha2use,
                                    group.by = groupby, split.by = splitby,
                                    split.plot = as.logical(opt$dodge)
                                    )+
                labs(title = "",y = x) +
                theme(legend.position = "none",
                # panel.spacing = unit(.05, "lines"),
                axis.title.x = element_text(size = 0),
                axis.title.y = element_text(size = 12),
                axis.text.y=element_text(size = 8)))
    pdf(file.path(output_dir,paste0("marker_gene_violin_plot.pdf", collapse = "_")),
                    width = 8, height = length(topn_markers2vis)*2)
    grid.arrange(grobs = gs, ncol=1)
    dev.off()
    png(file.path(output_dir,paste0("marker_gene_violin_plot.png", collapse = "_")),
                    width = 800,height = length(topn_markers2vis)*120)
    grid.arrange(grobs = gs, ncol=1)
    dev.off()
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}