sub_vis <- subparsers$add_parser("visualize", help = "visualize the specified features.")
sub_vis$add_argument("-l", "--markers",type ="character",
        help="the file of marker genes table to be visulized.")
sub_vis$add_argument("-n", "--topn", type="integer", default = 1,
        help = "the number of top markers for each cluster to visualizse.[default: %(default)s] ")
sub_vis$add_argument("-c", "--topby", type = "character", default = "avg_log2FC",
        help="the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.[default: %(default)s] ")
sub_vis$add_argument("-x", "--extraGene", type = "character", default = NULL,
        help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.[default: %(default)s] ")
sub_vis$add_argument("-d","--diffGene",type="character",default=NULL,
        help = "[OPTIONAL]The diff gene list of interest to visualize specified by the user.")
sub_vis$add_argument("-m", "--vismethod",type= "character",default="vlnplot,featureplot",
        help = "the visulization methods for the marker genes. Choices can be: heatmap,diff_heatmap,ridgeplot,vlnplot,dotplot,featureplot.[default: %(default)s] ")
sub_vis$add_argument("-g", "--groupby", type = "character", default = "clusters",
        help = "[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize marker genes.[default: %(default)s] ")
sub_vis$add_argument("--dodge", type = "character", default = "FALSE",
        help = "[OPTIONAL]visualize the feature between the contrast groups separately for each level in variable specified by --groupby.Only valid for violin plot![default: %(default)s] ")
sub_vis$add_argument("-y", "--splitby", type = "character", default = NULL,
        help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.[default: %(default)s] ")
sub_vis$add_argument("--ccolors", type = "character", default = "grey,red",
        help = "[OPTIONAL]the name of customized continious color palatte, recommandations:spectral, solar_extra,flame_light or color scale used to map the continious expression value for feature plot or dotplot,format:low_color,high_color.[default: %(default)s]")
sub_vis$add_argument("--vcolors", type = "character", default = "blindless",
        help = "[OPTIONAL]the name of customized discrete color used to map to the cell groupping variable.[default: %(default)s]Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12[default: %(default)s]")
sub_vis$add_argument("-s", "--pointsize", type = "double", default = 1,
        help = "[OPTIONAL]the point size in the plot.[default: %(default)s] ")
sub_vis$add_argument("--sample_ratio", type = "double", default = 0.6,
        help = "[OPTIONAL]the ratio of random subsample for each group when drawing heatmap.")
sub_vis$add_argument("--reduct", type = "character", default = "umap",
        help = "[OPTIONAL]the previous calculated reduction result used in the featureplot.[default: %(default)s] ")
sub_vis$add_argument("-a", "--alpha2use", type = "double", default = 0,
        help = "[OPTIONAL]the opacity of the pooints on the violin plot.[default: %(default)s] ")

# ====================== Subcmd: visualize the markers in different ways =========
args <- commandArgs(TRUE)
if ( "visualize"  %in% args ){
  opt<-intial_setting()
  if ( is.null( opt$markers ) & is.null(opt$extraGene) & is.null(opt$diffGene)){
      stop("NO marker genes is AVAILABLE!")
  }
 # read the specified assay and data slot in data object into memory
 # for visualization usually only "data" slot is loaded,which can save a lot of time and memory
 data_ob <- OESingleCell::ReadX(
              input = opt$input, informat = opt$informat,
              assays = assays,
              data.use = dataslots, # "data"
              reductions = opt$reduct, # only the used reduction results needed
              graphs = FALSE, # no graph object needed here
              verbose = FALSE)
  #get the subset of cells used  if necessay
  if ( !is.null(opt$predicate) ){
    desired_cells <- subset(OESingleCell::colData(data_ob), eval(parse(text=opt$predicate)))
    # same for Seurat, SingleCellExperiment, summarizedExperiment
    data_ob <- data_ob[, rownames(desired_cells)]
  }

  data_ob[["clusters"]] <- data_ob[[opt$groupby]]
  clusters <- Seurat::FetchData(data_ob, vars = "clusters")[[1]]
  data_ob[["clusters"]] <- factor(clusters, levels = sort(unique(as.numeric(clusters))))

  topn_markers <- data.frame()
  root_dir <- output_dir
  if ( !is.null(opt$markers) ){
      markers2vis <- read.table(opt$markers, sep="\t", header = T )
      topn_markers <- markers2vis %>%
          dplyr::group_by(cluster) %>%
          # dplyr::arrange(p_val, desc(avg_log2FC)) %>%
          dplyr::arrange(p_val,desc(avg_log2FC),desc(gene_diff)) %>%
          dplyr::top_n(opt$topn, wt = !!rlang::sym(opt$topby) ) %>%
          dplyr::arrange(cluster) %>%
          dplyr::mutate(folder_suffix = paste0("cluster",cluster)) %>%
          dplyr::select(cluster,gene,folder_suffix)
  }
  if ( !is.null(opt$extraGene) ){
    extra_gene <- readxl::read_xlsx(opt$extraGene)
    if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1] <- "extra"}
    formated_extra_gene <- as.data.frame(tidyr::gather(extra_gene, key = "cluster", value = "GENE"))
    match <- OESingleCell::CaseMatch(search = as.vector(formated_extra_gene$GENE), match = rownames(data_ob))
    filtered_gene <- formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
    if(length(filtered_gene)!=0){
        filtered_gene <- as.data.frame(filtered_gene)
        colnames(filtered_gene) <- "Gene"
        write.table(filtered_gene,file.path(root_dir,"genes_not_matched.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check genes_not_matched.xls for the genename.")}
    formated_extra_gene <- formated_extra_gene %>%
                          dplyr::filter(GENE %in% names(match)) %>%
                          dplyr::mutate(gene = match,folder_suffix = cluster) %>%
                          dplyr::select(cluster, gene,folder_suffix)
    topn_markers <- rbind(topn_markers, formated_extra_gene)
  }
  if ( !is.null(opt$diffGene) ) {
    diffGene_path <- normalizePath(opt$diffGene)
    diffGene_name <- gsub(".*/", "", diffGene_path)
    if (grepl("-diff-", diffGene_name)) {
      name <- gsub("-diff-.*", "", diffGene_name)
    }else if (grepl("-all_", diffGene_name)) {
      name <- gsub("-all_.*", "", diffGene_name)
    }

    markers2vis <- read.delim(opt$diffGene, sep = "\t", header = T, quote = "")
    markers2vis <- subset(markers2vis, !grepl("^(mt-|Rps|Rpl|MT-|RPS|RPL)", markers2vis$gene))
    up <- dplyr::filter(markers2vis, FoldChange > 1) %>%
      dplyr::arrange(desc(log2FoldChange)) %>%
      dplyr::top_n(opt$topn, log2FoldChange)
    down <- dplyr::filter(markers2vis, FoldChange < 1) %>%
      dplyr::arrange(log2FoldChange) %>%
      dplyr::top_n(as.numeric(paste0("-", opt$topn)), log2FoldChange)
    topn_markers <- rbind(up, down)
    write.table(topn_markers, file.path(root_dir, paste0("top", opt$topn, "_", name, "_genes.xls", collapse = "")), quote = F, row.names = FALSE, sep = "\t")
  }

  groupby <- opt$groupby
  splitby <- opt$splitby
  if (identical(groupby, splitby)) {
    warning("The variable specified by --splitby conflicted with the --groupby parameter, NULL will be used!")
    splitby <- NULL
  }

  if (as.logical(opt$dodge)) {
    if (is.null(splitby)) {
      stop("--dodge should not be set TRUE when --splitby set as NULL!")
    }
  }

  #determine the color schema for discrete and continious variables
  vcolors <- unlist(strsplit(opt$vcolors, ","))
  if ( length(vcolors) == 1 ){ # this means to use customized discrete color schema in OESingleCell package
    if ( is.null(OESingleCell::discrete_palette[[vcolors]]) ){
      info <- glue::glue("the customized discrete colors are: \n {paste(names(OESingleCell::discrete_palette), lengths(OESingleCell::discrete_palette), sep =': ', collapse = '\n')}.")
      warning("NO specified color schema found!")
      warning(info)
      stop()
    }
    discrete_colors <- OESingleCell::SelectColors(palette = vcolors)
  }else{ # this means to use the customized color schema from command line specified by user
    discrete_colors <- vcolors
  }
  ccolors <- unlist(strsplit(opt$ccolors, ","))
  if ( length(ccolors) == 1 ){ # this means to use customized discrete color schema in OESingleCell package
    if ( is.null(OESingleCell::continuous_palette[[ccolors]]) ){
      info <- glue::glue("the customized continious colors are: \n {paste(names(OESingleCell::continuous_palette), lengths(OESingleCell::continuous_palette), sep =': ', collapse = '\n')}.")
      warning("NO specified color schema found!")
      warning(info)
      stop()
    }
      continious_colors <- OESingleCell::SelectColors(palette = ccolors[1], is.discrete = FALSE)
  }else { # this means to use the customized color schema from command line specified by user
    continious_colors <- ccolors
  }

  if (is.null(opt$cluster_name)) {
    data_ob[["clusters"]] <- Seurat::Idents(data_ob)
  }else {
    # data_ob = Seurat::SetIdent( data_ob, value = opt$cluster_name)
    data_ob[[opt$cluster_name]] <- Seurat::Idents(data_ob)
  }
  clusters <- Seurat::FetchData(data_ob, vars = "clusters")[[1]]
  data_ob[["clusters"]] <- factor(clusters, levels = sort(unique(as.numeric(clusters))))


  if (is.null(opt$vismethod)) {
    print("NO marker gene visulization method provided,the default method vlnplot and featureplot will be used!")
    vismethods <- c("vlnplot", "featureplot")
  }else if (opt$vismethod == "all") {
    vismethods <- c("vlnplot", "featureplot", "ridgeplot", "dotplot", "heatmap")
  }else {
    vismethods <- unlist(strsplit(opt$vismethod, ","))
  }

cellmeta <- OESingleCell::colData(data_ob)
for ( vismethod in vismethods ){
  if ( vismethod == "vlnplot" ){
    # Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)
    for ( clusterx in unique(topn_markers$folder_suffix) ){
      topn <- topn_markers %>%
        dplyr::filter( folder_suffix == clusterx) %>%
        dplyr::select(cluster,gene,folder_suffix)
      topn_markers2vis <- as.vector(topn$gene)

      path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
      if ( file.exists( path4vis ) ){
          output_dir <- path4vis
      }else{
          output_dir <- path4vis
          dir.create(output_dir)
      }

      gs <- lapply(topn_markers2vis,
                   function(x) Seurat::VlnPlot(data_ob, features = x,
                                      cols = discrete_colors[1:length(unique(cellmeta[[groupby]]))],
                                      pt.size = opt$pointsize, alpha = opt$alpha2use,
                                      group.by = groupby,
                                      split.by = splitby,
                                      split.plot = as.logical(opt$dodge) )+
                  labs(title = "",y = x, x=NULL) +
                  theme(
                    legend.position="none",
                    legend.margin = margin(0,1,0,0, unit = "cm"),
                    plot.margin = margin(0, 0, 0, 0.1, unit = "cm"),
                      # unit(c(-0.1,2,-0.1,0.1), "null"),
                    # panel.spacing.x = unit(0.1,"null"),
                    axis.text = element_text(margin = grid::unit(0,"null")),
                    axis.title.x = element_text(size = 0),
                    axis.title.y = element_text(size = 12),
                    axis.text.x = element_text(size = 10, angle = 30),
                    axis.text.y=element_text(size = 8)))
      # gga = gridExtra::grid.arrange(grobs = gs, ncol=1, padding = unit(0, "line") )
      ggb <- do.call(ggpubr::ggarrange,
                     c(gs,list( ncol = 1, align = "v", common.legend = TRUE, legend = "right")))
      OESingleCell::save_ggplots(file.path(output_dir,"marker_gene_violin_plot"),
                                 plot = ggb,
                                 limitsize = FALSE,
                                 width = length(unique(cellmeta[[groupby]])),
                                 height = length(topn_markers2vis)*2)
    }
  }

  if ( vismethod == "featureplot" ){
      for ( clusterx in unique(topn_markers$folder_suffix) ){
          if ( !opt$reduct %in% Seurat::Reductions(data_ob) ){
              stop( "NO specified reduction found in the object!")
          }else{
              reduct <- opt$reduct
          }
          topn <- topn_markers %>%
                  dplyr::filter( folder_suffix == clusterx) %>%
                  dplyr::select(cluster,gene,folder_suffix)
          topn_markers2vis <- as.vector(topn$gene)

          path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
          if ( file.exists( path4vis ) ){
              output_dir <- path4vis
          }else{
              output_dir <- path4vis
              dir.create(output_dir)
          }

          suppressMessages({
            ggfeatures <- Seurat::FeaturePlot(data_ob, reduction= reduct,
                                              features = topn_markers2vis, split.by = splitby,
                                              keep.scale = "all",
                                              max.cutoff = 'q99',
                                              order = T,
                                              ncol = ifelse(length(topn_markers2vis) > 4,yes = 4,no=length(topn_markers2vis)),
                                              pt.size = opt$pointsize) &
                          ggplot2::scale_color_gradientn(colours = continious_colors) &
                          theme(legend.position = "none",
                                legend.margin = margin(0,0.5,0,0, unit = "cm"),
                                plot.margin = margin(0, 0, 0, 0.1, unit = "cm"),
                                plot.title = element_text(hjust = 0.5) )
          })
          ggm <- ggpubr::ggarrange(ggfeatures, common.legend = T, legend = "right")
          OESingleCell::save_ggplots(
              file.path(output_dir, "marker_gene_featureplot"),
              width =  4*ifelse(is.null(splitby),yes = ifelse(length(topn_markers2vis) > 4,
                                                                                            yes = 4,
                                                                                            no=length(topn_markers2vis)),
                                                                    no =length(unique(cellmeta[[splitby]]))),
              height = 3*ifelse(is.null(splitby),
                                yes = ceiling(length(topn_markers2vis)/4),
                                no  = length(topn_markers2vis)),
              plot = ggm, dpi = 1000 ,limitsize = F )
      }
  }

  if ( vismethod == "dotplot" ){
    # data_ob = Seurat::SetIdent( data_ob, value = groupby )
    SeuratObject::Idents(data_ob) <- groupby
    markers2vis4dotplot <- unique(as.vector(topn_markers$gene))
    ggdots <- Seurat::DotPlot(object = data_ob,
                              cols = continious_colors,
                              features = markers2vis4dotplot ) +
                    Seurat::RotatedAxis() +
                    ggplot2::coord_flip()
    OESingleCell::save_ggplots( file.path(root_dir,"marker_gene_dotplot"), plot = ggdots,
                                height = length(markers2vis4dotplot) * 0.2,
                                width = length(unique(cellmeta[[groupby]])) * 0.5,
                                dpi = 1000 ,limitsize = F )
  }

  if ( vismethod == "ridgeplot"){
      for ( clusterx in unique(topn_markers$folder_suffix) ){
        topn <- topn_markers %>% filter(folder_suffix == clusterx) %>% select(cluster, gene, folder_suffix)
        topn_markers2vis <- as.vector(topn$gene)

        path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
        if ( file.exists( path4vis ) ){
            output_dir <- path4vis
        }else{
            output_dir <- path4vis
            dir.create(output_dir)
        }
        # data_ob = Seurat::SetIdent( data_ob, value = groupby )
        Seurat::Idents(data_ob) <- groupby
        ggridge <- Seurat::RidgePlot(object = data_ob,
                                     features = topn_markers2vis,
                                     cols = discrete_colors[1:length(unique(cellmeta[[groupby]]))],
                                     ncol = 1)
        OESingleCell::save_ggplots(
          file.path(output_dir,"marker_gene_ridgeplot"),
          plot = ggridge, dpi = 1000 ,limitsize = F )
      }
    }
  }
  if ( vismethod == "heatmap" ) {
    markers2vis4heatmap <- unique(as.vector(topn_markers$gene))
    if (is.null(opt$sample_ratio)) {
      subseted_seurat <- data_ob
    }else {
      sampled_cellmeta <- data_ob@meta.data %>%
        rownames_to_column() %>%
        group_by(.dots = groupby) %>%
        sample_frac(size = opt$sample_ratio, replace = F) %>%
        column_to_rownames()
      subseted_seurat <- SubsetData(data_ob, cells = rownames(sampled_cellmeta))
    }
    subseted_seurat@meta.data[, groupby] <- as.factor(subseted_seurat@meta.data[, groupby])
    colors2use <- SelectColors(unique(subseted_seurat@meta.data[, groupby]), palette = 'customecol2', value = groupby)
    if (length(markers2vis4heatmap) > 135) {
      sz <- 4 - log(length(markers2vis4heatmap) / 100)
      heig <- 5 + log2(length(markers2vis4heatmap) / 10)
      wid <- 7.5
    }else if (length(markers2vis4heatmap) < 75) {
      sz <- 6 - log2(length(markers2vis4heatmap) / 80); heig <- 7; wid <- 7
    }else {
      sz <- 4 - log2(length(markers2vis4heatmap) / 120); heig <- 7; wid <- 7
    }

    subseted_seurat <- subseted_seurat %>% Seurat::ScaleData(assay = "SCT", features = markers2vis4heatmap)
    ggheat <- Seurat::DoHeatmap(object = subseted_seurat,
                                features = markers2vis4heatmap,
                                group.colors = colors2use,
                                group.by = groupby, group.bar = T, label = F) +
      theme(axis.text.y = element_text(size = sz, face = "bold"))
    # group.cex = 10, cex.row = 4,
    # slim.col.label = T, group.label.rot = F)
    p_heatmap <- ggheat + ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", order = 1),
                                          color = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1)))

    OESingleCell::save_ggplots(file.path(output_dir, "topmarker_gene_heatmap"),
                               plot = p_heatmap,
                               width = heig,
                               height = wid)
      # ggsave(file.path(root_dir,paste0("top","marker_gene_heatmap.pdf", collapse = "_")),height = heig, width = wid)
      # ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap.png", collapse = "_")),height = heig, width = wid, dpi = 1000 ,limitsize = F)
  }

  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}