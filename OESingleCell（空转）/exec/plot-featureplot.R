## usage help:
docstring <- "example1(Select the topn marker gene from findAllmarker results and Visiual):\\n\\n\\
scVis -i seurat.rds  -f rds --dataslots data --reduct umap -o results --assay SCT featureplot  -l  top10_markers_for_each_cluster.xls -g cluster  -y sampleid  \\n\\n\\
scVis -i seurat.rds  -f rds --dataslots data --reduct umap -o results --assay SCT featureplot  -x  Ly6c-hi.gene.list  -y sampleid  \\n\\n\\
scVis -i seurat.rds  -f rds --dataslots data --reduct umap -o results --assay SCT featureplot  -m  other_metadata_column.list -y sampleid \\n\\n\\
"
sub_features <- subparsers$add_parser("featureplot",
                                      description = docstring,
                                      formatter_class = 'argparse.RawTextHelpFormatter',
                                      #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                      argument_default = "True",
                                      help = "visualize the bunch of features in a featureplot.")
sub_features$add_argument("-l", "--markers", type = "character", default = NULL,
                          help = "the file of marker genes table to be visulized.")
sub_features$add_argument("-g", "--groupby", type = "character", default = "clusters",
                          help = "[OPTIONAL]The grouppinig variable in input markers  to visulize marker genes.[default: %(default)s] ")
sub_features$add_argument("-x", "--extraFeatures", type = "character", default = NULL,
                          help = "[OPTIONAL]The extra feature/gene list of interest to visualize specified by the user.")
sub_features$add_argument("--min.cutoff", type = "character", default = "q10",
                          help = "min.cutoff value in plot show ")
sub_features$add_argument("--max.cutoff", type = "character", default = "q90",
                          help = "max.cutoff value in plot show  ")

sub_features$add_argument("-y", "--splitby", type = "character", default = NULL,
                          help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
sub_features$add_argument("-c", "--colorscale", type = "character", default = "spectral",
                          help = "[OPTIONAL]the color scale used to map the gene expression value,format:low_color,high_color.[default: %(default)s]")
sub_features$add_argument("-s", "--pointsize", type = "double", default = 0.02,
                          help = "[OPTIONAL]the point size in the plot.")
sub_features$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)

# visulaize the feature of cell on the dimension reduction plot, the feature can be anything pulled by FetechData
args <- commandArgs(TRUE)
if ("featureplot" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "featureplot") {
    # step3 read the specified assay and data slot in data object into memory===========================================
    # for visualization usually only "data" slot is loaded,which can save a lot of time and memory
    suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                    informat = opt$informat,
                                                    assays = assays,
                                                    data.use =  dataslots, # "data","scale.data"
                                                    reductions = opt$reduct,
                                                     verbose = FALSE))
    SeuratObject::DefaultAssay(data_ob) <- assays
    #get subset of cells used col.name = col.name = or visualization if necessay
    if (!is.null(opt$predicate)) {
      df <- OESingleCell::colData(data_ob)
      desired_cells <- subset(df, eval(parse(text = opt$predicate)))
      data_ob <- subset(data_ob,cells=rownames(desired_cells))
      if (!is.null(Seurat::Images(data_ob))) {
        unuse_images<-Seurat::Images(data_ob)[ ! Seurat::Images(data_ob)  %in% (data_ob@meta.data$sampleid%>%unique)]
        if(length(unuse_images)>0){ data_ob@images[unuse_images]<-NULL}
        }
    }
    # get genelist
    if (is.null(opt$markers) & is.null(opt$extraGene) & is.null(opt$extraFeature)) {
        stop("NO marker genes/features is AVAILABLE!")
    }

    if (!is.null(opt$markers)) {
      futile.logger::flog.info("step1 prepare markers to be visualized") #================================================
      topn_markers <- read.table(opt$markers, sep = "\t", header = T) %>%
                      dplyr::arrange(!!!syms(c(opt$groupby))) %>%
                      dplyr::mutate(folder_suffix = paste0(opt$groupby, !!!syms(c(opt$groupby)))) %>%
                      dplyr::select(!!!syms(c(opt$groupby)), gene, folder_suffix)
    }
    if (!is.null(opt$extraFeature)) {
      futile.logger::flog.info("step2 add extra customized genes if existis") #==========================================
      extra_gene <- read.table(opt$extraFeature, sep = "\t", header = T)
      if(colnames(extra_gene)=="gene"){
        match <- OESingleCell::CaseMatch(search = as.vector(extra_gene$gene), match = rownames(data_ob))

        ## if have unmatch genes,save
        unmatch <- extra_gene$gene [ ! extra_gene$gene %in% names(match)]
        if (length(unmatch) > 0) {
          futile.logger::flog.warn("There are some gene name unmatch gene in query object,please check:", unmatch, capture = TRUE)
          unmatch %>% as.data.frame %>% readr::write_tsv("unuse_gene.list", col_names = FALSE)
          #quit()
        }
      }else{
      colnames(extra_gene) <- "gene"
      }
      topn_markers <- extra_gene

    }
    cellmeta <- OESingleCell::colData(data_ob)
    ccolors <- unlist(strsplit(opt$colorscale, ","))
    if (length(ccolors) == 1) { # this means to use customized discrete color schema in OESingleCell package
      # info <- glue::glue("the customized continious colors are: \n {paste(names(OESingleCell::continuous_palette),
      # lengths(OESingleCell::continuous_palette), sep =': ', collapse = '\n')}.")
      # warning("NO specified color schema found!")
      # warning(info)
      # stop()
      continious_colors <- OESingleCell::SelectColors(palette = ccolors, is.discrete = FALSE)
    }else { # this means to use the customized color schema from command line specified by user
      continious_colors <- ccolors
    }

    ##step4  ploting featureplot  using Seurat::Featureplot functions ==================================================
    #lapply(1:dim(topn_markers)[1], function(i) {
    parallel::mclapply(1:dim(topn_markers)[1], function(i) {
      pdf(NULL)
      x <- topn_markers[i, "gene"]
      if (!is.null(opt$markers)) {
        prefix = topn_markers[which(topn_markers$gene == x), "folder_suffix"]
        output_prefix = glue::glue("{output_dir}/markers_vis4{prefix}") }
      else { output_prefix = output_dir }
      if (!is.null(opt$extraFeature)) {
        data <- Seurat::FetchData(object = data_ob, vars = x)
      }else {
        data <- Seurat::FetchData(object = data_ob, vars = x, slot = dataslots)
      }

      ##for featureplot 目前版本seurat中FeaturePlot split.by使用时存在bug,无法统一分图标尺,需要人为修改
      min.feature.value <- Seurat::SetQuantile(opt$min.cutoff, data[, x])
      futile.logger::flog.info(glue::glue("min.feature.value:{min.feature.value}"))#==========================================
      max.feature.value <- Seurat::SetQuantile(opt$max.cutoff, data[, x])
      futile.logger::flog.info(glue::glue("max.feature.value:{max.feature.value}"))
      ###merge
      ggf1 <- Seurat::FeaturePlot(data_ob,
                                  slot = dataslots,
                                  reduction = opt$reduct,
                                  features = x,
                                  keep.scale = "all",
                                  min.cutoff = opt$min.cutoff,
                                  max.cutoff = opt$max.cutoff,
                                  #cols =continious_colors,
                                  #blend=TRUE,
                                  combine = TRUE,
                                  order = T,
                                  pt.size = opt$pointsize) &
        ggplot2::scale_color_gradientn(colours = continious_colors,
                                       limits = c(min.feature.value, max.feature.value))
      ggf1 <- ggf1 + theme(legend.position = "right")

      ##split
      ggf2 <- Seurat::FeaturePlot(data_ob,
                                  slot =  dataslots,
                                  reduction = opt$reduct,
                                  features = x,
                                  split.by = opt$splitby,
                                  keep.scale = "all",
                                  min.cutoff = opt$min.cutoff,
                                  max.cutoff = opt$max.cutoff,
                                  #cols =continious_colors,
                                  #blend=TRUE,
                                  combine = TRUE,
                                  order = T,
                                  pt.size = opt$pointsize) &
        ggplot2::scale_color_gradientn(colours = continious_colors,
                                       limits = c(min.feature.value, max.feature.value))
      ggf2 <- ggf2 + theme(legend.position = "right")
      ##for spatial featureplot

      if (!is.null(Seurat::Images(data_ob))) {

        crop <- as.logical(opt$crop)
        gsp <- OESingleCell::SpatialPlot(data_ob,
                                         slot =  dataslots,
                                         assay = assays,
                                         HE=opt$HE,
                                         features = x,
                                         combine = TRUE,
                                         alpha = c(0.5,1),
                                         min.cutoff = opt$min.cutoff,
                                         max.cutoff = opt$max.cutoff,
                                         cols = opt$colorscale,
                                         ncol = length(Seurat::Images(data_ob)),
                                         pt.size.factor = 0.8,
                                         crop = crop
        )
      }
      ##输出构图

      splitby_len <- length(unique(cellmeta[[opt$splitby]]))
      image_len <- ifelse( !is.null(Seurat::Images(data_ob)) ,length(Seurat::Images(data_ob)),0)
      max_colnum <- max(splitby_len, image_len)
      layout <- ifelse( !is.null(Seurat::Images(data_ob)),
                       glue::glue("A{strrep('#',max_colnum-1)}\n{strrep('B',splitby_len)}{strrep('#',max_colnum-splitby_len)}\n{strrep('C',image_len)}{strrep('#',max_colnum-image_len)}"),
                       glue::glue("A{strrep('#',max_colnum-1)}\n{strrep('B',splitby_len)}{strrep('#',max_colnum-splitby_len)}"))

       if ( !is.null(Seurat::Images(data_ob))) {
         ggf_merge <-patchwork::wrap_plots(A = ggf1, B = ggf2, C = gsp, heights = c(1, 1, 1.2), design = layout)
       }else{
         ggf_merge <-patchwork::wrap_plots(A = ggf1, B = ggf2, heights = c(1, 1), design = layout)
       }
      OESingleCell::save_ggplots(plot = ggf_merge,
                                 filename = glue::glue("{output_prefix}/{x}_featureplot"),
                                 width = 3* max_colnum + 1,
                                 height = 3 * ifelse((!is.null(Seurat::Images(data_ob))), 3, 2),
                                 dpi = 100,
                                 limitsize = F)
      dev.off()
    },mc.cores =10)
    #})
    ## output session information
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  } }
