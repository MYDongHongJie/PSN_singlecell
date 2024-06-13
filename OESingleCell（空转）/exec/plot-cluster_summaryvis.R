sub_clustervis = subparsers$add_parser("clustervis", help = "integration of multiple visualizion for cluster.")
sub_clustervis$add_argument("-l","--level", type = "character", default = "cell",
                            help = "the analysis object level: cell or spot.[default: %(default)s]")
sub_clustervis$add_argument("-s","--ptsize",type="double", default=0.5,
                            help='the point size in the dimplot.[default: %(default)s]')
sub_clustervis$add_argument("-r","--resolution",type="double", default=0.4,
                            help="resolution used in cell clustering.[default: %(default)s]")
sub_clustervis$add_argument('-t','--idents',type='character',default=NULL,
                            help='The column name in cell metadata used as identity of each cell combined with --cell2use.[default: %(default)s]')
sub_clustervis$add_argument("-b","--splitby",type="character",default="NULL",
                            help='visualize cells in seperate plot split by this groupping variable, only 2D supported.[default: %(default)s]')
sub_clustervis$add_argument("--from",type="character",default="NULL",
                            help='the original cluster id used in calibrate the clustering results.[default: %(default)s]')
sub_clustervis$add_argument("--to",type="character",default="NULL",
                            help='the adjusted cluster id used in calibrate the clustering results.[default: %(default)s]')
sub_clustervis$add_argument("-c","--groupby",type="character",default="clusters",
                            help='visualize cell metadata by coloring cells in different color according to cell clustering metadata.[default: %(default)s]')
sub_clustervis$add_argument("--sptsize",type="double",default="1.5",
                            help='the dimension reduction result to use in this run.[default: %(default)s]')
sub_clustervis$add_argument("--dims",type="integer",default="2",
                            help='the 2D/3D of the plot.Currently 3D plot is only support in plotly style.[default: %(default)s]')
sub_clustervis$add_argument("--propby",type="character",default="group",
                            help='calculate the proportion of cell in the specified groups.[default: %(default)s]')
sub_clustervis$add_argument("--plot",type="character",default="barplot",
                            help='the plot type used to summarize the cell clustering results.Options can be:barplot, pie, boxplot.[default: %(default)s]')
sub_clustervis$add_argument('--groups',type="character",default="clusters",
                            help="visualize cell metadata by coloring cells in different color according to cell clustering metadata.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "clustervis" %in% args ){
  opt<-intial_setting()
  ptsize <- as.numeric(opt$ptsize)
  reduct <- opt$reduct


  assays <- opt$assay
  seurat_ob <- OESingleCell::ReadX(input = opt$input,
                               informat = opt$informat,
                               assays = assays,
                               graphs = FALSE,
                               images = opt$image, # no graph object needed here
                               verbose = FALSE)
  Seurat::DefaultAssay(seurat_ob) <- assays

  if (is.null(seurat_ob@meta.data$clusters)) {
      seurat_ob <- Seurat::StashIdent(seurat_ob, save.name = "clusters")
    }else {
      seurat_ob <- Seurat::SetIdent(seurat_ob, value = "clusters")
    }
  if(opt$image==TRUE){
    if("slice" %in% names(seurat_ob@meta.data)==FALSE){
      seurat_ob$slice <- seurat_ob$sampleid
    }
  }
    # extract the preserved resolution setting
  if (is.null(opt$resolution)) {
    findcluster_record <- Command(seurat_ob, command = "FindClusters")
    resolution <- findcluster_record$resolution
  }else {
    resolution <- opt$resolution
  }

  if (is.null(opt$idents)) {
    print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
    ident2use <- "clusters"
  }else {
    ident2use <- opt$idents
  }

  if (!is.null(opt$splitby) & opt$splitby != "NULL") {
    facetby <- unlist(strsplit(opt$splitby, ",", perl = T))
  }else {
    facetby <- NULL
  }

    seurat_ob <- Seurat::StashIdent(seurat_ob, save.name = "clusters")
    cell_count_tsne <- table(Seurat::Idents(seurat_ob))
    tsne_cell_count_labels <- paste(paste(names(cell_count_tsne), cell_count_tsne, sep = "-"), " cells")

    if (!is.null(opt$from) & !is.null(opt$to)) {
      from_ident <- unlist(strsplit(opt$from, ",", perl = T))
      to_ident <- unlist(strsplit(opt$to, ",", perl = T))
      Seurat::Idents(seurat_ob) <- plyr::mapvalues(x = Seurat::Idents(seurat_ob),
                                           from = from_ident, to = to_ident)
    }

  if ( !is.null(opt$predicate) ){
    predicate <- opt$predicate
    df <- OESingleCell::colData(seurat_ob)
    desired_cells <- subset(df, eval(parse(text=predicate)))
    seurat_ob <- seurat_ob[, desired_cells]
  }


  cellmeta <- OESingleCell::colData(seurat_ob)

      if (is.null(opt$groupby)) {
        print("The groupping information is not specified. The sampleid will be used for plot")
        cellfeature2vis <- "clusters"
      }else {
        cellfeature2vis <- unlist(strsplit(opt$groupby, ",", perl = T))
      }
      for (groupfactor in cellfeature2vis) {
        output_dir <- file.path(output_dir, paste0("visualize_cluster_by_", groupfactor, collapse = ""))
        if (!file.exists(output_dir)) {
          dir.create(output_dir)
        }
        color_counter <- 0
        if (!is.null(facetby)) {
          for (facetbyx in facetby) {
            if (!is.null(facetbyx)) {
              nrow <- ceiling(length(unique(seurat_ob@meta.data[, facetbyx])) / 2)
            }
            nlevel <- length(unique(cellmeta[[groupfactor]]))
            color_counter <- nlevel
            groupvis_split <- Seurat::DimPlot(object = seurat_ob,
                                      dims = c(1, 2),
                                      reduction = reduct,
                                      pt.size = ptsize,
                                      # nrow = 1,#原本为ncol=1
                                      group.by = groupfactor,
                                      split.by = facetbyx) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
            groupvis_split <- groupvis_split + ggplot2::scale_colour_manual(values = OESingleCell::SelectColors(1:nlevel,palette='customecol2'),
                                                                   labels = sort(as.integer(unique(Seurat::Idents(seurat_ob)))))
            if (opt$image==TRUE) {

              gsp <- OESingleCell::SpatialPlot(seurat_ob,HE=F,
                            group.by = groupfactor,
                            ncol = 2,
                            alpha = c(0.1, 1),
                            pt.size.factor = as.numeric(opt$sptsize),
                            pt.alpha=FALSE,
                            stroke=0.2,combine=FALSE,
                            col= OESingleCell::SelectColors(1:nlevel,palette='customecol2')
                            )
              gsb <- do.call(ggpubr::ggarrange,
                     c(gsp, list(ncol = ifelse(length(gsp)<4,length(gsp),4),
                                 nrow = ceiling(length(gsp)/4),
                                 common.legend = TRUE,
                                 legend = "none",
                                 align = "none")))
              plot <- cowplot::plot_grid(groupvis_split, gsb,ncol = 1,nrow=2,align = "v")#rel_widths= c(1, length(gsp)))
              width <-  3 * length(gsp)#3 * length(gsp) + length(groupvis_split)
              height <- ceiling(length(gsp)/4)*8#ceiling(length(gsp)/4)*4
            }
            else {
              plot <- groupvis_split
              width<-8
              height<-6
            }
            OESingleCell::save_ggplots(file.path(output_dir,
                             paste0("splitby-",
                                    facetbyx,
                                    "_resolution",
                                    resolution,
                                    "_split_plot",
                                    collapse = "")),
                   plot = plot,
                   width = width,
                   height = height,
                   dpi=300,bg='white')


            nfacet <- length(unique(seurat_ob@meta.data[, facetbyx]))

            groupvis_all <- Seurat::DimPlot(object = seurat_ob,
                                    dims = c(1, 2),
                                    reduction = reduct,
                                    pt.size = ptsize,
                                    group.by = facetbyx) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
                          ggplot2::scale_colour_manual(values = unname(OESingleCell::SelectColors((nlevel + 1):(nlevel + nfacet),palette='customecol2')))##原来为+2
            OESingleCell::save_ggplots(file.path(output_dir, paste0("groupby-",
                                                                    facetbyx,
                                                                    "_resolution",
                                                                    resolution,
                                                                    "_contrast_plot",
                                                                    collapse = "")),
                                       plot = groupvis_all,
                                       width = 10,
                                       height = 12,
                                       dpi=300,bg='white')
            color_counter <- color_counter + nfacet

          }
        }else {
          if (as.numeric(opt$dims) == 2) {
            nlevel_col <- unique(seurat_ob@meta.data[,groupfactor])
            groupvis <- Seurat::DimPlot(object = seurat_ob, dims = c(1, 2), reduction = reduct,
                                pt.size = ptsize, group.by = groupfactor) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_colour_manual(values =OESingleCell::SelectColors(nlevel_col,palette='customecol2'))

                      OESingleCell::save_ggplots(file.path(output_dir, paste0("groupby-",
                                                                    groupfactor,
                                                                    "_resolution",
                                                                    resolution,
                                                                    "_contrast_plot",
                                                                    collapse = "")),
                                       plot = groupvis,
                                       width = 10,
                                       height = 12,
                                       dpi=300,bg='white')
          }else {
            suppressPackageStartupMessages(library("SeuratPlotly"))
            groupvis <- Seurat::DimPlotly3d(seurat_ob, reduction = reduct,
                                    pt_size = ptsize,
                                    grouping_var = groupfactor,
                                    plot_grid = F)
            htmlwidgets::saveWidget(groupvis, file = file.path(output_dir, paste0("groupby-", groupfactor, "_resolution", resolution, "_contrast_3D_plot.html", collapse = "")))
          }
        }
      }



    # summarize the cell type composition in each sample and group
      all_cell <- dim(seurat_ob)[1]
      DATA <- as.data.frame(seurat_ob@meta.data[, c(opt$groups, opt$propby)]) %>%
        dplyr::group_by(.dots = c(opt$groups, opt$propby)) %>%
        dplyr::summarize(count = dplyr::n()) %>%
        dplyr::mutate(freq = (count / all_cell) * 100)
      write.table(as.data.frame(DATA),
                  file.path(output_dir, file = "clust_cond_freq_info.xls"),
                  sep = "\t", col.names = T, row.names = F)
      clust_sum_all <- OESingleCell::PlotAbundances(seurat_ob, prop.by = opt$propby, group.by = opt$groups,
                                      split.by = NULL, method = opt$plot, ncol = ifelse(opt$plot == "pie", 4, 1),
                                      cols = unname(OESingleCell::SelectColors(1:length(unique(seurat_ob@meta.data[, opt$propby])),palette='customecol2'))
      )
      OESingleCell::save_ggplots(file.path(output_dir,
                       paste0("groupby-", opt$groups, "_resolution-", resolution, "_summary_plot", collapse = "")),
                                       plot = clust_sum_all,
                                       dpi=200,bg='white')




}