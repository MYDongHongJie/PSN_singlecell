# create the parsers for subcommand summarize
sub_summ <- subparsers$add_parser("summarize", help = "summarized statistics for clustering results")
sub_summ$add_argument( "-c", "--groupby", type = "character", default = "clusters",
             help = "[OPTIONAL]visualize cell metadata by coloring cells in different color according to cell clustering metadata.[default: %(default)s].")
sub_summ$add_argument( "-b", "--facetby", type = "character", default = NULL,
             help = "[OPTIONAL]visualize cells in seperate plot split by this groupping variable, only 2D supported.[default: %(default)s] ")
sub_summ$add_argument( "-s", "--ptsize", type = "double", default = 1.0,
             help = "[OPTIONAL]the point size in the plot.[default: %(default)s]")
sub_summ$add_argument( "-k", "--dims", type = "integer", default = 2,
             help = "[OPTIONAL]the 2D/3D of the plot.Currently 3D plot is only support in plotly style.[default: %(default)s]")
sub_summ$add_argument( "--reduct", type = "character", default = "umap",
             help = "[OPTIONAL]the previous calculated reduction result used in the scatter plot.[default: %(default)s]")
sub_summ$add_argument("--palette",type = "character",default= "blindless",
             help="the discrete color schema for each cell groups, blindless as default. Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20, paired:12, purplegray12:12.[default: %(default)s] ")
sub_summ$add_argument( "--dosummary", type = "character", default = "TRUE",
             help = "wether to plot the summary statistics of cells clustering.[default: %(default)s] ")

# =============== Subcmd: summarize, summarizing and visualizing the clustering analysis results =========
args <- commandArgs(TRUE)
if ( "summarize"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name ==  "summarize" ){
    data_ob <- OESingleCell::ReadX(input = opt$input,
                                   informat = opt$informat,
                                   assays = assays,
                                   reduction = opt$reduct,
                                   data.use = dataslots,  # data slot is enough
                                   verbose = F)
    groupfactors <- unlist(strsplit(opt$groupby, ",", perl = T))
    color_schema <- OESingleCell::SelectColors(palette = opt$palette)

    #===================================================================================================================
    ##
    #===================================================================================================================
    for (groupfactor in groupfactors ){
      output_dir <- glue::glue("{output_dir}/visualize_cluster_by_{groupfactor}")
      if ( !file.exists(output_dir) ){dir.create(output_dir)}
      color_counter <- 0
      nlevel <- length(unique(OESingleCell::colData(data_ob)[[groupfactor]]))
      if ( !is.null(opt$facetby) ){
        for ( facetbyx in unlist( strsplit( opt$facetby, "," ) )){
           if ( !is.null(facetbyx) ){
               nrow <- ceiling(length(unique(OESingleCell::colData(data_ob)[[facetbyx]]))/2)
           }
           #nlevel = length(unique(OESingleCell::colData(data_ob)[[groupfactor]]))
           color_counter <- nlevel
           data_ob <- Seurat::SetIdent(data_ob, value = groupfactor)
           groupvis_split <- Seurat::DimPlot(object = data_ob,
                                             dims = c(1,2),
                                             reduction = opt$reduct,
                                             pt.size = as.numeric(opt$ptsize),
                                             ncol = 2,
                                             group.by = groupfactor,
                                             split.by = facetbyx)+
               ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5))
           groupvis_split <- groupvis_split +
                             ggplot2::scale_colour_manual(values = color_schema[1:nlevel],
                                                          labels = sort(unique(Seurat::Idents(data_ob))))

           OESingleCell::save_ggplots(glue::glue("{output_dir}/splitby-{facetbyx}_split_plot"),
                                      limitsize = FALSE,
                                      plot = groupvis_split,
                                      width = 14,
                                      height = 5*nrow )

           nfacet <- length(unique(OESingleCell::colData(data_ob)[[facetbyx]]))
           groupvis_all <- Seurat::DimPlot(object = data_ob,
                                           dims = c(1, 2),
                                           reduction = opt$reduct,
                                           pt.size = as.numeric(opt$ptsize),
                                           group.by = facetbyx)+
                           ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
                           ggplot2::scale_colour_manual( values = color_schema[14:(14+nfacet)])
           OESingleCell::save_ggplots(glue::glue("{output_dir}/groupby-{facetbyx}_contrast_plot"),
                                      limitsize = FALSE,
                                      plot = groupvis_all,
                                      width = 7,
                                      height = 7 )
           color_counter <- color_counter + nfacet
       }
      }else{
        if ( as.numeric(opt$dims) == 2 ){
            groupvis <- Seurat::DimPlot(object = data_ob,
                                        dims = c(1, 2),
                                        reduction = opt$reduct,
                                        pt.size = as.numeric(opt$ptsize),
                                        group.by = groupfactor)+
                ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
                ggplot2::scale_colour_manual( values = color_schema[nlevel:nlevel*2])
            OESingleCell::save_ggplots(file.path(output_dir,paste0("groupby-",groupfactor, "_contrast_plot",collapse="")),
                                       plot = groupvis,
                                       limitsize = FALSE, width = 15, height = 10 )
        }else{
            suppressPackageStartupMessages(library("SeuratPlotly"))
            groupvis <- SeuratPlotly::DimPlotly3d(data_ob, reduction= opt$reduct,
                                                  pt_size = as.numeric(opt$ptsize),
                                                  grouping_var = groupfactor,
                                                  plot_grid = F)
            htmlwidgets::saveWidget(groupvis, file = file.path(output_dir,paste0("groupby-",groupfactor,"_contrast_3D_plot.html",collapse="")))
        }
      }
    }
    #===================================================================================================================
    ##
    #===================================================================================================================
    if ( as.logical(opt$dosummary) ){
          data_ob <- Seurat::SetIdent(data_ob, value = opt$groupby )
          DATA <- suppressMessages(OESingleCell::colData(data_ob)[,c("sampleid", opt$groupby)] %>%
                                   dplyr::group_by( .dots= c("sampleid", opt$groupby)) %>%
                                   dplyr::summarize(cell_number = dplyr::n()) %>%
                                   dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100))
          write.table(as.data.frame(DATA),
                      file.path(output_dir,file="clust_cond_freq_info.xls"),
                      sep="\t",
                      col.names=T,
                      row.names =F)
          # visulaize the summary statistics of cells clusters in each groupping factor
          if ( !is.null(opt$facetby) ){
              for (facetbyx in unlist( strsplit( opt$facetby, "," ) )){
                  clust_sum_all <- OESingleCell::PlotAbundances(data_ob,
                                                                prop.by = opt$groupby,
                                                                group.by = facetbyx,
                                                                method = "barplot",
                                                                cols= unname(OESingleCell::SelectColors(1:length(unique(Seurat::Idents(data_ob))),
                                                                                                       palette = opt$palette)))
                  OESingleCell::save_ggplots(file.path(output_dir,paste0("groupby-",facetbyx,"_summary_plot",collapse="")),
                                             dpi=1000,
                                             plot = clust_sum_all,
                                             width = 7,
                                             height = 7 )
               }
          }
          data_ob <- Seurat::SetIdent(data_ob, value = "sampleid")
          clust_sum_all2 <- OESingleCell::PlotAbundances(data_ob,
                                                         prop.by = "sampleid",
                                                         group.by = opt$groupby,
                                                         method = "barplot",
                                                         cols= unname(OESingleCell::SelectColors(1:length(unique(Seurat::Idents(data_ob))),
                                                                                                palette = opt$palette)))
          OESingleCell::save_ggplots(file.path(output_dir,paste0("groupby-",opt$groupby, "_summary_plot",collapse="")),
                                     dpi= 1000,
                                     plot = clust_sum_all2,
                                     width = 7,
                                     height = 7 )
      }
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}