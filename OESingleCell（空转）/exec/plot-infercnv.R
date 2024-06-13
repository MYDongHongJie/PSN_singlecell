sub_infercnv <- subparsers$add_parser("infercnv", help = "visualize the result of infercnv.")
sub_infercnv$add_argument( "--infercnv_output_path", "-l", type = "character",
             help = "[OPTIONAL]the infercnv result directory." )
sub_infercnv$add_argument("--groupby", "-g", type="character",
        help = "the number of top ranked markers for each cluster to visualizse.")
sub_infercnv$add_argument("--gorder", "-e", type = "character",
        help=paste0("[OPTIONAL]comma seperated list of levels for the specified parameter",
                    "'--groupby' to order the cells.The format can be gvariable1:l1,l2,l3;gvariable2:l1,l2,l3." ))
sub_infercnv$add_argument("--invert", "-v", type = "logical", default = FALSE,
        help = "reverse the selection conditions specified by all other conditions.[default: FALSE]")
sub_infercnv$add_argument("--colormapping", "-m", type = "character",
        help = paste0("The color mapping for groupping column of cells set by the parameters '--groupby'.",
                      "The exmaple format is variable1:colorschema1,variable2:colorschema2. ",
                      "The supported color schemas can be:blindless, col50, ditto, paired."))
sub_infercnv$add_argument("--vlnlegend", type = "logical", default =  FALSE,
        help = "whether to show the legend of the vlnplot.[default: FALSE]")
sub_infercnv$add_argument("--splitby","-y", type = "character",
        help = "the column used to split the vlnplot/featureplot")
sub_infercnv$add_argument("--vismethod",type = "character", default = "all",
        help = "plots to visualize. choose from all, heatmap, vlnplot, featureplot.[default: %(default)s]")

args <- commandArgs(TRUE)

if ( "infercnv"  %in% args ){
    opt<-intial_setting()
    if (opt$sub_name == "infercnv") {
      print(opt)
      #===================================================================================================================
      futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory")
      suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                      informat = opt$informat,
                                                      assays = assays,
                                                      data.use = dataslots,
                                                      verbose = F))
      ##subset
      if (!is.null(opt$predicate)) {
        futile.logger::flog.info(glue::glue("get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
        df <- OESingleCell::colData(data_ob)
        desired_cells <- subset(df, eval(parse(text = opt$predicate)))
        data_ob <- subset(data_ob,cells=rownames(desired_cells),invert = opt$invert)
        if (!is.null(Seurat::Images(data_ob))) {
          unuse_images<-Seurat::Images(data_ob)[ ! Seurat::Images(data_ob)  %in% (data_ob@meta.data$sampleid%>%unique)]
          if(length(unuse_images)>0){ data_ob@images[unuse_images]<-NULL}
          }
      }
      # set group order=================================================================================================
      if ( !is.null(opt$gorder)){
          group.order <- list()
          for ( x in  unlist(strsplit(opt$gorder, ";", perl =T))){
              m <- unlist(strsplit(x, ":|,", perl =T))
              group.order[[m[1]]] <- m[-1]
          }
      }else{
          group.order <- NULL
      }
      # set the groupby level===========================================================================================
      if ( !is.null( opt$groupby) ){
          cell.annos <- unlist(strsplit(opt$groupby, ",", perl =T))
          groupby <- unlist(strsplit(opt$groupby, ",", perl =T))[1]
      }else{
          stop("please provide the groupping of cells in the metadata of seurat object.")
      }

      # set the color schema of cell annotation bar=====================================================================

      if ( !is.null( opt$colormapping) ){
          group_colors <- list()
          for ( x in  unlist(strsplit(opt$colormapping, ",", perl =T))){
              m <- unlist(strsplit(x, ":", perl =T))
              group_colors[[m[1]]] <- m[-1]
          }
      }else{
          group_colors <- 1:length(cell.annos)
          names(group_colors) <- cell.annos
          group_colors <- as.list(group_colors)
      }

     if(opt$vismethod == "all"){
        vismethods <-   c("heatmap", "vlnplot", "featureplot")
     }else{
        vismethods <- unlist(strsplit(opt$vismethod, ","))
     }

      #readin the results/rds============================================================================================

      if ( !is.null(opt$output) ){
          data_ob <- OESingleCell::add_to_seurat(data_ob,
                                                 infercnv_output_path = output_dir,
                                                 nclones = opt$ncores)
          saveRDS( data_ob, file.path(output_dir, "cnv_seurat.rds"))
      }else{
          if ( !"CNV" %in% Seurat::Assays(data_ob) ){
             stop( "NO previous CNV assay results integration FOUND!")
          }
      }
      cnv_result <- data_ob@meta.data %>%
                    tibble::rownames_to_column(var = "barcode_inuse") %>%
                    dplyr::rename( "Barcode" = "orig.ident") %>%
                    dplyr::select(barcode_inuse,
                                  Barcode,
                                  sampleid,
                                  groupby ,
                                  clusters,
                                  group,
                                  cnv_level,
                                  cnv_group)
      write.table(cnv_result, file.path( output_dir, "cnv_result.xls"), quote = F,sep ="\t",row.names =F)

      if ( !is.null(opt$splitby)){
          if (opt$splitby %in% colnames(data_ob@meta.data)) {
              splitby <- opt$splitby
          } else {
              stop("can not find -splitby in the metadata")
          }
      }else{
          splitby <- NULL
      }
      # visulization====================================================================================================
      for ( vismethod in vismethods ){
          ## heatmap====================================================================================================
          if ( vismethod == "heatmap" ){
              gene_order_f <- data_ob[["CNV"]]@meta.features
              plot<- OESingleCell::CnvHeatmap(data_ob,
                                              group.by = groupby,
                                              group.order = group.order,
                                              features = rownames(gene_order_f),
                                              cell.annotation = cell.annos,
                                              gene.annotation = gene_order_f[,1,drop=F],
                                              row.split.by = NULL,
                                              col.split.by = "chr",
                                              group.colors = group_colors)
              OESingleCell::save_ggplots(plot = ggplotify::as.ggplot(plot),
                                         filename = file.path( output_dir, "cnv_heatmap"),
                                         width = 14,
                                         height = 8)
          }
          ## dittoPlot==================================================================================================
          if ( vismethod == "vlnplot" ){
              for (group_by in cell.annos){
                  nlevel <- length(unique(data_ob@meta.data[, group_by] ))
                  nlevel_splitby <- ifelse(is.null(splitby), 1, length(unique(data_ob@meta.data[, splitby] )))
                  vlnplot <- dittoSeq::dittoPlot(data_ob,
                                                 var= "cnv_level",
                                                 group.by=group_by,
                                                 plots =c("vlnplot", "boxplot"),
                                                 legend.show = opt$vlnlegend ,
                                                 split.by = splitby) +
                             ggplot2::scale_fill_manual(values = OESingleCell::SelectColors(object=NULL,
                                                                                            palette = group_colors[[group_by]],
                                                                                            value = group_by ,
                                                                                            n = nlevel))+
                          ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5))
                  OESingleCell::save_ggplots(filename=file.path(output_dir,
                                                                paste0("cnv_vlnplot_groupby_",
                                                                       group_by,
                                                                       ifelse(is.null(splitby),"",paste0("_splitby_",splitby)),
                                                                       collapse = "") ),
                                              plot = vlnplot ,
                                              height=6,
                                              width = nlevel*nlevel_splitby^0.4 + 2 )
              }
          }
          ## featureplot================================================================================================
          if ( vismethod == "featureplot" ){
              pointsize <- 0.8
              ggfeature <- Seurat::FeaturePlot(data_ob,
                                               features = "cnv_level",
                                               split.by = splitby,
                                               reduction= opt$reduct,
                                               cols = rev(OESingleCell::SelectColors(palette = "spectral",
                                                                                     is.discrete = FALSE)),
                                               pt.size = pointsize,
                                               order= T ) +
                           ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5))
              OESingleCell::save_ggplots(filename=file.path( output_dir,
                                                             paste0("cnv_level_on_",opt$reduct,
                                                                    ifelse(is.null(opt$splitby),
                                                                               "",
                                                                           paste0("_splitby_",splitby)))),
                                         plot=ggfeature,
                                         height=6,
                                         width = ifelse(is.null(opt$splitby),
                                                        7,
                                                        6*length(unique(data_ob@meta.data[,splitby] ))))
          }
      }
      ## output session information
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
    }
}
