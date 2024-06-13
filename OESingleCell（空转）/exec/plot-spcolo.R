docstring<- "example:\\n\\n\\
  scVis -i seurat.rds -f rds -o results --assay SCT spcolo  --misclist spotlight_results  --interest T_cells,B_cells   --colors yellow,blue \\n\\n\\
  scVis -i seurat.rds -f rds -o results --assay SCT spcolo  --interest gene1,gene2  --colors yellow,blue"
sub_spcolo<- subparsers$add_parser(
  "spcolo",
  description = docstring,
  formatter_class= 'argparse.RawTextHelpFormatter' ,
  #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Colocalization plot for two celltype or features ")
sub_spcolo$add_argument(
  "--misclist",
  type = "character",
  default = NULL,
  help = "[REQUIRED]the spotlight result's name in seurat's miscs,eg:spotlight_results,RCTD_results" )
sub_spcolo$add_argument(
  "--interest",
  type = "character",
  default = NULL,
  help = "[REQUIRED]the two interest features in the plot[default: %(default)s] " )
sub_spcolo$add_argument("--colors",
   type = "character",
   default = "yellow,blue",
   help = "the discrete color schema[default: %(default)s] ")
sub_spcolo$add_argument("--singlegene",#是否增加单个基因展示图
   type = "character",
   default = "FALSE",
   help = "if add the plot of singlegene?[default: %(default)s] ")
sub_spcolo$add_argument('--L_R_heatmap',
   type = "character",
   default = "FALSE",
   help = "do you need the heatmap ?[default: %(default)s] ")

args <- commandArgs(TRUE)
if ( "spcolo"  %in% args ){
  opt <- intial_setting()
  if (opt$sub_name == "spcolo") {
    #===================================================================================================================
    futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory")
    suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                    informat = opt$informat,
                                                    assays = assays,
                                                    data.use = dataslot,
                                                    verbose = F))
    ##remove unused sample image in data_ob@image
    if ('slice' %in% names(data_ob@meta.data) == FALSE) { data_ob$slice <- data_ob$sampleid }
    if (!is.null(opt$predicate)) {
      futile.logger::flog.info(glue::glue("get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
      df <- OESingleCell::colData(data_ob)
      desired_cells <- subset(df, eval(parse(text = opt$predicate)))
      data_ob <- subset(data_ob,cells=rownames(desired_cells))
      if (!is.null(Seurat::Images(data_ob))) {
        unuse_images<-Seurat::Images(data_ob)[ ! Seurat::Images(data_ob)  %in% (data_ob@meta.data$sampleid%>%unique)]
        if(length(unuse_images)>0){ data_ob@images[unuse_images]<-NULL}
        }
    }
    #===================================================================================================================
    futile.logger::flog.info(glue::glue("step2:add  feature informtaion into seurat object metadata"))

    interest<- unlist(stringr::str_split(opt$interest,","))
    opt$colors<- unlist(stringr::str_split(opt$colors,","))

    print(interest)
    if(!is.null(opt$misclist)){
      data_ob@meta.data <- data_ob@meta.data %>%
                           cbind(data_ob@misc[[opt$misclist]]  %>%
                                 tibble::column_to_rownames("barcodes"))
    }else{
      print("OK")
      data_ob@meta.data <- data_ob@meta.data %>%
                           cbind(Seurat::FetchData(data_ob,
                                                   vars=interest,
                                                   slot=dataslots))
    }
    #===================================================================================================================
    futile.logger::flog.info(glue::glue("step3: ploting "))
    images <- Seurat::Images(data_ob)
    #parallel::mclapply(images, function(slice) {
    for(slice in images ){
      img <- Seurat::GetImage(data_ob[[slice]], mode = "raw")
      img_grob <- grid::rasterGrob(img,
                                   interpolate = FALSE,
                                   width = grid::unit(1, "npc"),
                                   height = grid::unit(1, "npc"))

      metadata_ds <- data.frame(data_ob@meta.data)
      #colnames(metadata_ds) <- colnames(data_ob@meta.data)
      print(head(metadata_ds))
      metadata_ds <- metadata_ds %>%
                     tibble::rownames_to_column("barcode") %>%
                     dplyr::mutate(rsum = base::rowSums(.[, interest, drop = FALSE])) %>%
                     dplyr::filter(rsum != 0) %>%
                     dplyr::select("barcode") %>%
                     dplyr::left_join(metadata_ds %>%
                     tibble::rownames_to_column("barcode"), by = "barcode") %>%
                     tibble::column_to_rownames("barcode")

      spatial_coord <- data.frame(data_ob@images[[slice]]@coordinates) %>%
                       tibble::rownames_to_column("barcode") %>%
                       dplyr::mutate(imagerow_scaled = as.integer(imagerow) * data_ob@images[[slice]]@scale.factors$lowres,
                                     imagecol_scaled = as.integer(imagecol) * data_ob@images[[slice]]@scale.factors$lowres) %>%
                       dplyr::inner_join(metadata_ds %>%
                       tibble::rownames_to_column("barcode"), by = "barcode")
      if(!is.null(opt$misclist)){
        max1.cutoff <- max(spatial_coord[, interest[1]])
        spatial_coord[,interest[1]] <- spatial_coord[, interest[1]]/max1.cutoff
        max2.cutoff <- max(spatial_coord[, opt$interest[2]])
        spatial_coord[,interest[2]] <- spatial_coord[, interest[2]]/max2.cutoff
        min1.cutoff <- 0
        max1.cutoff <- 1
        min2.cutoff <- 0
        max2.cutoff <- 1
        max.cutoff <- 1
      }else{
        min1.cutoff <- 0
        min2.cutoff <- 0
        max1.cutoff <- max(spatial_coord[, interest[1]])
        max2.cutoff <- max(spatial_coord[, interest[2]])
        max.cutoff <- max(spatial_coord[, interest])
      }
      #是否加上底片
      if(opt$HE==TRUE){
      p_img <- ggplot2::annotation_custom(grob = img_grob,
                                                          xmin = 0,
                                                          xmax = ncol(img),
                                                          ymin = 0,
                                                          ymax = -nrow(img))

      }else{
        p_img <- NULL
      }

      if(opt$singlegene==TRUE){
        single_p <- list()
        for (i in 1:length(interest)) {
          singlecol <- colorRampPalette(c("white", opt$colors[i]))(10)
          single_p[[i]] <- OESingleCell::SpatialPlot(data_ob, cols = singlecol,
                                    slot = dataslots,
                                    assay = 'Spatial',
                                    HE = FALSE,
                                    features = interest[i],
                                    combine = FALSE,
                                    alpha = 1,
                                    min.cutoff = min1.cutoff,
                                    max.cutoff = max2.cutoff,
                                    ncol = 1, images = slice,
                                    pt.size.factor = 1.2)[[1]]+
                                    ggplot2::theme(panel.border = ggplot2::element_blank())+
                                    ggplot2::ggtitle(paste0(interest[i]," expression"))
        }
        single_p <- do.call(ggpubr::ggarrange,
                         c(single_p,
                           list(nrow = 1,
                                ncol =  2,
                                common.legend = FALSE,
                                legend = "right",
                                align = "hv")))
    }


      plot <- suppressMessages(ggplot2::ggplot() +
                               p_img +
                               ggplot2::geom_point(data = spatial_coord,
                                                   aes(x = imagecol_scaled,
                                                       y = imagerow_scaled,
                                                       size = get(interest[1]),
                                                       alpha = get(interest[1])) ,
                                                       color = opt$colors[1]) +
                               ggplot2::labs(size = c(interest[1]))+
                               ggplot2::scale_size(range=c(0 ,max1.cutoff/max.cutoff ) ,guide = guide_legend(order = 1) ) +
                               ggplot2::scale_alpha_continuous(range=c(0,max1.cutoff/max.cutoff))+
                               # coord_cartesian(expand = FALSE) +
                               ggnewscale::new_scale_color() + # Define scales before initiating a new one
                               ggnewscale::new_scale("size") +
                               ggnewscale::new_scale("fill")+
                               ggplot2::geom_point(data = spatial_coord,
                                                   aes(x = imagecol_scaled,
                                                       y = imagerow_scaled,
                                                       size = get(interest[2]),
                                                       alpha = get(interest[2])) ,
                                                    color = opt$colors[2]) +
                               ggplot2::labs(size = c(interest[2])) +
                               ggplot2::scale_size(range=c(0,max2.cutoff/max.cutoff ) ,guide = guide_legend(order = 2) ) +
                               ggplot2::scale_alpha_continuous(range=c(0,max2.cutoff/max.cutoff))+
                               ggplot2::scale_y_reverse() +
                               ggplot2::ylim(nrow(img),0) +
                               ggplot2::xlim(0, ncol(img)) +
                               #cowplot::theme_half_open(11, rel_small = 1) +
                               ggplot2::theme_void() +
                               ggplot2::coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
                               ggplot2::theme(legend.direction = "vertical", legend.box = "vertical",legend.key.size=unit(0.3,'inches'))+
                               guides(alpha = "none"))+
                               ggplot2::ggtitle(paste0(paste(interest,collapse = "_")," expression"))+
                               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      if(opt$singlegene==TRUE){
         plot <- cowplot::plot_grid(single_p, plot,rel_widths = c(10, 6), nrow = 2, ncol=1,align = "hv")
         width = 8 + sum(stringr::str_length(interest))*0.14
         height = 8
      }else{width = 4 + sum(stringr::str_length(interest))*0.14
         height = 4
      }

      OESingleCell::save_ggplots( glue::glue("{output_dir}/spcoloction_{interest[1]}_{interest[2]}_for_{slice}") ,
                                  plot = plot ,
                                  width = width  ,
                                  height = height,
                                  limitsize = FALSE,
                                  dpi = 300 )
    if(opt$L_R_heatmap==TRUE){
    # matrix <- Seurat::GetAssayData(data_ob,assay = "Spatial",slot="count")[interest,]%>%data.frame() %>% t()%>%data.frame()
    # matrix$diff_expression <- as.numeric(matrix[,1])-as.numeric(matrix[,2])
    # p_heatmap <- ggplot2::ggplot(matrix , ggplot2::aes(x = matrix[,1], y =matrix[,2])) +
    #   ggplot2::geom_raster(ggplot2::aes(fill = diff_expression), interpolate=TRUE) +
    #   ggplot2::scale_fill_gradient2(low=opt$colors[2], mid="#F4DFEB", high=opt$colors[1],
    #                    midpoint=0 , limits=range(matrix$diff_expression)) + ggplot2::theme_classic()+
    #   ggplot2::labs(x=paste0(interest[1]," expression"),y=paste0(interest[2],' expression'))+
    #   ggplot2::scale_x_continuous(expand = c(0, 0))+
    #   ggplot2::scale_y_continuous(expand = c(0, 0))
    # OESingleCell::save_ggplots( glue::glue("{output_dir}/spcoloction_{interest[1]}_{interest[2]}_for_{slice}_heatmap") ,
    #                       plot = p_heatmap ,
    #                       width = 6  ,
    #                       height = 6,
    #                       limitsize = FALSE,
    #                       dpi = 300 )
      p_heatmap <- Seurat::FeaturePlot(data_ob,feature=c(interest[1],interest[2]),blend=T,combine=F,cols=c(colors[1],colors[2]))[[4]]
      OESingleCell::save_ggplots( glue::glue("{output_dir}/spcoloction_{interest[1]}_{interest[2]}_for_{slice}_heatmap") ,
                          plot = p_heatmap ,
                          width = 6  ,
                          height = 6,
                          limitsize = FALSE,
                          dpi = 300 )



    }



    }
      #}, mc.cores = 4)
    ## save session informations
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}