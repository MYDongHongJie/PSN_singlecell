sub_secondcelltype = subparsers$add_parser("secondcelltype", help = "plot firstcelltype and secondcelltype")
sub_secondcelltype$add_argument("-r", "--ratio", type="double", default=0.5,
             help="the absolute ratio in firstcelltyoe.[default: clusters]")
sub_secondcelltype$add_argument("--spotsize", type="double", default=1,
             help = "the spotsize to plot .[default: %(default)s]")
sub_secondcelltype$add_argument("--resultname", type="character",
             help = "'spotlight_results' or 'RCTD_results' in misc.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "secondcelltype" %in% args ){
  opt<-intial_setting()
  library(Seurat)
  object_celltype <- readRDS(opt$input)
  if(opt$resultname %in% names(object_celltype@misc)){
  result <- opt$resultname
  celltype_data <- object_celltype@misc[[result]][,-1]
  rownames(celltype_data) <- colnames(object_celltype)
  }else{
      celltype_data <- object_celltype@meta.data[,object_celltype@misc$celltype]
  }

  f_s_type <- apply(celltype_data,1,function(x){rang=sort(x)
                                                 rang2=rang[dim(celltype_data)[2]] %>% names()
                                                 if(rang[dim(celltype_data)[2]]<0.5){rang3 <- rang[(dim(celltype_data)[2]-1)] %>% names()
                                                              }else{
                                                     rang3 <- rang2
                                                 }
                                               rang_all <- c(rang2,rang3)

                                                 })


  object_celltype$firstcelltype <- f_s_type[1,]
  object_celltype$secondcelltype <-  f_s_type[2,]
  colors <- OESingleCell::SelectColors(unique(object_celltype$secondcelltype))

  p1_1 <- list()
  for(i in Seurat::Images(object_celltype)){
      p_celltype <- OESingleCell::SpatialPlot(seurat_ob= object_celltype,group.by='firstcelltype',cols=colors,
                                           min.cutoff = 0.1,
                                           combine=FALSE,
                                           ncol= 2,
                                           alpha=c(0.1,0.5),
                                           pt.size.factor=opt$spotsize,
                                           pt.alpha =1,
                                           images=i)
      p1_1[[i]] <- p_celltype[[1]]
  }


  p1_all <- do.call(ggpubr::ggarrange, c(p1_1, list(ncol = length(Seurat::Images(object_celltype)),
                                                                  nrow= 1,
                                                                  common.legend = TRUE,
                                                                  legend = "none",
                                                                  align = "hv")))

  p2_1 <- list()
  for(i in Seurat::Images(object_celltype)){
      p2_celltype <- OESingleCell::SpatialPlot(seurat_ob= object_celltype,group.by='secondcelltype',cols=colors,
                                           min.cutoff = 0.1,
                                           combine=FALSE,
                                           ncol= 2,
                                           alpha=c(0.1,0.5),
                                           pt.size.factor=opt$spotsize,
                                           pt.alpha =1,
                                           images=i)
      p2_1[[i]] <- p2_celltype[[1]]
  }

  p2_all <- do.call(ggpubr::ggarrange, c(p2_1, list(ncol = length(Seurat::Images(object_celltype)),
                                                                  nrow= 1,
                                                                  common.legend = TRUE, legend = "none",
                                                                  align = "hv")))

  combined_plot = cowplot::plot_grid(p1_all, p2_all,nrow = 2, align = "hv")
  combined_plot2 <- cowplot::plot_grid(combined_plot,ggpubr::get_legend(p1_1[[1]]),ncol=2,nrow=1,rel_widths = c(length(Seurat::Images(object_celltype)),1))
   OESingleCell::save_ggplots(file.path(output_dir,'celltype'),
                                plot=combined_plot2,
                                width=length(Seurat::Images(object_celltype))*4+3,
                                height=8,
                                dpi=200,bg='white')
  saveRDS(object_celltype,file=file.path(output_dir,'seurat_with_first_secondtype.rds'))



#scVis -i urther_analysis_2021-12-21/1.肿瘤上皮/result/seurat.rds -o ./ secondcelltype --resultname RCTD_results



}