#######################################################################################################################
docstring<- " example1:\\n\\n\\
sctool  -i  st.rds   -f rds  -o results   -assay sct  --dataslot counts  scDC  --celltype  celltype  --conditions group   --method BCa  --niter 1000 \\n\\n\\
example2:\\n\\n\\
sctool   -o results  scDC  --cellmeta cellmeta.tsv  --celltype  celltype  --conditions group   --method BCa  --niter 1000 "
sub_scDC <- subparsers$add_parser("scDC",
                                       description = docstring,
                                       formatter_class= 'argparse.RawTextHelpFormatter' ,
                                       #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
                                       argument_default = "True" ,
                                       help="Using scDC to perform single cell differential composition analysis")
sub_scDC$add_argument("--cellmeta", type = "character",
                      help = " celltype meta files ,required  formart:tsv ." )
sub_scDC$add_argument("--celltype", type = "character",
                      help = "[REQUIRED] celltype column names in cell meta dataframe ." )
sub_scDC$add_argument("--conditions", type = "character",
                      help = "[REQUIRED] conditions column names in cell meta dataframe ." )
sub_scDC$add_argument("--method", type = "character", default = "BCa",choices=c("percentile", "BCa", "multinom"),
                      help = "[REQUIRED]. method  for calCI [default: %(default)s] " )
sub_scDC$add_argument("--niter", default = 1000,type = "integer",
                    help = " iteration  for resampling [default:\"%(default)s\"]")
sub_scDC$add_argument("--order_list", default = NULL,
                    help = " conditions order list  for line graph plot [default:\"%(default)s\"]")
sub_scDC$add_argument("--rerun", default = TRUE,
                    help = " rerun  scDC_noClustering: True or False [default:\"%(default)s\"]")
# ================ Subcmd: st_deconv, deconvolute the cell type composition for spatial transcriptomics or bulk transcriptomics using scRNA-seq data ========
args <- commandArgs(TRUE)

if ( "scDC"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "scDC"){
    futile.logger::flog.info("step1: import celltype meta data information ")#===================================================
    if(!is.null(opt$cellmeta)){
      suppressMessages(cellmeta <- readr::read_tsv(opt$cellmeta))
    }else if(!is.null(opt$input)){
      suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                      informat = opt$informat,
                                                      assays = assays,
                                                      data.use = dataslots,  # data slot is enough
                                                      verbose = F))
       cellmeta <- OESingleCell::colData(data_ob)
    }
    cellmeta <- cellmeta %>%
                dplyr::select(opt$celltype,opt$conditions)%>%
                dplyr::rename(group=opt$conditions)%>%
                dplyr::rename(celltype=opt$celltype)
    use_celltype<-cellmeta$celltype%>%unique%>%c()
    subject<- as.character(cellmeta$group)
    cellTypes <- as.character(cellmeta$celltype) # vector of cell types or cell clusters
    cond <- as.character(cellmeta$group) # vector for condition
    #####
    condition_vec <- c()
    for(s in unique(subject)[order(unique(subject))]){
      # get condition associated with this sample:
      condition <- subset(cellmeta, as.character(group) == s) %>% .$group %>% unique %>% as.character
      condition_vec <- c(condition_vec, rep(condition, length(unique(cellTypes))))
    }
    ### set colors
    suppressMessages(cols<- OESingleCell::SelectColors( palette ="blindless",value =  length(use_celltype)))

    futile.logger::flog.info("step2:run scDC ")#========================================================================
    if(file.exists(glue::glue("{output_dir}/by_{opt$conditions}/res_scDC_noClust.rds")) & opt$rerun==FALSE ){
      futile.logger::flog.info(glue::glue("Using the existed scDC result file:{output_dir}/by_{opt$conditions}/res_scDC_noClust.rds"))
      res_scDC_noClust<-readRDS(glue::glue("{output_dir}/by_{opt$conditions}/res_scDC_noClust.rds"))
    }else{
      res_scDC_noClust <- scDC::scDC_noClustering(cellTypes,
                                                  subject,
                                                  calCI = TRUE,
                                                  calCI_method = opt$method,
                                                  ncores = opt$ncores,
                                                  nboot = opt$niter)
      ##mkdirs
      if (! file.exists(glue::glue("{output_dir}/by_{opt$conditions}"))) {
        dir.create(glue::glue("{output_dir}/by_{opt$conditions}") ,recursive = T)
      }

      saveRDS(res_scDC_noClust,glue::glue("{output_dir}/by_{opt$conditions}/res_scDC_noClust.rds"))
    }
    futile.logger::flog.info("step3:barplot and density plot ")#========================================================
    p1<-scDC::barplotCI(res_scDC_noClust, condition_vec)
    p2<-scDC::densityCI(res_scDC_noClust, condition_vec)
    OESingleCell::save_ggplots(filename= glue::glue("{output_dir}/by_{opt$conditions}/1.柱状图/barplot_for_diff_conditions"),
                                    plot =ggplotify::as.ggplot(p1) ,
                                    width =2*length(use_celltype),
                                    limitsize = FALSE,
                                    height= 6 )
    OESingleCell::save_ggplots(filename= glue::glue("{output_dir}/by_{opt$conditions}/2.密度图/density_plot_for_diff_conditions"),
                                    plot =ggplotify::as.ggplot(p2) ,
                                    width = 2*length(use_celltype),
                                    limitsize = FALSE,
                                    height= 6 )
    futile.logger::flog.info("step4:heatmap plot ")#====================================================================

    if(is.null(opt$order_list)){ opt$order_list <- subject%>%unique()%>%unlist()}
    proportion_df <- res_scDC_noClust$results
    proportion_df$median <- apply(res_scDC_noClust$thetastar, 1, median)
    proportion_df$cond <-  factor(proportion_df$subject, levels = c(opt$order_list%>%strsplit(",")%>%unlist()))
    proportion_df$method <- factor(proportion_df$method, levels = c("BCa"))
    n_method <- length(unique(proportion_df$method))
    n_celltype <- length(unique(proportion_df$cellTypes))
    proportion_df<-proportion_df %>%dplyr::arrange()%>%dplyr::filter(method!="")
    proportion_df[is.na(proportion_df)] <- 0
    proportion_df$median <-proportion_df$median *100
    proportion_df$conf_low<-as.numeric(proportion_df$conf_low)*100
    proportion_df$conf_high<-as.numeric(proportion_df$conf_high)*100
    proportion_df%>%head
    merge_data<- proportion_df%>%
                 dplyr::select("cellTypes","subject","median")%>%
                 tidyr::spread("subject","median")%>%
                 tibble::column_to_rownames("cellTypes")

    merge_data<- merge_data %>% scale(center = TRUE, scale = TRUE)
    pdf(NULL) ##for removing  Rplot.pdf files
    p<-ComplexHeatmap::Heatmap(as.matrix(merge_data), # 输入数据为矩阵
                                    name = 'condition special celltype', # 热图图例名称
                                    col = rev(RColorBrewer::brewer.pal(11, "Spectral")), #配色方案
                                    border = FALSE, # 显示边框
                                    cluster_rows = T, #按照行进行聚类
                                    cluster_columns = F, #按照列进行聚类
                                    #column_km = 2, # 划分列聚类
                                    #row_km = 2, # 划分行聚类
                                    show_column_names = T, #显示列名
                                    show_row_names = T, #显示行名与否
                                    column_names_side ="bottom", #ifelse(opt$colcluster == "FALSE", "top", "bottom"), , #列名位置 top bottom
                                    row_names_rot = 0, ##行名旋转角度
                                    column_names_rot = 45, ##列名旋转角度
                                    row_names_side = "right",#ifelse(opt$rowcluster == "FALSE", "left", "right"), #行名位置 left right
                                    column_names_gp = grid::gpar(fontsize = 10), ##列字体大小
                                    row_names_gp = grid::gpar(fontsize = 10), ##行字体大小
                                    #alpha=0.5, #设置图片颜色透明度
                                    # width = ncol(gene_data)*unit(5, "mm"), # 格子的宽度
                                    # height = nrow(gene_data)*unit(2, "mm"),# 格子的高度
                                    top_annotation = NULL, # 添加左侧注释信息
                                    right_annotation = NULL, # 添加右侧注释信息
                                    left_annotation = NULL, # 添加左侧注释信息
                                    bottom_annotation = NULL, # 添加下方注释信息
                                    heatmap_legend_param = list(title = "",
                                                                #title_position = "topcenter", # 标题相对图例的位置 topcenter, topleft, leftcenter, lefttop.
                                                                #at=c(0,10), #图例范围
                                                                labels_gp = grid::gpar(fontsize = 10), ##图例字体
                                                                legend_direction = "vertical", ##图例放置
                                                                legend_height = unit(5, "cm") #图例长度
                                    ))

    OESingleCell::save_ggplots(filename= glue::glue("{output_dir}/by_{opt$conditions}/3.热图/heatmap_for_diff_conditions"),
                                    plot =ggplotify::as.ggplot(p) ,
                                    width = 6,
                                    height= 8 )
    futile.logger::flog.info("step5: line graph plot ")#====================================================================
    ## split by conditions
    for( i in 1:length(use_celltype)){
           p<- proportion_df %>%
               dplyr::filter(cellTypes==use_celltype[i]) %>%
               ggplot(aes(x = cond,
                          y = median,
                          group = 1)) +
               geom_point(colour = "black", size =1) +
               geom_line(color="grey", size=1) +
               geom_ribbon(aes(ymin =conf_low, ymax=conf_high, fill = cellTypes), alpha = 0.3)+
               theme_bw() +
               theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
               ggtitle(use_celltype[i]) +
               xlab("Conditions") +
               ylab("Percentage(Bootstrap)(%)")+
               guides(fill="none")+
               scale_fill_manual(values = cols[i]  )#+
         OESingleCell::save_ggplots(filename= glue::glue("{output_dir}/by_{opt$conditions}/4.折线图/split_by_celltype/{use_celltype[i]}_line_point"),
                                    plot = p ,
                                    width = 8,
                                    height=4 )
        }
    ## summary plot
     p<- proportion_df %>%
               ggplot(aes(x = cond,
                          y = median,
                          group = 1)) +
               geom_point(colour = "black", size =1) +
               geom_line(color="grey", size=1) +
               geom_ribbon(aes(ymin =conf_low, ymax=conf_high, fill = cellTypes), alpha = 0.3)+
               theme_bw() +
               theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
               ggtitle("") +
               xlab("Conditions") +
               ylab("Percentage(Bootstrap)(%)")+
               guides(fill="none")+
               scale_fill_manual(values = cols  )+facet_wrap(~cellTypes,ncol=4)
         OESingleCell::save_ggplots(filename= glue::glue("{output_dir}/by_{opt$conditions}/4.折线图/combined_line_point"),
                                    plot = p ,
                                    width = 16,
                                    height= length(use_celltype)/2 )
    # -
    readr::write_tsv(proportion_df, file=glue::glue('{output_dir}/by_{opt$conditions}/scDC_results.xls'))
    futile.logger::flog.info("step6: running fitGLM")#====================================================================
    library(broom.mixed)
    res_GLM <- scDC::fitGLM(res_scDC_noClust, condition_vec, subject_effect =FALSE, pairwise = FALSE)
    sum_GLM <- summary(res_GLM$pool_res_fixed)
    readr::write_tsv(sum_GLM, file=glue::glue('{output_dir}/by_{opt$conditions}/scDC_summary.xls'))
    ## save session information
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}