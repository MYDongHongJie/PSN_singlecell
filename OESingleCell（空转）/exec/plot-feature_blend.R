feature_blend <- function(data_ob,
                          sample,
                          genelist,
                          cols_set = c("lightgrey", "red", "darkblue"),
                          blend.threshold = 0.2,
                          image.alpha = 1 ,
                          pt.size.1 = 1.0,
                          pt.size.2 = 1.0,
                          crop=crop,
                          output_dir) {
    futile.logger::flog.info(glue::glue("1:{sample}-subset单个样本的seurat对象")) ##======================================
    data_ob_sub <- subset(data_ob, subset = sampleid == sample)
    if (!is.null(Seurat::Images(data_ob_sub))) {
        unuse_images <- Seurat::Images(data_ob_sub)[!Seurat::Images(data_ob_sub) %in% (data_ob_sub@meta.data$sampleid %>% unique)]
        if (length(unuse_images) > 0) { data_ob_sub@images[unuse_images] <- NULL }
    }
    futile.logger::flog.info(glue::glue("2:{sample}-绘制Dimplot图")) ##==================================================
    plot_1 <- Seurat::FeaturePlot(data_ob_sub,
                                  features = genelist,
                                  blend = TRUE,
                                  cols = cols_set,
                                  combine = T,
                                  blend.threshold = blend.threshold,
                                  pt.size = pt.size.1)
    futile.logger::flog.info(glue::glue("3:{sample}-绘制空间featplot blend图,并获取三个子图配色情况")) ##====================
    plot <- Seurat::FeaturePlot(data_ob_sub,
                                features = genelist,
                                blend = TRUE,
                                cols = cols_set,
                                combine = F,
                                blend.threshold = blend.threshold,
                                pt.size = pt.size.2)
    data_ob_sub <- Seurat::AddMetaData(data_ob_sub, metadata = cbind(plot[[1]]$data, plot[[2]]$data, plot[[3]]$data))
    # ggplot_build 函数是 R 语言中 ggplot2 包的一部分，它用于构建一个 ggplot 图形对象的详细信息，包括数据、
    # 图层、标尺、坐标轴、统计变换等等的细节。这个函数返回一个列表，其中包含了用于绘图的所有信息。
    get_cols <- function(plot_object, gene) {
        g <- ggplot2::ggplot_build(plot_object)
        a <- cbind(tibble::rownames_to_column(plot_object$data), g$data[[1]]) %>%
                dplyr::select(gene, "colour") %>%
                unique
        col_sum <- a$colour
        names(col_sum) <- a[, gene] %>% as.character
        return(col_sum)
    }

    genelist_merge <- paste(genelist[1], genelist[2], sep = "_")
    col_sum_1 <- get_cols(plot[[1]], genelist[1])
    col_sum_2 <- get_cols(plot[[2]], genelist[2])
    col_sum_3 <- get_cols(plot[[3]], genelist_merge)
    futile.logger::flog.info(glue::glue("4:{sample}-根据3中提取的配色，绘制空间Dimplot图片")) #===========================
    plot_2_1 <- Seurat::SpatialDimPlot(data_ob_sub,
                                       group.by = genelist[1],
                                       image.alpha = image.alpha,
                                       pt.size = pt.size.2,
                                       crop = crop,
                                       cols = col_sum_1) +
            ggtitle(genelist[1]) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5,size=16)) & Seurat::NoLegend()
    plot_2_2 <- Seurat::SpatialDimPlot(data_ob_sub,
                                       group.by = genelist[2],
                                       image.alpha = image.alpha,
                                       pt.size = pt.size.2,
                                       crop = crop,
                                       cols = col_sum_2) +
            ggtitle(genelist[2]) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5,size=16)) & Seurat::NoLegend()
    plot_2_3 <- Seurat::SpatialDimPlot(data_ob_sub,
                                       group.by = genelist_merge,
                                       image.alpha = image.alpha,
                                       pt.size = pt.size.2,
                                       crop = crop,
                                       cols = col_sum_3) +
            ggtitle(genelist_merge) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5,size=16)) & Seurat::NoLegend()
    plot_2 <- patchwork::wrap_plots(plot_2_1,
                                    plot_2_2,
                                    plot_2_3,
                                    plot[[4]],
                                    nrow = 1,
                                    widths = c(1, 1, 1, 0.9))
    futile.logger::flog.info(glue::glue("5:{sample}-合并并保存图片")) #====================================================
    plot_1_2 <- patchwork::wrap_plots(plot_1, plot_2, nrow = 2)
    OESingleCell::save_ggplots(glue::glue("{output_dir}/feature_blend_{genelist[1]}_{genelist[1]}_for_{sample}"),
                               plot = plot_1_2,
                               width = 16,
                               height = 9,
                               limitsize = FALSE,
                               dpi = 300)
}

docstring <- "example:\\n\\n\\
  scVis -i seurat.rds -f rds -o results --HE TRUE --assay SCT feature_blend  --misclist RCTD_results  --interest T_cells,B_cells   --colors lightgrey,red,darkblue --image.alpha 0.7   --crop TRUE --pt.size.1 0.1 --pt.size.2 1.0  \\n\\n\\
  scVis -i seurat.rds -f rds -o results --HE TRUE --assay SCT feature_blend  --interest gene1,gene2   --colors lightgrey,red,darkblue --image.alpha 0.7  --crop TRUE --pt.size.1 0.1 --pt.size.2 1.0  "
sub_feature_blend <- subparsers$add_parser(
        "feature_blend",
        description = docstring,
        formatter_class = 'argparse.RawTextHelpFormatter',
        #formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
        argument_default = "True",
        help = "绘制2种基因或者细胞类型feature blend图")
sub_feature_blend$add_argument(
        "--misclist",
        type = "character",
        default = NULL,
        help = "如果展示的是细胞类型鉴定情况，表示在seurat对象misc中，细胞类型鉴定结果存放的list名称，需要根据鉴定算法来调整，如：RCTD_results [default: %(default)s] ")
sub_feature_blend$add_argument(
        "--interest",
        type = "character",
        default = NULL,
        help = "2个感兴趣的基因或者细胞类型，以逗号分割[default: %(default)s] ")
sub_feature_blend$add_argument(
        "--colors",
        type = "character",
        default = "lightgrey,red,darkblue",
        help = "针对2个感兴趣的基因或者细胞类型，提供颜色，以逗号分割,如:lightgrey,#ff0000,#00ff00  [default: %(default)s] ")
sub_feature_blend$add_argument(
        "--image.alpha",
        type = "double",
        default = 1.0,
        help = "空间HE图片透明度，范围从0到1.[default: %(default)s] ")
sub_feature_blend$add_argument(
        "--crop",
        type = "character",
        default = FALSE,
        help = "空间HE图片是否进行裁剪.[default: %(default)s] ")
sub_feature_blend$add_argument(
        "--pt.size.1",
        type = "double",
        default = 0.1,
        help = "dimplot图中点大小设置[default: %(default)s] ")
sub_feature_blend$add_argument(
        "--pt.size.2",
        type = "double",
        default = 1.2,
        help = "空间dimplot图中点大小设置[default: %(default)s] ")
sub_feature_blend$add_argument(
        "--blend.threshold",
        type = "double",
        default = 0.2,
        help = "feature blend值设置，用于控制调色板范围，从弱信号到强信号的颜色阈值，范围从0到1[default: %(default)s] ")
args <- commandArgs(TRUE)
if ("feature_blend" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "feature_blend") {
        #===================================================================================================================
        futile.logger::flog.info("step1:导入seurat对象")
        suppressMessages(data_ob <- OESingleCell::ReadX(input = opt$input,
                                                        informat = opt$informat,
                                                        assays = assays,
                                                        data.use = dataslot,
                                                        verbose = F))
        futile.logger::flog.info("判断是否是细胞类型鉴定结果，根据情况进行数据整理")
        if (!is.null(opt$misclist)) {
            data_ob@meta.data <- data_ob@meta.data %>% cbind(data_ob@misc[[opt$misclist]] %>% tibble::column_to_rownames("barcodes"))
        }
        futile.logger::flog.info("判断是否需要进行数据抽提")
        if (!is.null(opt$predicate)) {
            futile.logger::flog.info(glue::glue("get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
            df <- OESingleCell::colData(data_ob)
            desired_cells <- subset(df, eval(parse(text = opt$predicate)))
            data_ob <- subset(data_ob, cells = rownames(desired_cells))
            if (!is.null(Seurat::Images(data_ob))) {
                unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique)]
                if (length(unuse_images) > 0) { data_ob@images[unuse_images] <- NULL }
            }
        }
        #===============================================================================================================
        images <- Seurat::Images(data_ob)
        if (opt$HE==FALSE) { image.alpha <- 0 }else{image.alpha<-opt$image.alpha}
        futile.logger::flog.info("step2:针对每个样本，绘制图片并保存")
        parallel::mclapply(images, function(sample) {
            tryCatch({
                feature_blend(data_ob = data_ob,
                              sample = sample,
                              genelist = unlist(stringr::str_split(opt$interest, ",")),
                              cols_set = unlist(stringr::str_split(opt$colors, ",")),
                              blend.threshold = opt$blend.threshold,
                              image.alpha = image.alpha,
                              pt.size.1 = opt$pt.size.1,
                              pt.size.2 = opt$pt.size.2,
                              crop=as.logical(opt$crop),
                              output_dir = output_dir)
            }, error = function(e) { futile.logger::flog.info(glue::glue("{sample}绘图失败，报错如下: {e}")); break })
        }, mc.cores = opt$ncores)
        futile.logger::flog.info("step3:输出运行日志")
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}