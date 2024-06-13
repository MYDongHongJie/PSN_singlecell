#' Creates a matrix for celltype interaction in spots
#'
#' @return  a matrix for celltype-celltype interactions summary for all spots
#'
#' @param seurat.obj  A singular preprocessed Seurat object
#' @param celltypes_spatial  celltypes from seurat_ob_ref
#' @sample  names prefix for output interaction csv
#' @outdir output directory
cell_interaction <- function(
  obj,
  celltypes_spatial,
  sample,
  group.by="clusters",
  min.prop=0.1,
  output_dir
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    if (sample != "combined_all_slice") {
        obj <- base::subset(obj, subset = sampleid == sample)
    }
    celltypes <- obj@meta.data %>% dplyr::select(celltypes_spatial)
    #celltypes["unknown"] <- 0
    #celltypes[which(rowSums(celltypes) == 0), "unknown"] <- 1
    clusters <- obj@meta.data %>% dplyr::select(group.by)
    colnames(clusters) <- "clusters"
    celltypes_out <- cbind(clusters, celltypes)
    interaction_matrix <- matrix(0, ncol = length(celltypes_spatial), nrow = length(celltypes_spatial))
    rownames(interaction_matrix) <- unique(celltypes_spatial)
    colnames(interaction_matrix) <- unique(celltypes_spatial)
    for (i in 1:nrow(celltypes)) {
        temp1 <- sort(colnames(celltypes[i, celltypes[i, ] > min.prop]), decreasing = F)
        if (length(temp1) > 0) {
            temp <- t(combn(temp1, 2))
            for (j in 1:nrow(temp)) {
                interaction_matrix[
                    as.character(temp[j, 1]),
                    as.character(temp[j, 2])
                ] <- interaction_matrix[
                    as.character(temp[j, 1]),
                    as.character(temp[j, 2])
                ] + 1
            }
            rm(temp)
        }
        rm(temp1)
    }

    write.table(
        celltypes_out,
        file = glue::glue("{output_dir}/../{sample}-celltype.xls"),
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )
    write.table(
        interaction_matrix,
        file = glue::glue("{output_dir}/{sample}-interactions.xls"),
        sep = "\t",
        col.names = NA,
        quote = FALSE
    )
    cell_list <- list(interaction_matrix, celltypes_out)
    return(cell_list)
}

#' plot circos for  celltype-celltype interaction
#'
#' @return   celltype's plot color for this celltype-celltype interaction  matrix
#' @param interaction_matrix a matrix for celltype-celltype interactions summary for all spots
#' @param color_pelette  color_pelette for all this project's celltype

plot_interaction <- function(
  slice,
  interaction_matrix,
  color_pelette,
  formart,
  output_dir
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    ## ======================================================================================================
    if (formart == "pdf") {
        pdf(
            file.path(
                output_dir,
                glue::glue("{slice}_celltype_interaction.{formart}")
            ),
            width = 12,
            height = 8
        )
    } else if (formart == "png") {
        png(
            file.path(
                output_dir,
                glue::glue("{slice}_celltype_interaction.{formart}")
            ),
            width = 12,
            height = 8,
            units = "in",
            res = 1000
        )
    }
    # =======================================================================================================
    layout(matrix(c(1, 2), nrow = 1), widths = c(0.8, 0.5))
    color_used <- color_pelette[colnames(interaction_matrix)]
    col <- matrix(
        rep(color_used, each = ncol(interaction_matrix), T),
        nrow = nrow(interaction_matrix),
        ncol = ncol(interaction_matrix)
    )
    rownames(col) <- rownames(interaction_matrix)
    colnames(col) <- colnames(interaction_matrix)
    circlize::chordDiagram(
        as.matrix(interaction_matrix),
        grid.col = color_used,
        col = col,
        annotationTrack = "grid"
    )
    lgd_points <- ComplexHeatmap::Legend(
        labels = names(color_used),
        ncol = ceiling(length(names(color_used)) / 40),
        size = 1,
        legend_gp = grid::gpar(fill = color_used),
        title_position = "topcenter",
        title = "Cell Types"
    )
    grid::pushViewport(grid::viewport(x = 0.8, y = 0.5))
    grid::grid.draw(lgd_points)
    grid::upViewport()
    dev.off()
    rm(color_used)
}
## plot heatmap
#' plot circos for  celltype-celltype interaction
#'
#' @return   celltype's plot color for this celltype-celltype interaction  matrix
#' @param interaction_matrix a matrix for celltype-celltype interactions summary for all spots
#' @param color_pelette  color_pelette for all this project's celltype

plot_heatmap <- function(
  celltype,
  sample_name,
  formart,
  prefix,
  palette,
  output_dir
) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    ## =======================================================================================================
    celltype <- celltype %>%
        tibble::rownames_to_column("barcodes") %>%
        dplyr::arrange(clusters)
    rownames(celltype) <- celltype[["barcodes"]]
    spot_info <- celltype %>%
        dplyr::select("barcodes", clusters) %>%
        tidyr::separate("barcodes", into = c("sampleid", "raw_barcode"), sep = "-", remove = FALSE)
    if (!all(unique(spot_info$sampleid) %in% sample_name)) {
        spot_info <- spot_info %>% dplyr::mutate(sampleid = as.character(sample_name))
    }
    celltype_new <- celltype %>% dplyr::select(!barcodes & !clusters)
    ##
    sampleid_name <- unique(spot_info$sampleid)
    sampleid.colors <- structure(
        names = sampleid_name,
        OESingleCell::SelectColors(1:(length(unique(spot_info$sampleid))+10), palette = palette)[11:(length(unique(spot_info$sampleid))+10)]
    )
    ##
    cluster_name <- unique(spot_info$clusters)
    clusters.colors <- structure(
        names = as.character(1:max(as.integer(cluster_name))),
        OESingleCell::SelectColors(1:max(as.integer(cluster_name)), palette = palette)
    )
    clusters.colors <- clusters.colors[as.character(cluster_name)]

    if (prefix == "combined_all_slice") {
        colAnn <- ComplexHeatmap::HeatmapAnnotation(
                df = data.frame(
                        sampleid = factor(spot_info$sampleid,levels=unique(spot_info$sampleid)),
                        clusters = as.factor(spot_info$clusters)
                ),
                col = list(
                        sampleid = sampleid.colors,
                        clusters = clusters.colors
                ),
                show_annotation_name = T,
                annotation_legend_param = list(
                        sampleid = list(ncol = 1, direction = "horizontal"),
                        clusters = list(ncol = 1, direction = "horizontal")
                )
        )
        sampleid <- NULL
    }else {
        colAnn <- ComplexHeatmap::HeatmapAnnotation(
                df = data.frame(
                        clusters = as.factor(spot_info$clusters)
                ),
                col = list(
                        clusters = clusters.colors
                ),
                show_annotation_name = T,
                annotation_legend_param = list(
                        clusters = list(ncol = 1, direction = "horizontal")
                )
        )
        sampleid <- names(sampleid.colors)
    }
    ## =======================================================================================================
    if (formart == "pdf") {
        pdf(
            file.path(output_dir, glue::glue("{prefix}-celltype_heatmap.{formart}")),
            width = 10,
            height = 8
        )
    } else if (formart == "png") {
        png(
            file.path(output_dir, glue::glue("{prefix}-celltype_heatmap.{formart}")),
            width = 10,
            height = 8,
            units = "in",
            res = 600
        )
    }

    ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
        as.matrix(t(celltype_new)),
        col = circlize::colorRamp2(c(0, 1), c("whitesmoke", "mediumvioletred")), # grDevices::colorRampPalette(c("#406AA8", "white", "#D91216"))(256),
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        column_title = sampleid, #列标题
        column_title_side = "top", #列标题位置
        show_row_names = T,
        row_names_rot = 0,
        row_names_side = "left",
        column_names_gp = grid::gpar(fontsize = 10), ## 列字体大小
        row_names_gp = grid::gpar(fontsize = 10), ## 行字体大小
        heatmap_legend_param = list(
            title = "cell_type",
            # title_position = "topcenter", # 图例标题位置
            at = c(0, 1), # 图例范围
            legend_direction = "vertical",
            legend_height = unit(2, "cm") # 图例长度
        ),
        top_annotation = colAnn
    ))
    dev.off()
}

## /results/seurat.rds
docstring <- "example:\\n\\n\\
  scVis -i seurat.rds -f rds -o results --assay Spatial spatialpieplot --pt.alpha TRUE  "
sub_spatialpieplot <- subparsers$add_parser(
    "spatialpieplot",
    description = docstring,
    formatter_class = "argparse.RawTextHelpFormatter",
    # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "plot for visium celltype parse results"
)
sub_spatialpieplot$add_argument(
    "--group.by",
    type = "character",
    default = "clusters",
    help = "[REQUIRED]cell barcodes group by in the plot[default: %(default)s] "
)
sub_spatialpieplot$add_argument(
    "--piescale",
    type = "double",
    default = 0.38,
    help = "[REQUIRED]the pie scale  in the plot[default: %(default)s] "
)
sub_spatialpieplot$add_argument(
    "--misclist",
    type = "character",
    default = "spotlight_results",
    help = "[REQUIRED]the spotlight result's name in seurat's miscs [default: %(default)s] "
)
sub_spatialpieplot$add_argument(
    "--pt.alpha",
    type = "character",
    default = "TRUE",
    help = "[OPTIONAL]the pie scale  in the plot[default: %(default)s] "
)
sub_spatialpieplot$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)
sub_spatialpieplot$add_argument(
    "--palette",
    type = "character",
    default = "customecol2",
    help = paste0(
        "the discrete color schema mapped to the cell annotations specified by --groupby.[default: %(default)s]",
        " Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
        " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"
    )
)
sub_spatialpieplot$add_argument(
    "--common.legend",
    default = 'FALSE',
    help = "whether to generate common legend for spatialplot.[default: %(default)s]"
)
args <- commandArgs(TRUE)
if ("spatialpieplot" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "spatialpieplot") {
    # =====================================================================================================================
    ####  load query object data and setting color palette
    suppressMessages(data_ob <- OESingleCell::ReadX(
        input = opt$input,
        informat = opt$informat,
        assays = assays[1],
        data.use = dataslots,
        verbose = F
    ))
    Seurat::DefaultAssay(data_ob) <- assays[1]
    ## remove used images
    if (!is.null(Seurat::Images(data_ob))) {
        unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique())]
        if (length(unuse_images) > 0) {
            data_ob@images[unuse_images] <- NULL
        }
    }

    if (is.factor(data_ob@meta.data[, opt$group.by])) {
        data_ob@meta.data[, opt$group.by] <- droplevels(data_ob@meta.data[, opt$group.by])
    }
    # data_ob[[opt$groupby]] <- factor(data_ob[[opt$groupby]], levels = sort(unique(data_ob[[opt$groupby]])))
    ## add  celltype informtaion into seurat object metadata
    # data_ob@meta.data <- data_ob@meta.data %>%
    #     cbind(data_ob@misc[[opt$misclist]] %>% tibble::column_to_rownames("barcodes"))
    images <- Seurat::Images(data_ob)

    ### get all celltype's name
    celltypes_spatial <- colnames(data_ob@misc[[opt$misclist]] %>% dplyr::select(!c("barcodes")))
    ## set solor palette
    # if (length(celltypes_spatial) <= 9) {
    #     palette <- "visium9"
    # }
    # if (length(celltypes_spatial) > 9 && length(celltypes_spatial) <= 22) {
    #     palette <- "ditto"
    # }
    # if (length(celltypes_spatial) > 22 && length(celltypes_spatial) <= 50) {
    #     palette <- "col50"
    # }
    # if (length(celltypes_spatial) > 50 && length(celltypes_spatial) <= 90) {
    #     palette <- "blindless"
    # }
    palette <- opt$palette
    # =====================================================================================================================
    # split  celltype's percentage for each celltype
    crop <- as.logical(opt$crop)
    parallel::mclapply(celltypes_spatial, function(celltype) {
        plot <- OESingleCell::SpatialPlot(
          seurat_ob = data_ob,
          features = celltype,
          min.cutoff = "q0",
          max.cutoff = "q90",
          combine = TRUE,
          cols = "spectral",
          ncol = length(images),
          crop = crop
          # alpha = 1,
          # pt.alpha = opt$pt.alpha
        )
        OESingleCell::save_ggplots(
            file.path(
                glue::glue("{output_dir}/1.celltype_spatialplot"),
                "split_celltype_results",
                glue::glue("{celltype}_percent_for_each_sample_plot")
            ),
            plot = plot,
            width = length(images) * 3 + sum(stringr::str_length(celltype)) * 0.14,
            height = 3,
            limitsize = FALSE,
            dpi = 600
        )
    }, mc.cores = 8)
    # =====================================================================================================================
    ## split celltype's percentage for each sample
    # for(slice in images ){
    parallel::mclapply(images, function(slice) {
        # plot_slice <-list()
        data <- Seurat::FetchData(object = data_ob, vars = c("sampleid", celltypes_spatial)) %>%
            dplyr::filter(sampleid == slice) %>%
            dplyr::select(!sampleid)
        min.cutoff <- Seurat::SetQuantile("q0", data[, celltypes_spatial])
        max.cutoff <- Seurat::SetQuantile("q90", data[, celltypes_spatial])
        plots <- list()
        for (i in 1:length(celltypes_spatial)) {
            p <- OESingleCell::SpatialPlot(
              seurat_ob = data_ob,
              features = celltypes_spatial[i],
              min.cutoff = min.cutoff,
              max.cutoff = max.cutoff,
              combine = FALSE,
              images = slice,
              cols = "spectral",
              crop = crop
              #alpha = 1,
              #pt.alpha = opt$pt.alpha
            )
            plots[[i]] <- p[[1]]
            # ggplot2::ggsave(plots[[i]] ,filename=glue::glue("{slice}.{celltypes_spatial[i]}.pdf"))
        }
        ncol <- ifelse(length(celltypes_spatial) > 3, 4, length(celltypes_spatial))
        nrow <- ceiling(length(celltypes_spatial) / ncol)
        plot_sample <- do.call(
            ggpubr::ggarrange,
            c(
                plots,
                list(
                    ncol = ncol,
                    nrow = nrow,
                    # common.legend = TRUE,
                    ## legend = "right",
                    align = "hv"
                )
            )
        )
        # plot_sample<-patchwork::wrap_plots(plots,ncol=1)
        OESingleCell::save_ggplots(
            glue::glue("{output_dir}/1.celltype_spatialplot/split_sample_results/celltype_percent_for_{slice}_plot"),
            plot = plot_sample,
            width = ncol * 5,
            height = nrow * 5,
            limitsize = FALSE,
            dpi = 300
        )
    }, mc.cores = 8)

    # =====================================================================================================================
    # scatterpieplot for each slice
    parallel::mclapply(images, function(slice) {
        OESingleCell::SpatialScatterPie(
            object = data_ob,
            misclist = opt$misclist,
            cell_types_interest = celltypes_spatial,
            assay = "Spatial",
            images = slice,
            pie.alpha = 1,
            pie.scale = opt$piescale,
            output_dir = glue::glue("{output_dir}/1.celltype_spatialplot"),
            cols = OESingleCell::SelectColors(celltypes_spatial, palette = palette)
        )
    }, mc.cores = 8)
    # =====================================================================================================================
    ## plot circos  for all slice's spots celltype-celltype interaction
    interaction_matrix <- cell_interaction(data_ob, celltypes_spatial, "combined_all_slice", opt$group.by,0.1, glue::glue("{output_dir}/2.celltype_interaction"))

    color_pelette <- OESingleCell::SelectColors(celltypes_spatial, palette = palette)
    plot_interaction("all", interaction_matrix[[1]], color_pelette, "png", glue::glue("{output_dir}/2.celltype_interaction"))
    plot_interaction("all", interaction_matrix[[1]], color_pelette, "pdf", glue::glue("{output_dir}/2.celltype_interaction"))
    plot_heatmap(interaction_matrix[[2]], sample_name = data_ob$sampleid, "pdf", "combined_all_slice", palette, glue::glue("{output_dir}/3.celltype_heatmap"))
    plot_heatmap(interaction_matrix[[2]], sample_name = data_ob$sampleid, "png", "combined_all_slice", palette, glue::glue("{output_dir}/3.celltype_heatmap"))
    ## plot circos  for each slice's spots celltype-celltype interaction
    parallel::mclapply(images, function(slice) {
        interaction_matrix <- cell_interaction(data_ob, celltypes_spatial, slice, opt$group.by, 0.1, glue::glue("{output_dir}/2.celltype_interaction"))
        plot_interaction(glue::glue("sample_{slice}"), interaction_matrix[[1]], color_pelette, "png", glue::glue("{output_dir}/2.celltype_interaction"))
        plot_interaction(glue::glue("sample_{slice}"), interaction_matrix[[1]], color_pelette, "pdf", glue::glue("{output_dir}/2.celltype_interaction"))
        plot_heatmap(interaction_matrix[[2]], sample_name = slice, "pdf", glue::glue("sample_{slice}"), palette, glue::glue("{output_dir}/3.celltype_heatmap"))
        plot_heatmap(interaction_matrix[[2]], sample_name = slice, "png", glue::glue("sample_{slice}"), palette, glue::glue("{output_dir}/3.celltype_heatmap"))
    }, mc.cores = 8)

    # ==================================================================================================================
    # spatial dimplot of top1_celltype and top2_celltype
    top_celltype <- c("top1_celltype", "top2_celltype")
    lapply(top_celltype, FUN = function(top_celltype) {
        p_top <- OESingleCell::SpatialPlot(
                seurat_ob = data_ob,
                group.by = top_celltype,
                combine = TRUE,
                cols = palette,
                crop = crop,
                ncol = length(images),
                image.alpha = 0.7,
                HE = FALSE,
                common.legend = as.logical(opt$common.legend)
        )
        OESingleCell::save_ggplots(
            glue::glue("{output_dir}/4.top_celltype_spatialplot/{top_celltype}_spatial_plot"),
            plot = p_top,
            width = length(images) * 6,
            height = 6,
            limitsize = FALSE,
            dpi = 300
        )
    })
    ############################### save top2 celltype information #######################################
    if( "rawbc" %in%  colnames(data_ob@meta.data) ){
        celltype <-  data_ob@meta.data %>% dplyr::select("rawbc","sampleid","group","clusters","top1_celltype","top2_celltype") %>% dplyr::rename("barcode" = "rawbc")
    }else{
        celltype <-  data_ob@meta.data %>% dplyr::select("orig.ident","sampleid","group","clusters","top1_celltype","top2_celltype") %>% dplyr::rename("barcode" = "orig.ident")
    }

    write.table(celltype,
        file=file.path(output_dir,"celltype_infor.csv"),
        quote=F,
        sep=",",
        row.names=F)
    
    #### copy readme file 
    if (grepl("SPOTlight", opt$misclist, ignore.case = T)) {
        if (!file.exists(file.path(output_dir, "SPOTlight细胞类型鉴定结果说明.docx"))) {
            file.copy("/public/dev_scRNA/oesinglecell3_test/document/SPOTlight细胞类型鉴定结果说明.docx", output_dir)
            message("SPOTlight细胞类型鉴定结果说明.doc 完成拷贝")
        }else {
            message("SPOTlight细胞类型鉴定结果说明.doc 已完成拷贝, 不再重复拷贝")
        }
    }else if (grepl("RCTD", opt$misclist, ignore.case = T)) {
        if (!file.exists(file.path(output_dir, "RCTD细胞类型鉴定结果说明.docx"))) {
            file.copy("/public/dev_scRNA/oesinglecell3_test/document/RCTD细胞类型鉴定结果说明.docx", output_dir)
            message("RCTD细胞类型鉴定结果说明.docx 完成拷贝")
        }else {
            message("RCTD细胞类型鉴定结果说明.docx 已完成拷贝, 不再重复拷贝")
        }
    }else if (grepl("spatialDWLS", opt$misclist, ignore.case = T)) {
        if (!file.exists(file.path(output_dir, "SpatialDWLS细胞类型鉴定分析说明文档.docx"))) {
            file.copy("/public/dev_scRNA/oesinglecell3_test/document/SpatialDWLS细胞类型鉴定分析说明文档.docx", output_dir)
            message("SpatialDWLS细胞类型鉴定分析说明文档.docx 完成拷贝")
        }else {
            message("SpatialDWLS细胞类型鉴定分析说明文档.docx 已完成拷贝, 不再重复拷贝")
        }
    }
    #  =====================================================================================================================
    ## save session information
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
}
}
