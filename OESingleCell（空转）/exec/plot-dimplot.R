## 绘制 面积 图
Plot_area=function(x, prop.by = "res.1", group.by = "sampleid", split.by = NULL, 
     ncol = NULL, cols = NULL) 
{
    df <- as.data.frame( x@meta.data[,c(prop.by, group.by)]) %>%
                    table() %>% as.data.frame()  %>% dplyr::rename(cell_number=Freq) %>% 
                    dplyr::arrange(get(group.by)) %>% 
                    dplyr::group_by(.dots= group.by) %>% 
                    dplyr::mutate(freq =(cell_number / sum(cell_number)) * 100) %>% 
                    as.data.frame()
    p =  ggplot(df, aes_string(x=group.by, y="freq", fill=prop.by,group=prop.by)) + 
  geom_area(size=0.4,colour="white")+
  labs(x = NULL, y = "Proportion [%]")+
  theme_bw()+
  scale_x_discrete(expand = c(0, 0.02)) +scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25))+
  scale_fill_manual(values=cols)+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        strip.background = element_rect(fill = NA, color = NA),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black", size=10), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
        geom_vline(aes_string(xintercept = group.by),linetype="dashed", size=0.2, colour="white")
}



# dimplot
#' This function takes in a seurat object and cell types of interest and returns a scatterpie plot with each spot
#' situated in its spatial location.
#'
#' @param data_ob: Object of class Seurat with the spatial data and the cell type proportions in the metadata.
#' @param group.by: misc list name for spotlight_results
#' @param split.by: Object of class numeric between 0-1 indicating the degree of transparency of the image.
#' @param reduct: Object of class character, name of the slice image to load as found in object@images, by default it
#' will grab the first one on the list.
#' @param ptsize: Object of class numeric between 0-1 indicating the degree of transparency of the scatterpie.
#' @param spointsize: Object of class numeric containing the size of the pie charts.
#' @param palette: object of class dataframe containing a color for each cell type, the first column is the cell type
#' and the second is the color assigned to it.
#' @param output_dir: object of class dataframe containing a color for each cell type, the first column is the cell type
#' @importFrom  Seurat ScaleFactors GetImage GetTissueCoordinates
#' @importFrom  scatterpie geom_scatterpie
#' @import cowplot
#' @export
#' @examples
seu_sub_dimplot <- function(
  data_ob,
  group.by,
  split.by = NULL,
  reduct,
  ptsize = 1,
  spointsize = 2,
  palette,
  resolution,
  label = FALSE,
  output_dir,
  crop = TRUE,
  common.legend = TRUE
) {
    cellmeta <- OESingleCell::colData(data_ob)
    if (!is.factor(data_ob@meta.data[, group.by])) {
        data_ob@meta.data[, group.by] <- factor(data_ob@meta.data[, group.by])
    }
    data_ob@meta.data[, split.by] <- factor(data_ob@meta.data[, split.by])
    Seurat::Idents(data_ob) <- group.by
    dim_1_range <- Seurat::Embeddings(data_ob, reduction = glue::glue("{reduct}"))[, 1] %>% range()
    dim_2_range <- Seurat::Embeddings(data_ob, reduction = glue::glue("{reduct}"))[, 2] %>% range()
    width_height_ratio <- (dim_1_range[2] - dim_1_range[1]) / (dim_2_range[2] - dim_2_range[1])
    width_height_ratio <- ifelse(width_height_ratio > 1.3, 1.3, width_height_ratio)
    reordered_cell_count_by_cluster <- table(Seurat::Idents(data_ob))
    cell_count_labels <- glue::glue("{names(reordered_cell_count_by_cluster)}") #-{reordered_cell_count_by_cluster} cells")

    if (is.null(split.by)) {
        futile.logger::flog.info("1.1.dimplot for all") # = ==============================================================
        colors2use <- OESingleCell::SelectColors(levels(Seurat::Idents(data_ob)), palette = palette)
        pdf(NULL)
        plot_dimplot <- Seurat::DimPlot(
            object = data_ob,
            reduction = reduct,
            label = label,
            pt.size = as.numeric(ptsize),
            group.by = group.by
        ) +
            # cols=colors2use)# +
            labs(fill = "") +
            theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
            # theme(plot.title = element_text(hjust = 0.5)) +
            scale_colour_manual(
                values = colors2use,
                breaks = levels(Seurat::Idents(data_ob)),
                labels = cell_count_labels
            ) & guides(colour = guide_legend(
          title = group.by,
          override.aes = list(size = 3),
          nrow = 15
          #ncol = ifelse(length(cell_count_labels) >= 10, ceiling(length(cell_count_labels) / 10), 1)
        ))
    } else {
      futile.logger::flog.info("1.1.dimplot for split")  #===========================================================

        colors2use <- OESingleCell::SelectColors(levels(Seurat::Idents(data_ob)), palette = palette)
        pdf(NULL)
        plot_dimplot <- Seurat::DimPlot(
            object = data_ob,
            reduction = reduct,
            label = label,
            ncol = length(unique(cellmeta[[split.by]])),
            pt.size = as.numeric(ptsize),
            group.by = group.by,
            split.by = split.by
        ) +
            # cols=colors2use)# +
            # theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
            labs(fill = "") +
            scale_colour_manual(
                values = colors2use,
                breaks = levels(Seurat::Idents(data_ob)),
                labels = levels(Seurat::Idents(data_ob))
            ) & guides(colour = guide_legend(
          title = group.by,
          override.aes = list(size = 3),
          nrow = 15
          #ncol = ifelse(length(cell_count_labels) >= 10, ceiling(length(cell_count_labels) / 10), 1)
        ))
    }

   futile.logger::flog.info("1.2.spatialdimplot")#===================================================================
    crop <- as.logical(crop)
    if (!is.null(Seurat::Images(data_ob))) {
        plot_gsp <- OESingleCell::SpatialPlot(
          data_ob,
          group.by = group.by,
          # alpha = c(0.1, 1),
          HE = FALSE,
          ncol = length(Seurat::Images(data_ob)),
          pt.size.factor = as.numeric(spointsize),
          #pt.alpha = FALSE,
          stroke = 0.2,
          combine = TRUE,
          cols = palette,
          crop = crop,
          common.legend = as.logical(common.legend)
        ) # &  theme(legend.text=element_text(size=10),legend.title=element_text(size=12))
    }
    futile.logger::flog.info("1.3.拼图")#============================================================================
    if (!is.null(Seurat::Images(data_ob))) {
        splitby_len <- ifelse(is.null(split.by), 1, length(unique(cellmeta[[split.by]])))
        image_len <- length(Seurat::Images(data_ob))
        max_colnum <- max(splitby_len, image_len)
        layout <- glue::glue("{strrep('B',splitby_len)}{strrep('#',max_colnum-splitby_len)}\n{strrep('C',image_len)}{strrep('#',max_colnum-image_len)}")
        plot_dimplot <- patchwork::wrap_plots(
            B = plot_dimplot,
            C = plot_gsp,
            # widths = c(splitby_len,image_len),
            heights = c(4.5, 4.5),
            # guides = "collect",
            design = layout
        )
        width <- 4.5 * max(splitby_len, image_len) +
            max(
                stringr::str_length(group.by),
                stringr::str_length(cell_count_labels) * ceiling(length(cell_count_labels) / 10)
            ) / 5
        height <- 4.5 * 2
    } else {
        width <- 4 * ifelse(width_height_ratio > 1.3, 1.3, width_height_ratio) +
            max(
                stringr::str_length(group.by),
                stringr::str_length(cell_count_labels) * ceiling(length(cell_count_labels) / 10)
            ) / 5
        #height <- ifelse(!is.null(split.by), 4 * length(unique(cellmeta[[split.by]])), 4 * 1)
        height <- 4
    }
    # redcut_name<-stringr::str_split(reduct,'_')[[1]][2]
    OESingleCell::save_ggplots(
        filename = ifelse(
            !is.null(split.by),
            glue::glue("{output_dir}/splitby-{split.by}_resolution{resolution}_split_plot"),
            glue::glue("{output_dir}/{reduct}_groupby-{group.by}_resolution{resolution}")
        ),
        plot = plot_dimplot,
        height = height,
        width = width,
        to.pdf = TRUE,
        to.png = TRUE,
        limitsize = FALSE
    )
    futile.logger::flog.info("2 contrast_plot")#  =================================================================
    if (!is.null(split.by)) {
        nlevel_col <- levels(data_ob@meta.data[,split.by])
        print(nlevel_col)
        groupvis <- Seurat::DimPlot(
            object = data_ob,
            dims = c(1, 2),
            reduction = reduct,
            pt.size = ptsize,
            group.by = split.by
        ) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::scale_colour_manual(values = OESingleCell::SelectColors(c(seq(1:10),nlevel_col), palette = palette)[11:(10+length(nlevel_col))])
        OESingleCell::save_ggplots(
            glue::glue("{output_dir}/groupby-{split.by}_resolution{resolution}_contrast_plot"),
            plot = groupvis,
            width = 8 * ifelse(width_height_ratio > 1.3, 1.3, width_height_ratio) + 2,
            height = 8,
            dpi = 300,
            bg = "white"
        )
    }
}

docstring <- " example1:\\n\\n\\
  scvis  -i  query.rds   -f rds  -o results  -d rds   --assay SCT  --dataslot data  clusterplot -g   --refcelltype  celltype --refassay  RNA "
sub_dimplot <- subparsers$add_parser("dimplot", help = "scatter plot with cell annotations")
sub_dimplot$add_argument(
    "-r",
    "--resolution",
    type = "double",
    default = 0.4,
    help = "resolution used in cell clustering.[default: %(default)s]"
)
sub_dimplot$add_argument(
    "-g",
    "--groupby",
    type = "character",
    default = "ident",
    help = paste0(
        "[OPTIONAL]The grouppinig variable in the metadata for groupping cells in the scatter plot.",
        "If not specified, %(default)s will be used as default."
    )
)
sub_dimplot$add_argument(
    "-y",
    "--splitby",
    type = "character",
    default = NULL,
    help = paste0(
        "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing",
        " the gene expression difference in different levels."
    )
)
# sub_dimplot$add_argument(
#     "--propby",
#     type = "character",
#     default = "group",
#     help = "calculate the proportion of cell in the specified groups.[default: %(default)s]"
# )
# sub_dimplot$add_argument(
#     "--groups",
#     type = "character",
#     default = "clusters",
#     help = "visualize cell metadata by coloring cells in different color according to cell clustering metadata.[default: %(default)s]"
# )
sub_dimplot$add_argument(
    "--palette",
    type = "character",
    default = "customecol2",
    help = paste0(
        "the discrete color schema mapped to the cell annotations specified by --groupby.[default: %(default)s]",
        " Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
        " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"
    )
)
sub_dimplot$add_argument(
    "-s",
    "--ptsize",
    type = "double",
    default = 0.01,
    help = "[OPTIONAL]the point size in the plot."
)
sub_dimplot$add_argument(
    "--spointsize",
    type = "double",
    default = 1.2,
    help = "[OPTIONAL]the point size in the spatialfeature plot."
)

sub_dimplot$add_argument(
    "--shuffle",
    type = "character",
    default = "FALSE",
    help = paste0(
        "[OPTIONAL]Whether to randomly shuffle the order of points. This can be useful for crowded plots ",
        "if points of interest are being buried.[default: %(default)s]."
    )
)
sub_dimplot$add_argument(
    "--label",
    type = "character",
    default = "FALSE",
    help = "Whether to label the clusters.[default: %(default)s]"
)
sub_dimplot$add_argument(
    "--plot",
    type = "character",
    default = "barplot",
    help = "the plot type used to summarize the cell clustering results.Options can be:barplot, pie, boxplot.[default: %(default)s]"
)
sub_dimplot$add_argument(
    "--highlight",
    type = "character",
    default = NULL,
    help = paste0(
        "the logical expression for selecting specific cells to be highlighted on the scatter plots.",
        "[default: %(default)s]"
    )
)
sub_dimplot$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)
sub_dimplot$add_argument(
    "--common.legend",
    default = 'TRUE',
    help = "whether to generate common legend for spatialplot.[default: %(default)s]"
)
args <- commandArgs(TRUE)
if ("dimplot" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "dimplot") {
        # ===================================================================================================================
        futile.logger::flog.info("step1:read the specified assay and data slot in data object into memory") # ===============
        # for visualization usually only "data" slot is loaded,which can save a lot of time and memory
        suppressMessages(data_ob <- OESingleCell::ReadX(
            input = opt$input,
            informat = opt$informat,
            assays = assays,
            data.use = dataslots, # "data","scale.data"
            reductions = opt$reduct,
            graphs = FALSE,
            images = FALSE, # no graph object needed here
            verbose = FALSE
        ))
        # ===================================================================================================================
        futile.logger::flog.info(glue::glue("step2:get subset of cells used col.name = col.name = or visualization if necessay:{opt$predicate}"))
        if (!is.null(opt$predicate)) {
            df <- OESingleCell::colData(data_ob)
            desired_cells <- subset(df, eval(parse(text = opt$predicate)))
            data_ob <- subset(data_ob, cells = rownames(desired_cells))
        }
        if (!is.null(Seurat::Images(data_ob))) {
            unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique())]
            if (length(unuse_images) > 0) {
                data_ob@images[unuse_images] <- NULL
            }
        }
        # ===================================================================================================================
        futile.logger::flog.info("step3: dimplot")
        group.bys <- unlist(strsplit(opt$groupby, ",", perl = T))
        if (!is.null(opt$splitby)) {
            split.bys <- unlist(strsplit(opt$splitby, ",", perl = T))
        } else {
            split.bys <- opt$splitby
        }
        ###
        for (group.by in group.bys) {
            output_dir <- glue::glue("{output_dir}/visualize_cluster_by_{group.by}")
            ## for总图
            futile.logger::flog.info(glue::glue("main  dimplot for {group.by} "))
            seu_sub_dimplot(
                data_ob = data_ob,
                group.by = group.by,
                reduct = opt$reduct,
                ptsize = opt$ptsize,
                spointsize = opt$spointsize,
                palette = opt$palette,
                label = FALSE,
                resolution = opt$resolution,
                output_dir = output_dir,
                crop = opt$crop,
                common.legend = opt$common.legend
            )
            ## for 拆分
            if (!is.null(split.bys)) {
                for (split.by in split.bys) {
                    futile.logger::flog.info(glue::glue("split  dimplot for {group.by} split by {split.by}"))
                    seu_sub_dimplot(
                        data_ob = data_ob,
                        group.by = group.by,
                        split.by = split.by,
                        reduct = opt$reduct,
                        ptsize = opt$ptsize,
                        spointsize = opt$spointsize,
                        palette = opt$palette,
                        label = FALSE,
                        resolution = opt$resolution,
                        output_dir = output_dir,
                        crop = opt$crop,
                        common.legend = opt$common.legend
                    )
                }
            }
        }
        # ===================================================================================================================
            futile.logger::flog.info("step4:  summarize the cell type composition in each sample and group") # SpatialColors(n = 100)====================
            clust_cond_freq <- as.data.frame(data_ob@meta.data[, c("sampleid", group.by)]) %>%
                dplyr::group_by(.dots = c("sampleid", group.by)) %>%
                dplyr::summarize(count = dplyr::n()) %>%
                dplyr::mutate(freq = (count / sum(count)) * 100)
            write.table(
                as.data.frame(clust_cond_freq),
                glue::glue("{output_dir}/clust_cond_freq_info.xls"),
                sep = "\t",
                col.names = T,
                row.names = F
            )
            pdf(NULL)
            if (!is.null(split.bys)) {
                for (split.by in split.bys) {
                    clust_sum_all <- OESingleCell::PlotAbundances(
                    data_ob,
                    prop.by = group.by,
                    group.by = split.by,
                    # split.by = 'sampleid',
                    method = opt$plot,
                    ncol = ifelse(opt$plot == "pie", 4, 1),
                    cols = unname(OESingleCell::SelectColors(1:length(unique(data_ob@meta.data[, group.by])), palette = opt$palette))
                    )
                    OESingleCell::save_ggplots(
                        glue::glue("{output_dir}/groupby-{split.by}_resolution-{opt$resolution}_summary_plot"),
                        plot = clust_sum_all,
                        width = ifelse(length(unique(data_ob@meta.data[, split.by]))<=1,3,length(unique(data_ob@meta.data[, split.by])) * 1+1),
                        height = 6,
                        dpi = 200
                    )
                }
            }
            Seurat::Idents(data_ob) <- "sampleid"
            nlevel_col <- levels(data_ob@meta.data[,"sampleid"])
            print(nlevel_col)
            clust_sum_all <- OESingleCell::PlotAbundances(
            data_ob,
            prop.by = "sampleid",
            group.by = group.by,
            # split.by = 'sampleid',
            method = opt$plot,
            ncol = ifelse(opt$plot == "pie", 4, 1),
            cols = unname(OESingleCell::SelectColors(c(seq(1:10),nlevel_col), palette = opt$palette)[11:(10+length(nlevel_col))])
            )
            OESingleCell::save_ggplots(
                glue::glue("{output_dir}/groupby-{group.by}_resolution-{opt$resolution}_summary_plot"),
                plot = clust_sum_all,
                width = ifelse(length(unique(data_ob@meta.data[, group.by]))<=1,3,length(unique(data_ob@meta.data[, group.by])) * 1),
                height = 6,
                dpi = 200
            )   
            futile.logger::flog.info("step5:  summarize the cell type composition in each sample and group by area plot") 
            clust_area_all = Plot_area(data_ob, group.by = group.by, 
                                        prop.by = "sampleid",
                                        cols= unname(OESingleCell::SelectColors(c(seq(1:10),nlevel_col), palette = opt$palette)[11:(10+length(nlevel_col))]))
            OESingleCell::save_ggplots(
                glue::glue("{output_dir}/groupby-{group.by}_resolution-{opt$resolution}_area_plot"),
                plot = clust_area_all,
                width = ifelse(length(unique(data_ob@meta.data[, group.by]))<=1,3,length(unique(data_ob@meta.data[, group.by])) * 1),
                height = 6,
                dpi = 200
            )  

            if (!is.null(split.bys)) {
                for (split.by in split.bys) {
                    clust_area_all = Plot_area(data_ob, group.by = split.by, 
                            prop.by = group.by,
                            cols= unname(OESingleCell::SelectColors(1:length(unique(data_ob@meta.data[, group.by])), palette = opt$palette))
                    )
                    OESingleCell::save_ggplots(
                        glue::glue("{output_dir}/groupby-{split.by}_resolution-{opt$resolution}_area_plot"),
                        plot = clust_area_all,
                        width = ifelse(length(unique(data_ob@meta.data[, split.by]))<=1,3,length(unique(data_ob@meta.data[, split.by])) * 1+1),
                        height = 6,
                        dpi = 200
                    )
                }
            }
            
        # ===================================================================================================================
        ## step5: output session information===============================================================================
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}
