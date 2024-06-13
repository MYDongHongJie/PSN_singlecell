sub_markervis <- subparsers$add_parser("markervis", help = "visualize for markers.")
sub_markervis$add_argument("-l","--markers", type = "character",
                            help = "the list file of marker genes to be visulized.")
sub_markervis$add_argument("--topn", "-n",
        type = "integer", default = 1,
        help = "the number of top markers for each cluster to visualizse.[default: %(default)s]"
    )
sub_markervis$add_argument("--HE",
        type = "character", default = "T",
        help = "切片背景"
    )
sub_markervis$add_argument("--topby", "-c",
        type = "character", default = "avg_log2FC",
        help = "the column used to pick top n marker gene to visulize.The option can be one of the column in the input marker genes table.[default: %(default)s]"
    )
sub_markervis$add_argument("--extraGene", "-x",
        type = "character", default = NULL,
        help = "The extra gene list of interest to visualize specified by the user."
    )
sub_markervis$add_argument("--diffGene", "-d",
        type = "character", default = NULL,
        help = "The diff gene list of interest to visualize specified by the user."
    )
sub_markervis$add_argument("--groupby", "-g",
        type = "character", default = "clusters",
        help = "The grouppinig variable in the metadata for separate the cells to visulize marker genes.[default: %(default)s]"
    )
sub_markervis$add_argument("--splitby", "-y",
        type = "character", default = NULL,
        help = "The variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels[default: %(default)s]"
    )

sub_markervis$add_argument("--pointsize", "-s",
        type = "double", default = 1,
        help = "the point size in the plot.[default: %(default)s]"
    )
sub_markervis$add_argument("--alpha2use", "-a",
        type = "double", default = 1,
        help = "the opacity of the pooints on the violin plot.[default: %(default)s]"
    )
sub_markervis$add_argument("--spatialpointsize", "-p",
        type = "double", default = 1.5,
        help = "the point size in the spatial  plot.[default: %(default)s]"
    )
sub_markervis$add_argument("--sample_ratio",
        type = "double", default = 0.6,
        help = "the ratio of random subsample for each group when drawing heatmap.[default: %(default)s]"
    )
sub_markervis$add_argument("--reduct",
        type = "character", default = "tsne",
        help = "the previous calculated reduction result used in the featureplot.[default: %(default)s]"
    )

sub_markervis$add_argument("--vismethod", "-m",
        type = "character", default = "vlnplot,featureplot",
        help = "the visulization methods for the marker genes of each cell cluster. he methods can be ridgeplot,vlnplot,dotplot,featureplot,heatmap,featurebygroup.[default: %(default)s]"
    )
sub_markervis$add_argument("--var2use", "-q",
        type = "character", default = NULL,
        help = "The column name in cell metadata used as identity of each cell combined with levels4var.[default: %(default)s]"
    )
sub_markervis$add_argument("--levels4var", "-u",
        type = "character", default = NULL,
        help = " subset of factor levels for the specified factor by --var2use.[default: %(default)s]"
    )
sub_markervis$add_argument(
    "--palette",
    type = "character",
    default = "customecol2",
    help = paste0(
        "the discrete color schema mapped to the cell annotations specified by --groupby.[default: %(default)s]",
        " Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,",
        " tableau20:20,Buen:17,UKBB:18,TF1:17,paired:12"
    )
)
sub_markervis$add_argument(
    "--plotlegend",
    type = "character",
    default = 'right',
    help = " can be 'right', 'left','none'.[default: %(default)s]"
)
sub_markervis$add_argument(
    "--crop",
    default = 'FALSE',
    help = "whether to crop in spatialplot, data from cytassist project should be 'TRUE'.[default: %(default)s]"
)


# =command line parameters setting=============================
args <- commandArgs(TRUE)
if ("markervis" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "markervis") {
    suppressWarnings({
        suppressPackageStartupMessages(library("Seurat"))
        suppressPackageStartupMessages(library("dplyr"))
        suppressPackageStartupMessages(library("RColorBrewer"))
        suppressPackageStartupMessages(library("gridExtra"))
        suppressPackageStartupMessages(library("ggplot2"))
        suppressPackageStartupMessages(library("tibble"))
        suppressPackageStartupMessages(library("cowplot"))
        suppressPackageStartupMessages(library("ggrastr"))
        suppressPackageStartupMessages(library("OESingleCell"))
    })
    # =================================================================================
    # parse the command line parameters
    # =================================================================================
    if (is.null(opt$input)) {
        stop("the seurat object is NOT AVAILABLE!")
    } else {
        seurat_ob <- readRDS(opt$input)
        if (seurat_ob@version < 3) {
            seurat_ob <- UpdateSeuratObject(seurat_ob) # make sure the seurat object match with the latest seurat package
        }

        # change the default assay for reduction if necessary
        if (!is.null(opt$assay)) {
            Seurat::DefaultAssay(seurat_ob) <- opt$assay
        } else {
            Seurat::DefaultAssay(seurat_ob) <- "RNA"
        }

        metadata <- seurat_ob@meta.data
        if (is.null(metadata$clusters)) {
            seurat_ob <- StashIdent(seurat_ob, save.name = "clusters")
        } else {
            # if it is the first time to run this script on a seurat object, the
            # clusters here is actually the sample index but cell cluster ids.
            # After running this script, the cluster id will be overwrite with
            # the actual cell cluster id.
            seurat_ob <- Seurat::SetIdent(seurat_ob, value = "clusters")
        }
    }

    # get the subset of cells used for visualization if necessay
    if (!is.null(opt$levels4var)) {
        if (is.null(opt$var2use)) {
            print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
            ident2use <- "clusters"
        } else {
            ident2use <- opt$var2use
        }
        cluster_list <- unlist(strsplit(opt$levels4var, ",", perl = T))
        subset_cell <- colnames(seurat_ob)[seurat_ob@meta.data[, ident2use] %in% cluster_list]
        seurat_ob <- subset(seurat_ob, cells = subset_cell)
        seurat_ob@meta.data[[ident2use]] <- factor(seurat_ob@meta.data[[ident2use]], levels = sort(unique(seurat_ob@meta.data[[ident2use]])))
    }

    # output directory setting
    if (is.null(opt$output)) {
        print("NO output directory specified,the current directory will be used!")
        root_dir <- getwd()
    } else {
        if (file.exists(opt$output)) {
            root_dir <- opt$output
        } else {
            root_dir <- opt$output
            dir.create(root_dir)
        }
    }

    if (is.null(opt$groupby)) {
        print("NO groupping variable AVAILABLE for cell groupping! The default cell clusters id will be used!")
        groupby <- "clusters"
    } else {
        groupby <- opt$groupby
    }
    print(groupby)

    if (!is.null(opt$splitby)) {
        splitby <- opt$splitby
        if (splitby == groupby) {
            stop("The variable specified by --splitby conflicted with the --groupby parameter, NULL will be used! ")
            splitby <- NULL
        }
    } else {
        splitby <- NULL
    }

    if (is.null(opt$markers) & is.null(opt$extraGene) & is.null(opt$diffGene)) {
        stop("NO marker genes is AVAILABLE!")
    }

    topn_markers <- data.frame()
    if (!is.null(opt$markers)) {
        markers2vis <- read.table(opt$markers, sep = "\t", header = T,quote="",check.names=F)
        p_name <- ifelse('p_val' %in% colnames(markers2vis), 'p_val', 'p-value')
        topn_markers <- markers2vis %>%
            dplyr::group_by(cluster) %>%
            dplyr::arrange(p_name, desc(avg_log2FC), desc(gene_diff)) %>%
            dplyr::top_n(opt$topn, .data[[opt$topby]]) %>%
            dplyr::arrange(cluster) %>%
            dplyr::mutate(folder_suffix = paste0("cluster", cluster)) %>%
            dplyr::select(cluster, gene, folder_suffix)
    }
    if (!is.null(opt$extraGene)) {
        extra_gene <- read.table(opt$extraGene, sep = "\t", header = T,quote="")
        if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene") {
            colnames(extra_gene)[1] <- "extra"
        }
        formated_extra_gene <- as.data.frame(tidyr::gather(extra_gene, key = "cluster", value = "GENE"))
        match <- CaseMatch(search = as.vector(formated_extra_gene$GENE), match = rownames(seurat_ob))
        filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
        if(length(filtered_gene)!=0){
        filtered_gene = as.data.frame(filtered_gene)
        colnames(filtered_gene) = "Gene"
        write.table(filtered_gene,file.path(root_dir,"genes_not_matched.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check genes_not_matched.xls for the genename.")
        }
        formated_extra_gene <- formated_extra_gene %>%
            dplyr::filter(GENE %in% names(match)) %>%
            dplyr::mutate(gene = match, folder_suffix = cluster) %>%
            dplyr::select(cluster, gene, folder_suffix)
        topn_markers <- rbind(topn_markers, formated_extra_gene)
    }

    if (!is.null(opt$diffGene)) {
        diffGene_path <- normalizePath(opt$diffGene)
        diffGene_name <- gsub(".*/", "", diffGene_path)
        if (grepl("-diff-", diffGene_name)) {
            name <- gsub("-diff-.*", "", diffGene_name)
        } else if (grepl("-all_", diffGene_name)) {
            name <- gsub("-all_.*", "", diffGene_name)
        }

        markers2vis <- read.delim(opt$diffGene, sep = "\t", header = T)
        up <- dplyr::filter(markers2vis, FoldChange > 1) %>%
            dplyr::arrange(desc(log2FoldChange)) %>%
            dplyr::top_n(opt$topn, log2FoldChange)
        down <-dplyr:: filter(markers2vis, FoldChange < 1) %>%
            dplyr::arrange(log2FoldChange) %>%
            dplyr::top_n(as.numeric(paste0("-", opt$topn)), log2FoldChange)
        topn_markers <- rbind(up, down)
        write.table(topn_markers, file.path(root_dir, paste0("top", opt$topn, "_", name, "_genes.xls", collapse = "")), quote = F, row.names = FALSE, sep = "\t")
    }


    if (is.null(opt$vismethod)) {
        print("NO marker gene visulization method provided,the default method vlnplot and featureplot will be used!")
        vismethods <- c("vlnplot", "featureplot")
    } else if (opt$vismethod == "all") {
        vismethods <- c("vlnplot", "featureplot", "ridgeplot", "dotplot", "heatmap", "violinEnsemble")
    } else {
        vismethods <- unlist(strsplit(opt$vismethod, ","))
    }

    crop <- as.logical(opt$crop)
    # =================================================================================
    # visualize the markers in different ways
    # =================================================================================
    for (vismethod in vismethods) {
        if (vismethod == "vlnplot") {
            # Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)

            for (clusterx in unique(topn_markers$folder_suffix)) {
                topn <- topn_markers %>%
                    dplyr::filter(folder_suffix == clusterx) %>%
                    dplyr::select(cluster, gene, folder_suffix)
                topn_markers2vis <- as.vector(topn$gene)

                path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
                if (file.exists(path4vis)) {
                    output_dir <- path4vis
                } else {
                    output_dir <- path4vis
                    dir.create(output_dir, recursive = TRUE)
                }

                colors2use <- OESingleCell::SelectColors(sort(unique(seurat_ob@meta.data[, groupby])),palette=opt$palette)
                gs <- lapply(topn_markers2vis, function(x) {
                    VlnPlot(seurat_ob,
                        features = x, cols = colors2use,
                        pt.size = opt$pointsize,
                        group.by = groupby, split.by = splitby
                    ) +
                        labs(title = "", y = x) +
                        theme(
                            legend.position = "none",
                            # panel.spacing = unit(.05, "lines"),
                            axis.title.x = element_text(size = 0),
                            axis.title.y = element_text(size = 12),
                            axis.text.y = element_text(size = 8)
                        )
                })
                i <- 1
                while (i <= length(gs)) {
                    gs_sub <- gs[i:(i + 49)]
                    gs_sub <- gs_sub[which(!sapply(gs_sub, is.null))]

                    pdf(file.path(output_dir, paste0("top", "marker_gene_violin_plot", ifelse(i == 1, "", i %/% 10 + 1), ".pdf", collapse = "_")),
                        width = 8, height = length(gs_sub)*2 + max(nchar(as.character(seurat_ob@meta.data[,groupby])))/2
                    )
                    grid.arrange(grobs = gs_sub, ncol = 1)
                    dev.off()
                    png(file.path(output_dir, paste0("top", "marker_gene_violin_plot", ifelse(i == 1, "", i %/% 10 + 1), ".png", collapse = "_")),
                        width = 8, height = length(gs_sub)*2 + max(nchar(as.character(seurat_ob@meta.data[,groupby])))/2, res = 96, units = "in"
                    )
                    grid.arrange(grobs = gs_sub, ncol = 1)
                    dev.off()
                    i <- i + 50
                }
            }
        }

        if (vismethod == "featureplot") {
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
            sc <- scale_colour_gradientn(colours = myPalette(100))
            suppressMessages({
                #lapply(unique(topn_markers$folder_suffix), function(clusterx) {
                for (clusterx in unique(topn_markers$folder_suffix)) {
                    if (!opt$reduct %in% names(Key(seurat_ob))) {
                        stop("NO specified reduction found in the object!")
                    } else {
                        reduct <- opt$reduct
                    }
                    topn <- topn_markers %>%
                        dplyr::filter(folder_suffix == clusterx) %>%
                        dplyr::select(cluster, gene, folder_suffix)
                    topn_markers2vis <- as.vector(topn$gene)

                    path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
                    if (file.exists(path4vis)) {
                        output_dir <- path4vis
                    } else {
                        output_dir <- path4vis
                        dir.create(output_dir, recursive = TRUE)
                    }
                      parallel::mclapply(topn_markers2vis,function(x) {
                          #lapply(topn_markers2vis, function(x) {
                          ggfeature <- Seurat::FeaturePlot(seurat_ob,
                                                           features = x,
                                                           cols = c("grey", "red"),
                                                           split.by = splitby,
                                                           reduction = reduct,
                                                           ncol = 1,
                                                           min.cutoff = NA,#最小值，
                                                           max.cutoff = NA, #最大值
                                                           pt.size = opt$pointsize) +
                            theme(plot.title = element_text(hjust = 0.5)) +
                            sc

                          if (opt$image == TRUE) {
                              if ( opt$HE == FALSE){
                                ggfeature_spatial <- OESingleCell::SpatialPlot(
                                    seurat_ob,
                                    assay = opt$assay,
                                    features = x,
                                    alpha = opt$alpha2use,
                                    pt.size.factor = opt$spatialpointsize,
                                    stroke = 0.1,
                                    min.cutoff = NA,
                                    max.cutoff = NA,
                                    combine = FALSE,
                                    HE = F,
                                    cols = "spectral",
                                    ncol = length(Images(seurat_ob)),
                                    crop = crop
                                )}else{
                                    ggfeature_spatial <- OESingleCell::SpatialPlot(
                                    seurat_ob,
                                    assay = opt$assay,
                                    features = x,
                                    alpha = opt$alpha2use,
                                    pt.size.factor = opt$spatialpointsize,
                                    stroke = 0.1,
                                    min.cutoff = NA,
                                    max.cutoff = NA,
                                    combine = FALSE,
                                    HE = T,
                                    cols = "spectral",
                                    ncol = length(Images(seurat_ob)),
                                    crop = crop
                                    )}

                              # SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
                              # gene_range <- range(FetchData(object = seurat_ob,slot="data",vars = x))
                              # min.cutoff <- gene_range[1]
                              # max.cutoff <- gene_range[2]
                              # ggfeature_spatial <- lapply(ggfeature_spatial, function(x) {
                              #     x + scale_fill_gradientn(
                              #             colors = SpatialColors(n = 100),
                              #             limits = c(min.cutoff, max.cutoff)
                              #         )+ theme(legend.title = element_text(face = "bold", size = 16))
                              # })

                              nrow <- ifelse(length(ggfeature_spatial) >= 4, ceiling(length(ggfeature_spatial) / 4), 1)
                              ncol <- ifelse(length(ggfeature_spatial) >= 4, 4, length(ggfeature_spatial))
                              ggfeature_spatial_all <- do.call(
                                ggpubr::ggarrange,
                                c(
                                  ggfeature_spatial,
                                  list(
                                    nrow = nrow,
                                    ncol = ncol,
                                        common.legend = TRUE,
                                        legend = opt$plotlegend,
                                         hjust = 0.5,
                                        align = "hv"
                                    )
                                )
                            )
                            plot <- plot_grid(ggfeature,
                                ggfeature_spatial_all,
                                rel_widths = c(1, ifelse(length(ggfeature_spatial)==1,1,2+0.5))
                            )
                             ggsave(file.path(output_dir, paste(clusterx, x, "marker_gene_Spatialfeatureplot.pdf", sep = "_")),
                            plot = plot,
                            limitsize = FALSE,
                            width = 4 * ncol + 4 * nrow+1,
                            height = 4 * nrow,bg='white')
                        ggsave(file.path(output_dir, paste(clusterx, x, "marker_gene_Spatialfeatureplot.png", sep = "_")),
                            plot = plot,
                            limitsize = FALSE,
                            width = 4 * ncol + 4 * nrow+1,
                            height = 4 * nrow,bg='white' )
                        } else {
                            plot <- ggfeature
                            OESingleCell::save_ggplots(
                              filename = file.path(output_dir, paste(clusterx, x, "marker_gene_featureplot", sep = "_")),
                              plot = plot,
                              height = 4,
                              width = 5,
                              to.pdf = TRUE,
                              to.png = TRUE,
                              limitsize = FALSE
                            )
                        }
                    } , mc.cores = opt$ncores)
                }
                #})
            })
        }

        if (vismethod == "dotplot") {
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
            sc <- scale_colour_gradientn(colours = myPalette(100))
            for (clusterx in unique(topn_markers$folder_suffix)) {
                topn <- topn_markers %>%
                    dplyr::filter(folder_suffix == clusterx) %>%
                    dplyr::select(cluster, gene, folder_suffix)
                topn_markers2vis <- as.vector(topn$gene)
                path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
                if (file.exists(path4vis)) {
                    output_dir <- path4vis
                } else {
                    output_dir <- path4vis
                    dir.create(output_dir)
                }
                seurat_ob <- SetIdent(seurat_ob, value = groupby)
                dot_plot <- DotPlot(object = seurat_ob, features = topn_markers2vis) + RotatedAxis() + sc
                dot_plot + guides(color = guide_colorbar(order = 1, title = "Average Expression"))
                ggsave(file.path(output_dir, paste0("top", "marker_gene_dotplot.pdf", collapse = "_")), width = (0.3 * length(topn_markers2vis) + 2.5 + max(nchar(names(table(seurat_ob@meta.data[, groupby])))) / 10))
                ggsave(file.path(output_dir, paste0("top", "marker_gene_dotplot.png", collapse = "_")), dpi = 300, limitsize = F, width = (0.3 * length(topn_markers2vis) + 2.5 + max(nchar(names(table(seurat_ob@meta.data[, groupby])))) / 10),bg='white')
            }
        }

        if (vismethod == "ridgeplot") {
            for (clusterx in unique(topn_markers$folder_suffix)) {
                topn <- topn_markers %>%
                    dplyr::filter(folder_suffix == clusterx) %>%
                    dplyr::select(cluster, gene, folder_suffix)
                topn_markers2vis <- as.vector(topn$gene)
                colors2use <- OESingleCell::SelectColors(sort(unique(seurat_ob@meta.data[, groupby])),palette=opt$palette)

                path4vis <- file.path(root_dir, paste0("markers_vis4", clusterx, collapse = ""))
                if (file.exists(path4vis)) {
                    output_dir <- path4vis
                } else {
                    output_dir <- path4vis
                    dir.create(output_dir)
                }
                seurat_ob <- SetIdent(seurat_ob, value = groupby)
                nrow <- ceiling(length(topn_markers2vis) / 2)
                RidgePlot(
                    object = seurat_ob, features = topn_markers2vis,
                    cols = colors2use, ncol = 2
                )
                ggsave(file.path(output_dir, paste0("top", "marker_gene_ridgeplot.pdf", collapse = "_")), width = ifelse(length(topn_markers2vis) == 1, yes = 6, no = 12), height = nrow * 5)
                ggsave(file.path(output_dir, paste0("top", "marker_gene_ridgeplot.png", collapse = "_")), dpi = 300, limitsize = F, width = ifelse(length(topn_markers2vis) == 1, yes = 6, no = 12), height = nrow * 5)
            }
        }

        if (vismethod == "violinensemble") {
            colors2use <- OESingleCell::SelectColors(sort(unique(seurat_ob@meta.data[, groupby])),palette=opt$palette)
            ggensemble <- ViolinEnsemble(
                object = seurat_ob, features = topn_markers2vis, cols = colors2use,
                group.by = groupby, show_point = T, pt_size = opt$pointsize
            )
            ggsave(file.path(output_dir, paste0("all_top", "marker_gene_enhanced_violin_plot.pdf", collapse = "_")))
            ggsave(file.path(output_dir, paste0("all_top", "marker_gene_enhanced_ciolin_plot.png", collapse = "_")),
                dpi = 300, limitsize = F
            )
            # for ( clusterx in unique(topn_markers$cluster) ){
            #     topn = topn_markers %>% filter( cluster == clusterx) %>% select(cluster,gene)
            #     topn_markers2vis = as.vector(topn$gene)
            #
            #     path4vis = file.path(root_dir,paste0("markers_vis4cluster",clusterx,collapse = ""))
            #     if ( file.exists( path4vis ) ){
            #         output_dir = path4vis
            #     }else{
            #         output_dir = path4vis
            #         dir.create(output_dir)
            #     }
            #     ggensemble = violinEnsemble( object =  seurat_ob, genes = topn_markers2vis,
            #                                 group_var = groupby,
            #                                 slot_use = slot2use,assay_use = "RNA")
            #     ggsave(file.path(output_dir,paste0( "top", "marker_gene_enhanced_violin_plot.pdf",collapse = "_" )))
            #     ggsave(file.path(output_dir,paste0( "top", "marker_gene_enhanced_ciolin_plot.png", collapse = "_")),
            #                     dpi = 300 ,limitsize = F)
            # }
        }

        if (vismethod == "heatmap") {
            markers2vis4heatmap <- unique(as.vector(topn_markers$gene))
            seurat_ob <- seurat_ob %>% Seurat::ScaleData(assay = "SCT", features = markers2vis4heatmap)
            # if (is.null(opt$sample_ratio)) {
            # subseted_seurat <- seurat_ob
            # } else {
            #     sampled_cellmeta <- seurat_ob@meta.data %>%
            #         rownames_to_column() %>%
            #         group_by(.dots = groupby) %>%
            #         sample_frac(size = opt$sample_ratio, replace = F) %>%
            #         column_to_rownames()
            #     subseted_seurat <- subset(seurat_ob, cells = rownames(sampled_cellmeta))
            # }
            seurat_ob@meta.data[, groupby] <- factor(seurat_ob@meta.data[, groupby],
                                                           levels = sort(unique(seurat_ob@meta.data[, groupby])))
            print(unique(sort(seurat_ob@meta.data[, groupby])))
            if (length(markers2vis4heatmap) > 135) {
                sz <- 4 - log(length(markers2vis4heatmap) / 100)
                heig <- 5 + log2(length(markers2vis4heatmap) / 10)
                wid <- 7.5
            }else if (length(markers2vis4heatmap) > 75) {
                sz <- 4 - log2(length(markers2vis4heatmap) / 120); heig <- 7; wid <- 7
            }else {
                sz <- 6 - log2(length(markers2vis4heatmap) / 80); heig <- 7; wid <- 7
            }
            #colors2use <- OESingleCell::SelectColors(1:length(unique(seurat_ob@meta.data[, groupby])),palette=opt$palette)
            colors2use <- OESingleCell::SelectColors(levels(seurat_ob@meta.data[, groupby]), palette = opt$palette)

            ggheat <- Seurat::DoHeatmap(
              object = seurat_ob,
              features = markers2vis4heatmap,
              group.colors = colors2use,
              group.by = groupby,
              group.bar = T,
              label = F,
              draw.lines = T
            ) +
              #  theme(axis.text.y = element_text(size = 4))
              # group.cex = 10, cex.row = 4,
            # slim.col.label = T, group.label.rot = F)
            #ggheat + guides(fill = guide_colorbar(title. position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
            #ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap.pdf", collapse = "_")))
            #ggsave(file.path(root_dir, paste0("top", "marker_gene_heatmap.png", collapse = "_")), dpi = 300, limitsize = F)
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = sz,face = "bold"))+
                ggplot2::guides(fill = ggplot2::guide_colorbar( title.position = "top", order = 1)  ) +#, color = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 1)))+
                ggplot2::scale_fill_gradientn(colors = c("#406AA8", "white", "#D91216"), na.value ="white")
            OESingleCell::save_ggplots(
                        file.path(opt$output, paste0("top","marker_gene_heatmap", collapse = "")) ,
                        plot = ggheat,
                        dpi = 300 ,
                        limitsize = F,
                        width = heig,
                        height = wid)

        }

        if (vismethod == "diff_heatmap") {
            markers2vis4heatmap <- unique(as.vector(topn_markers$gene))
            if (is.null(opt$sample_ratio)) {
                subseted_seurat <- seurat_ob
            } else {
                sampled_cellmeta <- seurat_ob@meta.data %>%
                    rownames_to_column() %>%
                    group_by(.dots = opt$groupby) %>%
                    sample_frac(size = opt$sample_ratio, replace = F) %>%
                    column_to_rownames()
                subseted_seurat <- subset(seurat_ob, cells = rownames(sampled_cellmeta))
            }
            subseted_seurat@meta.data[, opt$groupby] <- as.vector(subseted_seurat@meta.data[, opt$groupby])
            colors2use <- OESingleCell::SelectColors(unique(subseted_seurat@meta.data[, opt$groupby]),palette=opt$palette)
            print(opt$groupby)
            ggheat <- DoHeatmap(
                object = subseted_seurat,
                features = markers2vis4heatmap,
                group.colors = colors2use,
                group.by = opt$groupby, group.bar = T, label = F
            ) +
                theme(axis.text.y = element_text(size = 6))
            # group.cex = 10, cex.row = 4,
            # slim.col.label = T, group.label.rot = F)
            ggheat + guides(fill = guide_colorbar(title.position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
            ggsave(file.path(root_dir, paste0("top", opt$topn, "_", name, "_heatmap.pdf", collapse = "")))
            ggsave(file.path(root_dir, paste0("top", opt$topn, "_", name, "_heatmap.png", collapse = "")), dpi = 300, limitsize = F)
        }

        if (vismethod == "geneset") {
            path4vis <- file.path(root_dir, "geneset_visualization")
            if (file.exists(path4vis)) {
                output_dir <- path4vis
            } else {
                output_dir <- path4vis
                dir.create(output_dir)
            }
            seurat_ob <- SetIdent(seurat_ob, value = groupby)
            topn_markers2vis <- list()
            for (clusterx in unique(topn_markers$folder_suffix)) {
                topn_markers2vis[[clusterx]] <- subset(topn_markers, folder_suffix == clusterx)$gene
            }
            seurat_ob <- AddModuleScore(seurat_ob, features = topn_markers2vis, name = names(topn_markers2vis))
            colnames(seurat_ob@meta.data)[(dim(seurat_ob[[]])[2] - length(topn_markers2vis) + 1):dim(seurat_ob[[]])[2]] <- names(topn_markers2vis)
            metadata <- seurat_ob@meta.data
            metadata <- metadata[, c(groupby, names(topn_markers2vis))]
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
            sc <- scale_colour_gradientn(colours = myPalette(100))
            # vln=list();feat=list()
            for (geneset in names(topn_markers2vis)) {
                ## featureplot
                if (!opt$reduct %in% names(Key(seurat_ob))) {
                    stop("NO specified reduction found in the object!")
                } else {
                    reduct <- opt$reduct
                }

                ####
                ggfeature <- FeaturePlot(seurat_ob,
                                         features = geneset, cols = c("grey", "red"),
                                         split.by = splitby, reduction = reduct,
                                         ncol = 2, pt.size = opt$pointsize
                ) +
                  theme(plot.title = element_text(hjust = 0.5)) +
                  sc
                ggfeature_spatial <- OESingleCell::SpatialPlot(seurat_ob,
                                                               features = geneset,
                                                               alpha = opt$alpha2use,
                                                               pt.size.factor = opt$spatialpointsize,
                                                               combine = FALSE,
                                                               ncol = length(Images(seurat_ob)),
                                                               crop = crop
                )
                # SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
                # gene_range <- range(FetchData(object = seurat_ob, vars = geneset))
                # min.cutoff <- gene_range[1]
                # max.cutoff <- gene_range[2]
                # ggfeature_spatial <- lapply(ggfeature_spatial, function(x) {
                #     x
                #     + scale_fill_gradientn(
                #             colors = SpatialColors(n = 100),
                #             limits = c(min.cutoff, max.cutoff)
                #         )
                # })
                # + theme(legend.title = element_text(face = "bold", size = 16)))
                nrow <- ifelse(length(ggfeature_spatial) >= 4, ceiling(length(ggfeature_spatial) / 4), 1)
                ncol <- ifelse(length(ggfeature_spatial) >= 4, 4, length(ggfeature_spatial))
                ggfeature_spatial_all <- do.call(
                  ggpubr::ggarrange,
                  c(
                    ggfeature_spatial,
                    list(
                      nrow = nrow,
                      ncol = ncol,
                            common.legend = TRUE,
                            legend = 'right',
                            hjust = 0.5,
                            align = "hv"
                        )
                    )
                )
                plot <- plot_grid(ggfeature,
                    ggfeature_spatial_all,
                    rel_widths = c(1, ncol / nrow)
                )

                ggsave(file.path(output_dir, paste(geneset, "score_Spatialfeatureplot.pdf", sep = "_")),
                    plot = plot,
                    limitsize = FALSE,
                    width = 4 * ncol + 4 * nrow,
                    height = 4 * nrow
                )
                ggsave(file.path(output_dir, paste(geneset, "score_Spatialfeatureplot.png", sep = "_")),
                    plot = plot,
                    limitsize = FALSE,
                    width = 4 * ncol + 4 * nrow,
                    height = 4 * nrow
                )
                # feat[[geneset]]=FeaturePlot(seurat_ob,features = geneset,
                #                             cols = c("grey","red"),split.by = NULL, reduction = reduct, ncol = 2, pt.size = 0.4, order = T) +
                #   theme( plot.title = element_text(hjust = 0.5)) + sc
                # ggsave(file.path(output_dir,paste0(geneset,"_score_featureplot.pdf")),plot=feat[[geneset]])
                # ggsave(file.path(output_dir,paste0(geneset,"_score_featureplot.png")),plot=feat[[geneset]])
            }
        }

        if (vismethod == "velocity") {
            markers2vis4velo <- unique(as.vector(topn_markers$gene))
            emb <- Embeddings(seurat_ob, reduction = reduct)
            cell.dist <- as.dist(1 - armaCor(t(emb)))
            emat <- seurat_ob[["spliced"]]
            nmat <- seurat_ob[["unspliced"]]
            colors2use <- OESingleCell::SelectColors(sort(unique(seurat_ob@meta.data[, groupby])), palette = opt$palette)
            rvel.cd <- Tool(object = seurat_ob, slot = "RunVelocity")
            gene.relative.velocity.estimates(emat, nmat,
                                             deltaT = 1,
                                             kCells = 25, kGenes = 1, fit.quantile = 0.2,
                                             cell.emb = emb, cell.colors = colors2use, cell.dist = cell.dist,
                                             show.gene = gene, old.fit = rvel.cd, do.par = T
            )
        }
    }
        # ==============================================================================================================
        ## save session information
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}


