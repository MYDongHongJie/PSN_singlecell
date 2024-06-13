sub_cellcycle <- subparsers$add_parser("cellcycle", help = "infer the cell cycle.")
sub_cellcycle$add_argument("--splitby", "-b", type = "character", default = "sampleid",
                           help = "[OPTIONAL] The cell type annotation column name to use in seurat metadata.")
sub_cellcycle$add_argument("--species", "-s", type = "character",
                           help = "the species of sample")
sub_cellcycle$add_argument("--bartype", "-t", type = "character", default = "stack",
                           help = "[OPTIONAL]barplot type is stack or dodge.[default: %(default)s]")
sub_cellcycle$add_argument("--method", "-m", type = "character", default = 'seurat_3ph',
                           choices=c("seurat_3ph","dropbead","scran"),
                           help = "the method use to calculate cellcyle,can be 'seurat_3ph','dropbead','scran'.[default:%(default)s]")
sub_cellcycle$add_argument("--ident2use", "-q", type = "character", default = NULL,
                           help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell.[default: %(default)s]")
sub_cellcycle$add_argument("--which_cells", "-u", type = "character", default = NULL,
                           help = "[OPTIONAL]The subset of levels used for subtyping.[default: %(default)s]")
sub_cellcycle$add_argument("--proportion", "-v", type = "double", default = NULL,
                           help = "the proportion of cells used to subsample from each cluster.[default: %(default)s]")
sub_cellcycle$add_argument("--pt.sze", type = "character", default = 1,
                           help = "[OPTIONAL]the point size setting in the dimension reduction plot.[default: %(default)s]")
sub_cellcycle$add_argument("--cycle", "-k", type = "character",
                           help = "the cell cycle gene markers for human and mouse")
sub_cellcycle$add_argument("-c","--colorschema", type="character", default="customecol2",
                          help = "[OPTIONAL]the supported color schema currently:seurat,blindless,paried,colx22,jet,tableau20,tableau10medium,colorblind10,cyclic.[default: %(default)s]")
sub_cellcycle$add_argument("--reduct", type = "character", default = "tsne",
                           help = "the dimension reduction result used to visualize the cells.[default: %(default)s]")
sub_cellcycle$add_argument("--groupby", "-g", type = "character", default = "clusters",
                           help = "groupping cells on heatmap.[default: %(default)s]")

args <- commandArgs(TRUE)
if ('cellcycle' %in% args) {
  opt <- intial_setting()
  if(opt$sub_name == "cellcycle" ){
    cellcycle_name <- paste0(method, "_CellCycle")
    if (is.null(opt$species)) { stop("NO SPECIE SPECIFIED!") }else { spe <- tolower(opt$species) }
    if (!is.null(opt$splitby) & opt$splitby != "NULL") {
      facetby <- unlist(strsplit(opt$splitby, ",", perl = T))
    }else {
      facetby <- NULL
    }
    cellfeature2vis <- unlist(strsplit(opt$groupby, ",", perl = T))
    #===================================================================================================================
    #read in the 10X data in different format:  whatever format the input file has, all the input expression matrix
    #should be transformed to the seurat object finally.
    if (!is.null(opt$input)) {
      # the input is a seurat object which may contain more than one sample
        seurat_ob <- OESingleCell::ReadX(input = opt$input,
                                         informat = opt$informat,
                                         assays = assays,
                                         data.use = dataslots,
                                         verbose = F)
        findcluster_record <- Seurat::Command(seurat_ob, command = "FindClusters")
        resolution <- findcluster_record$resolution

        if (is.null(seurat_ob@meta.data$clusters)) {
          seurat_ob <- Seurat::StashIdent(seurat_ob, save.name = "clusters")
        }else {
          # if it is the first time to run this script on a seurat object, the clusters here is actually the sample index
          # but cell cluster ids. After running this script, the cluster id will be overwrite with the actual cell cluster id.
          seurat_ob <- Seurat::SetIdent(seurat_ob, value = "clusters")
        }

        #get the subset of cells used for visualization
        if (!is.null(opt$which_cells)) {
          cluster_list <- unlist(strsplit(opt$which_cells, ",", perl = T))
          seurat_ob <- subset(seurat_ob,
                              cells = OldWhichCells(seurat_ob,
                                                    subset.name = opt$ident2use,
                                                    accept.value = cluster_list))
        }
    }



    #===================================================================================================================
    ## Cell cycle phase with Dropbead:point towards cell cycle gene list supplied with Dropbead
    runcycle <- switch(method,
                       "dropbead" = {
                           suppressPackageStartupMessages(library("dropbead"))
                           if (is.null(opt$cycle)) {
                             stop("NO cell cycle marker gene file is AVAILABLE!")
                           }else {
                             cell_cycle_marker <- opt$cycle
                           }
                           # random subsample the cells in each cluster for heatmap visualization in case of big data
                           if (!is.null(opt$proportion)) {
                             sampled_cells <- seurat_ob@meta.data %>%
                                              dplyr::mutate(barcode = rownames(seurat_ob@meta.data)) %>%
                                              dplyr::group_by(clusters) %>%
                                              dplyr::sample_frac(size = opt$proportion, replace = F)
                             seurat_subset <- subset(seurat_ob, cells = sampled_cells$barcode)

                             ## Create a new single species object
                             SingleSS_ob <- new("SingleSpeciesSample",
                                                species1 = spe, ## Change to mouse if working with mouse data!
                                                cells = rownames(seurat_subset@meta.data),
                                                genes = rownames(Seurat::GetAssayData(seurat_subset, slot = "data")),
                                                dge = as.data.frame(as.matrix(Seurat::GetAssayData(seurat_subset, slot = "data"))))
                           }else {
                             ## Create a new single species object
                             SingleSS_ob <- new("SingleSpeciesSample",
                                                species1 = spe, ## Change to mouse if working with mouse data!
                                                cells = rownames(seurat_ob@meta.data),
                                                genes = rownames(Seurat::GetAssayData(seurat_ob, slot = "data")),
                                                dge = as.data.frame(as.matrix(Seurat::GetAssayData(seurat_ob, slot = "data"))))
                           }

                           phases <- assignCellCyclePhases(SingleSS_ob, gene.file = cell_cycle_marker, do.plot = F)
                           if (opt$groupby %in% names(seurat_ob@meta.data)) {
                             # order1 get original annotation
                             cells <- SingleSS_ob@cells[as.integer(rownames(phases))]
                             anno <- seurat_ob@meta.data[cells, opt$groupby]
                             # reorder
                             phases_reorder <- phases[order(anno),]
                             cells <- SingleSS_ob@cells[as.integer(rownames(phases_reorder))]
                             anno <- as.character(seurat_ob@meta.data[cells, opt$groupby])

                             annotation_col <- data.frame(clusters = anno, row.names = rownames(phases_reorder))
                             col <- OESingleCell::SelectColors(palette = opt$colorschema, n= length(unique(seurat_ob@meta.data[, opt$groupby])))
                             #names(col)=as.character(1:length(unique(seurat_ob@meta.data[,opt$groupby])))
                             names(col) <- unique(annotation_col$clusters)
                             annotation_colors <- list(clusters = col)
                             names(annotation_colors) <- opt$groupby
                             names(annotation_col) <- opt$groupby
                             # visualize on heatmap
                             ggheat<- pheatmap::pheatmap(as.matrix(t(phases_reorder)),
                                                cluster_row = F,
                                                cluster_col = F,
                                                fontsize = 20,
                                                color = colorRampPalette(c("#3794bf", "#FFFFFF", "#cc4140"))(100),
                                                show_colnames = F,
                                                annotation_col = annotation_col,
                                                annotation_colors = annotation_colors) %>%
                                       grid::grid.grabExpr() %>%
                                       ggplotify::as.ggplot()
                             OESingleCell::save_ggplots(glue::glue("cell_cycle_annotation_heatmap_orderedby_{opt$groupby}"),
                                                        width = 20,
                                                        height = 20)
                           } else {
                             warning("groupby not found in metadata")
                           }
                           pdf(file.path(output_dir, paste0("cell_cycle_annotation_heatmap.pdf")), width = 20, height = 20)
                           pheatmap::pheatmap(as.matrix(t(phases)),
                                              treeheight_row = 0,
                                              treeheight_col = 0,
                                              fontsize = 20,
                                              color = colorRampPalette(c("#3794bf", "#FFFFFF", "#cc4140"))(100),
                                              show_colnames = F)
                           dev.off()
                           sorted_phases <- phases %>% arrange(rownames(phases))
                           ## For each cell, check which phase has highest score and assign this to the cell
                           ## Code by akrun (https://stackoverflow.com/questions/36274867/getting-column-names-for-max-value-in-each-row)
                           j1 <- max.col(sorted_phases, "first")
                           value <- sorted_phases[cbind(1:nrow(sorted_phases), j1)]
                           cluster <- names(sorted_phases)[j1]
                           res <- data.frame(CellCycle = cluster)
                           res <- cbind(sorted_phases, res)
                           colnames(res) <- paste0(method, "_", colnames(res))
                           rownames(res) <- rownames(seurat_ob@meta.data)
                           ## Add the cell cycle phase to the Seurat object as Metadata
                           seurat_ob <- Seurat::AddMetaData(seurat_ob, metadata = res, col.name = colnames(res))
                           seurat_ob@meta.data[, cellcycle_name] <- factor(seurat_ob@meta.data[, cellcycle_name],
                                                                           levels = colnames(phases))
                       },
                       "scran" <- {
                           if (spe %in% c("human", "mouse")) {
                             ref.pairs <- readRDS(system.file("exdata", paste0(spe, "_cycle_markers.rds"), package = "scran"))
                             a <- read.table(paste0("/home/dongjiaoyang/database/genename2id/", spe, ".txt"), header = F)
                             a <- tibble::column_to_rownames(as.data.frame(a), "V1")
                             a[, 1] <- as.character(a[, 1])
                           } else stop("species not supported.")
                           ref <- list(
                             G1 = data.frame(first = a[ref.pairs$G1[, 1],], second = a[ref.pairs$G1[, 2],], stringsAsFactors = F),
                             S = data.frame(first = a[ref.pairs$S[, 1],], second = a[ref.pairs$S[, 2],], stringsAsFactors = F),
                             G2M = data.frame(first = a[ref.pairs$G2M[, 1],], second = a[ref.pairs$G2M[, 2],], stringsAsFactors = F)
                           )
                           #=================================================================================
                           # Classifying a new dataset:
                           #=================================================================================
                           if (length(intersect(rownames(seurat_ob), ref[[1]]$first)) == 0) stop("Gene names not matched. please check the species.")
                           sce_ob <- Seurat::as.SingleCellExperiment(seurat_ob) #to sce
                           assignments <- scran::cyclone(sce_ob, ref)
                           scores <- assignments$normalized.scores %>%
                             dplyr::mutate(barcodes = rownames(seurat_ob@meta.data), phase = assignments$phases) %>%
                             dplyr::select(barcodes, everything()) %>%
                             tibble::column_to_rownames("barcodes")
                           colnames(scores) <- paste0(method, "_", colnames(scores))
                           colnames(scores)[4] <- cellcycle_name
                           seurat_ob <- Seurat::AddMetaData(seurat_ob, metadata = scores, col.name = colnames(scores))
                       },
                       "seurat_3ph" = {

                           metadata <- seurat_ob@meta.data
                           if (!"Phase" %in% colnames(seurat_ob@meta.data)) {
                             genes.inuse <- rownames(Seurat::GetAssayData(seurat_ob, slot = "counts"))
                             s.genes <- Seurat::CaseMatch(search = Seurat::cc.genes$s.genes, match = genes.inuse)
                             g2m.genes <- Seurat::CaseMatch(search = Seurat::cc.genes$g2m.genes, match = genes.inuse)
                             seurat_ob <- Seurat::CellCycleScoring(object = seurat_ob,
                                                                   s.features = s.genes,
                                                                   g2m.features = g2m.genes,
                                                                   set.ident = F)
                           }
                           seurat_ob <- Seurat::SetIdent(seurat_ob, value = "Phase")
                           nlevel <- length(unique(seurat_ob@meta.data[, "Phase"]))
                           ggtsne <- Seurat::DimPlot(object = seurat_ob,
                                                     reduction = opt$reduct,
                                                     pt.size = ifelse(is.null(opt$pt.sze),0.5 ,as.numeric(opt$pt.sze))) +
                                     theme(plot.title = element_text(hjust = 0.5)) +
                                     scale_colour_manual(values = OESingleCell::SelectColors(palette = opt$colorschema, n= nlevel))
                           OESingleCell::save_ggplots(filename =file.path(output_dir, "Cellcycle.png"),
                                                      plot=ggtsne,
                                                      dpi = 1000,
                                                      limitsize = F)

                           cellcycle_results <- seurat_ob@meta.data %>%
                                                dplyr::rename("Barcode" = "orig.ident") %>%
                                                dplyr::select(Barcode, sampleid, clusters, group, S.Score, G2M.Score, Phase)
                           write.table(cellcycle_results,
                                       quote = F,
                                       file.path(output_dir, "Cellcycle_results.xls"),
                                       sep = "\t",
                                       row.names = F)
                       }
    )

    #===================================================================================================================
    # visualizataion
    if (method %in% c("dropbead", "scran")) {

      ggdim <- Seurat::DimPlot(object = seurat_ob,
                               dims = c(1, 2),
                               reduction = opt$reduct,
                               pt.size = ifelse(dim(seurat_ob)[2] < 500,1.5 ,as.numeric(opt$pt.sze)) ,
                               group.by = cellcycle_name) +
               ggplot2::scale_colour_manual(values = OESingleCell::SelectColors(palette = opt$colorschema,
                                                                                n= length(unique(seurat_ob@meta.data[, cellcycle_name]))))

      OESingleCell::save_ggplots(filename =glue::glue("{output_dir}/{opt$reduct}_resolution{resolution}_visualize_CellCycle_plot"),
                                 plot = ggdim,
                                 dpi=1000,
                                 limitsize = F)
      meta.data <-  seurat_ob@meta.data %>%
                    dplyr::mutate(cell_barcode = rownames(seurat_ob@meta.data)) %>%
                    dplyr::select(cell_barcode, everything())
      write.table(meta.data ,
                  file.path(output_dir, "cell_cycle_annotation_result.xls"),
                  col.names = T,
                  row.names = F,
                  sep = "\t",
                  quote = F)

      ## split visualizataion ###
      # output_dir = file.path(output_dir, paste0( "visualize_cluster_by_", cellcycle_name, collapse = ""))
      # if ( !file.exists(output_dir) ){
      #     dir.create(output_dir,recursive=T)
      # }
      if (!is.null(facetby)) {
        for (facetbyx in facetby) {
          seurat_ob <- Seurat::SetIdent(seurat_ob, value = cellcycle_name)
          DATA <- as.data.frame(meta.data[, c(facetbyx, cellcycle_name)]) %>%
                  table() %>%
                  as.data.frame() %>%
                  dplyr::rename(cell_number = Freq) %>%
                  dplyr::arrange(get(facetbyx)) %>%
                  dplyr::group_by(.dots = facetbyx) %>%
                  dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100) %>%
                  as.data.frame()
          write.table(DATA,
                      file.path(output_dir, file = paste0(facetbyx, "_clust_cond_freq_info.xls")),
                      sep = "\t",
                      col.names = T,
                      row.names = F,
                      quote = F)
          # clust_sum_all = PlotAbundances(seurat_ob, prop.by = cellcycle_name , group.by = facetbyx, method = "barplot",
          #     cols= CustomCol2(1:nlevel)
          # )
          sum_all <-  ggplot(DATA, aes_string(x = facetbyx, y = "freq", fill = cellcycle_name)) +
                      geom_bar(stat = "identity", position = opt$bartype) +
                      #geom_bar(stat ="identity",position ="stack") +
                      labs(x = " ", y = paste0("Proportion [%]")) +
                      theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            axis.line = element_line(color = "black")) +
                      scale_y_continuous(expand = c(0, 0)) + ## //这个可以去掉与X轴间隙
                      #scale_x_discrete(expand = c(0.5,0)) +
                      scale_fill_manual(values = COESingleCell::SelectColors(palette = opt$colorschema, n= nlevel)) +
                      guides(fill = guide_legend(title = "Cellcycle"))

          OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/groupby-{facetbyx}_resolution-{resolution}_summary_plot"),
                                     limitsize = FALSE,
                                     plot = sum_all)
          ## 降维聚类图
          seurat_ob <- Seurat::SetIdent(seurat_ob, value = cellcycle_name)
          nrow <- ceiling(length(unique(seurat_ob@meta.data[, facetbyx])) / 2)
          groupvis_split <- Seurat::DimPlot(object = seurat_ob,
                                            dims = c(1, 2),
                                            reduction = opt$reduct,
                                            pt.size = ifelse(dim(seurat_ob)[2] < 500,1.5 ,as.numeric(opt$pt.sze)) ,
                                            ncol = 2,
                                            group.by = cellcycle_name,
                                            split.by = facetbyx) +
                            theme(plot.title = element_text(hjust = 0.5))+
                            scale_colour_manual(values = OESingleCell::SelectColors(palette = opt$colorschema, n= nlevel))
          OESingleCell::save_ggplots(filename = glue::glue("{output_dir}/splitby-{facetbyx}_resolution{resolution}_split_plot"),
                                     limitsize = FALSE,
                                     plot = groupvis_split,
                                     height = 6 * nrow,
                                     width = 14)
        }
      }
    }
      ##  update or save seurat object
      if ( as.logical(opt$update) ){
        SeuratDisk::UpdateH5Seurat(file = opt$input,
                                   object = data_ob)
      }else{
        OESingleCell::SaveX(data_ob,
                            output = opt$output,
                            update = FALSE,
                            outformat = opt$outformat,
                            prefix = opt$prefix)
      }
      ## save session informations
      write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
      quit()
  }
}

