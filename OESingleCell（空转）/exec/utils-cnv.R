sub_cnv <- subparsers$add_parser("CNV", help = "Inferring CNV from Single-Cell RNA-Seq.")
sub_cnv$add_argument("--celltype", "-l", type = "character",
                     help = ("[Required] The cell type annotation column name to use in seurat metadata"))
sub_cnv$add_argument("--groupby", "-g", type = "character",
                     help = ("[Required] The groupping of cells in the metadata of seurat object to visualize"))
sub_cnv$add_argument("--colormapping", "-m", type = "character",
                     help = paste0("[OPTIONAL] The color mapping for groupping column of cells set by the ",
                                   "parameters '--groupby'.The exmaple format is variable1:colorschema1,",
                                   "variable2:colorschema2.The supported color schemas can be: blindless, col50, ",
                                   "ditto, paired, customecol2"))
sub_cnv$add_argument("--malignant", "-t", type = "character",
                     help = "[OPTIONAL] The cell type name to be assumed as cancer cells.")
sub_cnv$add_argument("--normalcells", "-r", type = "character",
                     help = "[OPTIONAL] A comma seperated list containing the classifications of normal cells.[default: %(default)s]")
sub_cnv$add_argument("--reduct", type = "character", default = "tsne",
                     help = "[OPTIONAL]The reduction used in the DimPlot.[default: %(default)s]")
sub_cnv$add_argument("--pointsize", type = "double", default = "0.8",
                     help = "[OPTIONAL]The pointsize used in the DimPlot.[default: %(default)s]")
sub_cnv$add_argument("--sample_ratio", "-s", type = "double", default = "0.8",
                     help = paste0("[OPTIONAL] The ratio of random subsample for each group. Only normal cells",
                                   " will be subseted.The number of all cells must be less than 60,000 cells, ",
                                   "or it will out of memory in step 4: measuring baselines.[default: %(default)s]"))
sub_cnv$add_argument("--ident2use", "-q", type = "character", default = NULL,
                     help = paste0("[OPTIONAL] The column name in cell metadata used as identity of each cell",
                                   " combined with which_cell.[default: %(default)s]"))
sub_cnv$add_argument("--which_cells", "-u", type = "character", default = NULL,
                     help = "[OPTIONAL] The subset of cluster ids used for analysis.[default: %(default)s]")

####inforcnv部分
sub_cnv$add_argument("--gene_order", type = "character",
                     help = paste0("[OPTIONAL]data file containing the positions of each gene along each ",
                                   "chromosome in the genome.The gene_order_file, contains chromosome, start, and ",
                                   "stop position for each gene, tab-delimited without header[default: %(default)s]"))
sub_cnv$add_argument("--mode", type = "character", default = "subclusters",
                     help = paste0("[OPTIONAL]Grouping level for image filtering or HMM predictions. ",
                                   "Options can be:samples,subclusters,cells.samples: fastest, but subclusters",
                                   "is ideal.subclusters: detect the subclusters  of tumor[default: %(default)s]"))
sub_cnv$add_argument("--clusting2use", type = "character", default = "ward.D2",
                     help = paste0("[OPTIONAL] the hierarchical clustering methods of cells to use.",
                                   "Valid choices are: 'ward.D', 'ward.D2', 'single', 'complete',",
                                   "'average', 'mcquitty', 'median', 'centroid'.[default: %(default)s]"))
sub_cnv$add_argument("--refexp", type = "character",
                     help = "[OPTIONAL]the reference (normal) cells expression count matrix for inferring cnv")
sub_cnv$add_argument("--gtf", "-z", type = "character",
                     help = "[OPTIONAL] the exact gtf file for the alignment for this project.[default: %(default)s]")
sub_cnv$add_argument("--cutoff", "-c", type = "double", default = 0.1,
                     help = paste0("[OPTIONAL] min average read counts per gene among reference cells.",
                                   "cutoff=1 works well for Smart-seq2, ",
                                   "and cutoff=0.1 works well for 10x Genomics"))
sub_cnv$add_argument("--doHMM", type = "logical", default = FALSE,
                     help = paste0("[OPTIONAL] wether to using HMM when predicting the CNV for each cell. If no",
                                   " CNV gene prediction needed, Set it False to save time.[default: %(default)s]"))
sub_cnv$add_argument("--pval", "-p", type = "double", default = 0.05,
                     help = "[OPTIONAL] max p-value for defining a significant tumor subcluster.[default: %(default)s]")

sub_cnv$add_argument("--method", type = "character",
                     help = "the mothed to infer CNV.Can be infercnv or copykat")
###input,informat,outformat ncores assays output_dir
args <- commandArgs(TRUE)
if ('CNV' %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "CNV") {
    ##step1 read the specified assay and data slot in data object into memory===========================================
    # seurat_ob <- OESingleCell::ReadX(input = opt$input,
    #                              informat = opt$informat,
    #                              assays = assays,
    #                              data.use = dataslots,  # data slot is enough
    #                              verbose = F)
    seurat_ob <- readRDS(opt$input)
    Seurat::DefaultAssay(seurat_ob) <- assays
    options(scipen=100)

    # step2 get the subset of cells used for visualization if necessay==================================================
    if (!is.null(opt$which_cells)) {
      cluster_list <- unlist(strsplit(opt$which_cells, ",", perl = T))
      seurat_ob <- seurat_ob[, which(x = factor(seurat_ob[[ident2use]][, 1]) %in% cluster_list)]
    }

   # step3  cnv reference using infercnv or copykat
    run_infer_cnv <- switch(
      opt$method,
      ##=======================================infercnv=================================================================
      "infercnv" = {
        message('step3: cnv reference using infercnv ')
        message("3.1: gene_order")
        if (!is.null(opt$gene_order)) {
          gene_order_f <- opt$gene_order
        }else if (!is.null(opt$gtf)) {
          gtf <- plyranges::read_gff(opt$gtf)
          gene.chr <- gtf %>%
                      plyranges::filter(type == "gene" & gene_name %in% rownames(seurat_ob)) %>%
                      as.data.frame() %>%
                      dplyr::select(gene_name, seqnames, start, end) %>%
                      dplyr::distinct(gene_name, .keep_all = T) %>%
                      dplyr::mutate(seqnames = paste0("chr", seqnames))
          gene_order_f <- glue::glue("{tempdir()}/gene_order_file.xls")
          print(gene_order_f)
          write.table(gene.chr, gene_order_f, col.names = F, row.names = F, sep = "\t", quote = F)
          print("gene_order file文件生成")
        }else {
          stop("NO gene coordination annotation file is not available!")
        }
        message("3.2: fetch the needed information")
        cellanno <- Seurat::FetchData(seurat_ob, vars = opt$celltype) %>%
                    tibble::rownames_to_column(var = "cellbarcode")
        print(head(cellanno))
        # if(!is.null(opt$normalcells)&!is.null(opt$malignant)){
        #   refcells <- unlist(strsplit(opt$normalcells, ",", perl = T))
        #   malignant <- unlist(strsplit(opt$malignant, ",", perl = T))
        # }
        #此处有逻辑问题
        if (!is.null(opt$normalcells)) {
          refcells <- unlist(strsplit(opt$normalcells, ",", perl = T))
          count_mat <- Seurat::GetAssayData(seurat_ob, "counts")#seurat_ob：全部细胞
        }else if (!is.null(opt$malignant)) {
          malignant <- unlist(strsplit(opt$malignant, ",", perl = T))
          print(malignant)
          refcells <- setdiff(unique(seurat_ob@meta.data[, opt$celltype]), malignant)
          count_mat <- Seurat::GetAssayData(seurat_ob, "counts")
          print(colnames(count_mat)[1:10])
        }else {
          print("NO reference normal cells or malignant cells are specified!
                                         The internal customized normal reference data will be used!")
          refexp <- OESingleCell::readRDSMC(opt$refexp, cores = 10)
          refcell_anno <- data.frame(cellbarcode = colnames(refexp), celltype = "normal")
          cellanno <- rbind(cellanno, refcell_anno)
          com.genes <- intersect(rownames(seurat_ob), rownames(refexp))
          count_mat <- cbind(Seurat::GetAssayData(seurat_ob, "counts")[com.genes,], refexp[com.genes,])
        }
        tempdir <- tempdir()
        cnv_celltyping <- file.path(tempdir, "cnv_celltype_group.xls")
        write.table(cellanno, cnv_celltyping, sep = "\t", col.names = F, row.names = F, quote = F)
        #gene_order_f <- file.path(tempdir, "gene_order_file.xls")
        #write.table(gene.chr, gene_order_f, col.names = F, row.names = F, sep = "\t", quote = F)
        print(cnv_celltyping)
        message('3.3.CreateInfercnvObject')
        infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = count_mat,
                                                       annotations_file = cnv_celltyping,
                                                       delim = "\t",
                                                       gene_order_file = gene_order_f,
                                                       ref_group_names = refcells)
        message('3.4.run infercnv')
        infercnv_obj <- infercnv::run(infercnv_obj,
                                      cutoff = opt$cutoff,
                                      analysis_mode = opt$mode,
                                      tumor_subcluster_pval = opt$pval,
                                      hclust_method = opt$clusting2use,
                                      out_dir = output_dir,
                                      num_threads = opt$ncores,
                                      cluster_by_groups = TRUE,
                                      denoise = TRUE,
                                      no_plot = T,
                                      no_prelim_plot = F,
                                      HMM = opt$doHMM)
                                      #tumor_subcluster_partition_method = 'leiden')
        message('3.5.plot heatmap')
        pdf(NULL) ##for removing  Rplot.pdf files
        plot<-ComplexHeatmap::Heatmap(t(as.matrix(infercnv_obj@expr.data)),
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                show_row_names = F,
                                show_column_names = F,
                                name = "CNV level",
                                use_raster = TRUE,
                                raster_quality = 4)
        ###save plot =======================================================================================================
        OESingleCell::save_ggplots(filename=glue::glue("{output_dir}/heatmap_infercnv"),
                                   plot = ggplotify::as.ggplot(plot),
                                   dpi=600)
      },
      #=================================================================================================================
      "copykat" = {
        message("step3: cnv reference using copykat")
        if (!is.null(opt$groupby)) {
          cell.annos <- unlist(strsplit(opt$groupby, ",", perl = T))
          groupby <- unlist(strsplit(opt$groupby, ",", perl = T))[1]
        }else {
          stop("please provide the groupping of cells in the metadata of seurat object.")
        }
        # set the color schema of cell annotation bar
        if (!is.null(opt$colormapping)) {
          group_colors <- list()
          for (x in  unlist(strsplit(opt$colormapping, ",", perl = T))) {
            m <- unlist(strsplit(x, ":", perl = T))
            group_colors[[m[1]]] <- m[-1]
          }
        }else {
          group_colors <- 1:(length(cell.annos) + 1)
          names(group_colors) <- c(cell.annos, "copykat.pred")
          group_colors <- as.list(group_colors)
          group_colors[[1]] <- "ditto"
          group_colors[[2]] <- 'blindless'
        }

        if (!is.null(opt$normalcells)) {
          normalcell <- unlist(strsplit(opt$normalcells, ",", perl = T))
          normal.cells <- rownames(seurat_ob@meta.data[which(seurat_ob@meta.data[, opt$celltype] %in% normalcell),])
        }else if (!is.null(opt$malignant)) {
          malignant <- unlist(strsplit(opt$malignant, ",", perl = T))
          malignant.cells <- rownames(seurat_ob@meta.data[which(seurat_ob@meta.data[, opt$celltype] %in% malignant),])
          normal.cells <- setdiff(Seurat::Cells(seurat_ob), malignant.cells)
        }else {
          print("NO normal cells or malignant cells are specified!")
          print("CopyKAT will predict normal cells automatically!")
          normal.cells <- ""
        }

        ## Needs to be less than 60,000 cells, or it will out of memory in step 4: measuring baselines...
        if (!is.null(opt$sample_ratio)) {
          if (normal.cells[1] != "") {
            sampled_cellmeta <- seurat_ob@meta.data[normal.cells,] %>%
                                tibble::rownames_to_column() %>%
                                dplyr::group_by(.dots = opt$celltype) %>%
                                dplyr::sample_frac(size = opt$sample_ratio, replace = F) %>%
                                tibble::column_to_rownames()
            subset_cells <- c(rownames(sampled_cellmeta), setdiff(Seurat::Cells(seurat_ob), normal.cells))
            seurat_ob <- subset(seurat_ob, cells = subset_cells)
          }else {
            sampled_cellmeta <- seurat_ob@meta.data %>%
                                tibble::rownames_to_column() %>%
                                dplyr::group_by(.dots = opt$celltype) %>%
                                dplyr::sample_frac(size = opt$sample_ratio, replace = F) %>%
                                column_to_rownames()
            seurat_ob <- subset(seurat_ob, cells = rownames(sampled_cellmeta))
          }
        }

        ##
        setwd(output_dir)
        copykat.res <- OESingleCell::copykat(rawmat = Seurat::GetAssayData(seurat_ob, slot = "counts"),
                                             id.type = "S",
                                             cell.line = "no",
                                             ngene.chr = 5,
                                             win.size = 25,
                                             KS.cut = 0.1,
                                             sam.name = "seurat",
                                             distance = "euclidean",
                                             norm.cell.names = normal.cells,
                                             n.cores = opt$ncores)
        CNA_mat <- data.frame(copykat.res$CNAmat)
        cell.pred <- data.frame(copykat.res$prediction)
        cell.pred <- cell.pred[, 2, drop = F]
        cell.pred <- cell.pred %>% dplyr::arrange(copykat.pred)

        meta.data <- seurat_ob@meta.data[rownames(cell.pred),]
        cell.anno <- cbind(meta.data, cell.pred)
        cell.anno <- cell.anno[, c(cell.annos, "copykat.pred")]
        cell.anno_data <- tibble::rownames_to_column(as.data.frame(cell.anno), "barcode")
        write.table(cell.anno_data, file.path(output_dir, "copykat_prediction.xls"), quote = F, row.names = F, sep = "\t")

        gene.anno <- as.data.frame(CNA_mat[, 1, drop = F])
        plot_data <- t(CNA_mat[, 4:ncol(CNA_mat)])

        color_schema_row <- list()
        color_schema_col <- list()

        for (x in colnames(cell.anno)) {
          nlevels <- length(unique(cell.anno[, x]))
          color_schema_row[[x]] <- OESingleCell::SelectColors(cell.anno,
                                                              value = x,
                                                              n = nlevels,
                                                              palette = group_colors[[x]])
        }
        color_schema_col[["chrom"]] <- OESingleCell::SelectColors(object = gene.anno,
                                                                  value = "chrom",
                                                                  n = length(unique(gene.anno[, "chrom"])),
                                                                  palette = "ditto")
        cellAnnotation <- ComplexHeatmap::HeatmapAnnotation(df = cell.anno,
                                                            col = color_schema_row,
                                                            which = "row",
                                                            show_annotation_name = T)

        geneAnnotation <- ComplexHeatmap::HeatmapAnnotation(df = gene.anno, name = "chr",
                                                            col = color_schema_col,
                                                            show_annotation_name = T)
        col.split.by <- "chrom"
        coldata_split <- CNA_mat[, col.split.by, drop = F]
        plot_data <- plot_data[rownames(cell.anno),]
        pdf(NULL) ##for removing  Rplot.pdf files
        plot<-ComplexHeatmap::Heatmap(t(scale(t(plot_data))),
                                name = "CNV level",
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                # clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
                                col = circlize::colorRamp2(c(-2, 0, 2), c("#406AA8", "white", "#D91216")),
                                # colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
                                # colorRamp2(c(-2, 0, 2), c("#053061", "white", "#67001F"))
                                row_order = rownames(cell.anno),
                                # row_split = rowdata_split,
                                row_gap = grid::unit(0, "mm"),
                                # column_order = colnames(plot_data),
                                column_split = coldata_split,
                                column_gap = grid::unit(0, "mm"),
                                column_title = unique(CNA_mat[, "chrom"]),
                                column_title_gp = grid::gpar(fontsize = 7),
                                column_title_side = "bottom",
                                border = T,
                                # top_annotation = geneAnnotation ,
                                left_annotation = cellAnnotation,
                                use_raster = TRUE,
                                raster_quality = 4,
                                show_heatmap_legend = TRUE,
                                show_row_names = F,
                                show_column_names = F)
        ##save ggplot
        OESingleCell::save_ggplots(plot = plot,
                                   filename=glue::glue("{output_dir}/copykat_heatmap"),
                                   width = 10,
                                   height = 10,
                                   dpi=1000)
        ## Dimplot======================================================================================================
        sub_ob <- subset(seurat_ob, cells = rownames(cell.anno))
        sub_ob@meta.data$copykat.pred <- cell.anno[rownames(sub_ob@meta.data), "copykat.pred"]
        sub_ob <- Seurat::SetIdent(sub_ob, value = "copykat.pred")
        ggtsne <- Seurat::DimPlot(object = sub_ob,
                                  reduction = opt$reduct,
                                  pt.size = opt$pointsize) +
                  ggplot2::scale_colour_manual(values = c("#7fc97f", "#beaed4"))
        ##save ggplot
        OESingleCell::save_ggplots(plot = ggtsne,
                                   filename=glue::glue("{output_dir}/copykat_prediction_tsne_plot"),
                                   width = length(unique(cellmeta[[opt$groupby]])),
                                   height = length(gset_list) * 2,
                                   dpi=1000)
      },
      "NO supported method supplied!"
    )
  }
}