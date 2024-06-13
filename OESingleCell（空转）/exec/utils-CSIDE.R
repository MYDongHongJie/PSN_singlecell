docstring <- " example1:\\n\\n\\
sctool -i st.rds -f rds  -o results -d rds  --assay Spatial  --dataslot counts  CSIDE --refobject sc.rds  --refcelltype celltype  --misclist RCTD_results "
sub_CSIDE <- subparsers$add_parser(
  "CSIDE",
  description = docstring,
  formatter_class = "argparse.RawTextHelpFormatter",
  # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Using CSIDE to get DE genes within cell type across multiple regions "
)
sub_CSIDE$add_argument(
  "--refobject",
  type = "character",
  help = "the customized reference expression matrix in seurat format."
)
sub_CSIDE$add_argument(
  "--refassay",
  type = "character",
  default = "RNA",
  help = "the default assay for reference data to use for the this run.[default: %(default)s]"
)
sub_CSIDE$add_argument(
  "--refcelltype",
  type = "character",
  help = "[REQUIRED]the cell type annotation column id in the reference expression matrix."
)
sub_CSIDE$add_argument(
  "--misclist",
  type = "character",
  help = "[REQUIRED]the matrix name of the celltype annotation or other functional matrix in the misclist.[default: %(default)s]."
)
sub_CSIDE$add_argument(
  "--contrast",
  type = "character",
  default = "all",
  help = "sampleid and groupby informations. the format is: 'sampleid:clusters:treatment:control01:...'.[default: %(default)s]."
)
sub_CSIDE$add_argument(
  "--mode",
  type = "character",
  default = "nonparam",
  help = "CSIDE mode for each sample, it should be choosed from 'discrete','continuous','regions' or 'nonparam'.[default: %(default)s]."
)
sub_CSIDE$add_argument(
  "--remove_ct",
  type = "character",
  default = NULL,
  help = "remove celltypes which has not converged on any genes. (noninformative).[default: %(default)s]."
)
sub_CSIDE$add_argument(
  "--thr_gene",
  type = "double",
  default = 5e-5,
  help = "minimum average normalized expression required for selecting genes.[default: %(default)s]."
)
sub_CSIDE$add_argument(
  "--thr_ct",
  type = "integer",
  default = 125,
  help = "min occurrence of number of cells for each cell type to be used as aggregated by aggregate_cell_types.[default: %(default)s]."
)
sub_CSIDE$add_argument(
  "--thr_logfc",
  type = "double",
  default = 0.4,
  help = "the natural log fold change cutoff for differential expression.[default: %(default)s]."
)
# sub_CSIDE$add_argument(
#   "--thr_weight",
#   type = "double",
#   default = NULL,
#   help = "the threshold of total normalized weights across all cell types in per pixel. Default 0.99 for doublet_mode or 0.8 for full_mode.[default: %(default)s]."
# )
sub_CSIDE$add_argument(
  "--fdr",
  type = "double",
  default = 0.01,
  help = "false discovery rate for hypothesis testing.[default: %(default)s]."
)

args <- commandArgs(TRUE)
if ("CSIDE" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "CSIDE") {
    ################################# prepare for object data ##########################################################
    futile.logger::flog.info("step1: load query and reference object data")
    #load query object data=============================================================================================
    suppressMessages(data_ob <- OESingleCell::ReadX(
      input = opt$input,
      informat = opt$informat,
      assays = assays[1],
      data.use = dataslots,
      verbose = F
    ))
    Seurat::DefaultAssay(data_ob) <- assays[1]

    #load reference object data=========================================================================================
    suppressMessages(ref_ob <- OESingleCell::ReadX(
      input = opt$refobject,
      informat = "rds",
      assays = opt$refassay,
      data.use = "data,counts",
      verbose = F
    ))

    #prepare for reference object data==================================================================================
    futile.logger::flog.info("step2: prepare for CSIDE object data")
    if (!"nCount_RNA" %in% colnames(ref_ob@meta.data)) {
      nCount_RNA <- SeuratObject::colSums(x = ref_ob, slot = "counts") %>%
        as.data.frame %>%
        dplyr::rename(nCount_RNA = '.')
      ref_ob <- Seurat::AddMetaData(ref_ob, metadata = nCount_RNA)
    }
    counts_ref <- as.data.frame(as.matrix(ref_ob@assays[["RNA"]]@counts))
    cell_types_ref <- ref_ob@meta.data[[opt$refcelltype]]
    cell_types_ref <- stringi::stri_replace_all_fixed(cell_types_ref, c(".","/","-"," "), "_", vectorize_all=FALSE)
    names(cell_types_ref) <- rownames(ref_ob@meta.data)
    cell_types_ref <- as.factor(cell_types_ref)
    nUMI_ref <- ref_ob@meta.data$nCount_RNA
    names(nUMI_ref) <- rownames(ref_ob@meta.data)
    reference <- spacexr::Reference(counts_ref, cell_types_ref, nUMI_ref)

    ################################# run pipeline for each sample #####################################################
    contrast <- stringr::str_split(opt$contrast, ":")[[1]]
    if (contrast[1] != "all") { data_ob <- subset(data_ob, sampleid == contrast[1]) }
    for (sampleid in unique(data_ob@meta.data$sampleid)) {
      #get subset of cells if necessary=================================================================================
      data_ob_sub <- subset(data_ob, sampleid == sampleid)
      if (!is.na(contrast[2])) {
        desired_cells <- subset(data_ob_sub@meta.data, data_ob_sub@meta.data[[contrast[2]]] %in% contrast[c(-1,-2)])
        data_ob_sub <- subset(data_ob_sub, cells = rownames(desired_cells))
      }
      data_ob_sub@meta.data[[contrast[2]]] <- factor(data_ob_sub@meta.data[[contrast[2]]],
                                                     levels = sort(unique(data_ob_sub@meta.data[[contrast[2]]])))
      #data_ob_sub@misc[[opt$misclist]] <- data_ob_sub@misc[[opt$misclist]] %>% dplyr::filter(barcodes %in% rownames(data_ob_sub@meta.data))
      if (!is.null(Seurat::Images(data_ob_sub))) {
        unuse_images <- Seurat::Images(data_ob_sub)[!Seurat::Images(data_ob_sub) %in% (data_ob_sub@meta.data$sampleid %>% unique())]
        if (length(unuse_images) > 0) {
          data_ob_sub@images[unuse_images] <- NULL
        }
      }

      #prepare for puck object data=====================================================================================
      counts <- as.data.frame(data_ob_sub@assays[["Spatial"]]@counts)
      coords <- do.call(rbind, lapply(Seurat::Images(data_ob_sub), function(slice) {
        spatial_coord <- data.frame(data_ob_sub[[slice]]@coordinates)
        colnames(spatial_coord) <- c("in_tissue", "y", "x", "pxl_col_in_fullres", "pxl_row_in_fullres")
        spatial_coord$y <- -1 * (spatial_coord$y) + 78
        return(spatial_coord)
      })) %>% dplyr::select(c(x, y))
      puck <- spacexr::SpatialRNA(coords, counts)

      #creat CSIDE object data and import weights=======================================================================
      futile.logger::flog.info("step3: creat CSIDE object data and import weights")
      myRCTD <- spacexr::create.RCTD(puck, reference)
      weights <- data_ob_sub@misc[[opt$misclist]] %>%
                 dplyr::filter(barcodes %in% rownames(data_ob_sub@meta.data)) %>%
                 tibble::column_to_rownames(var = "barcodes")
      colnames(weights) <- stringi::stri_replace_all_fixed(colnames(weights), c(".","/","-"," "), "_", vectorize_all=FALSE)
      weights <- sweep(
        weights %>% as.data.frame(),
        1,
        rowSums(weights %>% as.data.frame()),
        "/"
      )
      myRCTD <- spacexr::import_weights(myRCTD, weights)

      #run CSIDE========================================================================================================
      barcodes <- colnames(myRCTD@spatialRNA@counts)
      cell_types <- colnames(myRCTD@results$weights)
      if (!is.null(opt$remove_ct)) {
        cell_types <- cell_types[!cell_types %in% stringr::str_split(opt$remove_ct, ",")[[1]]]
      }
      myRCTD@config$max_cores <- opt$ncores
      if (!opt$mode %in% c("discrete", "continuous", "regions", "nonparam")) {
        stop("No CSIDE mode found, please select one mode from 'discrete', 'continuous','regions' or 'nonparam'.")
      }

      futile.logger::flog.info("step4: run CSIDE")
      if (opt$mode %in% c("discrete", "continuous")) {
        treatment <- subset(data_ob_sub@meta.data, data_ob_sub@meta.data[[contrast[2]]] == contrast[3])
        if (opt$mode == "discrete") {
          explanatory.variable <- as.integer(rownames(myRCTD@spatialRNA@coords) %in% rownames(treatment))
          names(explanatory.variable) <- rownames(myRCTD@spatialRNA@coords)
        }else if (opt$mode == "continuous") {
          continuous_coords <- myRCTD@spatialRNA@coords[rownames(treatment),]
          explanatory.variable <- spacexr::exvar.point.density(myRCTD, barcodes, continuous_coords, radius = 500)
        }
        myRCTD <- spacexr::run.CSIDE.single(myRCTD,
                                             explanatory.variable,
                                             cell_types = cell_types,
                                             doublet_mode = FALSE,
                                             gene_threshold = opt$thr_gene,
                                             cell_type_threshold = opt$thr_ct,
                                             fdr = opt$fdr,
                                             log_fc_thresh = opt$thr_logfc,
                                             weight_threshold = NULL)
        plot_puck <- spacexr::plot_puck_continuous(myRCTD@spatialRNA,
                                                   names(explanatory.variable),
                                                   explanatory.variable,
                                                   ylimit = c(0, 1),
                                                   title = 'plot of explanatory variable')
      }


      if (opt$mode == "regions") {
        region_list <- list()
        for (i in 3:length(contrast)) {
          region <- subset(data_ob_sub@meta.data, data_ob_sub@meta.data[[contrast[2]]] == contrast[i])
          region_list[[i-2]] <- rownames(region)
        }
        myRCTD <- spacexr::run.CSIDE.regions(myRCTD,
                                              region_list,
                                              cell_types = cell_types,
                                              doublet_mode = FALSE,
                                              gene_threshold = opt$thr_gene,
                                              cell_type_threshold = opt$thr_ct,
                                              fdr = opt$fdr,
                                              log_fc_thresh = opt$thr_logfc,
                                              weight_threshold = NULL)
        class_num <- rep(0, length(barcodes))
        names(class_num) <- barcodes
        plot_puck <- spacexr::plot_class(myRCTD@spatialRNA,
                                         barcodes,
                                         factor(class_num),
                                         title = 'plot of regions')
      }
      if (opt$mode == "nonparam") {
        myRCTD <- spacexr::run.CSIDE.nonparam(myRCTD,
                                               df = 15,
                                               cell_types = cell_types,
                                               doublet_mode = FALSE,
                                               gene_threshold = opt$thr_gene,
                                               cell_type_threshold = opt$thr_ct,
                                               fdr = opt$fdr,
                                               weight_threshold = NULL)
        X <- spacexr::build.designmatrix.nonparam(myRCTD, df = 15)
        plot_puck <- spacexr::plot_puck_continuous(myRCTD@spatialRNA,
                                                   myRCTD@internal_vars_de$barcodes,
                                                   X[, 3],
                                                   ylimit = c(0, 1),
                                                   title = 'plot of third basis function')
      }
      saveRDS(myRCTD, glue::glue("{output_dir}/myRCTD_{sampleid}.rds"))
      OESingleCell::save_ggplots(glue::glue("{output_dir}/explanatory_variable.plot"),
                   plot = plot_puck,
                   limitsize = FALSE,
                   height = 6,
                   dpi = 600
      )

      for (i in names(myRCTD@de_results$sig_gene_list)) {
          write.table(myRCTD@de_results$sig_gene_list[[i]] %>% tibble::rownames_to_column(var = "gene"),
                      glue::glue("{output_dir}/{i}_sig_gene_list.xls"),
                      sep = "\t",
                      quote = FALSE,
                      row.names = FALSE)
      }
      for (i in names(myRCTD@de_results$all_gene_list)) {
          write.table(myRCTD@de_results$all_gene_list[[i]] %>% tibble::rownames_to_column(var = "gene"),
                      glue::glue("{output_dir}/{i}_all_gene_list.xls"),
                      sep = "\t",
                      quote = FALSE,
                      row.names = FALSE)
      }
    }
    if (!file.exists(file.path(output_dir, "C-SIDE差异分析说明文档.docx"))) {
        file.copy("/public/dev_scRNA/oesinglecell3_test/document/C-SIDE差异分析说明文档.docx",
                  file.path(output_dir, "C-SIDE差异分析说明文档.docx")) }
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}