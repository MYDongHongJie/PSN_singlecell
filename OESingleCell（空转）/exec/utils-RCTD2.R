##seurat object to RCTD reference object for singlecell transcriptome data
#' @param ref_ob  the seurat object
#' @param ref_assay   the assay for the seurat object
#' @param ref_cell_type the celltype name in the seurat object metadata
#' @return Returns a \code{\linkS4class{SpatialRNA}} object containing the coordinates and counts
#' from the input files
#' @export
seurat_to_reference <- function(ref_ob, ref_assay, ref_cell_type) {
  if (!"nCount_RNA" %in% colnames(ref_ob@meta.data)) {
    nCount_RNA <- SeuratObject::colSums(x = ref_ob, slot = "counts") %>%
      as.data.frame %>%
      dplyr::rename(nCount_RNA = '.')
    ref_ob <- Seurat::AddMetaData(ref_ob, metadata = nCount_RNA)
  }
  counts_ref <- OESingleCell::GetAssayData(ref_ob, assay = ref_assay, slot = "counts")
  cell_types_ref <- OESingleCell::colData(ref_ob) %>%
    dplyr::select(ref_cell_type) %>%
    dplyr::mutate(cell_type = stringi::stri_replace_all_fixed(!!!rlang::syms(ref_cell_type), "/", "_") %>% as.factor()) %>%
    .$cell_type
  cell_types_ref <- setNames(cell_types_ref, OESingleCell::colData(ref_ob) %>% rownames)
  nUMI_ref <- OESingleCell::colData(ref_ob) %>%
    dplyr::select("nCount_RNA") %>%
    .$nCount_RNA
  nUMI_ref <- setNames(nUMI_ref, OESingleCell::colData(ref_ob) %>% rownames)
  reference <- spacexr::Reference(counts_ref, cell_types_ref, nUMI_ref)
  print("BB")
  return(reference)
}

##seurat object to RCTD object for spatial transcriptome data
#' @param st_ob  the seurat object
#' @param st_assay   the assay for the seurat object
#' @return Returns a \code{\linkS4class{SpatialRNA}} object containing the coordinates and counts
#' from the input files
#' @export
seurat_to_SpatialRNA <- function(st_ob, st_assay) {
  images <- Seurat::Images(st_ob)
  coords <- do.call(rbind, lapply(images, function(slice) {
    spatial_coord <- data.frame(st_ob[[slice]]@coordinates)
    colnames(spatial_coord) <- c("in_tissue", "y", "x", "pxl_col_in_fullres", "pxl_row_in_fullres")
    spatial_coord$y <- -1 * (spatial_coord$y) + 78
    spatial_coord <- spatial_coord %>% dplyr::select(c(x, y))
    return(spatial_coord)
  }))
  counts <- Seurat::GetAssayData(st_ob, assay = st_assay, slot = "counts")
  puck <- spacexr::SpatialRNA(coords, counts)
  return(puck)
}

#' Tidy RCTD results
#' @param rctd an RCTD object containing the results of the RCTD algorithm.documentation for more information on
#' interpreting the content of the RCTD object
#' @param rctd_mode the mode of RCTD analysis, such as in 'multi' mode,'doublet' mode, or 'full' mode
#' @param conf whether only return the confident weights when user choose 'multi' mode
#' @return a data.frame of cell type weights for every pixel
#' @export
rctd4weight <- function(rctd,
                        rctd_mode = c("multi", "full", "doublet"),
                        conf = TRUE) {
  rctd_mode <- match.arg(rctd_mode)
  if (!is.null(rctd) & !is.null(rctd_mode)) {
    spot_bc <- colnames(rctd@spatialRNA@counts)
    celltypes2use <- colnames(rctd@cell_type_info$info[[1]])
    if (rctd_mode == "multi") {
      if (conf) {
        nspot <- length(spot_bc)
        ncelltypes <- length(celltypes2use)
        rctd_weights_mtx <- matrix(rep(0, length(nspot * ncelltypes)),
                                   nrow = nspot, ncol = ncelltypes)
        rownames(rctd_weights_mtx) <- spot_bc
        colnames(rctd_weights_mtx) <- celltypes2use
        for (spot in seq_len(nrow(rctd_weights_mtx))) {
          rctd_weights_mtx[spot, match(rctd@results[[spot]]$cell_type_list[rctd@results[[spot]]$conf_list],
                                       colnames(rctd_weights_mtx))] <- rctd@results[[spot]]$sub_weights[rctd@results[[spot]]$conf_list]
        }
        rctd_weights <- data.frame(rctd_weights_mtx[, apply(rctd_weights_mtx, 2,
                                                            sum) > 0])
      }else {
        nspot <- length(spot_bc)
        ncelltypes <- length(celltypes2use)
        rctd_weights_mtx <- matrix(rep(0, length(nspot * ncelltypes)),
                                   nrow = nspot, ncol = ncelltypes)
        rownames(rctd_weights_mtx) <- spot_bc
        colnames(rctd_weights_mtx) <- celltypes2use
        for (spot in seq_len(nrow(rctd_weights_mtx))) {
          rctd_weights_mtx[spot, match(rctd@results[[spot]]$cell_type_list,
                                       colnames(rctd_weights_mtx))] <- rctd@results[[spot]]$sub_weights
        }
        rctd_weights <- data.frame(rctd_weights_mtx[, apply(rctd_weights_mtx, 2,
                                                            sum) > 0])
      }
    } else if (rctd_mode == "doublet") {
      doublet_results <- rctd@results$results_df
      doublet_weights <- as.data.frame(rctd@results$weights_doublet)
      rctd_doublet_weights <- data.frame(
        first_type = doublet_results@first_type,
        first_weight = doublet_weights@first_type,
        second_type = doublet_results@second_type,
        second_weight = doublet_weights@second_type,
        spot_class = doublet_results@spot_class)
      rownames(rctd_doublet_weights) <- rownames(doublet_results)
      doublet_celltype <- stringr::str_sort(union(doublet_results$first_celltype,
                                                  doublet_results$second_celltype))
      rctd_weights <- data.frame(ID = rownames(doublet_results))
      rctd_weights[, doublet_celltype] <- 0
      rownames(rctd_weights) <- rctd_weights$ID
      rctd_weights <- rctd_weights[, -1]
      for (ct in doublet_celltype) {
        rctd_weights[which(rctd_doublet_weights$first_type == ct), ct] <- rctd_doublet_weights[which(RCTD_doublet_weights$first_type == ct), "first_weight"]
        rctd_weights[which(rctd_doublet_weights$second_type == ct), ct] <- rctd_doublet_weights[which(rctd_doublet_weights$second_type == ct), "second_weight"]
      }
    } else if (rctd_mode == "full") {
      rctd_weights <- data.frame(rctd@results$weights)
    }
  }
  norm_weights <- sweep(
    rctd_weights %>% as.data.frame(),
    1,
    rowSums(rctd_weights %>% as.data.frame()),
    "/"
  )
  return(norm_weights)
}

#=======================================================================================================================
docstring <- " example1:\\n\\n\\
sctool  -i  st.rds   -f rds  -o results  -d rds   --assay Spatial  --dataslot counts  RCTD2 --mode full --refexp sc.rds  --refcelltype  celltype --refassay  RNA --multi_types 5"
sub_spacexr <- subparsers$add_parser(
  "RCTD2",
  description = docstring,
  formatter_class = "argparse.RawTextHelpFormatter",
  # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Using RCTD to refer Spaital transcriptome cell type"
)
sub_spacexr$add_argument(
  "--refexp",
  type = "character",
  help = "[REQUIRED]the customized reference expression matrix in seurat format, which may come from the quantification results of microarray, Bulk RNA-seq and scRNA-seq sequencing."
)
sub_spacexr$add_argument(
  "--refcelltype",
  type = "character",
  help = "[REQUIRED]the cell type annotation column id in the reference expression matrix."
)
sub_spacexr$add_argument(
  "--refassay",
  type = "character",
  default = "RNA",
  help = "[REQUIRED]the default assay for reference data to use for the this run.[default: %(default)s] "
)
sub_spacexr$add_argument(
  "-m",
  "--mode",
  default = "full",
  help = "mode: doublet, multi or full, for 10x visium data,default is full [default:\"%(default)s\"]"
)
sub_spacexr$add_argument(
  "--multi_types",
  type = "integer",
  default = NULL,
  help = "the number of celltype  for each spot when using the multi mode.[default: %(default)s] "
)

args <- commandArgs(TRUE)
if ("RCTD2" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "RCTD2") {
    # ==================================================================================================================
    futile.logger::flog.info("step1:load query object data")
    suppressMessages(data_ob <- OESingleCell::ReadX(
      input = opt$input,
      informat = opt$informat,
      assays = assays[1],
      data.use = dataslots,
      verbose = F
    ))
    Seurat::DefaultAssay(data_ob) <- assays[1]
    cellmeta <- slot(data_ob, "meta.data")
    # ==================================================================================================================
    futile.logger::flog.info("step2:load reference object data")
    suppressMessages(ref_ob <- OESingleCell::ReadX(
      input = opt$refexp,
      informat = "rds",
      assays = opt$refassay,
      data.use = "data,counts",
      verbose = F
    ))
    # the input is a seurat object which may contain more than one assay
    SeuratObject::DefaultAssay(ref_ob) <- opt$refassay
    ref_ob <- SeuratObject::SetIdent(ref_ob, value = opt$refcelltype)
    # ==================================================================================================================
    futile.logger::flog.info("step3:reformart data object")
    reference <- seurat_to_reference(ref_ob, opt$refassay, opt$refcelltype)
    puck <- seurat_to_SpatialRNA(data_ob, opt$assay)
    # ==================================================================================================================
    futile.logger::flog.info("step4:Creates an RCTD object from a scRNA-seq reference Seurat object and a SpatialRNA object")
    if (!file.exists(glue::glue("{output_dir}/myRCTD.rds"))) {
      RCTD_ob <- spacexr::create.RCTD(
        puck,
        reference,
        test_mode = FALSE,
        gene_cutoff = 1.25e-4,
        fc_cutoff = 0.5,
        gene_cutoff_reg = 2e-04,
        fc_cutoff_reg = 0.75,
        UMI_min = 1,
        UMI_max = 2e+06,
        counts_MIN = 10,
        UMI_min_sigma = 300,
        class_df = NULL,
        CELL_MIN_INSTANCE = 1,
        cell_type_names = NULL,
        MAX_MULTI_TYPES = opt$multi_types,
        keep_reference = F,
        #cell_type_info = NULL,
        CONFIDENCE_THRESHOLD = 10,
        DOUBLET_THRESHOLD = 25,
        max_cores = opt$ncores)
      # ==================================================================================================================
      futile.logger::flog.info("step5:Runs the RCTD pipeline on a 'RCTD' object")
      myRCTD <- spacexr::run.RCTD(RCTD_ob, doublet_mode = opt$mode)
      saveRDS(myRCTD, glue::glue("{output_dir}/myRCTD.rds"))
    }else {
      myRCTD <- readRDS(glue::glue("{output_dir}/myRCTD.rds"))
    }
    # ==================================================================================================================
    futile.logger::flog.info("step6:Tidy RCTD results, normalize the cell type proportions to sum to 1.")
    norm_weights <- rctd4weight(myRCTD, rctd_mode = opt$mode, conf = TRUE)
    out_celltype <- norm_weights %>%
      as.data.frame() %>%
      dplyr::select(sort(norm_weights %>% colnames)) %>%
      tibble::rownames_to_column(var = "barcodes") %>%
      data.frame()
    write.table(
      out_celltype,
      file = glue::glue("{output_dir}/all-celltype.xls"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE)
    data_ob <- subset(data_ob, cells = out_celltype$barcodes)
    futile.logger::flog.info("step7::save results to seurat object's mis part.") #=======================================
    #SeuratObject::Misc(data_ob, RCTD_results) <- out_celltype
    data_ob@misc$RCTD_results <- out_celltype
    ## 获得第一细胞类型及第二细胞类型
    maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
    metadata <- data_ob@misc$RCTD_results %>%
      tibble::column_to_rownames("barcodes") %>%
      dplyr::mutate(
        top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
        top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
      )
    data_ob <- Seurat::AddMetaData(data_ob, metadata = metadata)
    ## 保存seurat对象
    if (tolower(opt$outformat) == "h5seurat" & as.logical(opt$update)) {
      SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob)
    } else {
      saveRDS(
        data_ob,
        file = glue::glue("{output_dir}/seurat_RCTD.rds")
      )
    }
    # ==================================================================================================================
    ## save session information
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}
