# Title     : TODO
# Objective : TODO
# Created by: xiufeng.yang
# Created on: 2022/9/14
#######################################################################################################################
docstring <- "example1:\\n\\
sctool  -i  st.rds   -f rds  -o results  -d rds   --assay SCT --dataslot counts  CARD  --refmarker markerlist.tsv \\n\\
example2:\\n\\
sctool  -i  st.rds   -f rds  -o results  -d rds   --assay  SCT  --dataslot counts  CARD  --refexp sc.rds --refassay RNA --refcelltype  celltype --refsample sample  --refassay  RNA"
sub_CARD <- subparsers$add_parser(
  "CARD",
  description = docstring,
  formatter_class = "argparse.RawTextHelpFormatter",
  # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
  argument_default = "True",
  help = "Using CARD to refer Spaital transcriptome cell type"
)
sub_CARD$add_argument(
  "--reobj",
  type = "character",
  default = NULL,
  help = paste0("[REQUIRED]the customized reference expression matrix in seurat format, which may come from the",
                " quantification results of microarray, Bulk RNA-seq and scRNA-seq sequencing."))
sub_CARD$add_argument(
  "--refassay",
  type = "character",
  default = NULL,
  help = paste0("the ref assay in data object to use. When it comes to multimodal assay, this is the assay used to",
                "initialize the object, all the other assays will merged into it. "))
sub_CARD$add_argument(
  "--refcelltype",
  type = "character",
  help = "[REQUIRED]the cell type annotation column id in the reference expression matrix."
)
sub_CARD$add_argument(
  "--refsample",
  type = "character",
  default = "sampleid",
  help = "[REQUIRED]the  sample column id in the reference expression matrix."
)
sub_CARD$add_argument(
  "--refmarker",
  default = NULL,
  help = "the  refmarker  file  for reference data to use for the this run"
)
# === Subcmd: CARD, deconvolute the cell type composition for spatial transcriptomics  using scRNA-seq data ========
args <- commandArgs(TRUE)
if ("CARD" %in% args) {
  opt <- intial_setting()
  if (opt$sub_name == "CARD") {
    set.seed(100)
    # ====================================================================================================
    futile.logger::flog.info("step1:load query object data")
    suppressMessages(data_ob <- OESingleCell::ReadX(
      input = opt$input,
      informat = opt$informat,
      assays = assays[1],
      data.use = dataslots,
      verbose = F
    ))
    Seurat::DefaultAssay(data_ob) <- assays[1]
    cellmeta <- OESingleCell::colData(data_ob)
    if (!is.null(Seurat::Images(data_ob))) {
      unuse_images <- Seurat::Images(data_ob)[!Seurat::Images(data_ob) %in% (data_ob@meta.data$sampleid %>% unique)]
      if (length(unuse_images) > 0) { data_ob@images[unuse_images] <- NULL }
    }
    # ====================================================================================================
    if (!is.null(opt$refexp)) {
      futile.logger::flog.info("step2:load reference object data")
      suppressMessages(ref_ob <- OESingleCell::ReadX(
        input = opt$refexp,
        informat = "rds",
        assays = opt$refassay,
        data.use = "counts",
        verbose = F
      ))

      # the input is a seurat object which may contain more than one sample
      SeuratObject::DefaultAssay(ref_ob) <- opt$refassay
      ref_ob <- SeuratObject::SetIdent(ref_ob, value = opt$refcelltype)
      # ==============================================================================================================
      futile.logger::flog.info("step3:Run CARD  object from a scRNA-seq reference Seurat object and a SpatialRNA object")
      images <- Seurat::Images(data_ob)
      celltype_out <- list()
      for (slice in images) {
        st_ob_sub <- subset(data_ob, subset = sampleid == as.character(slice))
        if (!is.null(Seurat::Images(st_ob_sub))) {
          unuse_images <- Seurat::Images(st_ob_sub)[!Seurat::Images(st_ob_sub) %in% (st_ob_sub@meta.data$sampleid %>% unique)]
          if (length(unuse_images) > 0) { st_ob_sub@images[unuse_images] <- NULL }
        }
        spatial_count <- OESingleCell::GetAssayData(st_ob_sub, assay = assays[1], slot = "counts")
        spatial_location <- Seurat::GetTissueCoordinates(st_ob_sub, scale = "lowres", cols = c("imagerow", "imagecol"))
        colnames(spatial_location) <- c("x", "y")
        OESingleCell::GetAssayData(ref_ob, assay = opt$refassay, slot = "counts")
        sc_meta <- Seurat::FetchData(ref_ob, c(opt$refcelltype, opt$refsample))
        CARD_obj <- CARD::createCARDObject(sc_count = sc_count,
                                           sc_meta = sc_meta,
                                           spatial_count = spatial_count,
                                           spatial_location = spatial_location,
                                           ct.varname = opt$refcelltype,
                                           ct.select = opt$refselect,
                                           sample.varname = opt$refsample,
                                           minCountGene = 100,
                                           minCountSpot = 5)
        CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)
        celltype_out[[slice]] <- CARD_obj@Proportion_CARD %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "barcodes")
      }
    }
    #====================================================================================================================
    if (!is.null(opt$refmarker)) {
      futile.logger::flog.info("step2:load maker list data")
      markerList <- readr::read_tsv(opt$refmarker)
      markerList <- lapply(markerList, function(x) x[!is.na(x)])
      futile.logger::flog.info("step3:Run CARD  object from a scRNA-seq reference Seurat object and a SpatialRNA object")
      images <- Seurat::Images(data_ob)

      celltype_out <- list()
      for (slice in images) {
        st_ob_sub <- base::subset(data_ob, subset = sampleid == as.character(slice))
        if (!is.null(Seurat::Images(st_ob_sub))) {
          unuse_images <- Seurat::Images(st_ob_sub)[!Seurat::Images(st_ob_sub) %in% (st_ob_sub@meta.data$sampleid %>% unique)]
          if (length(unuse_images) > 0) { st_ob_sub@images[unuse_images] <- NULL }
        }
        spatial_count <- OESingleCell::GetAssayData(st_ob_sub, assay = assays[1], slot = "counts")
        spatial_location <- Seurat::GetTissueCoordinates(st_ob_sub, scale = "lowres", cols = c("imagerow", "imagecol"))
        colnames(spatial_location) <- c("x", "y")
        CARDfree_obj <- CARD::createCARDfreeObject(markerList = markerList,
                                                   spatial_count = spatial_count,
                                                   spatial_location = spatial_location,
                                                   minCountGene = 100,
                                                   minCountSpot = 5)
        CARDfree_obj <- CARD::CARD_refFree(CARDfree_object = CARDfree_obj)
        colnames(CARDfree_obj@Proportion_CARD) <- paste0("CT", colnames(CARDfree_obj@Proportion_CARD))
        celltype_out[[slice]] <- CARD_obj@Proportion_CARD %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "barcodes")
      }
    }

    # ==================================================================================================================
    futile.logger::flog.info("step5: save the cell type proportions in misc@CARD_results")
    out_celltype <- dplyr::bind_rows(celltype_out) %>%
      data.frame()
    write.table(
      out_celltype,
      file = glue::glue("{output_dir}/all-celltype.xls"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    data_ob <- subset(data_ob, cells = out_celltype$barcodes)
    ## save results to seurat object's mis part
    data_ob@misc$CARD_results <- out_celltype
    maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
    metadata <- data_ob@misc$CARD_results %>%
      tibble::column_to_rownames("barcodes") %>%
      dplyr::mutate(
        top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
        top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
      )
    data_ob <- Seurat::AddMetaData(data_ob, metadata = metadata)
    if (tolower(opt$outformat) == "h5seurat" & as.logical(opt$update)) {
      SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob)
    } else {
      OESingleCell::SaveX(
        data_ob,
        output = opt$output,
        update = FALSE,
        outformat = opt$outformat,
        prefix = opt$prefix
      )
    }
    #==============================================================================================================
    ## save session information
    write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
    quit()
  }
}
