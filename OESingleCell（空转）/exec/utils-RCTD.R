rename.Seurat <- function(
  object,
  ...
) {
    pd <- as.data.frame(object[[]])
    pd <- dplyr::rename(pd, ...)
    rownames(pd) <- colnames(object)
    slot(object, "meta.data") <- pd
    return(object)
}

# min weight to be considered a singlet as a function of nUMI
UMI_cutoff <- function(nUMI) {
    return(pmax(0.25, 2 - log(nUMI, 2) / 5))
}

#######################################################################################################################
docstring <- " example1:\\n\\n\\
sctool  -i  st.rds   -f rds  -o results  -d rds   --assay Spatial  --dataslot counts  RCTD --doublet_mode FALSE --refexp sc.rds  --refcelltype  celltype --refassay  RNA"
sub_RCTD <- subparsers$add_parser(
    "RCTD",
    description = docstring,
    formatter_class = "argparse.RawTextHelpFormatter",
    # formatter_class= '"lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=20,width=150)",
    argument_default = "True",
    help = "Using RCTD to refer Spaital transcriptome cell type"
)
sub_RCTD$add_argument(
    "--refexp",
    type = "character",
    help = "[REQUIRED]the customized reference expression matrix in seurat format, which may come from the quantification results of microarray, Bulk RNA-seq and scRNA-seq sequencing."
)
sub_RCTD$add_argument(
    "--refcelltype",
    type = "character",
    help = "[REQUIRED]the cell type annotation column id in the reference expression matrix."
)
sub_RCTD$add_argument(
    "--refassay",
    type = "character",
    default = "RNA",
    help = "[REQUIRED]the default assay for reference data to use for the this run.[default: %(default)s] "
)
sub_RCTD$add_argument(
    "-m",
    "--doublet_mode",
    default = "FALSE",
    help = "doublet_mode:TRUR or FALSE, for 10x visium data,default is FALSE [default:\"%(default)s\"]"
)
# ================ Subcmd: st_deconv, deconvolute the cell type composition for spatial transcriptomics or bulk transcriptomics using scRNA-seq data ========
args <- commandArgs(TRUE)
if ("RCTD" %in% args) {
    opt <- intial_setting()
    if (opt$sub_name == "RCTD") {
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
        # ====================================================================================================
        futile.logger::flog.info("step2:load reference object data")
        suppressMessages(ref_ob <- OESingleCell::ReadX(
            input = opt$refexp,
            informat = "rds",
            assays = opt$refassay,
            data.use = "data,counts",
            verbose = F
        ))
        if (!is.null(opt$refexp)) {
            # the input is a seurat object which may contain more than one sample
            SeuratObject::DefaultAssay(ref_ob) <- opt$refassay
            ref_ob <- SeuratObject::SetIdent(ref_ob, value = opt$refcelltype)
        }
        # ==============================================================================================================
        futile.logger::flog.info("step3:reformart reference data object")

        if( ! "nCount_RNA" %in% colnames(ref_ob@meta.data )  ){
            nCount_RNA <- SeuratObject::colSums(x = ref_ob, slot = "counts") %>%as.data.frame %>%dplyr::rename(nCount_RNA='.')
            ref_ob<-Seurat::AddMetaData(ref_ob,metadata = nCount_RNA)
        }

        reference <- rename.Seurat(ref_ob, liger_ident_coarse = as.character(opt$refcelltype)) %>% rename.Seurat(nUMI = nCount_RNA)
        reference@meta.data$liger_ident_coarse <- forcats::fct_inorder(reference@meta.data$liger_ident_coarse)
        puck <- RCTD::seurat_to_SpatialRNA(data_ob, opt$assay)

        # ==============================================================================================================
        futile.logger::flog.info("step4:Creates an RCTD object from a scRNA-seq reference Seurat object and a SpatialRNA object")
        RCTD_ob <- RCTD::create.RCTD(
            puck,
            reference,
            test_mode = FALSE,
            gene_cutoff = 1.25e-4,
            fc_cutoff = 0.5,
            gene_cutoff_reg = 2e-04,
            fc_cutoff_reg = 0.75,
            UMI_min = 1,
            UMI_max = 2e+06,
            class_df = NULL,
            CELL_MIN_INSTANCE = 1,
            max_cores = opt$ncores
        )
        # ===============================================================================================================
        futile.logger::flog.info("step5:Runs the RCTD pipeline on a 'RCTD' object")
        myRCTD <- RCTD::run.RCTD(RCTD_ob, doublet_mode = opt$doublet_mode)
        # saveRDS(myRCTD, file = glue::glue("{output_dir}/RCTD.rds"))
        # saveRDS(RCTD_ob, file = glue::glue("{output_dir}/raw_RCTD.rds"))
        # RCTD_ob<-readRDS(glue::glue("{output_dir}/raw_RCTD.rds"))
        # myRCTD<-readRDS(glue::glue("{output_dir}/RCTD.rds"))
        results <- myRCTD@results
        futile.logger::flog.info("normalize the cell type proportions to sum to 1.")
        norm_weights <- sweep(
            results$weights %>% as.data.frame(),
            1,
            rowSums(results$weights %>% as.data.frame()),
            "/"
        )
        out_celltype <- norm_weights %>%
            as.data.frame() %>%
            dplyr::select(sort(norm_weights%>%colnames))%>%
            # dplyr::filter_all(dplyr::any_vars(. > UMI_cutoff(RCTD_ob@spatialRNA@nUMI ))) %>%
           # round(., 2) %>%
            tibble::rownames_to_column(var = "barcodes") %>%
            data.frame()
        write.table(
            out_celltype,
            file = glue::glue("{output_dir}/all-celltype.xls"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        data_ob <- subset(data_ob, cells = out_celltype$barcodes)
        # SeuratObject::Misc(data_ob, slot = "RCTD_results") <- out_celltype
        ## save results to seurat object's mis part
        data_ob@misc$RCTD_results <- out_celltype
        maxn <- function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
        metadata <- data_ob@misc$RCTD_results %>%
            tibble::column_to_rownames("barcodes") %>%
            dplyr::mutate(
                top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
                top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
            )
        data_ob <- Seurat::AddMetaData(data_ob, metadata = metadata)
        # cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
        # spatialRNA <- myRCTD@spatialRNA
        # futile.logger::flog.info("make the plots")
        # # Plots all weights for each cell type as in full_mode. (saved as 'results/cell_type_weights.pdf')
        # RCTD::plot_weights(cell_type_names,
        #                    spatialRNA,
        #                    output_dir,
        #                    norm_weights)
        # # Plots the confident weights for each cell type as in full_mode (saved as 'results/cell_type_weights_unthreshold.pdf')
        # RCTD::plot_weights_unthreshold(cell_type_names,
        #                                spatialRNA,
        #                                output_dir,
        #                                norm_weights)
        # # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 'results/cell_type_occur.pdf')
        # RCTD::plot_cond_occur(cell_type_names,
        #                       output_dir,
        #                       norm_weights,
        #                       spatialRNA)
        # doublets module:双细胞模式#=========================================================================================
        # if (opt$doublet_mode) {
        #     # makes a map of all cell types, (saved as 'results/all_cell_types.pdf')
        #     RCTD::plot_all_cell_types(
        #         results$results_df,
        #         spatialRNA@coords,
        #         cell_type_names,
        #         output_dir
        #     )
        #
        #     # Plots the weights for each cell type as in doublet_mode. (saved as 'results/cell_type_weights_doublets.pdf')
        #     RCTD::plot_weights_doublet(
        #         cell_type_names,
        #         spatialRNA,
        #         output_dir,
        #         results$weights_doublet,
        #         results$results_df
        #     )
        #
        #     # obtain a dataframe of only doublets
        #     doublets <- results$results_df[results$results_df$spot_class == "doublet_certain", ]
        #     # Plots all doublets in space (saved as 'results/all_doublets.pdf')
        #     RCTD::plot_doublets(
        #         spatialRNA,
        #         doublets,
        #         output_dir,
        #         cell_type_names
        #     )
        #
        #     # Plots all doublets in space for each cell type (saved as 'results/all_doublets_type.pdf')
        #     RCTD::plot_doublets_type(
        #         spatialRNA,
        #         doublets,
        #         output_dir,
        #         cell_type_names
        #     )
        #
        #     # a table of frequency of doublet pairs
        #     doub_occur <- table(doublets$second_type, doublets$first_type)
        #
        #     # Plots a stacked bar plot of doublet ocurrences (saved as 'results/doublet_stacked_bar.pdf')
        #     RCTD::plot_doub_occur_stack(
        #         doub_occur,
        #         output_dir,
        #         cell_type_names
        #     )
        #
        #     # get a SpatialRNA object that has single cell types, each with a spatial coordinate and RNA counts.
        #     puck_d <- RCTD::get_decomposed_data(
        #         results$results_df,
        #         myRCTD@internal_vars$gene_list_reg,
        #         spatialRNA,
        #         results$weights_doublet,
        #         myRCTD@cell_type_info$renorm
        #     )
        #}
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
        # ==============================================================================================================
        ## save session information
        write_session_info(log_dir, sub_name = parser$parse_args()$sub_name)
        quit()
    }
}
