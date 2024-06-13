
# create the parsers for subcommand transfer
sub_transfer <- subparsers$add_parser("transfer", help = "integrate multiple data objects for cell annotation transfer for both same assays or alternative assays .")
sub_transfer$add_argument("--algrithms", type = "character", default = "tranfer",
          help = "the algrithms used integrate the different data objects, choices: transfer, umap_project.[default: %(default)s]")
sub_transfer$add_argument("--refdata", type = "character", default = NULL,
          help = "the reference data object with cell annotations.")
sub_transfer$add_argument("--nmethod", type = "character", default = "LogNormalize",
          help = "[OPTIONAL]Name of normalization method used: LogNormalize or SCT.[default: %(default)s]")
sub_transfer$add_argument("--refassay", type = "character", default = NULL,
          help = "Name of the Assay to use from reference data object.[default: %(default)s] ")
sub_transfer$add_argument("--reduct", type = "character", default = "pcaproject",
          help = "Dimensional reduction to perform when finding anchors.Options are: pcaproject: Project the PCA from the reference onto the query. when reference and query datasets are from scRNA-seq; cca: Run a CCA on the reference and query when they are from different assays.[default: %(default)s] ")
sub_transfer$add_argument("--refreduct", type = "character", default = NULL,
          help = "Name of dimensional reduction to use from the reference if running the pcaproject workflow. Optionally enables reuse of precomputed reference dimensional reduction. If NULL (default), use a PCA computed on the reference object.[default: %(default)s] ")
sub_transfer$add_argument("--project.query", type = "character", default = "FALSE",
          help = "[OPTIONAL] Project the PCA from the query dataset onto the reference. Used only in rare cases where the query dataset has a much larger cell number, but the reference dataset has a unique assay for transfer. In this case, the default features will be set to the variable features of the query object that are also present in the reference.[default: %(default)s] ")
sub_transfer$add_argument("--refid", type = "character", default = "idents",
          help = "the reference celltype annotation column id in reference data object.[default: %(default)s] ")
sub_transfer$add_argument("--ndims", type = "integer", default = 30,
          help = "[OPTIONAL]the number of dimensions to use from the reduction to specify the neighbor search space.[default: %(default)s]")
sub_transfer$add_argument("--kanchors", type = "integer", default = 5,
          help = "[OPTIOANL]the number of neighbors to use when finding anchors.[default: %(default)s]")
sub_transfer$add_argument("--kfilters", type = "integer", default = 200,
          help = "[OPTIONAL]the number of neighbors to use when filtering anchors. Set to NA to turn off filtering.[default: %(default)s]")


# ================ Subcmd: transfer, transfer cell annotation from reference data using data integration ========
args <- commandArgs(TRUE)
if ( "transfer"  %in% args ){
  opt<-intial_setting()
  # make sure data, scale.data are availiable for both query and reference data object
  # make sure standard preprocessing have been performed on each object
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays[1],
                                 data.use = dataslots,
                                 verbose = F)
  ref_ob <- OESingleCell::ReadX(input = opt$refdata,
                                informat = opt$informat,
                                assays = opt$refassay,
                                data.use = dataslots,
                                verbose = F)
  if ( !is.null(opt$predicate) ){
      df <- OESingleCell::colData(data_ob)
      desired_cells <- subset(df, eval(parse(text=opt$predicate)))
      data_ob <- data_ob[, rownames(desired_cells)]
  }

   data_ob <- Seurat::FindVariableFeatures(data_ob, selection.method = "vst")
   ref_ob  <- Seurat::FindVariableFeatures(ref_ob, selection.method = "vst")
   # select features that are repeatedly variable across datasets for integration
   integrate_features <- Seurat::SelectIntegrationFeatures(list(data_ob, ref_ob), nfeatures = 2000)
  # Find a set of anchors between a reference and query object.
  # These anchors can later be used to transfer data from the reference to query object
  anchors <- Seurat::FindTransferAnchors(reference = ref_ob,
                                         query = data_ob,
                                         features = integrate_features,
                                         normalization.method = opt$nmethod,
                                         reference.assay = opt$refassay,
                                         query.assay = assays[1],
                                         reduction = opt$reduct,
                                         reference.reduction = opt$refreduct,
                                         project.query = as.logical(opt$project.query),
                                         dims = 1:opt$ndims,
                                         k.anchor = opt$kanchors,
                                         k.filter = opt$kfilters )
  data_ob <- switch(opt$algrithms,
                    transfer = {
           predictions <- Seurat::TransferData(anchors, refdata = OESingleCell::colData(ref_ob)[[opt$refid]])
           data_ob <- OESingleCell::AddMetaData(data_ob, metadata = predictions)
           write.table(tibble::as_tibble(predictions, rownames = "cellbarcode"),
                       file.path(output_dir,"cell_type_label_transfer_result.xls"),
                       sep="\t",
                       col.names=T,
                       row.names=F)
         },
                    umap_project = {
           if ( !"umap" %in% OESingleCell::Reductions(ref_ob) ) {
             stop("NO Required UMAP embeddings FOUND in referece data!Please run the umap reduction first!")
           }
           data_ob <- Seurat::MapQuery(anchorset = anchors,
                                       reference = ref_ob,
                                       query = data_ob,
                                       refdata = list(celltype = opt$refid),
                                       reference.reduction = "pca",
                                       reduction.model =  "umap")
         },
                    symphony = {
         },
                    stop("NO specified cell annotation transfer method supported!")
  )

  if ( tolower(opt$outformat) == "h5seurat" & as.logical(opt$update) ){
    SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob )
  }else{
    OESingleCell::SaveX(data_ob,
                        output = opt$output,
                        update = FALSE,
                        outformat = opt$outformat,
                        prefix = opt$prefix)
  }

  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}