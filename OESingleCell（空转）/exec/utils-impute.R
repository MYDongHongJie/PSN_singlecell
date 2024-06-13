# create the parsers for subcommand impute
sub_impute = subparsers$add_parser("impute", help = "impute the missing values in count matrix.")
sub_impute$add_argument("--method", type = "character", default = "magic",
          help = "the method used to impute the sparse matrix.[default: %(default)s]")
sub_impute$add_argument("--slot", type = "character", default = "counts",
          help = "the name of the matrix to impute of the specified assay.[default: %(default)s] ")


# =============== Subcmd: impute, impute the missing value to solve the dropouts =========
args <- commandArgs(TRUE)
if ( "impute"  %in% args ){
  opt<-intial_setting()
  # read the specified assay and data slot in data object into memory
  data_ob = OESingleCell::ReadX(input = opt$input,
                                informat = opt$informat,
                                reductions = opt$reduct,
                                assays = assays[1], # only RNA assay is valid
                                data.use = opt$slot, # counts slot is enough
                                verbose = F)

  dims = min(opt$ndims, ncol(OESingleCell::Embeddings(data_ob, reduction = opt$reduct)))
  assayx_pre = SeuratObject::Assays(data_ob)
  # MAGIC, DCA, ALRA, Enhance
  data_ob <- switch(tolower(opt$method),
      "magic" = {
        OESingleCell::RunMagic(data_ob, assay = assays[1], n.jobs = opt$ncores, npca = ndims)
      },
      "dca" = {
        OESingleCell::RunDCA(data_ob, assay = assays[1], reduce_lr = ndims)
      },
      "alra" = {
        OESingleCell::RunALRA(data_ob, assay = assays[1] )
      },
      "enhance" = {
        OESingleCell::Enhance(data_ob, assay = assay[1])
      },
      "atac_impute" = {
        OESingleCell::RunATACImpute(data_ob, assay = assay[1],
                                    reduction = opt$reduct, dims = opt$ndim,
                                    threads = opt$ncores)
      },
      "NO specifed method SUPPORT!"
    )

  if ( tolower(opt$outformat) == "h5seurat" & as.logical(opt$update) ){
    assayx = setdiff(SeuratObject::Assays(data_ob), assayx_pre)
    SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob, assay = assayx )
  }else{
    OESingleCell::SaveX(data_ob, output = opt$output,update = FALSE,
                        outformat = opt$outformat, prefix = opt$prefix)
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}