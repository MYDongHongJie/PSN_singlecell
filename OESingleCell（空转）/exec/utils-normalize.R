# normalize the raw data before downstream data analysis
sub_normalize <- subparsers$add_parser("normalize", help = "normalize the raw data.")
sub_normalize$add_argument("-m", "--nmethod", type = "character", default = "sct",
              help = "the normalization method to choose from: sct, LogNormalize, scran, cr,clr for single cell data.[default: %(default)s] ")
sub_normalize$add_argument("-r", "--vars2regress", type = "character", default = "nCount_RNA,percent.mito",
                help = "comma separated list of sources of variation to regress out from the data. The other options can be percent.ribo, CC.Difference etc.[default: %(default)s]")
sub_normalize$add_argument( "-s", "--scale-factor", default = 10000, type = 'integer',
                help = "Sets the scale factor for cell-level normalization.[default: %(default)s]" )
sub_normalize$add_argument("--sct-method", default = "glmGamPoi", type = 'character',
                help = "the distribution type used to fit generalized linear models of each gene count against the sequencing depths, choices:poisson, glmGamPoi(faster).[default: %(default)s]")
sub_normalize$add_argument( "--margin", default = NULL, type = 'integer',
                help = "If performing CLR normalization, normalize across features (1) or cells (2).[default: %(default)s] " )
sub_normalize$add_argument("--model", type = "character", default = "linear",
          help = "Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'[default: %(default)s] ")
sub_normalize$add_argument("--use-umi", type = "character", default = "FALSE",
          help = "Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'[default: %(default)s] ")
sub_normalize$add_argument("--block-size", type = "integer", default = 1000,
          help = "Default size for number of features to scale at in a single computation. Increasing block.size may speed up calculations but at an additional memory cost.[default: %(default)s] ")

# =============== Subcmd: normalize, normalize the data ============
args <- commandArgs(TRUE)
if ( "normalize"  %in% args ){
  opt<-intial_setting()
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots, # counts only
                                 verbose =F )
  if ( !is.null(opt$vars2regress) ){
    vars2regress <- unlist(strsplit(opt$vars2regress, ",", perl =T))
  }else{
    vars2regress <- NULL
  }

  if ( tolower(opt$nmethod) == "sct" ){
    data_ob <- Seurat::SCTransform(data_ob,
                                   method = opt$sct_method,
                                   vars.to.regress = vars2regress,
                                   verbose = FALSE,
                                   return.only.var.genes = FALSE)
    SeuratDisk::UpdateH5Seurat(file = opt$input, data_ob, assay = "SCT", verbose = F)
    # OESingleCell::SaveX(data_ob, output = opt$input, assay = "SCT",
    #                   outformat = opt$outformat, update = T)
  }else{
    data_ob <- Seurat::NormalizeData(data_ob,
                                     normalization.method = opt$nmethod,
                                     scale.factor = opt$scale_factor,
                                     margin = opt$margin,
                                     block.size = opt$block_size,
                                     verbose = FALSE)
    data_ob <- Seurat::ScaleData(data_ob,
                                 assay = assays[1],
                                 vars.to.regress = vars2regress,
                                 split.by = NULL,
                                 model.use = opt$model,
                                 use.umi = as.logical(opt$use_umi),
                                 scale.max = 10,
                                 block.size = opt$block_size,
                                 min.cells.to.block = 3000,
                                 verbose = TRUE )
    OESingleCell::SaveX(data_ob,
                        output = opt$input,
                        assay = assays[1],
                        outformat = opt$outformat,
                        update = T)
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}