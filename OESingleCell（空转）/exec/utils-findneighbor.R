sub_findneighbor = subparsers$add_parser("findneighbor", help = "only find the clusters in the input data")
sub_findneighbor$add_argument( "-k", "--k_param", default = 20, type = 'integer',
    help = "Defines k for the k-nearest neighbor algorithm.[default: %(default)s] " )
sub_findneighbor$add_argument( "--nn_method", default = 'rann', type = 'character',
    help = "Method for nearest neighbor finding. Options include: rann (default), annoy(faster).[default: %(default)s] " )
sub_findneighbor$add_argument( "-a", "--dist_metric", type = 'character', default = "euclidean",
    help = "Distance metric. Options include: euclidean (default), cosine, manhattan, and hamming.[default: %(default)s] " )
sub_findneighbor$add_argument( "--graph_name", type = 'character', default = NULL,
    help = "Name of graph to use for the clustering algorithm.[default: %(default)s]" )
sub_findneighbor$add_argument( "--force_recalc", type = "character", default = "FALSE",
    help = "Force recalculation of SNN.[default: %(default)s] " )
sub_findneighbor$add_argument( "--reduction", type = 'character', default = "pca",
    help = "one or more reductions to use as input for building the (S)NN for mono-modal or multi-modual assays.[default: %(default)s] " )
# sub_findneighbor$add_argument( "-f", "--features", default = NULL, type = 'character',
#     help = "Comma-separated list of genes to use for building SNN. Alternatively, text file with one gene per line." )
# sub_findneighbor$add_argument( "--reduction", default = NULL, type = 'character',
#     help = "one or more Reduction  to use as input for building the SNN" )
sub_findneighbor$add_argument( "--dims", default = "1:15", type = 'character',
    help = "Dimensions of reduction to use as input. A comma-separated list of the dimensions to use in construction of the SNN graph.[default: %(default)s] " )


# =============== Subcmd: findneighbor, SNN or wSNN calculation =========
args <- commandArgs(TRUE)
if ( "findneighbor"  %in% args ){
  opt<-intial_setting()
  # read the specified assay and data slot in data object into memory
  data_ob = OESingleCell::ReadX(input = opt$input,
                                informat = opt$informat,
                                assays = assays,
                                data.use = dataslots, # data slot is enough
                                verbose = F)

  reductions = unlist(strsplit(opt$reduction, ",", perl = T))
  dims = lapply(unlist(strsplit(opt$dims, ",", perl = T)), function(x){eval( parse(text = x))})
  if ( length(assays) > 1 ){ # multimodal assays
    data_ob = Seurat::FindMultiModalNeighbors(data_ob,
                                              k.nn = opt$k_param,
                                              reduction.list = as.list(reductions),
                                              dims.list = dims,
                                              verbose = F)
  }else{ # monomodal assays
    data_ob = Seurat::FindNeighbors( data_ob,
                                   k.nn = opt$k_param,
                                   reduction = reductions[1],
                                   dims = dims[1],
                                   force.recalc = T, verbose = F)

  }

  if ( as.logical(opt$update) ){
    SeuratDisk::UpdateH5Seurat(file = opt$input,
                               object = data_ob,
                               reduction = reductions,
                               verbose = F)
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