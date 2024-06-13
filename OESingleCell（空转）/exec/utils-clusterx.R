sub_clusterx <- subparsers$add_parser("clusterx", help = "universal interface of clustering methods")
sub_clusterx$add_argument("--algorithm",
                          type="character",
                          default = "snn",
                          choices=c("snn","louvain","slm","leiden"),
                          help="the only supported clustering algorithm currently:snn,louvain,leiden. For other methods, use subcommand 'clusterx'.[default: %(default)s]" )
sub_clusterx$add_argument("-r", "--resolution",
                          type = "double",
                          default = 0.4,
                          help = "vaule used to set the resolution of cluster distiguish, use a value above(below)1.0 if you want to obtain a larger(smaller) number of communities.[default: %(default)s]")
sub_clusterx$add_argument( "--graphid",
                           type = "character",
                           default = "RNA_umap_0.4",
                           help = "the graph name used to find communities.[default: %(default)s]")
sub_clusterx$add_argument("--palette",
                          type = "character",
                          default = "blindless",
                          choices=c("blindless","cold","glasbey","ditto","alphabet","alphabet2","colx22","cyclic","tableau20", "paired", "purplegray12"),
                          help="the discrete color schema for each cell cluster, blindless as default. Choices:blindless:69,cold:32,glasbey:32,ditto:24,alphabet:26,alphabet2:26,colx22:22,cyclic:20,tableau20:20, paired:12, purplegray12:12[default: %(default)s]")
sub_clusterx$add_argument("-s", "--pointsize",
                          type = "double",
                          default = 0.5,
                          help = "[OPTIONAL]the point size in the plot.[default: %(default)s]")

# =============== Subcmd: clusterx, run data clustering =========
args <- commandArgs(TRUE)
if ( "clusterx"  %in% args ){
  opt<-intial_setting()
     # read the specified assay and data slot in data object into memory
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots,
                                 verbose = F)
  # get the subset of cells used for visualization if necessay
  # subset cells on the expression matrix using logical exprssion on features
  if ( !is.null(opt$predicate) ){
      df <- OESingleCell::colData(data_ob)
      desired_cells <- subset(df, eval(parse(text=opt$predicate)))
      data_ob <- data_ob[, rownames(desired_cells)]
  }

  if ( !is.null(opt$graphid) ){
    stopifnot( opt$graphid %in% SeuratObject::Graphs(data_ob) )
  }

  #Clustering with the reduction results using different clustering methods
  clustering.alg <- c("snn" = 1, "louvain" = 2, "slm" = 3, "leidn" = 4)
  message(glue::glue("Beginning to group the cells using algorithm {opt$algorithm}!" ) )
  res <-stringr::str_split(opt$resolution,",")[[1]]
  if(length( res)>1){res<- seq(res[1], res[2] , res[3]) }
  data_ob <- Seurat::FindClusters(object = data_ob,
                                  resolution = res,
                                  graph.name = opt$graphid ,
                                  algorithm = clustering.alg[opt$algorithm],
                                  verbose = F)
  message("Clustering Finished!")

  #visulize by clusttree
  plot<- clustree(data_oba@meta.data,
                  prefix = glue::glue("{opt$graphid}.res") ,
                  layout = "sugiyama")
                  # node_colour = "percent.mt",
                  # node_colour_aggr = "mean")

  OESingleCell::save_ggplots(glue::glue("{output_dir}/{opt$graphid}_{opt$algorithm}_clusttree_{opt$resolution}_plot"),
                           plot = plot,
                           width= 8 ,
                           height=length(res)*4,
                           dpi = 1000)
  if ( as.logical(opt$update) ){
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