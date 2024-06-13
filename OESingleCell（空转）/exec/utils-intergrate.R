# create the parsers for subcommand integrate
sub_integrate = subparsers$add_parser("integrate", help = "integrate multiple data across different conditions or sources.")
sub_integrate$add_argument("--method", type = 'character', default = "seurat",
          help = "the method used to integrate single cell data, choices:seurat, liger, harmony.[default: %(default)s] ")
sub_integrate$add_argument("--reduction", type = 'character', default = "cca",
          help = "cca,rpca,rlsi.[default: %(default)s] ")



# ================ Subcmd: integrate, integrate the data from different source and condditions ========
args <- commandArgs(TRUE)
if ( "intergrate"  %in% args ){
  opt<-intial_setting()
  data_ob = OESingleCell::ReadX(input = opt$input,
                                informat = opt$informat,
                                assays = assays[1],
                                data.use = dataslots,
                                verbose = F)

  switch(opt$method,
         seurat = {
           data_ob = OESingleCell::RunIntegrateX(data_ob, split.by = opt$splitby ,reduction='cca')
         },
         liger = {
           data_ob = ScaleData(data_ob, split.by = opt$splitby, do.center = FALSE)
           data_ob = SeuratWrappers::RunOptimizeALS(data_ob, k = 20, lambda = 5, split.by = opt$splitby)
           data_ob = SeuratWrappers::RunQuantileNorm(data_ob, split.by = opt$splitby)
         }
  )

  if ( tolower(opt$outformat) == "h5seurat" & as.logical(opt$update) ){
    SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob, assay = "integrated" )
  }else{
    OESingleCell::SaveX(data_ob, output = opt$output,update = FALSE,
                        outformat = opt$outformat, prefix = opt$prefix)
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}