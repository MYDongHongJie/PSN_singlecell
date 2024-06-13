sub_jackstraw = subparsers$add_parser("jackstraw", help = "the jackstraw plot")
sub_jackstraw$add_argument("--components", type = "integer", default = 20,
         help = "[OPTIONAL]the dimensions used for visualization or calibration.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "jackstraw" %in% args ){
  opt<-intial_setting()
  # read the specified assay and data slot in data object into memory
  # for visualization usually only "data" slot is loaded,which can save a lot of time and memory
  data_ob = OESingleCell::ReadX(
             input = opt$input, informat = opt$informat,
             assays = assays, data.use = dataslots, # "data","scale.data"
             reductions = opt$reduct, graphs = FALSE, images = FALSE, # no graph object needed here
             verbose = FALSE)

  #get subset of cells used col.name = col.name = or visualization if necessay
  if ( !is.null(opt$predicate) ){
    df = colData(data_ob)
    desired_cells= subset(df, eval( parse(text=opt$predicate)))
    data_ob = data_ob[, desired_cells]
  }

  components2use = min( opt$components, ncol(Embeddings(data_ob[[opt$reduct]])))
  data_ob <- Seurat::JackStraw(data_ob , num.replicate = 100,
                         reduction = opt$reduct, dims = components2use)
  data_ob <- Seurat::ScoreJackStraw(data_ob, dims = components2use, reduction = opt$reduct)
  ggjack = Seurat::JackStrawPlot(object = data_ob, dims = 1:components2use)
  OESingleCell::save_ggplots(file.path(output_dir,"jackstraw_plot"),plot=ggjack, dpi = 600)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}