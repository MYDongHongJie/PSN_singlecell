sub_elbow = subparsers$add_parser("elbow", help = "the elbow plot of components variance.")

args <- commandArgs(TRUE)
if ( "elbow" %in% args ){
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

  # to visulaize the metadata of each cell on the tsne plot
    ggelb = OESingleCell::ElbowPlot2(data_ob)
    OESingleCell::save_ggplots(file.path(output_dir,"Elbow_plot"),plot=ggelb, dpi = 600)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}