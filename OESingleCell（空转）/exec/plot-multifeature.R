args <- commandArgs(TRUE)
if ( "multifeature" %in% args ){
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

  if ( !is.null(opt$feature1) & !is.null(opt$feature2) ){
    ggscater = FeatureScatter(data_ob, feature1 = opt$feature1, slot = opt$slot,
                              feature2 = opt$feature2, group.by = opt$groupby )
    OESingleCell::save_ggplots(file.path(output_dir, paste0(feature1,"_VS_",feature2, "_visulazation")),
                               dpi = 600, plot=ggscater)
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}