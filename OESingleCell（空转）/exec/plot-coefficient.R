# coefficient heatmap
sub_coefficient = subparsers$add_parser("coefficient",
                                        help = "export the gene expression in different group sepecified by the user.")
sub_coefficient$add_argument("-g", "--groupby", type="character", default="clusters",
             help="[OPTIONAL]The grouppinig variable in the metadata for separate the cells to calculate the correlation.[default: clusters]")
# sub_coefficient$add_argument("--reduct", type="character", default="umap",
#              help = "[OPTIONAL]the reduction results used to embedding the coefficient heatmap.[default: %(default)s]")

# =============== Subcmd: coefficient heatmap =========

args <- commandArgs(TRUE)
if ( "coefficient" %in% args ){
  opt<-intial_setting()
  groupby = opt$groupby
  # read the specified assay and data slot in data object into memory
  data_ob = OESingleCell::ReadX(input = opt$input,
                                informat = opt$informat,
                                assays = assays,
                                data.use = dataslots, # "data"
                                reductions = opt$reduct, # only the used reduction results needed
                                graphs = FALSE, # no graph object needed here
                                images = FALSE, verbose = FALSE)

  # subset test data object
  if ( !is.null(opt$predicate) ){
      df = OESingleCell::colData(data_ob)
      desired_cells= subset(df, eval( parse(text=opt$predicate)))
      data_ob = data_ob[, rownames(desired_cells)]
  }

  # export the gene expression matrix
  groupby_data = vector()
  for (i in names(table(data_ob[[groupby]]))  ) {
      sub_ob = data_ob[, rownames(data_ob@meta.data[which(data_ob@meta.data$clusters %in% i),])]
      normalized_data = as.matrix(sub_ob[[opt$assay]]@data)
      meta.data = sub_ob@meta.data %>% tibble::rownames_to_column(var = "id")
      groupby_data = cbind(groupby_data,rowMeans(normalized_data))
  }
  colnames(groupby_data) = names(table(data_ob[[groupby]]))
  data = tibble::rownames_to_column(as.data.frame(groupby_data),var="GeneID")
  write.table(data, file.path(output_dir,paste0("normalized_data_groupby_",groupby,".xls")),quote = F, row.names = F, sep = "\t")

  # calculate the correlation and visulize by heatmap
  colnames(groupby_data) = gsub('^',paste0(groupby,"_"),colnames(groupby_data))
  matrix<-cor(groupby_data,method="pearson")
  wid<-5+1.5*log2(length(colnames(data)))
  hig<-5+1.5*log2(length(colnames(data)))
  coefficient = pheatmap::pheatmap(matrix,
                                   display_numbers = T,
                                   border_color = "white",
                                   scale = "none",
                                   fontsize_number=(10.0+0.0001*log2(length(colnames(data)))),
                                   number_format = "%.4f",
                                   fontsize_row = (10.0+0.0001*log2(length(colnames(data)))),
                                   fontsize_col = (10.0+0.0001*log2(length(colnames(data)))),
                                   number_color="black",
                                   angle_col=45)
  OESingleCell::save_ggplots(file.path(output_dir,"coefficient_heatmap"),
                          plot = coefficient, dpi = 1000 ,
                          height = hig, width = wid, limitsize = F)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}