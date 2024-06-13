# create the parsers for subcommand metacell
sub_metacell <- subparsers$add_parser("metacell", help = "produce cell count matrix of metacell")

args <- commandArgs(TRUE)
if ( "metacell"  %in% args ){
  opt<-intial_setting()
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = assays,
                                 data.use = dataslots)
  if (!is.null(opt$features)) {
    genes <- read.table(opt$features, sep = "\t", header = T)
    desired_features <- as.vector(genes$gene)
  }else {
    desired_features <- NULL
  }
  if (!is.null(opt$avgby)) {
    avgby <- opt$avgby
  }else {
    avgby <- "clusters"
  }
  # data_ob = SetIdent( data_ob, value = avgby)
  Seurat::Idents(data_ob) <- avgby
  avg_mtx <- Seurat::AverageExpression(data_ob, assays = opt$assay, features = desired_features, slot = opt$slot)
  write.table(avg_mtx, file.path(output_dir, paste0("average_expression_matrix_by_", avgby, ".xls")),
              sep = "\t", col.names = T, row.names = T)
  # output metadata either
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}