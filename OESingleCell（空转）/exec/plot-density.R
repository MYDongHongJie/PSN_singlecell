
sub_dense = subparsers$add_parser("density", help = "density plot with features in cell communities.")
sub_dense$add_argument("-l", "--markers",type ="character",
        help="the file of marker genes table to be visulized.")
sub_dense$add_argument("-n", "--topn", type="integer", default = 1,
        help = "the number of top markers for each cluster to visualizse.")
sub_dense$add_argument("-c", "--topby", type = "character", default = "avg_Log2FC",
        help="the column used to pick top n marker genes.The option can be one of the column in the input marker genes table.")
sub_dense$add_argument("-x", "--extraGene", type = "character", default = NULL,
        help = "[OPTIONAL]The extra gene list of interest to visualize specified by the user.")
sub_dense$add_argument("-g", "--groupby", type = "character", default = "clusters",
        help = "[OPTIONAL]The grouppinig variable in the metadata for separate the cells to visulize marker genes.")
sub_dense$add_argument("-y", "--splitby", type = "character", default = NULL,
        help = "[OPTIONAL]the variable in the metadata used to split the graph by the variable levels to comparing the gene expression difference in different levels.")
sub_dense$add_argument("--color-scale", type = "character", default = "grey,red",
        help = "[OPTIONAL]the color scale used to map the gene expression value,format:low_color,high_color.[default: %(default)s]")
sub_dense$add_argument("-s", "--pointsize", type = "double", default = 1,
        help = "[OPTIONAL]the point size in the plot.")
# sub_dense$add_argument("--reduct", type = "character", default = "umap",
#         help = "[OPTIONAL]the previous calculated reduction result used in the featureplot,.")
sub_dense$add_argument("--combine", type = "character", default = "FALSE",
        help = "whether to plot the gene pairs in the same plot.[default: %(default)s]")
args <- commandArgs(TRUE)
if ( "density" %in% args ){
  opt<-intial_setting()
  # TO DO
  print("NOT IMPLIMENTED YET!")
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}