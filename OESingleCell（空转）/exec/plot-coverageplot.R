sub_coverageplot <- subparsers$add_parser("coverageplot",
                                          help = "Visualize links along with DNA accessibility information by running CoveragePlot.")
sub_coverageplot$add_argument("--peaks", help = "the list of peak regions of interesting to be visulized in a file.")
args <- commandArgs(TRUE)
if ( "coverageplot" %in% args ){
  opt<-intial_setting()
  # CoveragePlot(bone, region = "chr1-40189344-40252549")
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}