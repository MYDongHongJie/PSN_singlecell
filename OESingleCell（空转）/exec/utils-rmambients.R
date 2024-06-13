# create the parsers for subcommand removeEmptyDrops
sub_removeEmpty = subparsers$add_parser("rmambients", help = "remove the empty droplet using DropletUitls")
sub_removeEmpty$add_argument("--fdr", type = "double", default = 0.25,
             help = "the significance level for a cell to be identified as an empty droplet.[default: %(default)s] ")
sub_removeEmpty$add_argument("--iters", type = "integer", default = 1000,
             help = "the iteration to be used in calculation.[default: %(default)s] ")

args <- commandArgs(TRUE)
if ( "rmambients"  %in% args ){
  opt<-intial_setting()
    data_ob = OESingleCell::ReadX(input = opt$input, informat = opt$informat,
                    assays = assays, data.use = dataslots)

    # remove the empty droplets using dropletuitls if necessary
    emptyRemoved = RemoveEmptyDrops(
                  fdrThreshold=as.numeric(opt$fdr),
                  emptyDropNIters= as.numeric(opt$iters),
                  raw_count = GetAssayData(data_ob, slot = "counts") )
    subset_se = data_ob[, colnames(emptyRemoved)]
    # TO DO
    # make a UMI rank plot
    OESingleCell::SaveX(data_ob, output = opt$input, assay = assays[1],
                        outformat = opt$outformat,
                        update = F)
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}