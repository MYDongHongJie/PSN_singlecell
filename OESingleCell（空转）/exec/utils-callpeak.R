sub_callpeak = subparsers$add_parser("callpeak", help = "call peak for each group of cell for scATAC-seq")
sub_callpeak$add_argument("--groupby", type = "character", default = "clusters",
                          help = "Grouping variable to use. If set, peaks will be called independently on each group of cells and then combined. Note that to call peaks using subsets of cells we first split the fragment file/s used, so using a grouping variable will require extra time to split the files and perform multiple MACS peak calls, and will store additional files on-disk that may be large. Note that we store split fragment files in the temp directory, and if the program is interrupted before completing these temporary files will not be removed. If NULL, peaks are called using all cells together (pseudobulk).[default: %(default)s]")
sub_callpeak$add_argument("--mac2", type = "character", default = NULL,
                          help = "Path to MACS program. If NULL, try to find MACS automatically.[default: %(default)s]")
sub_callpeak$add_argument("--broad", type = "character", default = "FALSE",
                          help = "whether to Call broad peaks (--broad parameter for MACS).[default: %(default)s]")
sub_callpeak$add_argument("--ident", type = "character", default = NULL,
                          help = "List of identities to include if grouping cells (only valid if also setting the group.by parameter). If NULL, peaks will be called for all cell identities.[default: %(default)s]")
sub_callpeak$add_argument("--combine", type = "character", default = NULL,
                          help = "Controls whether peak calls from different groups of cells are combined using GenomicRanges::reduce when calling peaks for different groups of cells (group.by parameter). If FALSE, a list of GRanges object will be returned. Note that metadata fields such as the p-value, q-value, and fold-change information for each peak will be lost if combining peaks.[default: %(default)s]")
sub_callpeak$add_argument("--chunk", type = "integer", default = 2000,
                          help = "Number of regions to load into memory at a time, per thread. Processing more regions at once can be faster but uses more memory.[default: %(default)s]")
# ================ Subcmd: callpeak, peak calling for the single cell Chromatin data object. ========
args <- commandArgs(TRUE)
if ( "callpeak"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name == "callpeak" ){
    idents = unlist( strsplit(opt$idents, ",") )
    # read the specified assay and data slot in data object into memory
    data_ob = OESingleCell::ReadX(input = opt$input,
                                  informat = opt$informat,
                                  assays = assays,
                                  outdir = tempdir(),
                                  data.use = dataslots, # data slot is enough
                                  verbose = F)

   # To use the peak calling functionality in Signac you will first need to install MACS2.
    peaks = Signac::CallPeaks( data_ob,
                               group.by = opt$groupby,
                               macs2.path = opt$mac2,
                               idents = idents,
                               broad = as.logical(opt$broad),
                               outdir = tempdir(),
                               combine.peaks = as.logical(opt$combine))

    # Construct a feature x cell matrix from a genomic fragments file
    peak_counts <- Signac::FeatureMatrix(fragments,
                                         features = peaks,
                                         process_n = opt$chunk,
                                         sep = c("-", "-"),
                                         verbose = FALSE )

    data_ob[["peaks"]] <- Seurat::CreateAssayObject(counts = peak_counts)

    OESingleCell::SaveX(data_ob,
                        output = output_dir,
                        outformat = opt$outformat,
                        prefix = opt$prefix,
                        update = FALSE )

    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}