sub_geneacitivity <- subparsers$add_parser("geneacitivity", help = "predicting gene activity from atac peak assay,CCANS and Links.")
sub_geneacitivity$add_argument("--window","-w", type = "integer", default = 500000,
                            help = "the genome sliding window to determine the distance threshold for calculating coaccess peaks.[default: %(default)s] .")
sub_geneacitivity$add_argument("--outassay","-s",type = "character", default = "FALSE",
                            help = " the assays name to store the geneAcitivities data .[default: %(default)s] .")
sub_geneacitivity$add_argument("--method","-m",type = "character", default = "signac",
                            help = "method to calculate the gene activity .[default: %(default)s] .")
sub_geneacitivity$add_argument("--run-ccan", type = "character", default = "FALSE",
                           help = "[OPTIONAL]wether to find cis-Co-accessibility Networks (CCANS).[default: %(default)s] ")
sub_geneacitivity$add_argument("--do-norm", type = "character", default = "TRUE",
                           help = "[OPTIONAL]wether to normalize the count matrix from cicero.[default: %(default)s] ")

# ================ Subcmd: runcicero, predicting gene activity. =======================================================
if ( "geneacitivity"  %in% args ){
  opt<-intial_setting()
  if(opt$sub_name =="geneactivity"){
    message("step1:read the specified assay and data slot in data object into memory")#===============================================
    data_ob <- OESingleCell::ReadX(input = opt$input,
                                   informat = opt$informat,
                                   assays = assays,
                                   data.use = dataslots, # data slot is enough
                                   verbose = F)

    # print("create geneActivities assay with counts")
    # data_ob[["geneActivities"]] <- Seurat::CreateAssayObject(counts = cicero_gene_activities)
    message(glue::glue("step2:create geneActivities assay with counts by using {opt$method}"))
    if(opt$method=="cicero"){
      data_ob <- OESingleCell::runCicero(object = data_ob,
                                         assay = assays[1],
                                         outassay=outassay,
                                         window = opt$window,
                                         out.dir = output_dir,
                                         run.ccan = as.logical(opt$run_ccan))
    }else if(opt$method=="signac"){
      #To create a gene activity matrix, we extract gene coordinates and extend them to include the 2 kb upstream region
      # (as promoter accessibility is often correlated with gene expression).We then count the number of fragments for each
      # cell that map to each of these regions, using the using the FeatureMatrix() function.
      gene.activities <- Signac::GeneActivity(data_ob)
      # add the gene activity matrix to the Seurat object as a new assay and normalize it
      data_ob[[outassay]] <- Seurat::CreateAssayObject(counts = gene.activities)
      data_ob <- Seurat::NormalizeData(object = data_ob,
                               assay = outassay,
                               normalization.method = 'LogNormalize',
                               scale.factor = median(data_ob$nCount_RNA)
      )
    }
    # if ( as.logical(opt$do_norm) ){
    #   print("add cicero_gene_activities(data) to geneActivities assay")
    #   data_ob <- Seurat::SetAssayData(object = data_ob,
    #                                   assay = opt$outassay,
    #                                   new.data = as(log2(Seurat::GetAssayData(data_ob, assay = "RNA", slot = "counts")*10^6 + 1 ) ,"dgCMatrix"),
    #                                   slot = 'data')
    #   DefaultAssay(data_ob) <- opt$outassay
    #   print("add scale data")
    #   data_ob <- Seurat::ScaleData(data_ob, assay =  opt$outassay)
    # }
    SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob,assay = opt$outassay )
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
  }
}