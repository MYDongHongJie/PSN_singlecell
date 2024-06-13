sub_removedbls <- subparsers$add_parser("rmdoublets", help = "remove doublets.")
sub_removedbls$add_argument("--dbl_rates", type = "double", default = 0.06,
                            help = "the expected doublets rate in each sample.")
sub_removedbls$add_argument("--method", type = "character", default = "scrublet",
                            help = "the method used to detect doublets, choices: scrublet, doubletfinder.[default: %(default)s]")
sub_removedbls$add_argument("--subset", type = "character", default = "TRUE",
                            help = "logical indication of wether to subset the data object for only singlets.[default: TRUE]")
sub_removedbls$add_argument("--ident", type = "character", default = "clusters",
                            help = "[OPTIONAL]ONLY AVAILABLE for doubletfinder,the name of cell group annotation.[default: %(default)s]")

# =============== Subcmd: rmdoublets, remove the doublets from the data matrix =========
args <- commandArgs(TRUE)
if ( "rmdoublets"  %in% args ){
  opt<-intial_setting()
  ####结束
  bimodality_coefficient <- function(x) {
    n <- length(x)
    S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
    G <- S*(sqrt(n*(n-1)))/(n-2)
    K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
    K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
    B <- ((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))
    return(B)
  }
  # read the specified assay and data slot in data object into memory
  data_ob <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 reductions = "pca",
                                 assays = assays, # only RNA assay is valid
                                data.use = "counts", # counts slot is enough
                                verbose = F)
  data_list <- OESingleCell::SplitObject(data_ob, split.by = "sampleid")

  is_doublets <- switch(opt$method,
    "scrublet" = {
        dbl_all <- future.apply::future_lapply(data_list, function(x){
          x <- OESingleCell::RunScrublet(x, doublet.rate = opt$dbl_rates)
          dbl_res <- Seurat::FetchData(x, var = "predicted_doublets")
        }, future.seed = 123)
      },
    "doubletfinder" = {
      cellmeta <- OESingleCell::colData(data_ob)
      if ( ! opt$ident %in% colnames(cellmeta) ) stop("NO specfied cell identity column FOUND in cell annotation!")
      dbl_all <- list()
      for ( x in 1:length(data_list) ){
        data_list[[x]] <- RunDoubletFinder(data_list[[x]], doublet.rate = opt$dbl_rates, identity = opt$ident)
        dbl_all[[x]] <- Seurat::FetchData(data_list[[x]], var = "predicted_doublets")
      }
      dbl_all
        # dbl_all <- future.apply::future_lapply(data_list, function(x){
      #   x <- RunDoubletFinder(x, doublet.rate= opt$dbl_rates, identity = opt$ident )
      #   dbl_res <- Seurat::FetchData(x, var = "predicted_doublets")
      # }, future.seed = 123)
    }
  )
  if ( as.logical(opt$subset) ){ # for only singlets
    data_ob <- OESingleCell::AddMetaData(data_ob, metadata = do.call(rbind, is_doublets)[[1]], col.name = "is_doublets")
    # append all other data not loaded before subsetting for h5seurat formated file.
    if ( opt$informat == "h5seurat" ) data_ob <- SeuratDisk::AppendData(opt$input, data_ob)
    data_ob <- subset(data_ob, subset = is_doublets == FALSE)
    OESingleCell::SaveX(data_ob,
                        output = opt$output,
                        update = FALSE,
                        outformat = opt$outformat,
                        prefix = opt$prefix)
  }else{
    if ( tolower(opt$outformat) == "h5seurat" & as.logical(opt$update) ){
      SeuratDisk::UpdateH5Seurat(file = opt$input,
                                 object = data_ob,
                                 reduction = unique(c(opt$reduct1, opt$reduct2)))
    }else{
      OESingleCell::SaveX(data_ob,
                          output = opt$output,
                          update = FALSE,
                          outformat = opt$outformat,
                          refix = opt$prefix)
    }
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}