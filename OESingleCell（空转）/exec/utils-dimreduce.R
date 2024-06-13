sub_dimreduce = subparsers$add_parser("dimreduce", help = "dimension reduction.")
sub_dimreduce$add_argument("-t", "--components", type = "integer", default = NULL,
                           help = "the appropriate number of dimensions for primary reduction. Recommendations:pca:30,mnn:10[default: %(default)s]")
sub_dimreduce$add_argument("--perplexity", "-p", type = "integer", default = 30,
                           help = "The value of the perplexity used for tSNE.[default: %(default)s] ")
sub_dimreduce$add_argument("--rerun", type = "character", default = "FALSE",
                           help = "[OPTIONAL]wether to rerun the reduction to overwrite the previous results.[default: %(default)s] ")
sub_dimreduce$add_argument("--reduct1", type = "character", default = NULL,
                           help = "[OPTIONAL]the primary reduction methods whose results can be used as input for secondary reduction. The supported methods can be one of the followings: ica, cca, pca, rpca, lsi, svd, mnn, harmony, tsne, flt-tsne, umap, diffusion.[default: %(default)s] ")
sub_dimreduce$add_argument("--reduct2", type = "character", default = "tsne",
                           help = "one or more of the scondary dimension reduction methods used for this run. the options can be tsne, flt-sne, umap, duffusion, swne. That is to say not all the primary redcution method can used as scondary one.[default: %(default)s] ")
sub_dimreduce$add_argument("--batchid", type = "character", default = NULL,
                           help = "[OPTIONAL]the batch information column name of sample in the metadata.[default: %(default)s] ")


# =============== Subcmd: dimreduce, run dimension reduction =========
args <- commandArgs(TRUE)
if ( "dimreduce"  %in% args ){
  opt<-intial_setting()
  if ( is.null(opt$reduct1 )){
    print("NO first level dimension reduction methods specified,the default PCA will be used!")
    reduct1 = "pca"
  }else{
    reduct1 = tolower(opt$reduct1)
  }

  if ( !is.null(opt$reduct2) ){
    reduct2 = tolower(unlist(strsplit(opt$reduct2,",",  perl = T)))
  }else{
    # if no secondary reduction method specified, use primary reduction instead
    reduct2 = reduct1
  }

  if ( !is.null(opt$batchid) ){
    batchid = opt$batchid
  }

  data_ob = OESingleCell::ReadX(input = opt$input,
                                informat = opt$informat,
                                assays = assays,
                                reduction = opt$reduct,
                                data.use = dataslots,  # data,scale.data
                                verbose = F)

  # get the subset of cells used for visualization if necessay
  # subset cells on the expression matrix using logical exprssion on features
  if ( !is.null(opt$predicate) ){
      df = OESingleCell::colData(data_ob)
      desired_cells= subset(df, eval( parse(text=opt$predicate)))
      data_ob = data_ob[, rownames(desired_cells)]
  }

  # check out wether the specified primary reduction has beed calculated
  is.rerun_reduct1 = FALSE; is.rerun_reduct2 = FALSE
  switch(as.character(opt$rerun),
         "reduct1" = {
           is.rerun_reduct1 = TRUE
         },
         "reduct2" = {
           is.rerun_reduct2 = TRUE
         },
         "all" = {
           is.rerun_reduct1 = TRUE
           is.rerun_reduct2 = TRUE
         } )

  if ( !is.null(reduct1 )){ # primary reduction is specified
    # if no previous reduction found, catch the error and rerun the
    # reduction without stop with error
    if ( !reduct1 %in% OESingleCell::Reductions(data_ob) | is.rerun_reduct1 ){
      print( "NO specified primary reduction found in the object! Rerun begins!")
      message(paste0("Beginning ", reduct1, " Dimension Reduction"))
      dim_outdir = file.path(output_dir,paste0(reduct1, "_Dimension_Reduction"))
      if ( !dir.exists(dim_outdir) ){
        dir.create(dim_outdir)
      }

      data_ob = OESingleCell::RunDimReduc(data_ob, reduct1 = reduct1,
                                          reduct2 = reduct1,
                                          feature.use = OESingleCell::VariableFeatures(data_ob),
                                          perplexity = perplexity,
                                          assay.use = OESingleCell::DefaultAssay(data_ob),
                                          batch.use = batchid, npcs.use = opt$components)
      reduct1_coord = Seurat::FetchData(data_ob,
                                        var = c("rawbc", paste0(Seurat::Key(data_ob)[reduct1], 1:2))) %>%
        dplyr::rename("Barcode" = "rawbc")
      write.table( reduct1_coord, file.path(dim_outdir, paste0(reduct1, "_Dimension_Reduction_coordination.csv")),
                   sep = ",", col.names = T, row.names = F, quote = F)
    }else{
      print("The specified reduction results have beed calculated!Skip it!")
    }
  }

  if ( opt$reduct1 %in% c("pca","cca","harmony","ica", "lsi") ){ #if the prevous reduction is primary reduction
      #check the components
      # find the optimal components for secondary reduction
      optimal_pc = tryCatch(
                expr = Seurat::Command(data_ob, command = glue::glue("FindNeighbors.{opt$assay}.{opt$reduct1}"), value = "dims"),
                error = function(...){ return(NULL) }
              )
      if ( is.null(optimal_pc) ){ #if previous optimal components not available
          print( "NO previous optimal components is AVAILABLE, the optimal components number be detected automatically.")
          elb = Seurat::ElbowPlot(data_ob, reduction = opt$reduct1)
          optimal_pc = OESingleCell::FindElbow(elb$data)
          # suppressWarnings({ Misc(data_ob, "optimal_pc") = optimal_pc})
      }
  }else{ # reduct1 is not primary reduction and previouly run without primary reduction
          # the exception is mnn, it's better to be specified from command line with relative low value
          optimal_pc = min(opt$components, ncol(OESingleCell::Embeddings(data_ob,reduction = opt$reduct1)))
  }

  # different reduction method from primary reduction specified or forced to rerun
  if ( !opt$reduct2 %in% OESingleCell::Reductions(data_ob) | (opt$reduct2 != reduct1 & is.rerun_reduct2) ){
      message(paste0("Beginning ", opt$reduct2, " Dimension Reduction"))
      output_dir = file.path(output_dir,paste0(opt$reduct2, "_Dimension_Reduction"))
      if ( !dir.exists(output_dir) ){
        dir.create(output_dir)
      }

    data_ob = OESingleCell::RunDimReduc(data_ob, reduct1 = reduct1, reduct2 = opt$reduct2,
                                        feature.use = VariableFeatures(data_ob),
                                        assay.use = DefaultAssay(data_ob),
                                        batch.use = batchid, npcs.use = optimal_pc)
    reduct1_coord = Seurat::FetchData(data_ob,
                                      var = c("rawbc", paste0(Key(data_ob)[opt$reduct2], 1:2))) %>%
      dplyr::rename("Barcode" = "rawbc")
      write.table( reduct1_coord, file.path(dim_outdir, paste0(opt$reduct2, "_Dimension_Reduction_coordination.csv")),
                   sep = ",", col.names = T, row.names = F, quote = F)
  }

  if ( as.logical(opt$update) ){
    SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob,
                                 reduction = unique(c(opt$reduct1, opt$reduct2)))
  }else{
    OESingleCell::SaveX(data_ob, output = opt$output,update = FALSE,
                        outformat = opt$outformat, prefix = opt$prefix)
  }
  write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
  quit()
}