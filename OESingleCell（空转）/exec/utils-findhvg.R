sub_findhvg = subparsers$add_parser("findhvg", help = "features extraction of highly variable features.")
sub_findhvg$add_argument("--select_method", type = "character", default = "vst",
          help = "methods to choose top variable features. Choices:vst, mean.var.plot, dispersion and markvariogram, moransi sepecifi for spatial data.[default: %(default)s]")
sub_findhvg$add_argument("--loess_span", type = "double", default = 0.3,
          help = "(vst only)Loess span parameter used when fitting the variance-mean relationship.[default: %(default)s] ")
sub_findhvg$add_argument("--nfeature", type = "integer", default = 2000,
          help = "Number of features to select as top variable features.[default: %(default)s]")
sub_findhvg$add_argument( "--hvgnumber", type="integer", default = 6,
          help="the number of top spatially variable features which were outputed.[default: %(default)s]")

args<-commandArgs(TRUE)

# =============== Subcmd: findhvg, find the variable features ============
if ( "findhvg" %in% args ){
  opt<-intial_setting()
  data_ob = OESingleCell::ReadX(input = opt$input, informat = opt$informat,
                  assays = assays, data.use = dataslots,  # counts is enough
                  verbose = F)

  if ( opt$select_method %in% c("vst", "dispsersion", "mean.var.plot") ){ ##单细胞
    data_ob = OESingleCell::FindVariableFeatures( data_ob,
                                                  assay = assays[1],
                                                  selection.method = opt$select_method,
                                                  loess.span = opt$loess_span,
                                                  nfeatures = opt$nfeature,
                                                  verbose = F)
  }else if ( opt$select_method %in% c("markvariogram", "moransi") ) { ##空间
      data_list = OESingleCell::SplitObject(data_ob,split.by="sampleid")
      seu_list_copy = data_list
      sct_var <- function(x){
		x = Seurat::SCTransform(x, assay = "Spatial", verbose = FALSE)
		x = Seurat::VariableFeatures(x)
		}
      var_gene = lapply(seu_list_copy,sct_var)
      rm(seu_list_copy)
      spatial_DEGs_list = future.apply::future_lapply(1:length(data_list), function(sidx) {
        datax = Seurat::FindSpatiallyVariableFeatures(data_list[[sidx]] ,
                                            assay = assays[1],###空间使用SCT
                                            slot = "scale.data",
                                            image = Seurat::Images(data_list[[sidx]])[sidx],
                                            selection.method = opt$select_method,
                                            r.metric = 5,
                                            nfeatures = opt$nfeature,
                                            verbose = TRUE,
                                            features = var_gene[[sidx]])
        spatial_DEGs = Seurat::SpatiallyVariableFeatures(datax, selection.method = tolower(opt$select_method) )
      }, future.seed = 2020)
      names(spatial_DEGs_list) = names(data_list)
      global_DEGs = stack(spatial_DEGs_list) %>% dplyr::rename( gene = values, slice = ind)
      write.table(global_DEGs,
                  file = file.path(output_dir,paste0("top", opt$nfeature, "_for_all_spots.xls", collapse = "")),
                  col.names =T,row.names = F,sep = "\t",quote = FALSE)
      tophvg <- function(geneset=NULL,genenumber=NULL){
        geneset$ord <- as.numeric(rownames(geneset))
		newgene <- geneset %>% dplyr::group_by(slice) %>% dplyr::top_n(n = genenumber, wt =dplyr::desc(ord))
		newgene <- data.frame(newgene[,1:2])
	  }
      if ( is.null(opt$hvgnumber)){
	    hvgnumber = 6
      }else{
	    hvgnumber = opt$hvgnumber
      }
	  hvg_topgene <- tophvg(geneset = global_DEGs,genenumber = hvgnumber)
	  write.table(hvg_topgene, file = file.path(output_dir,paste0("top_",hvgnumber,"_hvg", "_for_all_spots.xls",
                                                                  collapse = "")), col.names =T,row.names = F,
                                                                  sep = "\t",quote = FALSE)
  }else if ( opt$select_method %in% c("spark","sparkx") ) {  ##空间使用sparkx
      #suppressPackageStartupMessages( library("SPARK") )
      data_spark <- OESingleCell::RunSpark ( object=data_ob,
                                              assay = "Spatial",
                                              split.by = "sampleid",
                                              method = opt$select_method,
                                              fit.model = "poisson",
                                              only.positive = TRUE,
                                              jobs = future::nbrOfWorkers(),
                                              verbose = FALSE)
      data_spark$gene <- rownames(data_spark)
      write.table(data_spark,
                  file = file.path(output_dir,paste0("sparkx_hvg_for_all_spots.xls", collapse = "")),
                  col.names =T,row.names = F,sep = "\t",quote = FALSE)
      data_spark_filter <- data_spark %>% dplyr::group_by(sample) %>% dplyr::filter(as.numeric(pval)<0.05)
      write.table(data_spark_filter,
                  file = file.path(output_dir,paste0("sparkx_hvg_pval0.05_for_all_spots.xls", collapse = "")),
                  col.names =T,row.names = F,sep = "\t",quote = FALSE)
      if ( is.null(opt$hvgnumber)){
	    hvgnumber = 6
      }else{
	    hvgnumber = opt$hvgnumber
      }
      data_spark_filter_top <- data_spark_filter %>%
                               dplyr::top_n(n = hvgnumber, wt =dplyr::desc(pval))
      write.table(data_spark_filter_top,
                  file = file.path(output_dir,
                                   paste0("sparkx_top_",hvgnumber,"_hvg", "_for_all_spots.xls", collapse = "")),
                  col.names =T,
                  row.names = F,
                  sep = "\t",
                  quote = FALSE)
      ###plan(multicore)
      #data_list = OESingleCell::SplitObject(data_ob,split.by="sampleid")
      ##future.apply::future_lapply(1:length(data_list), function(sidx) {
      # for (sidx in 1:length(data_list)){
      #       OESingleCell::spark_analysis(data_list[[sidx]],
      #                                    sampleid=Seurat::Images(data_list[[sidx]])[sidx],
      #                                    outdir=output_dir,
      #                                    slot = "data",
      #                                    fit.model = "poisson",
      #                                    num_core = 1,
      #                                    method = opt$select_method)

      #}
      #, future.seed = 2020)
  }

  if ( tolower(opt$outformat) == "h5seurat" ){
    hfile =  SeuratDisk::h5Seurat$new(filename = opt$input, mode = 'r+')
    on.exit(expr = hfile$close_all())
    assay_group = hfile[[glue::glue("assays/{assays[1]}")]]
    # Write out variable features
    assay_group$link_delete("variable.features")
    SeuratDisk::WriteH5Group( OESingleCell::VariableFeatures(data_ob),
                              name = 'variable.features',
                              hgroup = assay_group,
                              verbose = F)
    # Write out meta.features
    # attrs = hdf5r::h5attr_names(assay_group[["meta.features"]])
    # if ( length(attrs) ){
    #   sapply(attrs, function(x) assay_group[["meta.features"]]$attr_delete(x) )
    # }
    # assay_group$link_delete("meta.features")
    # assayx = Seurat::GetAssay( data_ob, assay = assays[1] )
    # print("XXXX")
    # currently I can not fix the error yet:
    # Can't create dataset _index - already exists!
    # SeuratDisk::WriteH5Group( assayx[[]], name = 'meta.features',
    #                           hgroup = assay_group,
    #                           verbose = F )
  }else{
    OESingleCell::SaveX(data_ob, output = opt$input, outformat = opt$outformat)
  }
  sessionInfo()
  quit()
}


#sctool -i /public/dev_scRNA/yfang/project/HT2020-16980-human-breastcancer/result/cluster_seurat/singlecell_object.clustering_resolution0.4.rds -f rds -o /public/dev_scRNA/yfang/project/HT2020-16980-human-breastcancer/hvgtest/moransi/ -d seurat findhvg --select-method moransi --nfeature 2000 --hvgnumber 6