sub_convert = subparsers$add_parser("convertX", help = "convert h5ad and rds")
sub_convert$add_argument("-r", "--reduction", type = "character", default = 'tsne',
          help = "the reduction of data.Can be 'umap','tsne','umap,tsne'.[default: %(default)s]")
sub_convert$add_argument("-s", "--sample_col", type = "character", default = 'sampleid',
          help = "the sample matadata col.[default: %(default)s]")
##h5ad转rds
args <- commandArgs(TRUE)
if ( "convertX"  %in% args ){
  opt<-intial_setting()
  scanpy <- reticulate::import("scanpy")
  pandas <- reticulate::import("pandas")
  adata  <-  scanpy$read(opt$input)
  assays=opt$assay
  reduction <- opt$reduction

  if (assays=="RNA"){
      seurat_obj <- OESingleCell::h5ad_rds(adata=adata,reduction=reduction,assays=assays)
      saveRDS(seurat_obj,file=paste0(output_dir,'seurat_obj.rds'))
      #saveRDS()
      #picture

  }else if(assays=="Spatial"){
      sample_col <- opt$sample_col
      seurat_obj <- OESingleCell::h5ad_rds(adata=adata,reduction=reduction,assays=assays)
      sample_bar <- Seurat::FetchData(seurat_obj,vars=sample_col) %>% tibble::rownames_to_column()
      sample_bar_list <- split(sample_bar$rowname,sample_bar[,sample_col])
      for(i in 1:length(sample_bar_list)){

          image_obj <- OESingleCell::add_images(adata=adata,rowname=sample_bar_list[[i]],sample_col=sample_col)
          slice <- names(sample_bar_list[i])
          DefaultAssay(object = image_obj) <- "Spatial"
          seurat_obj[[slice]] <- image_obj
      }
      saveRDS(seurat_obj,file=paste0(output_dir,'seurat_obj.rds'))

  }else{print("can not")}




}
