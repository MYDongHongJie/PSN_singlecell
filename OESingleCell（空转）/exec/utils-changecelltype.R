# changecelltype
sub_changecelltype = subparsers$add_parser("changecelltype", help = "add annotation for genelist.")
sub_changecelltype$add_argument("-c","--celltype", type = "character", default = NULL,
       help = "The celltype file ")
sub_changecelltype$add_argument("--palette", type = "character", default = "customecol2",
       help = "The color scheme")
sub_changecelltype$add_argument("-r","--reduct", type = "character", default = "umap",
       help = "The reduction to visualize")
args <- commandArgs(TRUE)
if ( "changecelltype" %in% args){
    opt <- intial_setting()
    if(opt$sub_name == "changecelltype" ){
      data_ob = OESingleCell::ReadX(input = opt$input,
                                    informat = opt$informat,
                                    assays = assays,
                                    data.use = dataslots)
      tsv=read.table(opt$celltype,header=T,sep="\t",stringsAsFactors=F,colClasses="character")
      id=names(tsv)[1]
      tsv=tibble::column_to_rownames(tsv,id)
      a=list()
      for (i in rownames(tsv) ){
          a[[i]]=unlist(strsplit(tsv[i,],",",perl = T))
      }
      b=c();for (name in names(a)){ for (cluster in a[[name]]) b=append(b,cluster) }
      if (! identical(sort(as.integer(b)),sort(as.integer(levels(data_ob$clusters))))){
          stop("levels not equal")
      } else {
          to=c()
          for (name in names(a)){ for (cluster in a[[name]]) to[cluster]=name }
          data_ob@meta.data[,id] = plyr::mapvalues(x=data_ob@meta.data[,"clusters"],
              from=levels(data_ob@meta.data[,"clusters"]),
              to=to[levels(data_ob@meta.data[,"clusters"])])
          data_ob@meta.data[,id] = factor(data_ob@meta.data[,id] ,levels=names(a))
          #seurat_ob = SetIdent( seurat_ob, value = "clusters")
      }
      nlevel = length(unique(data_ob@meta.data[,id]))
      pointsize=0.5
      colors2use = OESingleCell::SelectColors(1:length(unique(OESingleCell::colData(data_ob)[[id]])),palette=opt$palette)
      gg = Seurat::DimPlot(object = data_ob,
                           reduction = opt$reduct,
                           pt.size = pointsize,
                           group.by=id)+
           ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
           ggplot2::scale_colour_manual( values = unname(colors2use))
      OESingleCell::save_ggplots(file.path(output_dir,paste0(id,".pdf")),
                                 gg,
                                 width = max(nchar(rownames(tsv)))/15+7,
                                 dpi = 1000,
                                 limitsize = F)

      simplified_meta = data_ob@meta.data %>%
                        dplyr::rename( "Barcode" = "orig.ident") %>%
                        dplyr::select( Barcode, sampleid, clusters,group,!!id)
      write.table(simplified_meta,
                  quote = F,
                  sep =",",
                  row.names = F,
                  file.path(output_dir,paste0(id,".metadata.csv",collapse = "")))
    if ( as.logical(opt$update) ){
      SeuratDisk::UpdateH5Seurat(file = opt$input, object = data_ob, verbose = FALSE )
    }else{
      OESingleCell::SaveX(data_ob, output = opt$output,update = FALSE,
                          outformat = opt$outformat, prefix = opt$prefix)
    }
    write_session_info(log_dir,sub_name = parser$parse_args()$sub_name )
    quit()
    }
}