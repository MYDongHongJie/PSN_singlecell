sub_pyscvelo <- subparsers$add_parser("pyscvelo", help = "do scvelo by python.")
sub_pyscvelo$add_argument("--loom_dir",type="character",default=NULL,
                         help = "the address of loom")
sub_pyscvelo$add_argument("--groupby",type="character",default = "clusters",
                         help = "color for plot and split by ','.[default: %(default)s]")
sub_pyscvelo$add_argument("--reduction",type="character",default = "umap",
                         help = "umap or tsne .[default: %(default)s]")
sub_pyscvelo$add_argument("--onlygeneplot",type="character",default = 'FALSE',
                         help = "only plot gene plot.TRUE or PLOT .[default: %(default)s]")
sub_pyscvelo$add_argument("--genelist",type="character",default = NULL,
                         help = "top 10 file for gene plot .[default: %(default)s]")
# sub_pyscvelo$add_argument("--metadata",type="character",
#                          help = "the metadata address")
args <- commandArgs(TRUE)
if ( "pyscvelo"  %in% args ){
  opt<-intial_setting()
  if (opt$sub_name == "pyscvelo") {
      ##===================================================================================================
      if (opt$onlygeneplot=="TRUE"){
          pyscvelo <-"./donscvelo.py"
          h5adfile <- opt$input
          cmd <- glue::glue('python {pyscvelo} --input={h5adfile} ',
                            '--output={output_dir} --onlygeneplot={opt$onlygeneplot} --genelist={opt$genelist}')
          futile.logger::flog.info(glue::glue("running:{cmd}"))
          system(cmd)
      }else{
          futile.logger::flog.info("step1 import seurat object(rds )  and convert to  h5seurat")
          library(SeuratDisk)
          object <- readRDS(opt$input)
          Seurat::DefaultAssay(object) <-  opt$assay
          SaveH5Seurat(object, filename = file.path(output_dir,'object.h5Seurat'))
          #h5ad
          Convert(file.path(output_dir,'object.h5Seurat'), dest = "h5ad")
          h5adfile <- file.path(output_dir,'object.h5ad')
          #metadata
          write.table(data.frame(object@meta.data),file=file.path(output_dir,'metadata.txt'),sep="\t",quote = F)
          groupby=unlist( strsplit( opt$groupby, ",", perl =T))
          #获取level
          level_data <- data.frame(clusters= paste(unique(object@meta.data$clusters),collapse=','))  
          for (i in groupby){
              if( ! is.null(levels(object@meta.data[,i])) ) {
                  if( length(levels(object@meta.data[,i])) == length(unique(object@meta.data[,i])) ) {  
                      level_data[,i]=paste(levels(object@meta.data[,i]),collapse=',')
                  }else{
                      level_data[,i]=paste(intersect(levels(object@meta.data[,i]),unique(object@meta.data[,i])),collapse=',')
                  }
              }else{
                  level_data[,i]=paste(sort(unique(object@meta.data[,i])),collapse=',')
              }
          }
          write.table(level_data,file=file.path(output_dir,'group_level.txt'),sep='\t',quote=F,row.names=F) 
          ##===================================================================================================
          pyscvelo <-"./donscvelo.py"
          metadata_dir=file.path(output_dir,'metadata.txt')
          group_level_dir=file.path(output_dir,'group_level.txt')
            if(!is.null(opt$genelist)){
                genelist=opt$genelist
            }else{
                genelist='None'
            }
          print(genelist)
          cmd <- glue::glue('python {pyscvelo} --input={h5adfile} --loom_dir={opt$loom_dir} --metadata={metadata_dir} ',
                            '--output={output_dir} --groupby={opt$groupby} --base={opt$reduction} --genelist={genelist} --order={group_level_dir}')
          futile.logger::flog.info(glue::glue("running:{cmd}"))
          system(cmd)
      }
      ##==================================================================================================
      ## save session informations
if(!file.exists(file.path(output_dir, "scVelo分析说明.docx"))){
    file.copy("/public/scRNA_works/Documents/scVelo分析说明.docx",
    file.path(output_dir, "scVelo分析说明.docx"))
    }    
      quit()
}
}


#  python donscvelo.py
#  --input=/public/scRNA/works/yfang/HT2020-14595_Report_2021-02-25/seuerat_withnewobj.h5ad
#  --loom_dir=/public/scRNA/works/yfang/HT2020-14595_Report_2021-02-25/
#  --metadata=/public/scRNA/works/yfang/HT2020-14595_Report_2021-02-25/seuerat_withnewobj_metadata.txt
#  --output=./scvelotest22
#  --groupby=clusters,area
