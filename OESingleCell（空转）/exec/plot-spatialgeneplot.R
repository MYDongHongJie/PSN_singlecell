sub_spgene_pt = subparsers$add_parser("spgene_pt", help = "plot for cellphonedb result.")
sub_spgene_pt$add_argument("-g", "--genes", type="character",
             help="the gene set will be ploted in slce")
sub_spgene_pt$add_argument("--min", type="character",default="q0",
             help="the min.cutoff for spatialplot.[default: %(default)s]")
sub_spgene_pt$add_argument("--max", type="character",default="q90",
             help="the max.cutoff for spatialplot.[default: %(default)s]")
sub_spgene_pt$add_argument("--method", type="character", default='naivemean',
             help = "the mathod to caculat the co-expression of gene set and can be 'naivemean','cumsum','addmodulescore'.[default: %(default)s]")
sub_spgene_pt$add_argument("--spt",type="double", default= 1.0,
             help = "the spot size.[default: %(default)s]")
sub_spgene_pt$add_argument("--slice",type="character", default=NULL,
             help = "select slices to plot and the default is all.[default: %(default)s]")

args <- commandArgs(TRUE)
if ( "spgene_pt" %in% args ){
  opt<-intial_setting()
  object <- OESingleCell::ReadX(input = opt$input,
                                 informat = opt$informat,
                                 assays = opt$assay,
                                 verbose = F,
                                 images=opt$image)
  Seurat::DefaultAssay(object) <- opt$assay
  savegene_plot <- function(object,input,filename,spotsize,slice,output,min,max){
       # print(slice)
       gsp <- OESingleCell::SpatialPlot(object, assay = opt$assay,
                                        features = input,
                                        min.cutoff = min,
                                        max.cutoff = max,
                                        alpha = c(0.1, 1),
                                        pt.size.factor = spotsize,
                                        pt.alpha = TRUE,
                                        stroke = 0.2, combine = FALSE, images = slice)
      gsb <-  do.call(ggpubr::ggarrange,
                         c(gsp,
                           list(nrow = ceiling(length(gsp)/2),
                                ncol =  2,
                                common.legend = TRUE,
                                legend = "right",
                                align = "hv")))
      width <- 8+ceiling(nchar(input)/10)
      height <- 4*ceiling(length(gsp)/2)
      OESingleCell::save_ggplots(file.path(output,filename),
                                plot=gsb,
                                width=width,
                                height=height,
                                dpi=200,bg='white')
      }
   print('s1')
   print('s1')
  genes_df <- read.table(opt$genes,sep="\t",quote="") %>% split(1:dim(.)[1] )
  genes<- lapply(genes_df,function(x){
      genes <- x[1,] %>% strsplit(split = ",") %>% unlist()

      genes <- intersect(genes,rownames(object))

      return(genes)


  })
   print('s2')
  if(!is.null(opt$slice)){
    slice <- opt$slice %>% strsplit(split = ",") %>% unlist()
      }else{slice <- Seurat::Images(object)}
   print(slice)

  for (i in 1:length(genes)){
      for(gene in genes[[i]]){
    filenames=paste0("geneset_",names(genes)[i],"_",gene,"_plot")
    savegene_plot(object,input=gene,filename=filenames,spotsize=opt$spt,slice=slice,output=output_dir,min=opt$min,max=opt$max)

  }}

  if(!is.null(opt$method)){
      plotgene <- switch(opt$method,
          'naivemean'={
                       for(i in 1:length(genes)){
                       gene_list <- Seurat::GetAssayData(object[genes[[i]],],slot = "counts") %>%
                       data.frame %>% t() %>%data.frame()
                       names(gene_list) <- genes[[i]]
                       gene_list$sum <-  rowSums(gene_list)
                       gene_list$Ratio <-  rowSums(gene_list[,genes[[i]]]> 0)/length(genes[[i]])
                        met_line <- paste(genes[[i]],collapse = "_")%>% gsub('-','.',.)
                       object[[met_line]] <- (gene_list$sum)*(gene_list$Ratio)
                       filenames=paste0('geneset_',names(genes)[i],"_gene_coexp_by_NaiveMean")
                       print('s3')
                       savegene_plot(object,
                                     input=met_line,
                                     filename=filenames,
                                     spotsize=opt$spt,
                                     slice=slice,
                                     output=output_dir,
                                     min=opt$min,
                                     max=opt$max)

                       }


                   },
          'cumsum' ={ for(i in 1:length(genes)){
                     met_line <- paste(genes[[i]],collapse = "_")%>% gsub('-','.',.)
                     gene_list <- Seurat::GetAssayData(object[genes[[i]],],slot = "counts") %>%
                          data.frame %>%
                          t() %>%
                          data.frame
                      names(gene_list) <- genes[[i]]
                      gene_list$sum <-  rowSums(gene_list)
                      object[[met_line]] <- gene_list$sum
                      filenames=paste0('geneset_',names(genes)[i],"_gene_coexp_by_CumSum")
                      print('s4')
                      savegene_plot(object,
                                input=met_line,
                                filename=filenames,
                                spotsize=opt$spt,
                                slice=slice,
                                output=output_dir,
                                min=opt$min,
                                max=opt$max)
          }
          },
          'addmodulescore'={for(i in 1:length(genes)){
                                      filenames=paste0('geneset_',names(genes)[1],"_addmodulescore")
                                      object=Seurat::AddModuleScore(object, features=list(genes[[i]]),name=paste0('geneset_',names(genes)[i]))
                                      savegene_plot(object,
                                          input=paste0('geneset_',names(genes)[i],'1'),
                                          filename=filenames,
                                          spotsize=opt$spt,
                                          slice=slice,
                                          output=output_dir,
                                          min=opt$min,
                                          max=opt$max)   }



          }



      )

  }


}


