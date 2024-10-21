
library(optparse)
option_list <- list(
  make_option(c('-r','--rds'),help="rds path"),
  make_option(c('-g','--genefile'),help="file of marker gene,two colnumn",default = NULL),
  make_option(c('-c','--column'),help="colname of cells label in meta.data",default = "celltype"),
  make_option(c('-o','--outdir'),help="outdir",default = "./")
)
opt <- parse_args( OptionParser(option_list=option_list))
suppressMessages({
library(grid)
library(Seurat)
library(ggplot2)
library(paletteer)
library(scales)
library(patchwork)
library(RColorBrewer)
})
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
cols = colorls$"NPG"
#cols=alpha(paletteer::paletteer_d('ggsci::category20c_d3'),0.9)
#cols <- c(cols,"#8FBCBBFF","#88C0D0FF","#81A1C1FF","#5E81ACFF","#2A363BFF","#019875FF","#99B898FF","#FECEA8FF","#FF847CFF","#E84A5FFF","#C0392BFF","#96281BFF","#96281BBF")

rds <- opt$rds
genefile <- opt$genefile
outdir = opt$outdir
if(!file.exists(outdir)){
       dir.create(outdir)
}
PRO <- readRDS(rds)
if(!(opt$column %in% colnames(PRO@meta.data))){
    print("Please check labels colnames,dont allow repeat cellname!")
    print(colnames(PRO@meta.data))
    quit()
}
DefaultAssay(PRO) <- "RNA"
markerdf <- read.delim(genefile,sep="\t")
markerdf= markerdf[,c("ABV","marker")]
#Idents(PRO) <- factor(as.character(PRO@meta.data[,opt$column]),levles <- markerdf[,"V1"])
print(unique(Idents(PRO)))

PtSize <- function(cellnum){
    if(cellnum > 10000){
        pt.size <- 0.1
    }else if(cellnum > 5000){
        pt.size <- 0.2
    }else if(cellnum > 2000){
        pt.size <- 0.3
    }else{
        pt.size <- 0.5
    }
   return(pt.size)
}

PlotSize <- function(leng){
     #width <- if(leng>1) 10 else 6
     width=10
     height <- ceiling(leng/2)*5
     if(leng==1){
          width <- 5
          height <- 5
      }
      if(leng==2){
          width <- 10
          height <- 5
      }
     plotsize <- list(width=width, height=height)
     return(plotsize)
}

GeneUnique <- function(alist){
    tmp = vector(mode = "list", length = length(alist))
    tmp[[1]] <- alist[[1]]
    #names(tmp) <- names(alist)[1]
    for (i in 1:length(alist)-1){
        other <- alist[0:i]
        other <- unique(unlist(other, recursive = FALSE))
        speci <- unique(alist[[i+1]])
        specileave <- speci[!speci %in% other]
        tmp[[i+1]] <-  specileave
    }
    names(tmp) <- names(alist)
    return(tmp)
}
# if(length(markerdf[,1]) != length(unique(markerdf[,2]))){
#     print("Please check celltype,dont allow repeat cellname!")
#     quit()
# }

pt.size <- PtSize(nrow(PRO@meta.data))

#list_genes=split(unlist(strsplit(markerdf$V2,",")), markerdf$V1)
list_genes <- split(markerdf[,2], markerdf[,1])
list_genes <- lapply(list_genes,function(x){CaseMatch(unlist(strsplit(x,",")),rownames(PRO))})
#list_genes <- lapply(list_genes,function(x){unlist(strsplit(x,","))[unlist(strsplit(x,",")) %in% rownames(PRO@assays$RNA["counts"])]})
unique_list_genes <- GeneUnique(list_genes)

order_unique_list_genes = unique_list_genes[names(table(PRO@meta.data[,opt$column]))]

p1 <- DotPlot(PRO,
           features=order_unique_list_genes,
           cols = c("grey", "red"),
           cluster.idents = F)+
  RotatedAxis()+
  theme(
    panel.border = element_rect(color="black"), #面板边框
    panel.spacing = unit(1, "mm"), #面板间距

    # 分面标题
    #strip.background = element_rect(color="red"),
    strip.text = element_text(angle = 45,vjust=0.6,size=8,margin=margin(b=2, unit="mm")),
    strip.placement = 'outlet', #
    # 坐标轴线
    axis.line = element_blank(),
  )+labs(x="", y="")
  q <- ggplotGrob(p1)
  lg <- linesGrob(x=unit(c(0,1),"npc"), y=unit(c(0,0)+0.2,"npc"),
                gp=gpar(col="black", lwd=4))

  grid.newpage(); #grid.draw(lg)
  for (k in grep("strip-t",q$layout$name)) {
    q$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
  }
    png(paste0(outdir,"/","dotplot.png"),width = 100*(8+ceiling(sum(lengths(unique_list_genes))/10)) )
    grid.draw(q)
    dev.off()
    pdf(paste0(outdir,"/","dotplot.pdf"),width=8+ceiling(sum(lengths(unique_list_genes))/10))
    grid.draw(q)
    dev.off() 

for(ctype in names(list_genes)){
    print(ctype)
    print(list_genes[[ctype]])
    plotsize <- PlotSize(length(list_genes[[ctype]]))
    print(plotsize)
    genes <- list_genes[[ctype]][list_genes[[ctype]] %in% rownames(PRO@assays$RNA["counts"])]
    if(length(genes)>0){
    if(length(genes)==1){
        p1 <- FeaturePlot(PRO, features = genes,cols = c("lightgrey", "red"),order=TRUE,ncol=1)
    }else{
        p1 <- FeaturePlot(PRO, features = genes,
            cols = c("lightgrey", "red"),order=TRUE,ncol=2)
    }
    ggsave(p1,file=paste0(outdir,"/",ctype,".umap.png"),width=plotsize[['width']],height=plotsize[['height']])
    ggsave(p1,file=paste0(outdir,"/",ctype,".umap.pdf"),width=plotsize[['width']],height=plotsize[['height']])
    }
}
