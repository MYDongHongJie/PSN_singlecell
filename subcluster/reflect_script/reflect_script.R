library(optparse)
library(Seurat)
library(clusterProfiler)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(RColorBrewer)
library(paletteer)
library(future)
options(future.globals.maxSize= 10*1024*1024^2)

colors = colorls$"NPG"
option_list <- list(
  make_option(c("-i","--mainrds"), help="大类rds/主rds"),
  make_option(c("-s","--subrds"),help="subrds,use , split"),
  make_option(c("-o","--out_dir"),help="out_dir",default="./"),
  make_option(c("-t","--type"),help="type eg. hsa mmu ..",default="hsa"),
  make_option(c("-p","--step"),help="need rds or all(rds,umap,marker)",default="rds"),
  make_option(c("-c","--cmpfile"),help="compare file",default="NULL"),
  make_option(c("-m","--marker"),help="each celltype marker for plot",type = "character",default ="NULL"),
  make_option(c("-r","--rmcluster"),help="remove cluster from result",type = "character",default ="NULL"),
  make_option(c("-a","--avg_log2FC"),help="threshold for group compare foldchange",default =0.1)
  )
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
mainrds <- opt$mainrds
subrds <- opt$subrds
out_dir <- opt$out_dir
step <- opt$step
type <- opt$type

seurat_exp_cluster_dir <- out_dir
if(!dir.exists(seurat_exp_cluster_dir)){
    dir.create(seurat_exp_cluster_dir,recursive = T)
}
scall <- readRDS(mainrds)
subrds = strsplit(subrds,split = ",")[[1]]
scuse <- scall
for (i in (1:length(subrds))){
  scsub <- readRDS(subrds[i])
  scuse$barcode <- rownames(scuse@meta.data)
  scsub$barcode <- rownames(scsub@meta.data)
  scuse$celltype <- as.character(scuse$celltype)
  scsub$celltype <- as.character(scsub$celltype)
  meda<- scsub@meta.data
  for (j in (1:length(unique(scsub$celltype)))){
  a <- subset(meda,celltype== unique(scsub$celltype)[j])
  scuse@meta.data[rownames(a),]$celltype <- unique(scsub$celltype)[j]
  print(unique(scuse$celltype))
}
}
 


	if(opt$marker != "NULL"){
	    cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s02/software/script/markerplot/gene.dot.R --rds ",paste0(out_dir,"/rename_seuratobj.rds"), " --genefile ",opt$marker," --column  celltype --outdir ",paste0(out_dir,"/marker"))
	    system(cmd)
	}
}

