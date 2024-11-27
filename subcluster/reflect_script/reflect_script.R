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

option_list <- list(
  make_option(c("-i","--mainrds"), help="大类rds/主rds"),
  make_option(c("-s","--subrds"),help="subrds,use , split"),
  make_option(c("-o","--out_dir"),help="out_dir",default="./"),
  make_option(c("-t","--type"),help="type eg. hsa mmu ..",default="hsa"),
  make_option(c("-p","--step"),help="need rds or all(rds,umap,marker)",default="rds"),
  make_option(c("-c","--cmpfile"),help="compare file",default=NULL),
  make_option(c("-m","--marker"),help="each celltype marker for plot",type = "character",default =NULL),
  make_option(c("-r","--rmcluster"),help="remove cluster from result",type = "character",default =NULL),
  make_option(c("-a","--avg_log2FC"),help="threshold for group compare foldchange",default =0.1),
	make_option(c("-C","--cover"),help="Does it cover comparative information for differential analysis",type = "logical", default = FALSE),
	make_option(c("-N","--topn"),help="The top number of KEGG pathways",type = "character", default = "20")
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
scuse$barcode <- rownames(scuse@meta.data)
scuse$celltype <- as.character(scuse$celltype)
celltype_levels=unique(scuse$celltype)


# 处理名字并获取最后一个下划线后的部分
scuse$barcode <- as.character(gsub(".*_(\\w+)[-_]\\d+$", "\\1", scuse$barcode))
scuse$barcode = paste0(scuse$sample,"_",scuse$barcode)

for (i in (1:length(subrds))){
  scsub <- readRDS(subrds[i])
  
  scsub$barcode <- rownames(scsub@meta.data)
  scsub$celltype <- as.character(scsub$celltype)
	#修改barcode
	scsub$barcode = as.character(gsub(".*_(\\w+)[-_]\\d+$", "\\1", scsub$barcode))
	scsub$barcode = paste0(scsub$sample,"_",scsub$barcode)

  meda<- scsub@meta.data
	celltype_levels = c(celltype_levels,unique(meda$celltype))
  for (ce in unique(meda$celltype)){
		temp = subset(meda,celltype == ce)
		scuse@meta.data[which(scuse@meta.data$barcode %in% temp$barcode),"celltype"] <- ce
	}
}

scuse$celltype = factor(scuse$celltype,levels = unique(celltype_levels))

if( !is.null(opt$rmcluster)){
	rmcluster = unlist(strsplit(opt$rmcluster,","))
	`%!in%` =Negate(`%in%`)
	scuse <- SetIdent(scuse, value = scuse$celltype)
	scuse = subset(scuse,celltype %!in% rmcluster)
}else{
	scuse <- SetIdent(scuse, value = scuse$celltype)
}
scuse$celltype = droplevels(scuse$celltype)

print("保存rds")
saveRDS(scuse,paste0(out_dir,"/rename_seuratobj.rds"))

print("rds保存完成")
if(step != "rds"){
	if (!is.null(opt$cmpfile)){
		cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/rename_clusterv3.R -r ",paste0(out_dir,"/rename_seuratobj.rds"), " -g celltype -t", opt$type, " -o ",out_dir,"-m ",opt$cmpfile," -a ",opt$avg_log2FC,"-C ",opt$cover,"-N ",opt$topn)
	}else{
		cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/rename_clusterv3.R -r ",paste0(out_dir,"/rename_seuratobj.rds"), " -g celltype -t", opt$type, " -o ",out_dir," -a ",opt$avg_log2FC,"-C ",opt$cover,"-N ",opt$topn)
	}
	print(cmd)
	system(cmd)
}


if(!is.null(opt$marker)){
		print("运行marker绘图")
		cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/gene.dot.R --rds ", paste0(out_dir,"/rename_seuratobj.rds"), " --genefile ",opt$marker," --column  celltype --outdir ",paste0(out_dir,"/marker"))
		system(cmd)
}


