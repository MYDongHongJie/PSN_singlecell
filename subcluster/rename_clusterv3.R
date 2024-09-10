#A scripts to add cell anno to seurat cluster
#.libPaths("/PERSONALBIO/work/singlecell/s00/software/R.library/library")
#.libPaths(c("/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/sc/lib/R/library","/PERSONALBIO/work/singlecell/s00/software/R.library/library"))
options(bitmapType='cairo')
#setwd("~/work/test/")
print("load package...")
library(optparse)
library(Seurat)
library(clusterProfiler)
library(reshape2)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(RColorBrewer)
library(paletteer)
library(future)
options(future.globals.maxSize= 1891289600)
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/plot.r")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/seurat.plot.R")
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/cloud.R")
colors = colorls$"NPG"

option_list <- list(
  make_option(c("-r", "--rds"), help="",default="All_sample_combined.rds"),
  make_option(c("-c","--cluster"),help="cluster_file",default=NULL),
  make_option(c("-o","--out"),help="out dir",default="cluster_rename_dir"),
  make_option(c("-t","--type"),help="type eg. hsa mmu ..",default="hsa"),
  make_option(c("-m","--cmpfile"),help="compare file",default=NULL),
  make_option(c("-n","--ncores"),help="each celltype marker for plot",type = "integer",default =10),
  make_option(c("-a","--avg_log2FC"),help="threshold for group compare foldchange",default =0.25),
  make_option(c("-l","--cloud"),help="produce cloud data",action = "store_true", default = FALSE),
	make_option(c("-g","--groupby"),help="Specify a column for direct analysis. If there is already a celltype column, change it to old_celltype",default = 'celltype'),
	make_option(c("-C","--cover"),help="Does it cover comparative information for differential analysis",action = "store_true", default = TRUE)
  )
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/enrichment2.r")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/PSN_pipeline/enrichment.r")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/gdiff.function.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/plot.r")
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
future::plan("multicore", workers = min(future::availableCores(), opt$ncores))

run_cluster<-function(immune.combined,seurat_exp_cluster_dir,type,idents,colors){
  if(!dir.exists(seurat_exp_cluster_dir)){
    dir.create(seurat_exp_cluster_dir,recursive = T)
  }
  DefaultAssay(immune.combined) <- "RNA"
  Idents(immune.combined)<-idents
  # immune.combined <- NormalizeData(immune.combined) 
  # immune.combined <- ScaleData(immune.combined, feature=rownames(immune.combined),verbose = FALSE)
  
  markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10)
  
  all_top10_markers=markers %>%  top_n(n = 50, wt = avg_log2FC) %>%dplyr::distinct(.,gene,.keep_all = T) %>% top_n(n = 10, wt = avg_log2FC)
  
  for( clust_num in  unique(Idents(immune.combined))){
    cluster_dir=file.path(seurat_exp_cluster_dir,'Each_celltype_marker',paste("cluster",clust_num,sep="_"))
    if(!file.exists(cluster_dir)){
      dir.create(cluster_dir,recursive = TRUE)
    }
    upcluster_dir_enrich=paste(cluster_dir,"enrichment/up",sep="/")
    downcluster_dir_enrich=paste(cluster_dir,"enrichment/down",sep="/")
		all_dir_enrich=paste(cluster_dir,"enrichment/all",sep="/")
    if(!file.exists(upcluster_dir_enrich)){dir.create(upcluster_dir_enrich,recursive =TRUE)}
    if(!file.exists(downcluster_dir_enrich)){dir.create(downcluster_dir_enrich,recursive =TRUE)}
		if(!file.exists(all_dir_enrich)){dir.create(all_dir_enrich,recursive =TRUE)}
    cluster_markers=subset(markers,cluster==clust_num)
    rownames(cluster_markers)<-cluster_markers$gene
    if(nrow(cluster_markers)>1){
      #genelist=cluster_markers$gene
      up =subset(cluster_markers,p_val < 0.05 & avg_log2FC > 0.25)
      down =subset(cluster_markers,p_val < 0.05 & avg_log2FC < -0.25)
      upgenelist=up$gene
      downgenelist=down$gene
			try(enrichment(species=type,outDir=all_dir_enrich,geneList=cluster_markers$gene))
      try(enrichment(species=type,outDir=upcluster_dir_enrich,geneList=upgenelist))
      try(enrichment(species=type,outDir=downcluster_dir_enrich,geneList=downgenelist))
      write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
      gc(TRUE)
    }
  }
  write.table(markers,paste(seurat_exp_cluster_dir,"allmarkers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
  return(markers)
}

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
DefaultAssay(seurat_obj) <- "RNA"

if (seurat_obj@version >=5){
		try(seurat_obj <- JoinLayers(seurat_obj))
}

if (is.null(seurat_obj@assays$RNA["data"])) {
  seurat_obj <- NormalizeData(seurat_obj)
}

  
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
seurat_obj@active.ident <- seurat_obj$seurat_clusters
seurat_obj@meta.data$sample <- factor(seurat_obj@meta.data$sample,levels=unique(seurat_obj@meta.data$sample))
out_dir=opt$out
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

#load anno data
#seurat_anno <- read.table(opt$cluster,sep="\t",header=F,col.names = c("cluster","anno"))
if (is.null(opt$cluster) && is.null(opt$groupby)){
	print("no anno file and no groupby")
	quit()
}else if (!is.null(opt$cluster) ){
  seurat_data <- read.delim(opt$cluster)
	if (!all(c("ABV","cluster") %in% colnames(seurat_data))){
	print("ABV and cluster must be in the first two columns of the anno file")
	quit()
	}else {
		if (!"CellType" %in% colnames(seurat_data) ) {print("CellType need to be in the anno file")}
		selected_cols <- c("cluster","ABV")
		seurat_anno <- seurat_data[,selected_cols]
		colnames(seurat_anno) <- c("cluster","anno")
		seurat_anno <-seurat_anno %>% as_tibble() %>% 
  separate_rows(cluster, sep = ",")
		if(length(seurat_anno$cluster)!=length(unique(seurat_anno$cluster))){return(message("错误，cluster重复"))}
		new.cluster.ids <-as.vector(seurat_anno$anno)
		names(new.cluster.ids) <-seurat_anno$cluster
		seurat_obj <- subset(seurat_obj,idents=as.vector(seurat_anno$cluster))
		#seurat_obj@active.ident <- seurat_obj$seurat_clusters
		Idents(seurat_obj) = "seurat_clusters"
		seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
		seurat_obj$celltype <- Idents(seurat_obj)
		if (!opt$cover){
			saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
		}
		if( c("marker") %in%  colnames(seurat_data)){
			seurat_data$marker = paste0(seurat_data$marker,",")
		}
		write.table(seurat_data,file.path(out_dir,"anno.xls"),row.names = F,quote=F,sep="\t")
		# cmd <- paste0("cp ",opt$cluster," ",out_dir)
		# system(cmd)
	}
}else{
	if(opt$groupby !='celltype' && 'celltype' %in% colnames(seurat_obj@meta.data)){
		seurat_obj@meta.data$celltype = seurat_obj@meta.data$old_celltype
		seurat_obj@meta.data$celltype = seurat_obj@meta.data[,opt$groupby]
	}else if( !opt$groupby %in% colnames(seurat_obj@meta.data)){
		quit('opt$groupby in not in colnames(seurat_obj@meta.data')
	}else {
		 seurat_obj@meta.data$celltype = seurat_obj@meta.data[,opt$groupby]
	}
}


Idents(seurat_obj) <- seurat_obj@meta.data$celltype

markers <- run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(out_dir,"/Marker"),type = opt$type, idents = "celltype",colors)
Clustering = file.path(out_dir,"Clustering")
Seurat.Plot(seurat_obj,colors=colors,seurat_exp_cluster_dir=Clustering,markers=markers)

#sample_list=unique(seurat_obj$sample)

#run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(out_dir,"/Diff_Cluster"),type = opt$type, idents = "celltype")
cluster_overviwe_dir = file.path(out_dir,"Clustering/cluster_overview")
write.table(table(seurat_obj@meta.data$sample,seurat_obj@meta.data$celltype),file=paste(cluster_overviwe_dir,"Cluster_sample_percent.csv",sep = "/"),sep=",",quote=FALSE,col.names=NA)
# p1 = plot.clusters.group(data = seurat_obj,clusters =  "celltype", xlab = "Cluster number", log =TRUE, group = "group",legend.title = "Group",widths = c(3,1),color = colors)
# ggsave(p1,filename = paste(cluster_overviwe_dir,"Cluster_group_percent.pdf",sep = "/"))
# ggsave(p1,filename = paste(cluster_overviwe_dir,"Cluster_group_percent.png",sep = "/"))
# p2 = plot.clusters.group(data = seurat_obj,clusters =  "celltype", xlab = "Cluster number", log =TRUE, group = "sample",legend.title = "Sample",widths = c(3,1),color = colors)
# ggsave(p2,filename = paste(cluster_overviwe_dir,"Cluster_sample_percent.pdf",sep = "/"))
# ggsave(p2,filename = paste(cluster_overviwe_dir,"Cluster_sample_percent.png",sep = "/"))


#saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
#groupDiffAuto(seurat_obj,opt$out,"celltype",opt$type)
if(!is.null(opt$cmpfile)){
		if(file.exists(opt$cmpfile)){
			groupDiffSpeci(seurat_obj,opt$out,"celltype",opt$cmpfile,opt$type,opt$avg_log2FC)
			if (opt$cover){
				cmpfile = read.delim(opt$cmpfile,header = T,sep = "\t")
				cmpfile$diff = paste0(cmpfile[,1],"/",cmpfile[,2])
				Misc(object = seurat_obj, slot = "cmplist") = cmpfile[["diff"]]
				saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
			}
		}else{
			cmpfile = unlist(strsplit(opt$cmpfile,","))
			groupDiffSpeci_Auto(seurat_obj,opt$out,"celltype",cmpfile,opt$type,opt$avg_log2FC)
			if (opt$cover){ 
				Misc(object = seurat_obj, slot = "cmplist") =cmpfile
				saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
			}
		}
}else if (!is.null(Misc(object = seurat_obj, slot = "cmplist"))) {
	 cmpfile = Misc(object = seurat_obj, slot = "cmplist")
	 groupDiffSpeci_Auto(seurat_obj,opt$out,"celltype",cmpfile,opt$type,opt$avg_log2FC)
}else{
    groupDiffAuto(seurat_obj,opt$out,"celltype",opt$type,opt$avg_log2FC)
}


#system(paste0("cp /PERSONALBIO/work/singlecell/s00/software/script/README/添加细胞标签结果反馈说明.p* ",out_dir))
if (!is.null(opt$cluster)){
	if( c("marker") %in%  colnames(seurat_data)){
    cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/gene.dot.R --rds ",paste0(out_dir,"/rename_seuratobj.rds"), " --genefile ",opt$cluster," --column  celltype --outdir ",paste0(out_dir,"/Ident_CellType_Markers"))
    system(cmd)    
	}
}

if(opt$cloud){
    Cloud(out_dir,seurat_obj,markers)
}
