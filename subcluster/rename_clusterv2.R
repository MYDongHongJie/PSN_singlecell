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
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/seurat.plot.R")
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/cloud.R")
colors = colorls$"NPG"

option_list <- list(
  make_option(c("-r", "--rds"), help="",default="All_sample_combined.rds"),
  make_option(c("-c","--cluster"),help="cluster_file",default="cluster.txt"),
  make_option(c("-o","--out"),help="out dir",default="cluster_rename_dir"),
  make_option(c("-t","--type"),help="type eg. hsa mmu ..",default="hsa"),
  make_option(c("-p","--cmpfile"),help="compare file",default="NULL"),
  make_option(c("-m","--marker"),help="each celltype marker for plot",type = "character",default ="NULL"),
  make_option(c("-a","--avg_log2FC"),help="threshold for group compare foldchange",default =0.1),
  make_option(c("-l","--cloud"),help="produce cloud data",action = "store_true", default = FALSE)
  )
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/enrichment2.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/stdpipeV4/enrichment.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/gdiff.function.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/plot.r")
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

run_cluster<-function(immune.combined,seurat_exp_cluster_dir,type,idents,colors){
  if(!dir.exists(seurat_exp_cluster_dir)){
    dir.create(seurat_exp_cluster_dir,recursive = T)
  }
  DefaultAssay(immune.combined) <- "RNA"
  Idents(immune.combined)<-idents
  immune.combined <- NormalizeData(immune.combined) 
  immune.combined <- ScaleData(immune.combined, feature=rownames(immune.combined),verbose = FALSE)
  
  markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10)
  
  all_top10_markers=markers %>%  top_n(n = 50, wt = avg_log2FC) %>%dplyr::distinct(.,gene,.keep_all = T) %>% top_n(n = 10, wt = avg_log2FC)
  
  height= ceiling(length(all_top10_markers$gene)/5)*3 
  for( clust_num in  unique(Idents(immune.combined))){
    cluster_dir=paste(seurat_exp_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
    if(!file.exists(cluster_dir)){
      dir.create(cluster_dir)
    }
    upcluster_dir_enrich=paste(cluster_dir,"enrichment/up",sep="/")
    downcluster_dir_enrich=paste(cluster_dir,"enrichment/down",sep="/")
    if(!file.exists(upcluster_dir_enrich)){dir.create(upcluster_dir_enrich,recursive =TRUE)}
    if(!file.exists(downcluster_dir_enrich)){dir.create(downcluster_dir_enrich,recursive =TRUE)}
    cluster_markers=subset(markers,cluster==clust_num)
    rownames(cluster_markers)<-cluster_markers$gene
    if(nrow(cluster_markers)>1){
      #genelist=cluster_markers$gene
      up =subset(cluster_markers,p_val < 0.05 & avg_log2FC > 0.25)
      down =subset(cluster_markers,p_val < 0.05 & avg_log2FC < -0.25)
      upgenelist=up$gene
      downgenelist=down$gene
      try(enrichment(species=type,outDir=upcluster_dir_enrich,geneList=upgenelist))
      try(enrichment(species=type,outDir=downcluster_dir_enrich,geneList=downgenelist))
      write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
      top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)
      VlnPlot(immune.combined, features = top10_markers$gene,pt.size = 0 ,ncol=5,cols=colors)
      
      ggsave(paste(cluster_dir,"top10_vilion.pdf",sep="/"),width =20,height = height,bg='#ffffff')
      ggsave(paste(cluster_dir,"top10_vilion.png",sep="/"),width =20,height = height,bg='#ffffff')
      
      FeaturePlot(immune.combined, features = top10_markers$gene, min.cutoff = "q9",ncol=5,order=T,cols=c("lightgrey", "red"))
      ggsave(paste(cluster_dir,"top10_umap.pdf",sep="/"),width = 20,height = height,bg='#ffffff')
      ggsave(paste(cluster_dir,"top10_umap.png",sep="/"),width = 20,height = height,bg='#ffffff')
      DotPlot(immune.combined, features = rev(unique(top10_markers$gene)))
      ggsave(paste(cluster_dir,"top10_dotplot.pdf",sep="/"),width = 15,height = height,bg='#ffffff')
      ggsave(paste(cluster_dir,"top10_dotplot.png",sep="/"),width = 15,height = height,bg='#ffffff')
      gc(TRUE)
    }
  }
  write.table(markers,paste(seurat_exp_cluster_dir,"allmarkers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
  return(markers)
}

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj@active.ident <- seurat_obj$seurat_clusters
seurat_obj@meta.data$sample <- factor(seurat_obj@meta.data$sample,levels=unique(seurat_obj@meta.data$sample))
out_dir=opt$out
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

#load anno data
seurat_anno <- read.table(opt$cluster,sep="\t",header=F,col.names = c("cluster","anno"))
seurat_anno <-seurat_anno %>% as_tibble() %>% 
  separate_rows(cluster, sep = ",")
if(length(seurat_anno$cluster)!=length(unique(seurat_anno$cluster))){return(message("错误，cluster重复"))}
seurat_anno

new.cluster.ids <-as.vector(seurat_anno$anno)
names(new.cluster.ids) <-seurat_anno$cluster

print(seurat_anno$cluster)
print(unique(seurat_obj@active.ident))
seurat_obj <- subset(seurat_obj,idents=as.vector(seurat_anno$cluster))
#seurat_obj@active.ident <- seurat_obj$seurat_clusters
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$celltype <- Idents(seurat_obj)

saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
markers <- run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(out_dir,"/Diff_Cluster"),type = opt$type, idents = "celltype",colors)
Seurat.Plot(seurat_obj,colors=colors,seurat_exp_cluster_dir=out_dir,markers=markers)

sample_list=unique(seurat_obj$sample)

#run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(out_dir,"/Diff_Cluster"),type = opt$type, idents = "celltype")

p1 = plot.clusters.group(data = seurat_obj,clusters =  "celltype", xlab = "Cluster number", log =TRUE, group = "group",legend.title = "Group",widths = c(3,1),color = colors)
ggsave(p1,filename = paste(out_dir,"Cluster_group_percent.pdf",sep = "/"))
ggsave(p1,filename = paste(out_dir,"Cluster_group_percent.png",sep = "/"))
p2 = plot.clusters.group(data = seurat_obj,clusters =  "celltype", xlab = "Cluster number", log =TRUE, group = "sample",legend.title = "Sample",widths = c(3,1),color = colors)
ggsave(p2,filename = paste(out_dir,"Cluster_sample_percent.pdf",sep = "/"))
ggsave(p2,filename = paste(out_dir,"Cluster_sample_percent.png",sep = "/"))
DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- seurat_obj@meta.data$celltype
#saveRDS(seurat_obj,file = paste0(out_dir,"/rename_seuratobj.rds"))
#groupDiffAuto(seurat_obj,opt$out,"celltype",opt$type)
if(opt$cmpfile!= "NULL"){
    groupDiffSpeci(seurat_obj,opt$out,"celltype",opt$cmpfile,opt$type,opt$avg_log2FC)
}else{
    groupDiffAuto(seurat_obj,opt$out,"celltype",opt$type,opt$avg_log2FC)
}
cmd <- paste0("cp ",opt$cluster," ",out_dir)
system(cmd)
write.table(table(seurat_obj@meta.data$sample,seurat_obj@meta.data$celltype),file=paste(out_dir,"Cluster_sample_percent.csv",sep = "/"),sep=",",quote=FALSE,col.names=NA)
system(paste0("cp /PERSONALBIO/work/singlecell/s00/software/script/README/添加细胞标签结果反馈说明.p* ",out_dir))

if(opt$marker != "NULL"){
    cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/gene.dot.R --rds ",paste0(out_dir,"/rename_seuratobj.rds"), " --genefile ",opt$marker," --column  celltype --outdir ",paste0(out_dir,"/marker"))
    system(cmd)    
}
if(opt$cloud){
    Cloud(out_dir,seurat_obj,markers)
}
