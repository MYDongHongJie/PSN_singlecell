#run subcluster
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
library(scales)
library(harmony)
library(tibble)
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/enrichment2.r")
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/cloud.R")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/PSN_pipeline/enrichment.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/plot.r")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/recluster_seurat_plot.R")

colors = colorls$"NPG"
option_list <- list(
  make_option(c("-r", "--rds"), help="",default="All_sample_combined.rds"),
  make_option(c("-c","--cluster"),help="cluster name,split by ,"),
  make_option(c("-i","--idents"),help="colname of celltype in meta.data ,",default="celltype"),
  make_option(c("-s","--res"),help="find cluster resolution",default=0.4),
  make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory",default = "Reclust"),
  make_option(c("-t", "--type"), help="species type",default = "hsa"),
  make_option(c("-b","--batch"),help="remove batch effect methods",default="harmony")
)

opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

print("load seurat data")
seurat_obj <- readRDS(opt$rds)
DefaultAssay(seurat_obj) <- "RNA"
if (seurat_obj@version >=5){
		try(seurat_obj <- JoinLayers(seurat_obj))
}
if (!is.null(opt$cluster)){
	seurat_obj@meta.data$celltype <- seurat_obj@meta.data[,opt$idents]
	cluster <- unlist(strsplit(opt$cluster,split=","))
	print("subset data")
	seurat_obj <-subset(seurat_obj,celltype %in% cluster)
}



ifnb.list <- SplitObject(seurat_obj, split.by = "sample")
#ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
#features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
#ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
sample_list <- unique(seurat_obj@meta.data$sample)
if(opt$batch=="harmony"){
    seurat_obj<-merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)])%>% JoinLayers()%>%NormalizeData()
    VariableFeatures(seurat_obj) <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
    #seurat_obj<-subset(seurat_obj,Double_status=="Singlet")
    seurat_obj<-ScaleData(seurat_obj,feature=rownames(seurat_obj)) %>% RunPCA(verbose = FALSE) %>% RunHarmony( group.by.vars = "sample")
    seurat_obj<-RunUMAP(seurat_obj, reduction = "harmony", dims = 1:20)
    seurat_obj<-RunTSNE(seurat_obj, reduction = "harmony", dims = 1:20)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = opt$res)    
}else if(opt$batch=="cca"){
    ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
    ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
    immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",anchor.features = features)
    seurat_obj <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
    DefaultAssay(seurat_obj) <- "integrated"
    seurat_obj <- RunPCA(seurat_obj,  features = VariableFeatures(object = seurat_obj) ,verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
    seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:20)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:20)
    seurat_obj <- FindClusters(seurat_obj,resolution = opt$res)
}else{
		seurat_obj <- seurat_obj%>%NormalizeData()
    seurat_obj <- ScaleData(seurat_obj,feature=rownames(seurat_obj), verbose = FALSE)
    seurat_obj <- FindVariableFeatures(object = seurat_obj,selection.method = 'vst', nfeatures = 3000)
    seurat_obj <- RunPCA(seurat_obj,  features = VariableFeatures(object = seurat_obj) ,verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
    seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:20)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:20)
    seurat_obj <- FindClusters(seurat_obj,resolution = opt$res)
}
file_out=opt$workdir
#plot fig
if(dir.exists(file_out)){
  print("dir exists")
}else{
  dir.create(file_out)
}
cluster_summary=dcast(as.data.frame(table(data.frame("cluster"=Idents(seurat_obj),"sample"=seurat_obj$sample))),sample~cluster)
write.table(cluster_summary,paste(file_out,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)
#sample_group <- table(seurat_obj@meta.data$sample,seurat_obj@meta.data$group)
#write.table(sample_group,file=paste(file_out,"sample.group.csv"),sep=",",quote=F,col.names=NA)

DefaultAssay(seurat_obj) <- "RNA"
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
# avg_per_cluster=AverageExpression(seurat_obj,"RNA")$RNA
# colnames(avg_per_cluster)=gsub("RNA.","",colnames(avg_per_cluster),perl=T)
# write.table(avg_per_cluster,paste(file_out,"avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)


#find marker
print("Do find markers")
file_out_marker <- paste0(file_out,"/Marker")
#plot fig
if(dir.exists(file_out_marker)){
  print("dir exists")
}else{
  dir.create(file_out_marker)
}
#DefaultAssay(seurat_obj) <- "RNA"
#seurat_obj <- NormalizeData(seurat_obj)
#seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

run_cluster<-function(immune.combined,seurat_exp_cluster_dir,type,idents,colors){
  if(!dir.exists(seurat_exp_cluster_dir)){
    dir.create(seurat_exp_cluster_dir,recursive = T)
  }
  DefaultAssay(immune.combined) <- "RNA"
  Idents(immune.combined)<-idents
  # immune.combined <- NormalizeData(immune.combined) 
  # immune.combined <- ScaleData(immune.combined, feature=rownames(immune.combined),verbose = FALSE)
  
  markers <- FindAllMarkers(immune.combined, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10)
  
  #all_top10_markers=markers %>%  top_n(n = 50, wt = avg_log2FC) %>%dplyr::distinct(.,gene,.keep_all = T) %>% top_n(n = 10, wt = avg_log2FC)
  
	return(markers)
}
markers <- run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(file_out,"/Marker"),type = opt$type, idents = "seurat_clusters",colors)

# markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(file_out_marker,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

#Seurat.Plot(seurat_obj,colors=colors,seurat_exp_cluster_dir=file_out,markers=markers)

saveRDS(object = seurat_obj,file = paste0(file_out,"/sub.rds"))
if (file.exists(file.path(file_out,"sub.rds"))){
	rds = file.path(file_out,"sub.rds")
	cmd = glue::glue('Rscript /PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/loupe.R -i {rds} -o {file_out} -n loupe_from_seurat')
	system(cmd)
}
