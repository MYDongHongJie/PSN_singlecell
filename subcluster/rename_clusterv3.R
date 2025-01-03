#A scripts to add cell anno to seurat cluster
#.libPaths("/PERSONALBIO/work/singlecell/s00/software/R.library/library")
#.libPaths(c("/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/sc/lib/R/library","/PERSONALBIO/work/singlecell/s00/software/R.library/library"))
options(bitmapType='cairo')
#setwd("~/work/test/")
print("load package...")
library(optparse)
library(Seurat)
#library(clusterProfiler)
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
	make_option(c("-C","--cover"),help="Does it cover comparative information for differential analysis",type = "logical", default = FALSE)
	#make_option(c("-N","--topn"),help="The top number of KEGG pathways",type = "character", default = "20")
  )
#source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/enrichment2.r")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/PSN_pipeline/enrichment.r")
source("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/gdiff.function.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/plot.r")
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

species <- opt$type
#支持的物种列表，看看放在那里比较合适
species_mapping <- list(
  "hsa" = "org.Hs.eg.db", # human
  "mmu" = "org.Mm.eg.db", # mouse
  "rno" = "org.Rn.eg.db", # rat
  "dme" = "org.Dm.eg.db", # fruit_fly
  "dre" = "org.Dr.eg.db", # zebrafish
  "ath" = "org.At.tair.db", # Arabidopsis
  "sce" = "org.Sc.sgd.db", # yeast
  "cel" = "org.Ce.eg.db", # C.elegans
  "bta" = "org.Bt.eg.db", # Bovine
  "mcc" = "org.Mmu.eg.db", # monkey
  "cfa" = "org.Cf.eg.db", # dog
  "ssc" = "org.Ss.eg.db", # pig
  "gga" = "org.Gg.eg.db", # chicken
  "xla" = "org.Xl.eg.db", # frog
  "ptr" = "org.Pt.eg.db", # chimpanzee
  "aga" = "org.Ag.eg.db"  # mosquito
)

if (species %in% names(species_mapping)) {
  # 获取包名
  species_package <- species_mapping[[species]]

  # 加载对应的物种数据库
  suppressMessages(library(species_package, character.only = TRUE))

  # 加载数据库对象
  species.org <- get(species_package)
}else{
	species.org =NULL
}






future::plan("multicore", workers = min(future::availableCores(), opt$ncores))

run_cluster<-function(immune.combined,seurat_exp_cluster_dir,type,idents,colors,topn){
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
		all_dir_enrich=paste(cluster_dir,"enrichment/all",sep="/")
    if(!file.exists(upcluster_dir_enrich)){dir.create(upcluster_dir_enrich,recursive =TRUE)}
		if(!file.exists(all_dir_enrich)){dir.create(all_dir_enrich,recursive =TRUE)}
    cluster_markers=subset(markers,cluster==clust_num)
    rownames(cluster_markers)<-cluster_markers$gene
    if(nrow(cluster_markers)>1){
      #genelist=cluster_markers$gene
      up =subset(cluster_markers,p_val < 0.05 & avg_log2FC > 0.25)
      upgenelist=up$gene
			try(enrichment(species=type,outDir=all_dir_enrich,geneList=cluster_markers$gene))		
      try(enrichment(species=type,outDir=upcluster_dir_enrich,geneList=upgenelist))
			if (!is.null(species.org)) {
				cluster_markers$description= unname(mapIds(x = species.org,
							keys = rownames(cluster_markers),
							column = "GENENAME",
							keytype = "SYMBOL",
							multiVals = "first"))}

			
      write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
			if (!is.null(species.org)) {
				for (enrichment_res in c(all_dir_enrich,upcluster_dir_enrich)) {
					kegg_data = perform_KEGG_enrichment(species.org, cluster_markers$gene, upgenelist, downdiff_genes=NULL, enrichment_res, species)
					if(is.null(kegg_data)){
						print(paste0(enrichment_res,"kegg file is empty"))
					}
				}
			}
      gc(TRUE)
    }
  }
	allmarker_file = paste(seurat_exp_cluster_dir,"allmarkers.xls",sep="/")
	if (!is.null(species.org)) {
				markers$description= unname(mapIds(x = species.org,
							keys = rownames(markers),
							column = "GENENAME",
							keytype = "SYMBOL",
							multiVals = "first"))}
  write.table(markers,allmarker_file,sep="\t",quote=F,row.names=T,col.names=NA)

	#system(glue::glue("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/PlotKEGGNet.R -m  {allmarker_file} -o {seurat_exp_cluster_dir} -s {type} -t marker --topn {topn} " ))
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
		###比对有没有错误的细胞类型
		allCellType <- read.delim("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/celltype_all.xls")
		if(!all(seurat_anno$ABV %in% allCellType$ABV)){
			print("ABV in anno file must be in celltype_all.xls")
			no_match_name = seurat_anno$ABV[!seurat_anno$ABV %in% allCellType$ABV]
			no_match_celltype_name = as.data.frame(no_match_name)
			write.table(no_match_name,file.path(out_dir,"no_match_celltype_name.xls"),quote=F,row.names=F,sep="\t")
		}		
		colnames(seurat_anno) <- c("cluster","anno")
		seurat_anno <-seurat_anno %>% as_tibble() %>% 
  separate_rows(cluster, sep = ",")
		if(length(seurat_anno$cluster)!=length(unique(seurat_anno$cluster))){return(message("错误，cluster重复"))}
		new.cluster.ids <-as.vector(seurat_anno$anno)
		names(new.cluster.ids) <-seurat_anno$cluster
		Idents(seurat_obj) = "seurat_clusters"
		seurat_obj <- subset(seurat_obj,idents=as.vector(seurat_anno$cluster))
		#seurat_obj@active.ident <- seurat_obj$seurat_clusters
		
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

markers <- run_cluster(seurat_obj,seurat_exp_cluster_dir =paste0(out_dir,"/Marker"),type = opt$type, idents = "celltype",colors,topn=opt$topn)
Clustering = file.path(out_dir,"Clustering")
Seurat.Plot(seurat_obj,colors=colors,seurat_exp_cluster_dir=Clustering,markers=markers)


cluster_overviwe_dir = file.path(out_dir,"Clustering/cluster_overview")
table_data <- table(seurat_obj@meta.data$sample, seurat_obj@meta.data$celltype)

# 转换为数据框以便操作
df <- as.data.frame.matrix(table_data)
 
# 计算每个元素占总和的比例
proportion_df <- df / sum(df)
write.table(df,file=paste(cluster_overviwe_dir,"Cluster_sample_statistics.csv",sep = "/"),sep=",",quote=FALSE,col.names=NA)
write.table(proportion_df,file=paste(cluster_overviwe_dir,"Cluster_sample_percent.csv",sep = "/"),sep=",",quote=FALSE,col.names=NA)
print("star diffexp Analysis")
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
print("差异分析结束")

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



