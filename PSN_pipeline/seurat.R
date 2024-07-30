.libPaths("/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/scanpy/lib/R/library")
options(bitmapType='cairo')

source("/PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/enrichment.r")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/color/color.R")
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/plot.cellpercent.R")
source("/PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/seurat.plot.R")
source("/PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/findmarker.R")
# source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/loupe.R")
source("/PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/doublet.batch.R")

source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/gdiff.function.r")
SCRIPT <- "/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/stdpipeV5" #TODO

colors = colorls$"NPG"

ParseArg <- function(){
    library("optparse")
    option_list <- list(
        make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="20"),
        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default="False"),
        make_option(c("-t","--type"),help="organism short name,default %default",default="hsa"),
        make_option(c("-g","--gather"),help="remove batch effect methods",default="harmony"),
        make_option(c("-d","--double"),help="remove double",default="False",action="store_true"),
        make_option(c("-i","--immune"),help="whether is immune data",default="FALSE",action="store_true"),
        make_option(c("-l","--avg_log2FC"),help="log2 fold change value cutoff",default="0.25")
    )
    opt_parser=OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)
    return(opt)
}

LoadLibray <- function(){
    suppressMessages(library(future))
    suppressMessages(options(future.globals.maxSize= 1891289600))
    suppressMessages(library(tibble))
    suppressMessages(library(dplyr))
    suppressMessages(library(Seurat))
    suppressMessages(library(cowplot))
    suppressMessages(library(ggplot2))
    suppressMessages(library(DoubletFinder))
    suppressMessages(library(harmony))
    suppressMessages(library(dittoSeq))
    suppressMessages(library(scales))
    suppressMessages(library(data.table))
    suppressMessages(library(jsonlite))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(BiocGenerics))
    suppressMessages(library(Cairo))
		suppressMessages(library(ragg))
		suppressMessages(library(tidyverse))
}

#merge multi sample cellranger result to one seurat object
MergeData <- function(project_dir,sample){
    sample_list=c(sample$sample)
    sample_group=c(sample$group)
    ifnb.list=list()
    for(x in  c(1:nrow(sample))){
        if(opt$immune){
            datadir=paste0(project_dir,"/02_cellranger/",sample[x,][,1],"/count/sample_filtered_feature_bc_matrix")}
        else{
            datadir=paste0(project_dir,"/02_cellranger/",sample[x,][,1],"/filtered_feature_bc_matrix")
        }
        print(datadir)
        single.data<-Read10X(datadir)
        single.ob<-CreateSeuratObject(counts = single.data, project = "scRNA", min.cells = 3, min.features = 200)
        for (i in c(1:length(sample[x,]))) {
            single.ob[[colnames(sample)[i]]]<-sample[x,][,i]
        }
        ifnb.list= c(ifnb.list,list(single.ob))
    }

    if(length(sample_list)>1){
        single.ob=merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_group)}else{
        single.ob=ifnb.list[[1]]
    }
    return(single.ob) 
}

QcFilter <- function(single.ob,project_dir,species,Feature=c("nFeature_RNA", "nCount_RNA", "percent.mt")){ #TODO 增加筛选条件，不要写死
    cell_num_raw=table(single.ob$sample)
    seurat_qc_dir=paste(project_dir,"03_QC_Filter",sep="/")
    if(!file.exists(seurat_qc_dir)){
        dir.create(seurat_qc_dir)
    }
		cellqc_raw_list = list()
		cellqc_filter_list = list()
    if(file.exists(paste0(SCRIPT,"/mtgene/",species,".txt"))){
				subset_condition = glue::glue("nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt <{opt$mt_cutoff} & nCount_RNA<50000")
        mtdf <- read.table(paste0(SCRIPT,"/mtgene/",species,".txt"))
				
        single.ob[["percent.mt"]] <- PercentageFeatureSet(single.ob, features = mtdf$V1)
				subset_condition1 <- eval(parse(text = subset_condition), envir = single.ob@meta.data)
    }else{
				Feature =  Feature[Feature != "percent.mt"]
				subset_condition = 'nFeature_RNA > 400 & nFeature_RNA < 7000 &  nCount_RNA<50000'
    }
		#绘制原始的数据结果图片
		for (i in Feature){
				cellqc_raw_list[[i]] = VlnPlot(single.ob, features = i, ncol = 1,pt.size = 0,group.by="sample",cols=colors)+
        geom_boxplot(width=.2,col="black",fill="white") +theme(plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))
		}
		cellqc.raw <- do.call(plot_grid, c(plot = lapply(cellqc_raw_list, function(p) p + theme(legend.position = 'none')), ncol = length(Feature)))
		#细胞过滤
		subset_condition1 <- eval(parse(text = subset_condition), envir = single.ob@meta.data)
		single.ob <- single.ob[,subset_condition1]
		#绘制过滤后的结果图片
		for (i in Feature){
				cellqc_filter_list[[i]] = VlnPlot(single.ob, features = i, ncol = 1,pt.size = 0,group.by="sample",cols=colors)+
        geom_boxplot(width=.2,col="black",fill="white") +theme(plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))
		}
		cellqc.filter <- do.call(plot_grid, c(plot = lapply(cellqc_filter_list, function(p) p + theme(legend.position = 'none')), ncol = length(Feature)))
    ##output experiment raw information
    out_gene_umi=as.matrix(single.ob@meta.data)
    row_name=c("Cell",rownames(out_gene_umi))
    out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
    out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
    write.table(out_gene_umi,file=paste(seurat_qc_dir,"Basicinfo_nGene_nUMI_mito.csv",sep="/"),sep=",",quote=F,row.names=F,col.names=F)    
    #plot_grid(cellqc.raw,cellqc.filter,nrow=2)
    cellqc.raw/cellqc.filter
    width <- 8+ceiling(length(unique(single.ob$sample))/3)
    ggsave(paste(seurat_qc_dir,"cells_qc_filter.pdf",sep="/"),width = width,height = 7,bg='white')
    ggsave(paste(seurat_qc_dir,"cells_qc_filter.png",sep="/"),width = width,height = 7,bg='white')

    cell_num_filter=table(single.ob$sample)

    filter_summary=cbind(cell_num_raw,cell_num_filter)
    filter_summary=data.frame("cell_num_raw"=filter_summary[,1],"cell_num_filter"=filter_summary[,2])
    filter_summary[nrow(filter_summary)+1,]=apply(filter_summary,2,sum)
    rownames(filter_summary)[nrow(filter_summary)]="total"
    if(length(sample_list)==1){
        rownames(filter_summary)[1]=sample_list[1]
    }

    filter_summary$percent=apply(filter_summary,1,function(x){round(x[2]*100/x[1],2)})
    write.table(filter_summary,paste(seurat_qc_dir,"cells_filter_stat.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    return(single.ob)
}


Cluster <- function(single.ob,sample_list,rmdouble,method,species,project_dir){
    diff_cluster_dir=paste(project_dir,"05_Marker",sep="/")
    seurat_exp_cluster_dir=paste(project_dir,"04_Clustering",sep="/")
    rds_dir=paste(project_dir,"07_Rds",sep="/")
    if(!file.exists(diff_cluster_dir)){dir.create(diff_cluster_dir)}
    if(!file.exists(seurat_exp_cluster_dir)){dir.create(seurat_exp_cluster_dir)}
		
		cluster_overviwe = file.path(seurat_exp_cluster_dir,'2.cluster_overview')
		if(!file.exists(cluster_overviwe)){dir.create(cluster_overviwe)}

    if(!file.exists(rds_dir)){dir.create(rds_dir)}
    ifnb.list <- SplitObject(single.ob, split.by = "sample")
    #ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
    #features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
    #ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
    
    #calculate double
    if(rmdouble){
        ifnb.list <- DoubletFinder(sample_list,ifnb.list,diff_cluster_dir) 
    }
    #remove batch effect
    immune.combined <- BatchEffect(method,species,ifnb.list,sample_list,seurat_exp_cluster_dir)

	cluster_summary=dcast(as.data.frame(table(data.frame("cluster"=Idents(immune.combined),"sample"=immune.combined$sample))),sample~cluster)
	write.table(cluster_summary,paste(cluster_overviwe,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)
  samplegroupdf <- table(immune.combined@meta.data$sample,immune.combined@meta.data$group)
  write.table(samplegroupdf,file=file.path(cluster_overviwe,"sample.group.csv"),sep=",",quote=F,col.names=NA)
  
  #ElbowPlot
	pc_use=20  #可修改，默认20
  Elbowplot = ElbowPlot(immune.combined)$data %>% ggplot() +geom_point(aes(x=dims,y=stdev)) + 
	#geom_vline(xintercept=pc_use,color='darkred')+
	 ggtitle('Elbow plot-quantitative approach')	+xlab("PCs")+ ylab("Standard Deviation")+
	theme(axis.text.x = element_text(size=14,face="bold"),
		axis.text.y = element_text(size=14,face="bold"),
		axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
		plot.title=element_text(hjust=0.5))  
	preprocess=file.path(seurat_exp_cluster_dir,'1.preprocess')
	if(!file.exists(preprocess)){dir.create(preprocess)}

  ggsave(plot=Elbowplot,filename=file.path(preprocess,'pca.ElbowPlot.png'),dpi=1000,width=10,bg='white')
	ggsave(plot=Elbowplot,filename=file.path(preprocess,'pca.ElbowPlot.pdf'),width=10)
  
	#高变基因绘图
	assay <- ifelse(method == 'CCA', 'SCT', 'RNA')
	vargene_plot <- VariableFeaturePlot(immune.combined,assay=assay,cols=c('#999999','#007ACC')) 
	vargene_top10<- head(VariableFeatures(immune.combined,assay=assay), 10)
	if ( all(nchar(vargene_top10) < 10)){
		vargene_plot <- LabelPoints(plot = vargene_plot, points = vargene_top10, repel = TRUE)
	}

	vargene_plot =vargene_plot+theme_bw()+theme(
		panel.grid=element_blank(),
		axis.text.x = element_text(size=14,face="bold"),
		axis.text.y = element_text(size=14,face="bold"),
		axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),aspect.ratio = 1/1) 

	ggsave(plot=vargene_plot,filename=file.path(preprocess,'VariableFeatures_distribution.png'),dpi=1000,bg='white',width=8)
	ggsave(plot=vargene_plot,filename=file.path(preprocess,'VariableFeatures_distribution.pdf'),width=8)

	DefaultAssay(immune.combined) <- "RNA"
  immune.combined <- NormalizeData(immune.combined)


  saveRDS(object=immune.combined,paste(rds_dir,"All_sample_combined.rds",sep="/"))

  markers <- run_cluster(immune.combined,diff_cluster_dir,species,"seurat_clusters",colors,avg_log2FC=opt$avg_log2FC)

    #画图
  Seurat.Plot(immune.combined,colors=colors,seurat_exp_cluster_dir=seurat_exp_cluster_dir,markers=markers,avg_log2FC=opt$avg_log2FC)
	
	Topmarker_plot(immune.combined,markers,diff_cluster_dir,type,'seurat_clusters',colors,avg_log2FC=opt$avg_log2FC)

	if("celltype" %in% colnames(immune.combined@meta.data)){
		collapseby ='celltype'
	}else{
		collapseby ='seurat_clusters'
	}
	topn = 30 #TODO 之后要加到opt里
	Topn_heatmap_marker(immune.combined=immune.combined,markers=markers,collapseby=collapseby,Topn=topn,seurat_exp_cluster_dir=diff_cluster_dir,colors=colors)

	##marker绘图
  return(immune.combined)
}


project_dir=getwd()
if(!file.exists(project_dir)){
   dir.create(project_dir)
}


#read sample info
sample=read.table("sample_info.txt",sep="\t",header=T)
sample_list=c(sample$sample)
sample_group=c(sample$group)

opt <- ParseArg()
LoadLibray()
species <- opt$type
rmdouble <- opt$double
method <- opt$gather
single.ob <- MergeData(project_dir,sample)
single.ob <- QcFilter(single.ob,project_dir,species)
immune.combined <- Cluster(single.ob,sample_list,rmdouble,method,species,project_dir)

#seurat_diff_cluster_dir=paste(project_dir,"06_Diff_Group",sep="/")
rds_dir=paste(project_dir,"07_Rds",sep="/")
#if(!file.exists(seurat_diff_cluster_dir)){
#   dir.create(seurat_diff_cluster_dir)
#}
#groupDiffAuto(immune.combined,seurat_diff_cluster_dir,"seurat_clusters",species,opt$avg_log2FC)
###
# 添加代码以运行LoupeR函数
#LoupeR(immune.combined,rds_dir)
#immune.combined[["RNA2"]] <- as(object = immune.combined[["RNA"]], Class = "Assay")
#immune.combined[["RNA"]] = immune.combined[["RNA2"]]
#immune.combined[["RNA2"]] =NULL
#saveRDS(immune.combined,file.path(rds_dir,"All_sample_combined_v4.rds"))

if (file.exists(file.path(rds_dir,"All_sample_combined.rds"))){
	rds = file.path(rds_dir,"All_sample_combined.rds")
	cmd = glue::glue('Rscript /PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/loupe.R -i {rds} -o {rds_dir} -n loupe_from_seurat')
	system(cmd)
}



