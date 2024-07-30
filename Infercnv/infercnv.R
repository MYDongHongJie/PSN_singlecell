suppressMessages({
library(argparser)
library(Seurat)
library(infercnv)
library(ggplot2)
library(AnnoProbe)
library(future)
library(phylogram)
library(gridExtra)
library(grid)
library(dendextend)
library(ggthemes)
library(miscTools)
#library(tidyverse)
})
argv <- arg_parser('debug:zhong.chen@personalbio.cn')
argv <- add_argument(argv,"--rdspath",help="the path of rds")
argv <- add_argument(argv,"--ident",help="use ident",default = "seurat_clusters")
argv <- add_argument(argv,"--cutoff",help="infercnv's cutoff",default="0.1")
argv <- add_argument(argv,"--ref_celltype",help="reference of celltype",default="0,1")
argv <- add_argument(argv,"--subset_celltype",help="interesting celltype",default="F")
argv <- add_argument(argv,"--species",help="species",default="mouse")
argv <- add_argument(argv,"--denoise_adj",help="denoise is TRUE or FALSE",default = TRUE)
argv <- add_argument(argv,"--HMM_adj",help="HMM if TRUE or FALSE",default = TRUE)
argv <- add_argument(argv,"--display_cluster",help="display CNV cluster",default = 'F')
argv <- add_argument(argv,"--celltype_to_cluster",help="which celltype_to_cluster",default ='celltype')
argv <- add_argument(argv,"--outdir",help="the path of output",default="./infercnv")
argv <- add_argument(argv,"--analysis_mode",help = "the analysis mode",default = "samples")
argv <- add_argument(argv,"--cluster_by_groups",help = "cluster by groups",default = TRUE)
argv <- add_argument(argv,"--plot",help="proceed with subsequent drawing",default=TRUE)
argv <- add_argument(argv,"--HMM_type",help="HMM Predicted Level Classification,eg i3 or i6",default = "i6")
argv <- add_argument(argv,"--tree",help="HMM Predicted Level Classification,eg i3 or i6",default = FALSE)
argv <- parse_args(argv)

cluster_by_groups=argv$cluster_by_groups
analysis_mode= argv$analysis_mode
rdspath = argv$rdspath
cutoff = argv$cutoff
ident <- argv$ident
display_cluster = argv$display_cluster
celltype_to_cluster = argv$celltype_to_cluster
ref_celltype = argv$ref_celltype
subset_celltype = argv$subset_celltype
outdir = argv$outdir
species = argv$species
denoise_adj = argv$denoise_adj
HMM_adj <- argv$HMM_adj
cutoff <- as.numeric(cutoff)
ref_celltype <- unlist(strsplit(ref_celltype,split=","))
celltype_to_cluster<- unlist(strsplit(celltype_to_cluster,split=","))
dir.create(outdir)

sc <- readRDS(rdspath)
Idents(sc) <- ident
DefaultAssay(sc)="RNA"
print(subset_celltype)
if (subset_celltype == "F"){
print("Don't subset celltype")
} else {
subset_celltype <- unlist(strsplit(subset_celltype,split =","))
sc <- subset(sc,idents = subset_celltype)
}
if(display_cluster == "T"){
sc@meta.data$celltype2<-as.character(sc@meta.data$celltype)
sc@meta.data$celltype2[which(sc@meta.data$celltype2==celltype_to_cluster)]=as.character(sc@meta.data$seurat_clusters)[which(sc@meta.data$celltype2==celltype_to_cluster)]
sc@meta.data$celltype2<-as.factor(sc@meta.data$celltype2)
Idents(sc)<-sc@meta.data$celltype2
} else{
print('display celltype') }
sc@meta.data$cellname  <- sc@active.ident
#
sc <- subset(sc,downsample=200)
#

matrix<-as.matrix(sc@assays$RNA@counts)
groupinfo= data.frame(cellId = colnames(matrix),cellType= sc@meta.data$cellname)

geneInfor=annoGene(rownames(matrix),"SYMBOL",species)

lev <- unique(geneInfor$chr)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]


geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
levels(geneInfor$chr) <- lev

matrix =matrix [rownames(matrix) %in% geneInfor[,1],]
matrix =matrix [match(geneInfor[,1], rownames(matrix)),]

write.table(matrix ,file =paste0(outdir,'/expFile.txt'),sep = '\t',quote = F)
write.table(groupinfo,file = paste0(outdir,'/metaFiles.txt'),sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file = paste0(outdir,'/geneFile.txt'),sep = '\t',quote = F,col.names = F,row.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(outdir,"/expFile.txt"),
                                    annotations_file=paste0(outdir,"/metaFiles.txt"),
                                    delim="\t",
                                    gene_order_file= paste0(outdir,"/geneFile.txt"),
                                    ref_group_names= ref_celltype)
future::plan("multicore",workers=10)
if(argv$tree){
	print("将进行肿瘤进化树的前置infercnv分析，注意该分析的展示形式将基于subcluster，而非samples")
	
	infercnv_obj = infercnv::run(infercnv_obj,cutoff=1,
															out_dir=paste0(outdir,'/infercnv_out/'),
															cluster_by_groups=FALSE,
															plot_steps=FALSE,
															scale_data=T,
															denoise=T,
															noise_filter=0.12,
															analysis_mode='subclusters',
															HMM_type=argv$HMM_type,
															HMM=T,
															num_threads=10,
															tumor_subcluster_partition_method="random_trees",
															no_plot=TRUE,no_prelim_plot=TRUE,output_format =NA)
}else{
	infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=cutoff,
                             out_dir=paste0(outdir,'/infercnv_out/'),
                             cluster_by_groups=cluster_by_groups,
                             denoise=denoise_adj,
                             HMM=HMM_adj,
                             no_prelim_plot = T,
                             num_threads=10,
														 analysis_mode=analysis_mode,
														 HMM_type= argv$HMM_type,
                             output_format = "pdf")
	library(RColorBrewer)
	infercnv::plot_cnv(infercnv_obj, #�������õ���infercnv����
                   out_dir = paste0(outdir,'/infercnv_out/'),
									 #cluster_by_groups=TRUE,
                   plot_chr_scale = T, #��Ⱦɫ��ȫ����Ĭ��ֻ�����������õ��ģ�����
                   output_filename = "better_plot",output_format = "pdf", #����Ϊpdf�ļ�
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)))
	dir.create(paste0(outdir,'/Result/'))
	saveRDS(infercnv_obj,paste0(outdir,'/Result/',"inferCNV.rds"))
	file.copy(paste0(outdir,'/infercnv_out/',"infercnv.pdf"),paste0(outdir,'/Result/'))
	file.copy(paste0(outdir,'/infercnv_out/',"better_plot.pdf"),paste0(outdir,'/Result/'))

	cmd=(glue::glue("convert  -verbose -density 500 -trim  {outdir}/Result/infercnv.pdf  -quality 500  -flatten  {outdir}/Result/infercnv.png"))
	system(cmd)
	cmd=(glue::glue("convert  -verbose -density 500 -trim  {outdir}/Result/better_plot.pdf  -quality 500  -flatten  {outdir}/Result/better_plot.png"))
	system(cmd)							 
}


if(argv$plot){
    cmd <- paste(
        "Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/infercnv/infercnv_vis.r",
        "--rds", rdspath,
        "--lbj", outdir,
        "--groupby", ident,
        "--output infercnv_vis",
        "--species", species,
        "-m all"
    )
    system(cmd)
}

if(argv$tree){
    cmd <- paste(
        "python /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/Tumor_Evolutionary_Tree/nith.py",
        "-i", paste0(outdir, '/infercnv_out/'),
        "-o", "Tumor_Evolutionary_Tree"
    )
    system(cmd)

    cp_cmd <- paste(
        "cp", paste0(outdir, "/infercnv_out/17_HMM_pred*.pred_cnv_*.dat"),
        "Tumor_Evolutionary_Tree"
    )
    system(cp_cmd)
		readme_cmd <- "cp /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/Tumor_Evolutionary_Tree/肿瘤进化结果说明.doc  Tumor_Evolutionary_Tree"
		system(readme_cmd)
}






