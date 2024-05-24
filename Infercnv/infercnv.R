#����λ�����ݵĴ���
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
argv <- add_argument(argv,"--ref_celltype",help="reference of celltype")
argv <- add_argument(argv,"--subset_celltype",help="interesting celltype",default="F")
argv <- add_argument(argv,"--species",help="species",default="mouse")
argv <- add_argument(argv,"--denoise_adj",help="denoise is TRUE or FALSE",default = TRUE)
argv <- add_argument(argv,"--HMM_adj",help="HMM if TRUE or FALSE",default = TRUE)
argv <- add_argument(argv,"--outdir",help="the path of output",default="./infercnv")
argv <- add_argument(argv,"--plot",help="proceed with subsequent drawing",default="TRUE")
argv <- parse_args(argv)

rdspath = argv$rdspath
cutoff = argv$cutoff
ident <- argv$ident

ref_celltype = argv$ref_celltype
subset_celltype = argv$subset_celltype
outdir = argv$outdir
species = argv$species
denoise_adj = argv$denoise_adj
HMM_adj <- argv$HMM_adj
cutoff <- as.numeric(cutoff)
ref_celltype <- unlist(strsplit(ref_celltype,split=","))
dir.create(outdir)

sc <- readRDS(rdspath)
Idents(sc) <- ident
sc@meta.data$cellname  <- sc@active.ident

print(subset_celltype)
if (subset_celltype == "F"){
print("Don't subset celltype")
} else {
subset_celltype <- unlist(strsplit(subset_celltype,split =","))
sc <- subset(sc,idents = subset_celltype)
}
#
#
#
sc =sc %>% SCTransform()
sc <- subset(sc,downsample=200)
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
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=cutoff,
                             out_dir=paste0(outdir,'/infercnv_out/'),
                             cluster_by_groups=TRUE,
                             denoise=denoise_adj,
                             HMM=HMM_adj,
                             no_prelim_plot = T,
                             num_threads=10,
                             output_format = "pdf")

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=cutoff,
                             out_dir=paste0(outdir,'/infercnv_out/'),
                             cluster_by_groups=TRUE,
                             denoise=denoise_adj,
                             HMM=HMM_adj,
                             num_threads=10,
                             no_prelim_plot = T,
                             output_format = "png")

library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj, #�������õ���infercnv����
                   out_dir = paste0(outdir,'/infercnv_out/'),
                   plot_chr_scale = T, #��Ⱦɫ��ȫ����Ĭ��ֻ�����������õ��ģ�����
                   output_filename = "better_plot",output_format = "pdf", #����Ϊpdf�ļ�
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)))

infercnv::plot_cnv(infercnv_obj, #�������õ���infercnv����
                   plot_chr_scale = T, #��Ⱦɫ��ȫ����Ĭ��ֻ�����������õ��ģ�����
                   out_dir = paste0(outdir,'/infercnv_out/'),
                   output_filename = "better_plot",output_format = "png", #����Ϊpdf�ļ�
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)))


dir.create(paste0(outdir,'/Result/'))
saveRDS(infercnv_obj,paste0(outdir,'/Result/',"inferCNV.rds"))
file.copy(paste0(outdir,'/infercnv_out/',"infercnv.pdf"),paste0(outdir,'/Result/'))
file.copy(paste0(outdir,'/infercnv_out/',"infercnv.png"),paste0(outdir,'/Result/'))
file.copy(paste0(outdir,'/infercnv_out/',"better_plot.png"),paste0(outdir,'/Result/'))
file.copy(paste0(outdir,'/infercnv_out/',"better_plot.pdf"),paste0(outdir,'/Result/'))

if(argv$plot){
    cmd <- paste("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/infercnv/infercnv_vis.r --rds ",rdspath, " --lbj ",outdir," --groupby " ident,"--output infercnv_vis --species ",species,"-m all ")
    system(cmd)    
}









