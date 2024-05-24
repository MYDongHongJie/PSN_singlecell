#!/usr/bin/env R
# coding: utf-8
# author: donghongjie
# time: 2024/04/30
#	email: 979853020@qq.com

library(infercnv)

library(dplyr)
library(optparse)
library(ggplot2)
library(AnnoProbe)


#参数列表
option_list <- list(
  make_option(c("-i", "--rds"), help="rds of inferCNV",default=NULL),
  make_option(c("-l","--lbj"),help="The result path of infercnv",default="./"),
  make_option(c("-o","--output"),help="out dir",default="./"),
  make_option(c("-g","--groupby"),help="The groups to be displayed should be consistent with infercnv",default="celltype"),
  make_option(c("-s","--species"),help="Species information, consistent with infercnv",default="NULL"),
  make_option(c("-p","--platte"),help="Color palette names in color palette scripts",type = "character",default ="NPG"),
  make_option(c("-m","--map"),help="Drawing method, default to all, optional for vlnplot and heatmap,segmentation by ,",default ="all")
  )
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
##颜色的选择
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
group_palette = colorls[[opt$platte]]

#前期数据准备

if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir)
    }
}
output_dir = normalizePath(output_dir )



infercnv_obj = readRDS(file.path(opt$lbj,"infercnv_out/run.final.infercnv_obj"))
seurat_obj = readRDS(opt$rds)
expr <- infercnv_obj@expr.data
#判断一下infercnv结果的Barcode是否全在rds里面
mapping_re = setdiff(colnames(expr),rownames(seurat_obj@meta.data))
if (length(mapping_re) !=0) {
	 print('The barcode of rds in the infercnv result does not fully correspond to the parameter rds')
	 quit()
}

ref <- infercnv_obj@reference_grouped_cell_indices
ref <- ref[[1]]
obser <- infercnv_obj@observation_grouped_cell_indices
obser_merge <- obser[[1]]
for (i in 2:length(obser)){
	obser_merge = c(obser_merge,obser[[i]])
}

anno.df=data.frame(
  CB=c(colnames(expr)[ref],colnames(expr)[obser_merge]),
  class=c(rep("reference",length(ref)),rep("observation",length(obser_merge)))
)


species = opt$species
gn <- rownames(expr)
geneFile <-annoGene(rownames(expr),"SYMBOL",species)

sub_geneFile <-  geneFile[which(geneFile$SYMBOL%in%intersect(gn,geneFile$SYMBOL)),]
sub_geneFile = sub_geneFile
sub_geneFile <- sub_geneFile[!duplicated(geneFile$SYMBOL),]
expr=expr[intersect(gn,geneFile$SYMBOL),]

###加上groupby的参数
groupby = opt$groupby
metadata = seurat_obj@meta.data %>% mutate(CB =rownames(.))
group_data = metadata %>% select(CB,!!rlang::sym(groupby))



#聚类,提取结果
set.seed(99999)


kmeans_df <- data.frame(kmeans_class=NA)
i=0
for (celltype in unique(metadata[[groupby]])){
	print(celltype)
	bar_need = metadata[which(metadata[[groupby]]==celltype),][["CB"]]
	expr_temp = expr[,colnames(expr)%in%bar_need]
	kmeans_result_temp <- kmeans(t(expr_temp), 2)
	kmeans_df_temp <- data.frame(kmeans_class=kmeans_result_temp$cluster)
	if(i!=0){
		kmeans_df_temp$kmeans_class = kmeans_df_temp$kmeans_class+i
	}
	kmeans_df=rbind(kmeans_df,kmeans_df_temp)
	i=i+2
}

kmeans_df$CB=rownames(kmeans_df)
kmeans_df = kmeans_df[-1,]
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df=kmeans_df %>% inner_join(group_data,by="CB")
kmeans_df_s=arrange(kmeans_df,!!rlang::sym(groupby),kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释

if (opt$map=="all"){
	vismethods=c("vlnplot","heatmap")
}else{
	vismethods = unlist(strsplit(opt$map,","))
}
color_g=group_palette[1:length(unique(kmeans_df_s[[groupby]]))]
for ( vismethod in vismethods ){
		if (vismethod == 'vlnplot'){
		###vlnplot
		expr2=expr-1
		expr2=expr2 ^ 2
		CNV_score=as.data.frame(colMeans(expr2))
		colnames(CNV_score)="CNV_score"

		CNV_score$CB=rownames(CNV_score)
		kmeans_df_s$CB=rownames(kmeans_df_s)
		CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
		CNV_score$class = factor(CNV_score$class,levels=c('observation','reference'))
		
		#CNV_score 
		p=ggplot(CNV_score, aes(x=!!rlang::sym(groupby),y=CNV_score))+
			
			geom_violin(aes(fill = !!rlang::sym(groupby),color=!!rlang::sym(groupby)), trim = FALSE)+
			
			geom_boxplot(aes(fill=!!rlang::sym(groupby)),notch = F,width=0.3)+
			
			labs(x="", y="CNV_score")+
			#theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
		scale_fill_manual(values=color_g)+ 
		scale_color_manual(values=color_g)+ 
		theme(axis.text=element_text(colour='black',size=10),axis.text.x=element_text(angle = 45, hjust = 1)
		,panel.background=element_blank(),axis.line = element_line( colour = "black"))

		nlevel=length(unique(CNV_score[,groupby]))

		ggsave(file.path( output_dir,paste0("CNV_vlnboxplot_groupby_",groupby,".pdf",collapse = "") ), 
											plot = p ,
											width = nlevel*2^0.5/2,limitsize = FALSE,bg="white")
		ggsave(file.path( output_dir,paste0("CNV_vlnboxplot_groupby_",groupby,".png",collapse = "") ), 
											plot = p ,dpi=1000,
											width = nlevel*2^0.5/2,limitsize = FALSE,bg="white")
		rownames(CNV_score)=CNV_score$CB
		CNV_score$CB = NULL
		write.table(CNV_score, file = file.path(output_dir,"CNV_score.xls"), quote = FALSE, sep = '\t', row.names = T, col.names = NA)
	}
	if (vismethod == 'heatmap'){
		#定义热图的注释，及配色
		.libPaths("/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/stdpipeV3/lib/R/library/")
		library(tidyverse)
		library(ComplexHeatmap)
		library(circlize)

		top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:length(unique(geneFile$chr)),labels_gp = gpar(cex = 1.2)))

		
		names(color_g) = unique(kmeans_df_s[[groupby]])

		kmeans_df_s$kmeans_class=NULL
		kmeans_df_s$CB=NULL
		left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("observation"="red","reference" = "blue"),Idents=color_g))
		#options(digits = 2)
		bounds=round(c((min(expr)-0.1),median(expr),max(expr)+0.1),2)
		pdf(file.path(output_dir,"inferCNV_heatmap_exprdata.pdf"),width = 12,height = 8)
		ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
										col = colorRamp2(bounds, c("#377EB8","white","#E41A1C")),
										cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
										column_split = factor(sub_geneFile$chr, paste("chr",1:length(unique(geneFile$chr)),sep = "")), 
										column_gap = unit(2, "mm"),
										use_raster =FALSE,
										heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=bounds,legend_height = unit(3, "cm")),
										top_annotation = top_anno,left_annotation = left_anno, #添加注释
										row_title = NULL,column_title = NULL)
		draw(ht, heatmap_legend_side = "right")
		dev.off()
		
		cmd=(glue::glue("convert  -verbose -density 500 -trim  {output_dir}/inferCNV_heatmap_exprdata.pdf  -quality 100  -flatten  {output_dir}/inferCNV_heatmap_exprdata.png"))
		system(cmd)
	}
}