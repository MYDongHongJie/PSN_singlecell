library(Seurat)
library(dplyr)
library(ggplot2)
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
colors = colorls$"NPG"
option_list <- list(
  make_option(c("-i", "--rds"), help="input of rds file",default=NULL),
  make_option(c("-o","--out"),help="out dir",default="ITHscore"),
  make_option(c("-g","--groupby"),help="groupby",default="group,sample,celltype,seurat_clusters"),
  make_option(c("-p","--pcs"),help="The num of pc", type="integer",default=20)
)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

data = readRDS(opt$rds)
oudir =opt$out
if (!file.exists(oudir)){dir.create(oudir,recursive = T)}
groupby = unlist(strsplit(opt$groupby,","))


data=Seurat::FindVariableFeatures(data, loess.span = 0.3,
                              clip.max = "auto", mean.function = "FastExpMean",
                              dispersion.function = "FastLogVMR", num.bin = 20,
                              nfeature = 3000, binning.method = "equal_width" )

data = data %>% ScaleData()%>% RunPCA(npcs = opt$pcs)


data_d <- dist(data@reductions[["pca"]]@cell.embeddings, method = "euclidean")
distance_matrix <- as.matrix(data_d)
rm(data_d)

# 计算平均距离
ITH = rowSums(distance_matrix)/(ncol(distance_matrix)-1)
ITH_score =as.data.frame(ITH) 

data <- AddMetaData(data,metadata=ITH_score)
metadata = data@meta.data
#构建分组信息
metadata$Barcode = rownames(metadata)

metadata_select = metadata %>% select(Barcode,all_of(groupby),ITH)


#colnames(metadata_select)[1] <- "Barcode"
write.table(metadata_select,file.path(oudir,"ITHscore.xls"),sep="\t",quote=F,row.names=F,col.names = T)


#绘图
for (ident in groupby){
	Idents(data) = ident
	nfeatureplot = FeaturePlot(data, reduction = "umap", label = T,repel=TRUE,features='ITH',cols=c("lightgrey", "red"))
	ggsave(plot=nfeatureplot,filename=file.path(oudir,paste0("ITHscore_featureplot_labelby_",ident,".pdf")))
	ggsave(plot=nfeatureplot,filename=file.path(oudir,paste0("ITHscore_featureplot_labelby_",ident,".png")))
}


for (parameter in groupby){
	len =length(unique(metadata_select[[parameter]]))
	colors_use = colors[1:len]
	p=ggplot(metadata_select, # 用来画图的数据集名称
				aes_string(x=parameter, # 数据集的Species列作为横坐标
						y="ITH",  # 数据集的value列作为纵坐标
						color = parameter, # 根据数据集的Species列设置箱线图的边框颜色
						fill = parameter))+ # 根据数据集的Species列设置散点图和箱线图的填充颜色
		geom_boxplot(width=.4, alpha=.2)+ # 将箱线图的宽度设置为0.4，缩窄了一些，不透明度设置为0.2，填充颜色变得很浅
		scale_color_manual(values = colors_use)+ # 自定义边框颜色
		scale_fill_manual(values = colors_use)+ # 自定义填充颜色，可以和边框颜色不一样，我觉得配套的比较好看，就用了一样的颜色
		theme_test()+ # 把背景颜色去掉并添加边框，可以修改为其他背景，比如theme_classic
		labs(y="ITH Score", x="")+ # 修改横纵坐标名称
		theme(strip.background = element_blank(), #把头上的灰框去掉
					axis.text = element_text(color="black", size=10),
					axis.line = element_blank())
	if (parameter != "seurat_clusters"){
    p=p+theme(axis.text.x =element_text(angle = 45, hjust = 1))
  }
   
	ggsave(plot=p,filename=file.path(oudir,paste0("ITHscore_groupby",parameter,".pdf")),width=0.1*len+5.8)
	ggsave(plot=p,filename=file.path(oudir,paste0("ITHscore_groupby",parameter,".png")),width=0.1*len+5.8)
}

#分组绘图
if ( all(c("sample","group") %in% groupby)){
	len =length(unique(metadata_select[["sample"]]))
	colors_use = colors[1:len]
	p=ggplot(metadata_select, # 用来画图的数据集名称
				aes_string(x="group", # 数据集的Species列作为横坐标
						y="ITH",  # 数据集的value列作为纵坐标
						color = "sample", # 根据数据集的Species列设置箱线图的边框颜色
						fill = "sample"))+ # 根据数据集的Species列设置散点图和箱线图的填充颜色
		geom_boxplot(width=.4, alpha=.2)+ # 将箱线图的宽度设置为0.4，缩窄了一些，不透明度设置为0.2，填充颜色变得很浅
		scale_color_manual(values = colors_use)+ # 自定义边框颜色
		scale_fill_manual(values = colors_use)+ # 自定义填充颜色，可以和边框颜色不一样，我觉得配套的比较好看，就用了一样的颜色
		theme_test()+ # 把背景颜色去掉并添加边框，可以修改为其他背景，比如theme_classic
		labs(y="ITH Score", x="")+ # 修改横纵坐标名称
		theme(strip.background = element_blank(), #把头上的灰框去掉
					axis.text = element_text(color="black", size=10),
					axis.line = element_blank())
	if (parameter != "seurat_clusters"){
    p=p+theme(axis.text.x =element_text(angle = 45, hjust = 1))
  }
	ggsave(plot=p,filename=file.path(oudir,"ITHscore_groupby_sample_spiltby_group.pdf"),width=0.1*len+4.8)
	ggsave(plot=p,filename=file.path(oudir,"ITHscore_groupby_sample_spiltby_group.png"),width=0.1*len+4.8)
}

for (parameter1 in c("sample","group")){
	print(parameter1)
	len1 = length(unique(metadata_select[[parameter1]]))
	if((parameter1 %in% groupby)){
		for (parameter2 in setdiff(groupby,c("sample","group"))){
		len2 =length(unique(metadata_select[[parameter2]]))
    colors_use = colors[1:len2]
    p=ggplot(metadata_select, # 用来画图的数据集名称
				aes_string(x=parameter1, # 数据集的Species列作为横坐标
						y="ITH",  # 数据集的value列作为纵坐标
						color = parameter2, # 根据数据集的Species列设置箱线图的边框颜色
						fill = parameter2))+ # 根据数据集的Species列设置散点图和箱线图的填充颜色
			geom_boxplot(width=.6, alpha=.2)+ # 将箱线图的宽度设置为0.4，缩窄了一些，不透明度设置为0.2，填充颜色变得很浅
			scale_color_manual(values = colors_use)+ # 自定义边框颜色
			scale_fill_manual(values = colors_use)+ # 自定义填充颜色，可以和边框颜色不一样，我觉得配套的比较好看，就用了一样的颜色
			theme_test()+ # 把背景颜色去掉并添加边框，可以修改为其他背景，比如theme_classic
			labs(y="ITH Score", x="")+ # 修改横纵坐标名称
			theme(strip.background = element_blank(), #把头上的灰框去掉
						axis.text = element_text(color="black", size=10),
						axis.line = element_blank())
		if (parameter1 != "seurat_clusters"){
			p=p+theme(axis.text.x =element_text(angle = 45, hjust = 1))
		}
		ggsave(plot=p,filename=file.path(oudir,paste0("ITHscore_groupby_",parameter2,"_spiltby_",parameter1,".pdf")),width=0.8*len1+1.3*len2-4)
		ggsave(plot=p,filename=file.path(oudir,paste0("ITHscore_groupby_",parameter2,"_spiltby_",parameter1,".png")),width=0.8*len1+1.3*len2-4)
		}
	}
}
setwd(setwd(oudir))
cmd = "cp /PERSONALBIO/work/singlecell/s04/Test/donghongjie/ITHscore/ITH_Score结果说明.docx ./"
system(cmd)
sink("Software version.txt")
sessionInfo()
sink()

