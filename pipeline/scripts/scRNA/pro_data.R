suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("magrittr"))

rm(list=ls())


# 分析需求:
# 1.数据下载链接：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893，其中需要对比的是IPF (n = 12) controls (n = 10)这两组；
# 2.看S1PR1这个在正常人和IPF病人的定位和表达量，说明，S1PR1在病人和正常人都在肺血管内皮高表达。需要两组的t SNE图，说明一下S1PR1在肺血管内皮上；
# 3.对比一下IPF和control两组S1PR1在内皮上的表达量。

# 关联合同号:
# 评估
# 可行性:
# 可执行
# 不可行原因:
# 分析方案:
# 数据下载链接：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893
# 其中需要对比的是IPF (n = 12) controls (n = 10)这两组

# 1. 绘制S1PR1基因分别在正常人和IPF病人的featureplot图。【说明：S1PR1在病人和正常人都在肺血管内皮高表达。需要两组的t SNE图，说明一下S1PR1在肺血管内皮上】

# 2. 绘制S1PR1基因在IPF和control两组内皮细胞的小提琴对比图


# seurat_ob <- readRDSMC("/public/scRNA/works/liuhongyan/Project/scRNA/HT2021-23293-h/GSE135893_ILD_annotated_fullsize.rds")

# ## PF 肺纤维 ##
# ## 肺纤维化 (PF) 是一种慢性肺病
# ## non-fibrotic control 非纤维化作为 control
# ## 其中需要对比的是IPF (n = 12) controls (n = 10)这两组

# #### 获取关注的 22 个样本的表达量信息 ######
# target <- subset(seurat_ob, subset=Diagnosis %in%  c("Control", "IPF") )

# saveRDSMC(target,"./target_seurat.rds")


# #### 进行标准化，归一化处理 ####
# all <- seurat_ob 
# seurat_ob <- target 

# opt=list()
# opt$normmeth="LogNormalize"
# opt$ncores=10
# if ( is.null(opt$ncores) ) {
#     nCores = 10
# }else{
#     nCores = opt$ncores
# }
# # #####################################################################
# #====== Setup the Seurat Object (Required)=========================
# # ####################################################################
# options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
# plan("multicore", workers = min(detectCores(), nCores)) 

# seurat_ob <- NormalizeData(object = seurat_ob,
#                 normalization.method = opt$normmeth,scale.factor = 10000)


# opt$vars2regress="nCount_RNA,percent.mt"
# vars2regress = unlist(strsplit(opt$vars2regress, ",", perl =T))


# vars2regress = unique(vars2regress)
# seurat_ob = FindVariableFeatures(object= seurat_ob, loess.span = 0.3,
#                     clip.max = "auto", mean.function = "FastExpMean",
#                     dispersion.function = "FastLogVMR", num.bin = 20,
#                     nfeature = 4000, binning.method = "equal_width" )

# # regress out all the specified the varibales
# seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob),
#                     vars.to.regress = vars2regress, verbose = T )

# output_dir="./"

# saveRDSMC(seurat_ob, file.path(dirname(output_dir), "scaledata_seurat.rds"), threads = nCores)






#### 带有 p 值的小提琴图 ####
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("magrittr"))

rm(list=ls())
setwd("/public/scRNA/works/liuhongyan/Project/scRNA/HT2021-23293-h/Further_analysis_2022-01-14/")
# seurat_ob <- readRDSMC("/public/scRNA/works/liuhongyan/Project/scRNA/HT2021-23293-h/target_seurat.rds")

# seurat_ob@meta.data <- seurat_ob@meta.data %>% dplyr::rename(sampleid=Sample_Name)
# saveRDSMC(seurat_ob,"sampleid_target_seurat.rds")

output_dir=paste0("violin_plot/")
dir.create(output_dir,recursive=T)

# seurat_ob <- subset(seurat_ob , subset=population %in% "Endothelial")

# saveRDSMC(seurat_ob,"Endothelial_target_seurat.rds")

seurat_ob <- readRDSMC("/public/scRNA/works/liuhongyan/Project/scRNA/HT2021-23293-h/Endothelial_target_seurat.rds")
# > dim(seurat_ob)
# [1] 26377  7236
### 筛选出 S1PR1 均表达的细胞进行比较差异 
seurat_ob <- subset(seurat_ob, S1PR1>0 )
# > dim(seurat_ob)
# [1] 26377  5166

suppressPackageStartupMessages( library("dittoSeq"))
seurat_ob@meta.data <- seurat_ob@meta.data %>% dplyr::rename(sampleid=Sample_Name)
saveRDSMC(seurat_ob,"S1PR1_Endothelial_target_seurat.rds")

gene="S1PR1"
data <- as.data.frame(seurat_ob@assays$SCT@data)[gene,]
data <- as.data.frame(t(data))

DATA <- cbind(seurat_ob@meta.data, data) %>% dplyr::select(gene,"sampleid","Diagnosis","population","celltype")


data <- DATA %>% tibble::rownames_to_column("Barcode")

write.table(data,file.path(output_dir,paste0(gene,"_expression_data.xls")),sep="\t",row.names=F,quote=F)


vlnplot = dittoPlot(seurat_ob, var= "cnv_level", group.by=group_by, plots =c("vlnplot", "boxplot"), legend.show = F ,split.by = facetby,split.ncol=2) +
                                scale_fill_manual(values = SelectColors(object=NULL,palette = group_colors[[group_by]], value =group_by ,n = nlevel))+
                                theme( plot.title = element_text(hjust = 10),axis.title.x = element_text(size =20),axis.title.y = element_text(size = 20),axis.text.x=element_text(size = 20),axis.text.y=element_text(size = 10))


library(ggpubr) #stat_compare_means() + 
library(ggsignif)

## 绿和浅绿
# col=c("#7fc97f","#beaed4")
# col=c("#17AC52","#24AE9B")
# col=c("#7fc97f","#009E73")
# col=c("#7fc97f","#1b9e77")
col=c("#17AC52","#7fc97f")

## https://zhuanlan.zhihu.com/p/385170087
# 一般我们转换的标准如下


# [0-0.001]-->***

# [0.001-0.01]-->**

# [0.01-0.05]-->*

# [0.05-0.1]-->.

# [0.1-1]-->NS(not significant) 不显著

# 当然这个也可以自己去定义，把所有>0.05的都定义成不显著。

groupby="Diagnosis"
comparisons_list=list(c("Control","IPF"))
p <- ggplot(data,aes_string(x=groupby,y=gene,fill=groupby))+
            # stat_boxplot(geom = "errorbar",width = 0.5,lwd=0.5) + 
            geom_violin()+
            geom_boxplot(width=0.05,position=position_dodge(0.9))+
            labs(x="",y="Expression", title=gene) +
            #stat_compare_means(method = "t.test") + 
            #stat_compare_means(label.y =ceiling(max(data[[gene]]))) + 
            
            #stat_compare_means(label.y =round(max(data[[gene]]),2)+0.05 ) + 
            theme_classic()+ 
            # 1.07326767790565e-11
            # annotate("text",x=1,y=3.1,label=expression("p = 1.07e-11"),
            #     size=5,color='black') + 
            ## , map_signif_level = T
            geom_signif(comparisons =comparisons_list, test = wilcox.test, annotations = c("****")) + 
        #    geom_signif( comparisons =comparisons_list, test = wilcox.test ) + 
            # geom_signif(annotations = c(1.07e-11), map_signif_level = T,y_position = c(3.2), xmin = c(1), xmax = c(1.5), vjust = -1) + 

            theme(panel.grid = element_blank(),legend.position = "none",
                plot.title = element_text(hjust = 0.5))+
            scale_fill_manual(values=col)

ggsave(file.path(output_dir,paste0(gene,"_violin_plot.pdf")),plot=p,height=6,width=6)
ggsave(file.path(output_dir,paste0(gene,"_violin_plot.png")),plot=p,dpi = 1000,height=6,width=6)



##### 总的 rds 中 S1PR1 #####
all <-  readRDSMC("/public/scRNA/works/liuhongyan/Project/scRNA/HT2021-23293-h/target_seurat.rds")
all_seurat_ob <- subset(all, S1PR1>0 )


all_seurat_ob@meta.data <- all_seurat_ob@meta.data %>% dplyr::rename(sampleid=Sample_Name)

saveRDSMC(all_seurat_ob,"S1PR1_all_target_seurat.rds")

gene="S1PR1"
data_all <- as.data.frame(all_seurat_ob@assays$SCT@data)[gene,]
data_all  <- as.data.frame(t(data_all ))

DATA_all <- cbind(all_seurat_ob@meta.data, data_all ) %>% dplyr::select(gene,"sampleid","Diagnosis","population","celltype")


data_all  <- DATA_all %>% tibble::rownames_to_column("Barcode")

write.table(data,file.path(output_dir,paste0(gene,"_all_expression_data.xls")),sep="\t",row.names=F,quote=F)




groupby="Diagnosis"
comparisons_list=list(c("Control","IPF"))
p=list()
col=c("#386cb0","#aed4ff")
## p=2.89815213223565e-08
p[[1]] <- ggplot(data_all,aes_string(x=groupby,y=gene,fill=groupby))+
            # stat_boxplot(geom = "errorbar",width = 0.5,lwd=0.5) + 
            geom_violin()+
            geom_boxplot(width=0.05,position=position_dodge(0.9))+
            labs(x="",y="S1PR1 gene expression in lung", title="") +
            #stat_compare_means(method = "t.test") + 
            #stat_compare_means(label.y =ceiling(max(data[[gene]]))) + 
            
            #stat_compare_means(label.y =round(max(data[[gene]]),2)+0.05 ) + 
            theme_classic()+ 
            # 1.07326767790565e-11
            # annotate("text",x=1,y=3.1,label=expression("p = 1.07e-11"),
            #     size=5,color='black') + 
            ## , map_signif_level = T
            geom_signif(comparisons =comparisons_list, test = wilcox.test, annotations = c("****"), size=0.8,textsize = 6 ) + 
        #    geom_signif( comparisons =comparisons_list, test = wilcox.test ) + 
            # geom_signif(annotations = c(1.07e-11), map_signif_level = T,y_position = c(3.2), xmin = c(1), xmax = c(1.5), vjust = -1) + 

            theme(panel.grid = element_blank(),legend.position = "none",
                plot.title = element_text(hjust = 0.5))+
              scale_fill_manual(values=col) + theme(axis.text.x=element_text(size=12),
                       axis.title.x=element_text(size=14),
                       axis.text.y=element_text(size=12),
                       axis.title.y=element_text(size=14),
                       axis.ticks =element_line(size=0.8,colour="black"),
                       axis.ticks.length =unit(0.3, "lines"),
                       #axis.ticks.margin =unit(1.2, "cm"),
                       axis.line = element_line(size=0.8,colour="black"))


col=c("#17AC52","#7fc97f")
## p = 1.07326767790565e-11
p[[2]] <- ggplot(data,aes_string(x=groupby,y=gene,fill=groupby))+
            # stat_boxplot(geom = "errorbar",width = 0.5,lwd=0.5) + 
            geom_violin()+
            geom_boxplot(width=0.05,position=position_dodge(0.9))+
            labs(x="",y="S1PR1 gene expression in Endothelial", title="") +
            #stat_compare_means(method = "t.test") + 
            #stat_compare_means(label.y =ceiling(max(data[[gene]]))) + 
            
            #stat_compare_means(label.y =round(max(data[[gene]]),2)+0.05 ) + 
            theme_classic()+ 
            # 1.07326767790565e-11
            # annotate("text",x=1,y=3.1,label=expression("p = 1.07e-11"),
            #     size=5,color='black') + 
            ## , map_signif_level = T
            geom_signif(comparisons =comparisons_list, test = wilcox.test, annotations = c("****"), size=0.8,textsize = 6) + 
        #    geom_signif( comparisons =comparisons_list, test = wilcox.test ) + 
            # geom_signif(annotations = c(1.07e-11), map_signif_level = T,y_position = c(3.2), xmin = c(1), xmax = c(1.5), vjust = -1) + 

            theme(panel.grid = element_blank(),legend.position = "none",
                plot.title = element_text(hjust = 0.5))+
            scale_fill_manual(values=col) + theme(axis.text.x=element_text(size=12),
                       axis.title.x=element_text(size=14),
                       axis.text.y=element_text(size=12),
                       axis.title.y=element_text(size=14),
                       axis.ticks =element_line(size=0.8,colour="black"),
                       axis.ticks.length =unit(0.3, "lines"),
                       #axis.ticks.margin =unit(1.2, "cm"),
                       axis.line = element_line(size=0.8,colour="black"))


suppressPackageStartupMessages(library("gridExtra"))
pdf(file.path(output_dir,paste0(gene,"_violin_plot.pdf")),
    width = 8, height = 8)
grid.arrange(grobs = p, nrow=1)
dev.off()
png(file.path(output_dir,paste0(gene,"_violin_plot.png")),
    width = 8, height = 8, res = 96, units = "in")
grid.arrange(grobs = p, nrow=1)
dev.off()

# ggsave(file.path(output_dir,paste0(gene,"_violin_plot.pdf")),plot=p,height=6,width=6)
# ggsave(file.path(output_dir,paste0(gene,"_violin_plot.png")),plot=p,dpi = 1000,height=6,width=6)