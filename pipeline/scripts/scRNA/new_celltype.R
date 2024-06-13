suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("Seurat")) 
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("future")) ## 可以学习下
suppressPackageStartupMessages(library("dplyr"))

seurat_ob <- readRDS("/public/scRNA/works/liuhongyan/Project/scRNA/done_report/HT2020-11162-3-human/Further_analysis_2022-03-09/1.picture/2.UMAP/singlecell_object.clustering_resolution0.1.rds")


# Cluster1：TC
# Cluster2、3、8：TFC
# Cluster4：FC
# Cluster5：EC
# Cluster6：MC
# Cluster7：BC
# Cluster9：CC
seurat_ob@meta.data$new_celltype = ""
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(1)),"new_celltype"] = "TC"
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(2,3,8)),"new_celltype"] = "TFC"
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(4)),"new_celltype"] = "FC"
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(5)),"new_celltype"] = "EC"
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(6)),"new_celltype"] = "MC"
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(7)),"new_celltype"] = "BC"
seurat_ob@meta.data[which(seurat_ob@meta.data$clusters %in% c(9)),"new_celltype"] = "CC"

## c("#2CA02C","#D52122","#9467BD","#8C564B","#E377C2","#BCBD22","#17BECF")
seurat_ob@meta.data$new_celltype=factor(seurat_ob@meta.data$new_celltype, levels=c("TC", "TFC","FC","EC","MC","BC","CC"))
saveRDS(seurat_ob,"new_celltype_seurat.rds")

## visualization 
toid="new_celltype"
nlevel = length(unique(seurat_ob@meta.data[,toid]))
opt=list();opt$reduct="umap"
output_dir="new_celltype"
dir.create(output_dir,recursive=T)
if (dim(seurat_ob)[2] < 500){
  pointsize=1.5
}else{
  pointsize=0.5
}

usecol=c("#2CA02C","#D52122","#9467BD","#8C564B","#E377C2","#BCBD22","#17BECF")
ggtsne = DimPlot(object = seurat_ob, reduction = opt$reduct,pt.size = pointsize,group.by=toid,label=T, label.size=9) +
    ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5), legend.position="none", title=element_text(size=24), axis.text=element_text(size=22), axis.ticks.length =unit(0.4, "lines")  ) + ggplot2::scale_colour_manual( values = usecol)
ggplot2::ggsave(file.path(output_dir,paste0(toid,".pdf")), ggtsne)
ggplot2::ggsave(file.path(output_dir,paste0(toid,".png")), ggtsne)

simplified_meta = seurat_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "orig.ident") %>%
                               dplyr::select( Barcode, sampleid, clusters,group,!!toid)
write.table(simplified_meta, quote = F,sep =",",row.names = F,
             file.path(output_dir,paste0(toid,".metadata.csv",collapse = "")))


### 柱状图
opt$groupby="new_celltype"
seurat_ob = SetIdent( seurat_ob, value = opt$groupby )
DATA <- as.data.frame( seurat_ob@meta.data[,c("sampleid", opt$groupby)] ) %>%
    group_by( .dots= c("sampleid", opt$groupby)) %>%
    dplyr::summarize(cell_number = n()) %>%
    mutate(freq = (cell_number / sum(cell_number)) * 100)


write.table(as.data.frame(DATA),
        file.path(output_dir,file="clust_cond_freq_info.xls"),
        sep="\t",col.names=T, row.names =F)

seurat_ob = SetIdent(seurat_ob, value = "sampleid")
clust_sum_all2 = PlotAbundances(seurat_ob, prop.by = opt$groupby, group.by ="sampleid" , method = "barplot",
cols= usecol) + theme( title=element_text(size=24), axis.text=element_text(size=22),legend.title= element_blank(),  legend.text= element_text(size=22) , axis.ticks.length =unit(0.4, "lines") , plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), axis.text.x=element_text(angle=0,hjust=0.5) )

ggsave(file.path(output_dir,paste0("groupby-",opt$groupby,"_summary_plot.pdf",collapse="")),plot=clust_sum_all2)
ggsave(file.path(output_dir,paste0("groupby-",opt$groupby,"_summary_plot.png",collapse="")),dpi=1000, plot = clust_sum_all2)


## mix plot 
facetbyx="sampleid"
groupvis_all = DimPlot(object = seurat_ob, dims = c(1,2),reduction = opt$reduct,
                   pt.size = pointsize, group.by = facetbyx)+
                       theme( plot.title = element_text(hjust = 0.5), title=element_text(size=24), axis.text=element_text(size=22),legend.text= element_text(size=22), axis.ticks.length =unit(0.4, "lines") ) +
                       scale_colour_manual( values = c("#1F7784", "#FF7F0E", "#2CA02C", "#D52122") )
                   ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_contrast_plot.pdf",collapse="")),
                   limitsize = FALSE, plot = groupvis_all, height = 7, width = 7)
                   ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_contrast_plot.png",collapse="")),
                   limitsize = FALSE, plot = groupvis_all, width = 7, height = 7 )

