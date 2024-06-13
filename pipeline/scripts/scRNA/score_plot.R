library(Seurat)
library(OESingleCell)

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
get_colors <- function(object, groupby){
                if(paste0(groupby,"_col") %in% colnames(colData(object))){
                    groupby_col = paste0(groupby,"_col")
                    nlevel_list = levels(factor(colData(object)[,groupby]))
                    tmp_df <- unique(colData(object)[c(groupby, groupby_col)])
                    groupby_pal <- as.vector(tmp_df[,groupby_col])
                    names(groupby_pal) <-  as.vector(tmp_df[,groupby])
                    groupby_pal = as.list(groupby_pal)
                    user_color_pal = unlist(groupby_pal[nlevel_list])
                }else if(groupby =="clusters"){
                    nlevel <- sort(unique(colData(object)[,groupby]))
                    user_color_pal = CustomCol2(nlevel)
                }else {
                    nlevel = length(unique(colData(object)[,groupby]))
                    user_color_pal = CustomCol2(1:nlevel)
                }
                return(user_color_pal)
              }


data_ob<-readRDS("../houxu-20220926/cca_integrat_umap.clustering_resolution0.4_group_regions.rds") 

nCount_Spatial<-FeaturePlot(data_ob,features="nCount_Spatial",pt.size=0.5,min.cutoff = 3000)
ggsave("nCount_Spatial.png",width=7.6) 
ggsave("nCount_Spatial.pdf",width=7.6)

nCount_Spatial<-FeaturePlot(data_ob,features="nFeature_Spatial",pt.size=0.5,min.cutoff = 2000)
ggsave("nFeature_Spatial.png",width=7.6)
ggsave("nFeature_Spatial.pdf",width=7.6)
##############################################################################################
library(dplyr)
DATA <- OESingleCell::FetchData(data_ob, vars= c("sampleid", "regions") ) %>%
            dplyr::group_by( .dots= c("sampleid", "regions")) %>%
            dplyr::summarize(count = n()) %>%
            dplyr::mutate(freq = (count / sum(count)) * 100)
write.table(as.data.frame(DATA),
              file.path(file="clust_cond_freq_info.xls"),
              sep="\t",col.names=T, row.names =F)

data_ob@meta.data$regions=as.vector(data_ob@meta.data$regions)
data_ob@meta.data$regions[which(is.na(data_ob@meta.data$regions))]="unknown" 
data_ob = SetIdent(data_ob, value = "regions")

clust_sum_all = OESingleCell::PlotAbundances(data_ob, prop.by = "regions" , group.by = "sampleid", method = "barplot",
                                 cols= OESingleCell::SelectColors(levels(data_ob))
)
ggsave(file.path(paste0("groupby-sampleid_regions_summary_plot.pdf",collapse="")),plot=clust_sum_all, height = 10)
ggsave(file.path(paste0("groupby-sampleid_regions_summary_plot.png",collapse="")),dpi=1000, plot = clust_sum_all, height = 10)
##############################################################################################
#score <- read.delim("/public/scRNA_works/works/guokaiqi/project/stRNA/HT2021-14842-6_human/houxu-20221114/addmodulescore_cca_SCT/gene_module_scores.xls",sep="\t",row.names=1)
score <- read.delim("/public/scRNA_works/works/guokaiqi/project/stRNA/HT2021-14842-6_human/houxu-20230510/addmodulescore_cca_SCT/gene_module_scores.xls",sep="\t",row.names=1)
data_ob <- AddMetaData(data_ob, metadata = score)
#geneset = c("Hepatocytes", "Hepaticstellatecell", "NK_cells", "B_cells", "T_cells","Myeloid_cells", "Endothelial_cells" ,"Fibroblasts")
geneset = c("Hepatocytes_geneset", "Malignancy", "T.cells", "B.cells", "kuffer.cell","Hepatic.stellate.cell", "Bile.duct" ,"Myeloid.cells","Neutrophils")
metadata = data_ob@meta.data[,c("sampleid","clusters","regions","group_regions", geneset)]
metadata$clusters = factor(metadata$clusters)
vln=list()
for (gg in geneset){
vln[[gg]] = ggplot(metadata, aes_string(x="clusters", y=gg, fill="clusters",color="clusters" ))+ geom_violin()+
            labs(x="Clusters",y=paste0(gg,"_score")) + facet_wrap(~sampleid, ncol = 4) +
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            scale_fill_manual(values=get_colors(data_ob,"clusters")) + scale_color_manual(values=get_colors(data_ob,"clusters"))
          ggsave(file.path("vlnplot",paste0(gg,"_score_vlnplot.pdf")),plot=vln[[gg]],width=12,height=4)
          ggsave(file.path("vlnplot",paste0(gg,"_score_vlnplot.png")),plot=vln[[gg]],width=12,height=4)
}

library(ggridges)
ridge=list()
for (gg in geneset){
# ridge[[gg]] = ggplot(metadata, aes_string(y="clusters", x=gg, fill="clusters",color="clusters" ))+ geom_density_ridges()+
#             labs(y="Clusters",x=paste0(gg,"_score")) + facet_wrap(~sampleid, ncol = 4,scales = "free_x") +
#             theme(panel.grid.major =element_blank(),
#                   panel.grid.minor = element_blank(),
#                   panel.background = element_blank(),
#                   axis.line = element_line(color = "black")) +
#             scale_fill_manual(values=get_colors(data_ob,"clusters")) + scale_color_manual(values=get_colors(data_ob,"clusters"))
#           ggsave(file.path("vlnplot",paste0(gg,"_score_ridgeplot.pdf")),plot=ridge[[gg]],width=20,height=4)
#           ggsave(file.path("vlnplot",paste0(gg,"_score_ridgeplot.png")),plot=ridge[[gg]],width=20,height=4)

# ridge[[paste0(gg,"_line")]] = ggplot(metadata)+ geom_density_ridges(fill="white", aes_string(y=1, x=gg ))+
#             labs(y=gg,x=paste0(gg,"_score")) + facet_wrap(~sampleid, ncol = 4) +
#             theme(panel.grid.major =element_blank(),
#                   panel.grid.minor = element_blank(),
#                   panel.background = element_blank(),
#                   axis.line = element_line(color = "black")) + scale_fill_manual(values="white")
#             #scale_fill_manual(values=get_colors(data_ob,"clusters")) + scale_color_manual(values=get_colors(data_ob,"clusters"))
#           ggsave(file.path("vlnplot",paste0(paste0(gg,"_line"),"_score_ridgeplot.pdf")),plot=ridge[[paste0(gg,"_line")]],width=20,height=2)
#           ggsave(file.path("vlnplot",paste0(paste0(gg,"_line"),"_score_ridgeplot.png")),plot=ridge[[paste0(gg,"_line")]],width=20,height=2)
ridge[[paste0(gg,"_line")]] = ggplot(metadata, aes_string(y="clusters", x=gg, fill="clusters",color="clusters" ))+ geom_density_ridges(scale = 1.2,rel_min_height = 0)+
            labs(y="Clusters",x=paste0("score")) + facet_wrap(~sampleid, ncol = 4,scales = "free_x") +
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            scale_fill_manual(values=get_colors(data_ob,"clusters")) + scale_color_manual(values=get_colors(data_ob,"clusters"))
          ggsave(file.path("vlnplot",paste0(paste0(gg,"_line"),"_score_ridgeplot.pdf")),plot=ridge[[paste0(gg,"_line")]],width=20,height=8)
          ggsave(file.path("vlnplot",paste0(paste0(gg,"_line"),"_score_ridgeplot.png")),plot=ridge[[paste0(gg,"_line")]],width=20,height=8)

}

metadata
library(reshape2)
metadata_melt=melt(metadata[,c("sampleid","clusters", geneset)], id = c("sampleid", "clusters"), variable.name = "geneset", value.name ="score")
metadata_melt$clusters=factor(metadata_melt$clusters)

ridgeplot = ggplot(metadata_melt, aes_string(y="clusters", x="score", fill="clusters",color="clusters" ))+ geom_density_ridges(scale = 1.2,rel_min_height = 0)+
            labs(y="Clusters",x=paste0("score")) + facet_wrap(~sampleid, ncol = 4,scales = "free_x") +
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            scale_fill_manual(values=get_colors(data_ob,"clusters")) + scale_color_manual(values=get_colors(data_ob,"clusters"))
          ggsave(file.path("vlnplot",paste0("score_ridgeplot.pdf")),plot=ridgeplot,width=20,height=8)
          ggsave(file.path("vlnplot",paste0("score_ridgeplot.png")),plot=ridgeplot,width=20,height=8)

ridgeplot2 = ggplot(metadata_melt, aes_string(y="geneset", x="score", fill="geneset",color="geneset" ))+ geom_density_ridges(fill="white",scale = 0.9,rel_min_height = 0)+
            labs(y="geneset",x=paste0("score")) + facet_wrap(~sampleid, ncol = 4,scales = "free_x") +
            theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            scale_color_manual(values=CustomCol2(30:38)) + scale_fill_manual(values="white")
          ggsave(file.path("vlnplot",paste0("score_ridgeplot_geneset.pdf")),plot=ridgeplot2,width=20,height=6)
          ggsave(file.path("vlnplot",paste0("score_ridgeplot_geneset.png")),plot=ridgeplot2,width=20,height=6)
####################################
library(pheatmap)
library(RColorBrewer)
score <- read.delim("/public/scRNA_works/works/guokaiqi/project/stRNA/HT2021-14842-6_human/houxu-20221114/addmodulescore_cca_SCT/gene_module_scores.xls",sep="\t",row.names=1)
data_ob <- AddMetaData(data_ob, metadata = score)
geneset = c("Hepatocytes", "Hepaticstellatecell", "NK_cells", "B_cells", "T_cells","Myeloid_cells", "Endothelial_cells" ,"Fibroblasts")
#geneset = c("Hepatocytes_geneset", "Malignancy", "T.cells", "B.cells", "kuffer.cell","Hepatic.stellate.cell", "Bile.duct" ,"Myeloid.cells","Neutrophils")
metadata = data_ob@meta.data[,c("sampleid","group","clusters","regions","group_regions", geneset)]
# metadata=data_ob@meta.data
metadata$regions_merge=gsub("[0-9]","",metadata$regions)
metadata$clusters_regions = paste0("C",metadata$clusters,"_",metadata$regions_merge)
metadata=metadata[which(!is.na(metadata$regions_merge)),]   
table(metadata$regions_merge)
# BD   CV   FB   HL   PL 
#104  144  412  517 2009
metadata$id = rownames(metadata)
cols=c("Greens","Reds","Purples","Blues")
names(cols)=unique(metadata$sampleid)

samples=list()
for (s in unique(metadata$sampleid)){
metadata_sub = subset(metadata,sampleid == s)
collapsed_group = metadata_sub %>% group_by(.dots = "clusters_regions") %>% do(data.frame(cellid = paste0(.$id,collapse = ",")))
collapsed_group = collapsed_group[!is.na(collapsed_group$clusters_regions),]
collapsed_count = vector()
for ( cells in collapsed_group$cellid ){
    samplex = unlist(strsplit(cells, ",", perl =T))
    collapsed_count= cbind(collapsed_count,colMeans( metadata_sub[samplex,geneset,drop=F] ))
}
collapsed_count = as.matrix( collapsed_count )
collapsed_group = as.data.frame(collapsed_group)
colnames(collapsed_count) = as.matrix(collapsed_group[,1])
samples[[s]]=collapsed_count
# palette <- colorRampPalette(c("#406AA8", "white", "#D91216"))(n=299)
palette <- colorRampPalette(brewer.pal(9, cols[s]))(100)
metadata_sub$clusters = factor(metadata_sub$clusters)
annotation_col = unique(metadata_sub[,c("clusters_regions","clusters","regions_merge")])
rownames(annotation_col) = annotation_col$clusters_regions
annotation_col = annotation_col[,-1]
    annotation_colors_row1 = CustomCol2(unique(sort(annotation_col[,1])))
    names(annotation_colors_row1) = unique(sort(annotation_col[,1]))
    annotation_colors_row2 = CustomCol2(15:(14+as.numeric(length(unique(annotation_col[,2])))))
    names(annotation_colors_row2) = unique(sort(annotation_col[,2]))

    annotation_colors_row = list(annotation_colors_row1,annotation_colors_row2)
    names(annotation_colors_row) = colnames(annotation_col)[1:2]
    annotation_colors = annotation_colors_row

p = pheatmap(collapsed_count,
        color = palette,
        cex=1,
        border=F,
        angle_col=0,
        lwd=1.1,
        cellheight=40, cellwidth=40,
        annotation_col = annotation_col,
        cluster_rows= F,annotation_colors = annotation_colors,
        cluster_cols=F )
ggsave(file.path("heatmap",paste0(s,"_heatmap.pdf")), plot= p,width = 18,height = 10)
ggsave(file.path("heatmap",paste0(s,"_heatmap.png")), plot= p, dpi=1000,width = 18,height = 10)

collapsed_count = tibble::rownames_to_column(as.data.frame(collapsed_count, var = "geneset"))
write.table(collapsed_count, file.path("heatmap",paste0(s,"_heatmap_data.xls")),quote = F, row.names = F, sep = "\t")
}
######################################
x1 <- a*b/(a^2+1)
y1 <- (a^2*b)/(a^2+1)
A <- c(10,4)

(x1+y1)/2
######################################
# 导入 ggplot2 和 ggpmisc 包
library(ggplot2)
library(ggpmisc)
library(ggalt)
# 创建一个随机散点图并绘制斜率为1的参考直线和投影线
set.seed(123)
df <- data.frame(x = runif(20), y = runif(20))
ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  geom_xspline(aes(x = ifelse(x > y, y, x), y = ifelse(x > y, y, x)),
               color = "red", size = 1, alpha = 0.5)
metadata2 = subset(metadata,regions_merge %in% c("HL","PL"))
metadata2 = subset(metadata,regions_merge %in% c("BD","CV"))

for (g in geneset){
   meta = metadata2 %>% select(sampleid,regions_merge,!!g) %>% group_by(sampleid,regions_merge) %>% summarise(median=median(.data[[g]]))
   df = reshape2::dcast(meta, sampleid~regions_merge, value.var = "median")

   df$dist <- abs(df$BD - df$CV)/sqrt(2)
   # 计算直线的截距和斜率
   a <- 0
   b <- 1
   # 计算垂线的长度和坐标
   df$dist <- abs((a+b)*df$BD - (a-b)*df$CV)/sqrt(2)
   df$x0 <- df$BD + (df$CV - df$BD)/2
   df$y0 <- df$CV - (df$CV - df$BD)/2
   p = ggplot(df, aes(x =BD, y = CV,color = sampleid)) +
       geom_point(size = 3) +
       geom_abline(intercept = 0, slope = 1, linetype = "dotted") + coord_fixed(ratio = 1, xlim = c(-0.5,0.5), ylim =  c(-0.5,0.5)) + 
       geom_segment(aes(x = BD, y = CV, xend = x0, yend = y0, color = sampleid))+  #linetype ="dashed"
       theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) 
      #  geom_line(aes(x = ifelse(BD > CV, CV, BD), y = ifelse(BD > CV, CV, BD)),color = "red", size = 1, alpha = 0.5)
      #  geom_xspline(aes(x = ifelse(x > y, y, x), y = ifelse(x > y, y, x)),
                    # color = "red", size = 1, alpha = 0.5)
  ggsave(paste0(g,"_point.png"),width=8)
  ggsave(paste0(g,"_point.pdf"),width=8)
# # 计算每个点到斜率为 1 的直线的垂线长度和坐标
# a <- -1/sqrt(2) # 直线的斜率
# b <- 1/sqrt(2)  # 直线的截距

# df$dist <- abs(a*df$x + b*df$y)/sqrt(a^2 + b^2)
# df$x0 <- (b*(df$y - a*df$x))/(a^2 + b^2) + df$x
# df$y0 <- (-a*(df$y - a*df$x))/(a^2 + b^2) + df$y
}
median()
segments
#######################
seurat_ob<-readRDS("v3_rds/data_ob_v3.rds")

seurat_ob@meta.data$regions_merge = gsub("[0-9]","",seurat_ob@meta.data$regions)

seurat_ob@meta.data$group_regions_merge = seurat_ob@meta.data$group_regions
seurat_ob@meta.data$sampleid_regions = paste0(seurat_ob@meta.data$sampleid,"_",seurat_ob@meta.data$regions)
seurat_ob@meta.data$sampleid_regions=gsub("_NA","",seurat_ob@meta.data$sampleid_regions)
seurat_ob@meta.data$group_regions = paste0(seurat_ob@meta.data$group,"_",seurat_ob@meta.data$regions)
seurat_ob@meta.data$group_regions=gsub("_NA","",seurat_ob@meta.data$group_regions)

# > head(seurat_ob@meta.data,2)
#                               orig.ident nCount_SCT nFeature_SCT nCount_Spatial
# HC_1_AAACAAGTATCTCCCA AAACAAGTATCTCCCA-1       7993         3026           8410
# HC_1_AAACACCAATAACTGC AAACACCAATAACTGC-1       8366         3656          12654
#                       nFeature_Spatial slice sampleid      slide slide_area
# HC_1_AAACAAGTATCTCCCA             3028  HC_1     HC_1 V11Y10-354         A1
# HC_1_AAACACCAATAACTGC             3823  HC_1     HC_1 V11Y10-354         A1
#                                             fastq                     image
# HC_1_AAACAAGTATCTCCCA raw_data/GZ_25/GZ_25_fastqs raw_data/images/GZ_25.tif
# HC_1_AAACACCAATAACTGC raw_data/GZ_25/GZ_25_fastqs raw_data/images/GZ_25.tif
#                       species group batchid integrated_snn_res.0.8
# HC_1_AAACAAGTATCTCCCA   human    HC       1                      8
# HC_1_AAACACCAATAACTGC   human    HC       1                      2
#                       seurat_clusters integrated_snn_res.0.4
# HC_1_AAACAAGTATCTCCCA               4                      4
# HC_1_AAACACCAATAACTGC               0                      0
#                       integrated.umap.res0.4 clusters regions group_regions
# HC_1_AAACAAGTATCTCCCA                      5        5    <NA>            HC
# HC_1_AAACACCAATAACTGC                      1        1    <NA>            HC
#                       regions_col regions_merge group_regions_merge
# HC_1_AAACAAGTATCTCCCA        <NA>          <NA>                  HC
# HC_1_AAACACCAATAACTGC        <NA>          <NA>                  HC
#                       sampleid_regions
# HC_1_AAACAAGTATCTCCCA             HC_1
# HC_1_AAACACCAATAACTGC             HC_1



