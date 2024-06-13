library(tibble)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

data_ob = OESingleCell::ReadX(input = "seurat.h5seurat", informat = 'h5seurat',verbose = F) 
data<-as.data.frame(data_ob@reductions$umap@cell.embeddings) %>% rownames_to_column(var = "barcodes")
anno<-data_ob@meta.data[,c('new_celltype','group')] %>% rownames_to_column(var = "barcodes")
data=merge(data,anno,by="barcodes")

#绘图
CustomCol2 <- function(n){
  my_palette=c(
    "#86CBEC","#9E9BCA","#673C98","#64A41A","#6498CA","#42A998","#004000","#D5FED1","#E58024","#D22A81","#E5AA04","#C42C2B","#C35E1A","#A2731B","#89CAE8","#716CAC","#673D93","#67A224",  
    "#6698CB","#479F91","#329F2A","#281A6F","#216DA1","#65C4A2","#9D9F0",
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
} 

#设置等高线
bwidth = .05*c(sum(abs(range(data$UMAP_1))),sum(abs(range(data$UMAP_2))))

if(reduct="tsne"){
    gg_tsne = ggplot(data,aes(tSNE_1,tSNE_2,colour=new_celltype)) +
        stat_density2d(aes(x=tSNE_1,y=tSNE_2,colour=new_celltype),alpha=0.5,h=bwidth, contour_var = "ndensity") +
        geom_point(size=0.7,stroke=0.1) +
        theme_test()+
        scale_colour_manual(values=CustomCol2(1:length(unique(data$new_celltype)))) +
        guides(fill=FALSE) + NoAxes() +
        geom_segment(aes(x = min(data$tSNE_1)-5, y = min(data$tSNE_2), xend = min(data$tSNE_1) + 5, yend = min(data$tSNE_2) ), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+   #设置坐标轴x
        geom_segment(aes(x = min(data$tSNE_1)-5, y = min(data$tSNE_2), xend = min(data$tSNE_1) - 5, yend = min(data$tSNE_2) + 10), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +   #设置坐标轴y
        annotate("text", x = min(data$tSNE_1), y = min(data$tSNE_2) -2, label = "tSNE_1", color="black",size = 3, fontface="bold" ) +   
        annotate("text", x = min(data$tSNE_1) - 7, y = min(data$tSNE_2)+5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) +
        theme(legend.title = element_blank(), legend.key=element_rect(fill='white'), legend.text = element_text(size=20), panel.border = element_blank()) +  
        guides(color = guide_legend(override.aes = list(size=5)))  
}else{
    gg_tsne = ggplot(data,aes(UMAP_1,UMAP_2,colour=new_celltype)) +
        stat_density2d(aes(x=UMAP_1,y=UMAP_2,colour=new_celltype),alpha=0.5,h=bwidth, contour_var = "ndensity") +
        geom_point(size=0.7,stroke=0.1) +
        theme_test()+
        scale_colour_manual(values=CustomCol2(1:length(unique(data$new_celltype)))) +
        guides(fill=FALSE) + NoAxes() +
        geom_segment(aes(x = min(data$UMAP_1)-1, y = min(data$UMAP_2), xend = min(data$UMAP_1) , yend = min(data$UMAP_2) ), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+   #设置坐标轴x
        geom_segment(aes(x = min(data$UMAP_1)-1, y = min(data$UMAP_2), xend = min(data$UMAP_1) - 1, yend = min(data$UMAP_2) + 1 ), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +   #设置坐标轴y
        annotate("text", x = min(data$UMAP_1)-0.5, y = min(data$UMAP_2)-0.25, label = "UMAP_1", color="black",size = 3, fontface="bold" ) +
        annotate("text", x = min(data$UMAP_1)-1.25, y = min(data$UMAP_2) + 0.5, label = "UMAP_2",color="black",size = 3, fontface="bold" ,angle=90) +
        theme(legend.title = element_blank(), legend.key=element_rect(fill='white'), legend.text = element_text(size=20), panel.border = element_blank()) +
        guides(color = guide_legend(override.aes = list(size=5)))
}
#ggsave("new_celltype.pdf",p,width=10,height=10)
#ggsave("new_celltype.png",p,width=10,height=10)

#如果不想要图例，设置标签位置
cell_type_med <- data %>%  group_by(new_celltype) %>%  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)  )
p = gg_tsne +
    geom_text_repel(aes(label=new_celltype), data = cell_type_med, point.padding=unit(0.5, "lines"), colour="black",size=3) +  theme(legend.position = "none")
    
ggsave("new_celltype.pdf",p,width=7,height=7)
ggsave("new_celltype.png",p,width=7,height=7)



