library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)


#数据地址:https://github.com/weiiioyo/singlecell-heatmap/blob/main/data/
data_ob = OESingleCell::ReadX(input = 'new_celltype_20221130.h5seurat', informat = 'h5seurat', verbose = F)
#gene=read.table("gene.txt",header=T,sep="\t",quote='',stringsAsFactors = F)
topn_markers = read.delim("genelist.txt",header=T)
formated_extra_gene = as.data.frame(tidyr::gather(topn_markers,key = "cluster",value = "GENE"))
formated_extra_gene = formated_extra_gene[-12,]
#去除CD1C
formated_extra_gene = formated_extra_gene[-5,]
formated_extra_gene
topn_markers = data.frame()

match = OESingleCell::CaseMatch(search = as.vector(formated_extra_gene$GENE),match = rownames(data_ob))
formated_extra_gene = formated_extra_gene %>%
                          dplyr::filter(GENE %in% names(match)) %>%
                          dplyr::mutate(gene = match,folder_suffix = cluster) %>%
                          dplyr::select(cluster, gene,folder_suffix)
topn_markers = rbind(topn_markers, formated_extra_gene)
head(topn_markers)
#从Seurat对象中提取细胞注释以及基因表达量

vln.dat=FetchData(data_ob,c(topn_markers$gene,"new_celltype"))
vln.dat.melt = data.frame()
vln.dat.melt=vln.dat %>% 
            reshape2::melt(,topn_markers$gene) %>%
            rename("Gene"="variable") %>%
            group_by(new_celltype,Gene) %>%
            mutate(fillcolor=mean(value))


# 小提琴图的填充颜色
pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,'YlOrRd'))(100)

# 堆积小提琴图
p1 = ggplot(vln.dat.melt,aes(x=new_celltype,y=value,fill=fillcolor))+
  # 把小提琴图的外缘轮廓去除
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  # 分面
  facet_grid(Gene~.) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 45),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank( ),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "right"
        
        ) 



p1 = ggplot(vln.dat.melt,aes(x=new_celltype,y=value,fill=fillcolor))+
  # 把小提琴图的外缘轮廓去除
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  # 分面
  facet_grid(Gene~.)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "left"
        
        )+theme_dark()

p1 = ggplot(vln.dat.melt,aes(x=new_celltype,y=value,fill=fillcolor))+
  # 把小提琴图的外缘轮廓去除
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  # 分面
  facet_grid(Gene~.)+
  theme(panel.background = element_rect(fill = 'grey90'),
        panel.grid = element_blank(),
        #panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "left"
        
        )

a=topn_markers%>%count(cluster)
a$yend=cumsum(a$n)
a$ystart=a$yend-a$n
a$label_position=a$ystart+(a$yend-a$ystart)/2

p2=ggplot(a,aes(x=1,y=ystart,color=cluster))+
  geom_segment(aes(xend=1,yend=yend),size=a$n)+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(position = "right",
                     breaks = a$label_position,
                     labels = a$cluster,expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  facet_grid(cluster~.,scales = 'free')+
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
pp = p1+p2+plot_spacer()+plot_layout(ncol = 2, widths  = c(3, .03),heights = c(3,.03))




ggsave(plot=pp,"vln_merge.png",dpi = 300)

ggsave(plot=pp,"vln_merge.pdf",dpi = 300)

