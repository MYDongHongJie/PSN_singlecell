library(Seurat)
#准备数据
data_ob = OESingleCell::ReadX(input = "seurat.h5seurat", informat = 'h5seurat',verbose = F)
gene <- read.delim("genelist.txt")
gene=as.vector(CaseMatch(search=gene$gene,match=rownames(data_ob)))

data=as.data.frame(t(as.data.frame(data_ob@assays$RNA@data[gene,])))

#添加要绘图得分组信息，这里添加得clusters
for (i in unique(data_ob@meta.data$clusters)){
  group = rownames(data_ob@meta.data[which(data_ob@meta.data$clusters==i),])
  data[which(rownames(data) %in% group),'group']=i
}

#绘图      
library(ggsignif) #添加统计检验
library(gghalves)
#library(ggdist)
library(ggplot2) #绘图
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

# 如果加P值得话加上下面这个

#         geom_signif(comparisons = list(c("Are","control")), #指定比较对象
#               test = "wilcox.test", #指定检验方法
#               size = 0.4, #指定标记中线条的尺寸
#               textsize = 2.6, #指定标记中文字部分的大小
#               vjust = -0.05, #指定标记中文字部分与横线之间的距离指定标记中文字部分与横线之间的距离
#               step_increase = 0.1, #指定每根线条距离高低
#               map_signif_level =T )  
                                    
#outlier.size = 0,outlier.alpha =0
for (i in colnames(data)[-dim(data)[2]]){
     p<-ggplot(data,aes_string(x="group",y=i,fill="group",color="group"))+
          scale_fill_manual(values = CustomCol2(1:length(unique(data$group))))+
          scale_colour_manual(values = CustomCol2(1:length(unique(data$group)))) + 
          geom_half_point(position=position_nudge(x=-0.35,y=0),size =1, shape =19,range_scale = 0.5,alpha=0.5) + 
          geom_boxplot(position = position_nudge(x = 0),width=0.1,outlier.shape = NA,alpha=0.5)+
          geom_half_violin(position=position_nudge(x=0.1,y=0),
                        side='R',adjust=1.2,trim=F,color=NA,alpha=0.5) +           
          theme_classic()+
          ggtitle(i)+ labs(y="Gene expression") + #设置主标题


                 theme(
                    axis.title.y = element_text(size = 13), #调整Y轴标题字体大小

                    plot.title = element_text(hjust = 0.5),
                axis.title.x=element_blank())  +
           ggplot2::theme(legend.position="none")                                                     
     ggsave(paste0(i,".pdf"),p,width=6,height=5)
     ggsave(paste0(i,".png"),p,width=6,height=5)                        
}
