
# module load /home/xfyang/modulefiles/OESingleCell/beta


library(tidyr)
library(dplyr)
library(tidyverse)
library(readr)
library(ggplot2)
library(ggalluvial)
library(latex2exp)
suppressWarnings(suppressMessages(library("OESingleCell")))
# suppressPackageStartupMessages(library("R.matlab"))
suppressWarnings(suppressMessages(library("magrittr")))
# suppressPackageStartupMessages(library("future"))
# suppressPackageStartupMessages(library("Seurat"))
suppressWarnings(suppressMessages(library("tidyverse")))
# suppressPackageStartupMessages(library("reshape2"))
# suppressPackageStartupMessages(library("scales"))
# suppressPackageStartupMessages(library('circlize'))
# suppressPackageStartupMessages(library('seriation'))
# suppressPackageStartupMessages(library("viridis"))
# suppressPackageStartupMessages(library("patchwork"))
# suppressPackageStartupMessages(library("networkD3"))
# suppressPackageStartupMessages(library("rbokeh"))
# suppressPackageStartupMessages(library("webshot"))
suppressWarnings(suppressMessages(library("argparse")))
suppressWarnings(suppressMessages(library("futile.logger")))
suppressWarnings(suppressMessages(library("funr")))
data_ob = OESingleCell::ReadX(input = '/public/scRNA_works/works/chenhaoruo/project/DOE20223173-DOE20221324/DOE20223173-b1-mouse/20230524/result/sub_merge/seurat.h5seurat', informat = 'h5seurat',verbose = F)

total=data_ob@meta.data[,c("rawbc","clusters","new_clusters")]
# total= read.table("/public/scRNA_works/works/chenhaoruo/project/DOE20223173-DOE20221324/DOE20223173-b1-mouse/20230524/merge_meta.csv",sep=",",header=T)
# total$Correlation=1
plotdata=total
library(tidyverse)
library(ggplot2)
library(cols4all)
library(reshape2)
# df <- plotdata %>%
# make_long(clusters, new_clusters)

# df$node <- factor(df$node,levels = c(plotdata$new_clusters %>% unique()%>% rev(),
# plotdata$clusters %>% unique() %>% rev()))
    cluster_ids <- total[,"new_clusters"]
    sample_ids <- total[, "clusters"]
    counts <- table(cluster_ids, sample_ids)
  df <- melt(t(round(t(counts)/colSums(counts) * 100, 2)), 
      varnames = c("clusters", "new_clusters"), value.name = "freq")

library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")

library(ggraph)
library(ggforce)
library(tidyr)

dotdata=data.frame(
  x = c(rep("Pre",9),rep("Mer",9),rep("Nai",9)),
  y = c(seq(0, 50.4, length=9),seq(0, 50.4, length=9),seq(0, 50.4, length=9))
)
df$group=""
df$group[grep("^Pre", df$clusters)] <- "Pre"
df$group[grep("^Nai", df$clusters)] <- "Nai"
df$group2 <- sub(".*_(\\d+)", "\\1", df$clusters)
df$group3 <- sub(".*_(\\d+)", "\\1", df$new_clusters)
df$group2 = as.numeric(df$group2)
df$group3 = as.numeric(df$group3)

df_nai=df[df$group=="Nai" & df$freq!=0,]
df_pre=df[df$group=="Pre" & df$freq!=0,]

line_nai=data.frame(
  x = rep("Mer",length(df_nai$group2)),
  y = seq(0, 50.4, length=9)[df_nai$group3],
  xend = rep("Nai",length(df_nai$group2)),
  yend = seq(0, 50.4, length=9)[df_nai$group2],
  freq = df_nai$freq
)

line_pre=data.frame(
  x = rep("Mer",length(df_pre$group2)),
  y = seq(0, 50.4, length=9)[df_pre$group3],
  xend = rep("Pre",length(df_pre$group2)),
  yend = seq(0, 50.4, length=9)[df_pre$group2],
  freq = df_pre$freq
)

dotdata$label=c(c(1:9),c(1:9),c(1:9))
dotdata$x=factor(dotdata$x,levels = c("Pre","Mer","Nai"))
dotdata$temp=paste0(dotdata$x,"_",dotdata$label)

cor=read.table("/public/scRNA_works/works/chenhaoruo/project/DOE20223173-DOE20221324/DOE20223173-b1-mouse/20230524/result/cor2/correlation_total.xls",sep="\t",header=T)[,c(1,2,3)]
cor$group2 <- sub(".*_(\\d+)", "\\1", cor$Name1)
cor$group3 <- sub(".*_(\\d+)", "\\1", cor$Name2)

for(i in dotdata$temp){
  if(i %in% cor$Name1){
    dotdata$Correlation[dotdata$temp==i]=max(cor[cor$Name1==i,"Correlation"])
  }else if (i %in% cor$Name2){
    dotdata$Correlation[dotdata$temp==i]=max(cor[cor$Name2==i,"Correlation"])
  }
}
dotdata$Correlation=as.numeric(dotdata$Correlation) %>% round(.,2)
dotdata[dotdata$x=="Mer","Correlation"]=""

dotplot = ggplot(dotdata, aes(x = x, y = y ,color = x )) +
  geom_point(aes(size = 25) ) +
  geom_diagonal(data=line_nai,aes(x = x,y = y, xend = xend, yend = yend, size = (freq/30) , alpha= (freq*1.5)),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_diagonal(data=line_pre,aes(x = x,y = y, xend = xend, yend = yend, size = (freq/30) , alpha= (freq*1.5)),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_point(aes(size = 25) ) +
  # scale_fill_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_color_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_size_continuous(range = c(1, 10)) +
  theme_bw()+theme(   axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid=element_blank(), # 去网格线
                      panel.border=element_blank(), # 去边框
                      legend.position = "right",
                      legend.text = element_text(size=12),
                      legend.title= element_text(size= 12),
                    axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1)) +
  guides(size=FALSE) +
  geom_text(aes(label = label) , size = 4, color = "black", fontface = "bold", show.legend = FALSE) +
  geom_text(aes(label = Correlation) , size = 4, color = "black", fontface = "bold", show.legend = FALSE , nudge_y = 3 , nudge_x = 0.05) 

ggsave(dotplot,file="freq_sankey.png",width = 30,height = 10,units = "cm",dpi = 300)
ggsave(dotplot,file="freq_sankey.pdf",width = 30,height = 10,units = "cm",dpi = 300)



for(i in rownames(df)){
  df[i,"Correlation"]=max(cor[cor$Name1==df[i,"clusters"] & cor$group3==df[i,"new_clusters"],"Correlation"])
}

df_nai=df[df$group=="Nai" & df$freq!=0,]
df_pre=df[df$group=="Pre" & df$freq!=0,]

line_nai=data.frame(
  x = rep("Mer",length(df_nai$group2)),
  y = seq(0, 50.4, length=9)[df_nai$group3],
  xend = rep("Nai",length(df_nai$group2)),
  yend = seq(0, 50.4, length=9)[df_nai$group2],
  freq = df_nai$freq,
  cor=df_nai$Correlation,
  log_cor= -log(df_nai$Correlation)
)

line_pre=data.frame(
  x = rep("Mer",length(df_pre$group2)),
  y = seq(0, 50.4, length=9)[df_pre$group3],
  xend = rep("Pre",length(df_pre$group2)),
  yend = seq(0, 50.4, length=9)[df_pre$group2],
  freq = df_pre$freq,
  cor=df_pre$Correlation,
  log_cor= -log(df_pre$Correlation)
)


dotplot = ggplot(dotdata, aes(x = x, y = y ,color = x )) +
  geom_point(aes(size = 25) ) +
  geom_diagonal(data=line_nai,aes(x = x,y = y, xend = xend, yend = yend, size = (1/log_cor)*0.05 , alpha= (1/log_cor)*2),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_diagonal(data=line_pre,aes(x = x,y = y, xend = xend, yend = yend, size = (1/log_cor)*0.05 , alpha= (1/log_cor)*2),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_point(aes(size = 25) ) +
  # scale_fill_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_color_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_size_continuous(range = c(1, 10)) +
  theme_bw()+theme(   axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid=element_blank(), # 去网格线
                      panel.border=element_blank(), # 去边框
                      legend.position = "right",
                      legend.text = element_text(size=12),
                      legend.title= element_text(size= 12),
                    axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1)) +
  guides(size=FALSE) +
  geom_text(aes(label = label) , size = 4, color = "black", fontface = "bold", show.legend = FALSE) +
  geom_text(aes(label = Correlation) , size = 4, color = "black", fontface = "bold", show.legend = FALSE , nudge_y = 3 , nudge_x = 0.05) 

ggsave(dotplot,file="cor_sankey.png",width = 30,height = 10,units = "cm",dpi = 300)
ggsave(dotplot,file="cor_sankey.pdf",width = 30,height = 10,units = "cm",dpi = 300)



top2_line_nai=line_nai %>% group_by(yend) %>% top_n(2,cor)
top2_line_pre=line_pre %>% group_by(yend) %>% top_n(2,cor)

dotplot = ggplot(dotdata, aes(x = x, y = y ,color = x )) +
  geom_point(aes(size = 25) ) +
  geom_diagonal(data=top2_line_nai,aes(x = x,y = y, xend = xend, yend = yend, size = (1/log_cor)*0.05 , alpha= (1/log_cor)*2),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_diagonal(data=top2_line_pre,aes(x = x,y = y, xend = xend, yend = yend, size = (1/log_cor)*0.05 , alpha= (1/log_cor)*2),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_point(aes(size = 25) ) +
  # scale_fill_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_color_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_size_continuous(range = c(1, 10)) +
  theme_bw()+theme(   axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid=element_blank(), # 去网格线
                      panel.border=element_blank(), # 去边框
                      legend.position = "right",
                      legend.text = element_text(size=12),
                      legend.title= element_text(size= 12),
                    axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1)) +
  guides(size=FALSE) +
  geom_text(aes(label = label) , size = 4, color = "black", fontface = "bold", show.legend = FALSE) +
  geom_text(aes(label = Correlation) , size = 4, color = "black", fontface = "bold", show.legend = FALSE , nudge_y = 3 , nudge_x = 0.05) 

ggsave(dotplot,file="top2_cor_sankey.png",width = 30,height = 10,units = "cm",dpi = 300)
ggsave(dotplot,file="top2_cor_sankey.pdf",width = 30,height = 10,units = "cm",dpi = 300)



top2_line_nai=line_nai %>% group_by(yend) %>% top_n(2,cor)
top2_line_pre=line_pre %>% group_by(yend) %>% top_n(2,cor)

top2_line_nai_1=line_nai %>% group_by(y) %>% top_n(2,cor)
top2_line_pre_1=line_pre %>% group_by(y) %>% top_n(2,cor)

top2_1= rbind(top2_line_nai,top2_line_nai_1)
top2_2= rbind(top2_line_pre,top2_line_pre_1)

dotplot = ggplot(dotdata, aes(x = x, y = y ,color = x )) +
  geom_point(aes(size = 25) ) +
  geom_diagonal(data=top2_1,aes(x = x,y = y, xend = xend, yend = yend, size = (1/log_cor)*0.05 , alpha= (1/log_cor)*2),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_diagonal(data=top2_2,aes(x = x,y = y, xend = xend, yend = yend, size = (1/log_cor)*0.05 , alpha= (1/log_cor)*2),color = "black" ,strength = 0.75, show.legend = FALSE) +
  geom_point(aes(size = 25) ) +
  # scale_fill_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_color_manual(values = c("Pre" = "red", "Mer" = "blue", "Nai" = "green")) +
  # scale_size_continuous(range = c(1, 10)) +
  theme_bw()+theme(   axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid=element_blank(), # 去网格线
                      panel.border=element_blank(), # 去边框
                      legend.position = "right",
                      legend.text = element_text(size=12),
                      legend.title= element_text(size= 12),
                    axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1)) +
  guides(size=FALSE) +
  geom_text(aes(label = label) , size = 4, color = "black", fontface = "bold", show.legend = FALSE) +
  geom_text(aes(label = Correlation) , size = 4, color = "black", fontface = "bold", show.legend = FALSE , nudge_y = 3 , nudge_x = 0.05) 

ggsave(dotplot,file="top2_cor_sankey.png",width = 30,height = 10,units = "cm",dpi = 300)
ggsave(dotplot,file="top2_cor_sankey.pdf",width = 30,height = 10,units = "cm",dpi = 300)


