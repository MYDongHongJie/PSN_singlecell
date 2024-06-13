# module purge && module  load OESingleCell/2.0.0
library("ggplot2")
library("dplyr")
library(RColorBrewer)
library("optparse") 


data=read.delim("/public/scRNA_works/works/guokaiqi/project/scRNA/Mayan_Xuhanfu_Can_MLC/DOE202210395_YC_DZ/houxu-20230427/terms.xls",sep="\t")
# name=gsub(".xls","",opt$input)
# data <- data[which(data["ListHits"]>2), ]

#data$pro_pval = ifelse( data$up_or_down=="Up",as.numeric(paste0("-",data$pval)),data$pval) 
## data$pro_Enrichment_score = ifelse( data$up_or_down=="Up",as.numeric(paste0("-",data$Enrichment_score)),data$Enrichment_score) 
#data$is.just = ifelse( data$up_or_down=="Up",1,0) 
# tata<-function(data){
# # data$just = ifelse( data$up_or_down=="Up",1,0)
# data <- data %>% mutate(pro_pval = -log(pval,10) )
# data = data %>% group_by(celltype) %>% arrange(data,-ListHits)
# data$term <- factor(data$term, levels=rev(unique(data$term)) )
# return(data)
# }
# data=tata(data)
data <- data %>% mutate(pro_pval = -log(pval,10) )
data$term <- factor(data$term, levels=rev(unique(data$term)))
data$celltype <- factor(data$celltype, levels=as.character(unique(data$celltype)))

# data <- data %>%
  # group_by(celltype) %>%
  # mutate(colorp = scales::col_numeric(pro_pval, palette = c("Reds")))
library("gridExtra")
ditto = c("#E69F00","#56B4E9","#009E73", "#F0E442", "#0072B2" ,"#D55E00" ,"#CC79A7", "#666666" ,"#AD7700", "#1C91D4", 
          "#007756", "#D5C711", "#005685", "#A04700" ,"#B14380" ,"#4D4D4D","#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71" )
gs=list()
for (i in 1:6){
plot_data=subset(data, celltype == levels(data$celltype)[i])
gs[[i]] = ggplot(plot_data, aes(x=term,y=ListHits)) +coord_flip() + 
				geom_col(aes(fill=pro_pval),  width=0.6) + #scale_fill_manual(values = ditto[1:8])+ #fill=pro_pval,
                  scale_fill_gradient(low = "gray", high = ditto[i]) + labs(x="", y= "",fill="-log10 Pval")+
                theme_bw()  + theme_classic()  + scale_y_continuous(expand = c(0,0))+ #ylim(0,60)+
				theme(panel.background = element_blank(), panel.grid =element_blank(),panel.border = element_blank()) + 
				theme( axis.ticks.x=element_line(), axis.line.x = element_line(), axis.line.y= element_line(),axis.ticks.y=element_blank(),
					  axis.text.y=element_text( size=11,color="black"),axis.text.x=element_text( size=10,color="black")) #+  #, axis.title.x=element_text(size=8)
if (i==6){ gs[[i]] = gs[[i]] + labs(x="", y= "ListHits")#+
            # theme(axis.ticks.x=element_line(), axis.line.x = element_line(),
                  # axis.text.x=element_text( size=10,color="black")) #+  #, axis.title.x=element_text(size=8)
 }

}
  pdf("term_plot.pdf", width = 12, height = 10)
  # grid.arrange(grobs = gs, ncol=1,align_axes = TRUE)
  cowplot::plot_grid (gs[[1]],gs[[2]],gs[[3]],gs[[4]],gs[[5]],gs[[6]], ncol=1, align = "v")
  dev.off()
  
  png("term_plot.png", width = 1280, height = 1000)
  # grid.arrange(grobs = gs, ncol=1,align_axes = TRUE)
  cowplot::plot_grid (gs[[1]],gs[[2]],gs[[3]],gs[[4]],gs[[5]],gs[[6]], ncol=1, align = "v")
  dev.off()
  
  
  
		    # ggsave(paste0("term_plot.pdf"),plot = pp,height=5,width=15)
		    # ggsave(paste0("term_plot.png"),plot = pp,height=5,width=15,dpi=1000) 

# # 对数据进行预处理，以便正确地设置颜色和渐变
# data <- data %>%
  # group_by(celltype) %>%
  # mutate(color = scales::col_numeric(p_value, palette = c("gray","#E69F00")))

# # 绘制一个基于数量值的条形图，纵坐标为 term 名字，颜色根据分组设置并且颜色的渐变表示 p 值
# ggplot(my_data, aes(x = value, y = term, fill = group)) +
  # geom_col() +
  # scale_fill_manual(values = unique(my_data$color)) +
  # scale_fill_gradient(low = "white", high = "red") +
  # theme_classic()
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())