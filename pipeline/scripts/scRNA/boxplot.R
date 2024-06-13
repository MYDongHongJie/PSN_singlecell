# module purge && module  load OESingleCell/2.0.0
library("ggplot2")
library("dplyr")
library(RColorBrewer)

data <- read.delim("/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-18966-m/Further_analysis_2022-03-07/term_data.xls",sep="\t")

#data$pro_pval = ifelse( data$up_or_down=="Up",as.numeric(paste0("-",data$pval)),data$pval) 
## data$pro_Enrichment_score = ifelse( data$up_or_down=="Up",as.numeric(paste0("-",data$Enrichment_score)),data$Enrichment_score) 
#data$is.just = ifelse( data$up_or_down=="Up",1,0) 
data$just = ifelse( data$up_or_down=="Up",1,0)
data <- data %>% mutate(pro_pval = -log(pval,10) )
data$pro2_pval = ifelse( data$up_or_down=="Up",data$pro_pval,as.numeric(paste0("-",data$pro_pval)))

data$term <- factor(data$term, levels=rev(unique(data$term)) )
pp = ggplot(data, aes(term,pro2_pval)) +
				geom_col(aes(fill=Enrichment_score),  width=0.6) +
                #scale_fill_distiller(palette = "Spectral") + ## https://cloud.tencent.com/developer/article/1092853, fill 填充渐变颜色
                scale_fill_distiller(palette = "RdBu") +
                #scale_fill_distiller(palette = "Set1" ) + 
				coord_flip() +
				labs(x="", y=expression('-log'[10]*' Pvalue')) +
				theme_minimal() +
                ## hjust = data$just, 
				geom_text( aes(x= term, y=0, label = term),hjust=data$just, size = 5)+
				theme(axis.text.y=element_blank()) +
				theme(panel.grid =element_blank(), axis.line.x = element_line( size = 0.5), axis.ticks.x=element_line(size=0.5), legend.position = "top",  legend.title=element_text(size=14), legend.text=element_text(size=12),axis.text.x=element_text(size=12), axis.title.x=element_text(size=14)) + 
                scale_y_continuous(limits=c(-15,15), breaks=seq(-15,15,5),label=c("15","10","5","0","5","10","15"))
                # + ylim(-15,15) 

		    ggsave("term_plot.pdf",plot = pp,height=5)
		    ggsave("term_plot.png",plot = pp,height=5,dpi=1000) 

