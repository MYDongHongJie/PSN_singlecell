# Title     : DONE
# Objective : DONE
# Created by: kun.sun
# Created on: 2023/10/9
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("ggplot2"))
option_list <- list(
    make_option( c("--input_path", "-i" ), type = "character",
                 help = "qualimap输出结果文件"),
    make_option( c("--palette", "-p" ), type = "character",
                 help = "选择指定调色板"),
    make_option( c("--outdir", "-o" ), type = "character",
                 help = "输出结果文件夹")
    )
opt <- parse_args(OptionParser(option_list=option_list))
output_dir <- opt$outdir
if (! file.exists(output_dir)) { dir.create(output_dir, recursive = T) }
#汇总合并
sample_name<-list.dirs(path=opt$input_path,full.names=F,recursive=F)
plot_data <- data.frame()
for (i in sample_name) {
  data<-read.table(paste0(opt$input_path,"/",i,"/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"),sep="\t",header=F)
  colnames(data) <- c("Transcript_position","Counts")
  Add_data <- data.frame(data,pos=(data$Counts/sum(data$Counts))*100,Sample=rep(i,length(data$Counts)))
  plot_data <- rbind(plot_data,Add_data)
}
#绘图
p <- ggplot(plot_data,mapping=aes(x=Transcript_position,y=pos,colour=Sample))+
  geom_line(size=1.1)+
  theme_bw()+
  scale_y_continuous(breaks=seq(0,round(max(plot_data$pos),1)+0.5,0.2),limits=c(0,round(max(plot_data$pos),1)+0.5))+
  scale_x_continuous(breaks=seq(0,100,10))+
  xlab("Transcript Position(%)")+
  ylab("Percentage total cumulative mapped-read depth")+
  ggtitle("Coverage Profile Along Genes (total)")+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12),plot.title = element_text(hjust = 0.5))+
  scale_colour_manual( values = OESingleCell::SelectColors(unique(plot_data$Sample),palette = opt$palette))

ggsave(plot=p,filename=glue::glue("{output_dir}/qualimap_plot.png"),width=10,height=5)
ggsave(plot=p,filename=glue::glue("{output_dir}/qualimap_plot.pdf"),width=10,height=5)
