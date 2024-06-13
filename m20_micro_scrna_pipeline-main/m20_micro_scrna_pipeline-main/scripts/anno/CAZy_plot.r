#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("optparse"))

option_list = list(
    make_option( c( "--input","-i"), type = "character",
        help = "the result of CAZy blast file."
    ),
    make_option( c("--class","-c"),type="character",
        help = "the path of class_definition.txt file"
    ),
    make_option( c("--famly","-f"),type="character",
        help = "the path of family_information.txt file"
    ),
    make_option( c("--output","-o"),type="character",
        help = "the path of family_information.txt file"
    )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if ( ! is.null(opt$output) ){
    output_dir <-  opt$output
    if (! file.exists(output_dir)) { dir.create(output_dir, recursive = T) }
    output_dir<-normalizePath(output_dir)
}

if ( ! is.null(opt$input)){
    inputfile <- file.path(opt$input)
}

if( ! is.null(opt$famly)){
    familyfile <- file.path(opt$famly)
}

if( ! is.null(opt$class)){
    classfile <- file.path(opt$class)
}
data<-read.table(file=inputfile,sep="\t",header=T)
data_new<-data %>% separate(Description,sep="\\|",into=c("Subject","Familyid"))
refdata<-read.table(file=familyfile,head=T,sep="\t")
data_all<-left_join(data_new,refdata,by="Familyid")
write.table(data_all,file.path(output_dir,"CAZy_annodata.xls"),sep="\t",quote=F,row.names=F)
plot_data<-read.table(file=classfile,sep="\t",header=T)
class_num<-as.data.frame(table(data_all$Class))
colnames(class_num)<-c("Class","Genes_Count")
plot<-merge(plot_data,class_num,by="Class")
write.table(plot,file.path(output_dir,"CAZy_anno_plotdata.xls"),sep="\t",row.names=F,quote=F)
p=ggplot(data=plot,aes(x=Definition,y=Genes_Count,fill=Definition))+geom_bar(stat="identity", width = 0.5, position = position_dodge(0.7))+
guides(fill="none")+
theme(panel.background=element_rect(fill='transparent',color ="gray"),  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,color = "black",size=9))+
theme(panel.grid =element_blank())+
geom_text(mapping = aes(label = Genes_Count),size=3.5,vjust=-0.2)+
xlab("Class definition") + ylab("Genes number")
ggsave(filename=glue::glue("{output_dir}/CAZy_class_plot.pdf"),width=8,height=7,plot=p)
ggsave(filename=glue::glue("{output_dir}/CAZy_class_plot.png"),width=8,height=7,plot=p)