#!/usr/bin/env Rscript
#this script is a substitution of the former function of the Run_DEG.R, which 
#is aimed to visulize the results of the differential expression in the last
#step.

##################################################################################
##diff_stat_barplotS#########
diff_stat_barplot<-function(diff_stat,type, species){
  library(ggplot2)
  library(stringr)
  data<- read.table(diff_stat, header=T, sep="\t", quote="", stringsAsFactors=FALSE,check.names=F);
  for(i in 1:length(data$Control)){
    if( ! is.na(str_match(data$Case[i], "\\(")[1] == "(") ) {
      data$Control[i]<-substr(data$Control[i],str_locate_all(data$Control[i],"\\(")[[1]]+1,str_locate_all(data$Control[i],"\\)")[[1]]-1)
    }
    if( ! is.na(str_match(data$Case[i], "\\(")[1] == "(") ) {
      data$Case[i]<-substr(data$Case[i],str_locate_all(data$Case[i],"\\(")[[1]]+1,str_locate_all(data$Case[i],"\\)")[[1]]-1)
    }
  }
  stat<-matrix(nrow=2*length(data$Control),ncol=3)
  for(i in 1:length(data$Control)){
    stat[(2*i-1),]<-cbind(paste0(data$Case[i],"-vs-",data$Control[i]),"Up",data$Up_diff[i])
    stat[(2*i),]<-cbind(paste0(data$Case[i],"-vs-",data$Control[i]),"Down",data$Down_diff[i])
  }
  data2<-as.data.frame(stat)
  if ( !is.null(species) ){
    pic_title = paste0("Statistic of Differently Expressed ",type)
  }else{
    pic_title =  paste0("Statistic of Differently Expressed ",type," of ",species )
  }
  colnames(data2)<-c("Group","Type","Gene_number")
  write.table(data2,paste0(type,"_diff_stat_barplot.xls"),quote=F,row.names=F,col.names=T,sep="\t",na="")
  data2 <- read.table(paste0(type,"_diff_stat_barplot.xls"), header=T, sep="\t", quote="", check.names=F);
  uniq<-data2$Group[!duplicated(data2$Group)]
  data2$Group <- factor(data2$Group, levels=uniq)
  data2$Type <- factor(data2$Type, levels=c( "Up","Down"));
  p=ggplot(data=data2,aes(x=Group, y=Gene_number, fill=Type, width=0.7, space=0))+
    geom_bar(stat="identity", position="dodge")+
    geom_text(aes(label=data2$Gene_number),hjust=0.5,vjust=-0.5,size=2.7,position = position_dodge(0.7))+
    xlab("")+ylab(paste0("Differently Expressed ",type," number"))+
    labs(title= pic_title) +
    theme(plot.title = element_text(hjust = 0.5,size=13)) +
    theme(panel.border=element_rect(fill=NA,colour="black"),
          panel.background = element_rect(fill="transparent",colour=NA),
          plot.background = element_rect(fill="transparent",colour=NA))+
    theme(axis.text.y=element_text(size=10,color="black")) +
    theme(axis.text.x=element_text(angle = 45, hjust=1, vjust = 1, size=10,color="black")) +
    theme(legend.text=element_text(size=10))+
    theme(panel.grid =element_blank())
  ggsave(paste0(type,"_diff_stat_barplot.pdf"),height=6, width=7, plot=p)
  ggsave(paste0(type,"_diff_stat_barplot.png"),type="cairo-png", height=6, width=7, plot=p)
  file.remove(paste0(type,"_diff_stat_barplot.xls"))
}
##################################################################################
## heatmap_plot function define####
heatmap_plot<-function(condition,expression,group,FC,Pvalue,FDR, species ){
  #palette <-colorRampPalette(c("green", "black", "red"))(n=299)
  case_s_names<-unlist(strsplit(as.character(group$case),split=','))
  control_s_names<-unlist(strsplit(as.character(group$control),split=','))
  samples<-c(as.character(control_s_names),as.character(case_s_names))
  if(!is.null(Pvalue)){
    f <- c(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-diff-","pval-",Pvalue,"-FC-",FC,".",opt$type,".xls",sep=""))
  }
  if(!is.null(FDR)){
    f <- c(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-diff-","padj-",FDR,"-FC-",FC,".",opt$type,".xls",sep=""))
  }
  f_rn <- rownames(read.delim(f, header=T, sep="\t",row.names=1))
  diff.expression <- expression[match(f_rn,expression[,1]),samples]
  rownames(diff.expression)=f_rn
  ind <- apply(diff.expression, 1, mean) > 0
  diff.expression <- diff.expression[ind, ]
  
  if(!is.null(Pvalue)){
    picname=c(paste0(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-heatmap-", "pval-",Pvalue,"-FC-",FC, ".",opt$type))
    type=c("pvalue")
    thresh=Pvalue
  }
  if(!is.null(FDR)){
    picname=c(paste0( condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-heatmap-","padj-",FDR,"-FC-",FC, ".",opt$type))
    type<-c("padj")
    thresh=FDR
  }
  #check if the analysis is dualRNA-seq,if true,the title shold be specified for each speices
  if ( !is.null(species) ){
     pic_title = c(condition,"_",paste0(as.character(group$case_name),"-vs-",as.character(group$control_name),":",type,"<",thresh,"&& |log2FC|>",log2(FC))," for ",species )
  }else{
     pic_title = c(condition,"_",paste0(as.character(group$case_name),"-vs-",as.character(group$control_name),":",type,"<",thresh,"&& |log2FC|>",log2(FC)))
  }
  case_rep=rep(as.character(group$case_name), length(case_s_names))
  control_rep=rep(as.character(group$control_name), length(control_s_names))
  anno=c(control_rep,case_rep)
  annotation_col = data.frame(condition =factor(anno))
  rownames(annotation_col)=samples
  breaksList = seq(0, 5, by = 0.1)
  data<-log2(diff.expression+0.0001)
  if(length(f_rn)<30){
    xx<-pheatmap(data,annotation_col = annotation_col,
                 scale="row",
                 main=pic_title,
                 color =colorRampPalette(c("blue", "white", "red"))(256),
                 show_rownames=T)
  }else{
    xx<-pheatmap(data,annotation_col = annotation_col,
                 scale="row",
                 main=pic_title,
                 color =colorRampPalette(c("blue", "white", "red"))(256),
                 show_rownames=F)
  }
  save_pheatmap <- function(x, filename) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    #pdf
    pdf(paste0(filename, ".pdf"))
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
    #png
    png(paste0(filename, ".png"), height=2500, width=2500, res=300)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  save_pheatmap(xx, picname)
  if(file.exists("Rplots.pdf")){ file.remove("Rplots.pdf") }
  
  #convert -verbose -density 150  -trim *.pdf  -quality 100  -flatten -sharpen 0x1.0 *.tiff  1>&2>/dev/null
}
##################################################################################
## volcano_plot function define####
volcano_plot<-function(condition,group,FC,Pvalue,FDR, species ){
  fv <-function (d, p,  n1, n2, alpha,f1,name,title,cex = 0.4, pch = 19,
                 cols = c("#838B8B", "red", "green","#838B8B"),
                 ltys = c(1, 3), use_legend = TRUE, ...){   tt=f1
  f1 <- abs(d) <1
  col <- rep(cols[1], length(d))
  col[!f1 & p < alpha & d> 0] <- cols[2]
  col[!f1 & p < alpha & d< 0] <- cols[3]
  col[!f1 & p >=alpha ] <- cols[4]
  xlab = expression(paste(log[2], " fold change"))
  if(name=="pvalue"){ylab = expression(paste("-", log[10],"pvalue"))}
  if(name=="padj"){ylab = expression(paste("-", log[10],"padj"))}
  plot(d, -log10(p), cex = cex, pch = pch, xlab = xlab, ylab = ylab, col = col, main=title,...)
  abline(h=-log10(alpha), v=-f1, col = "gray60")
  abline(h=-log10(alpha), v=f1, col = "gray60")
  legend("topleft", c("Normal","Up" ,"Down"), pch = pch,
         col = cols, inset = 0.025,cex=0.8,,bg = "white")
  }
  data <- read.delim(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-all.",opt$type,".xls",sep=""),header = T, row.names = 1,sep="\t",quote="")
  #picname <-c(paste("Volcano-","-vs-",".pdf",sep=""))
  if(!is.null(Pvalue)){
    data<-data[which(data$pval!=""),]
    P<-data$pval
    picname <-c(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-volcano","-pval-",Pvalue,"-FC-",FC, ".",opt$type,sep=""))
    type=c("pvalue")
    thresh=Pvalue
  }
  if(!is.null(FDR)){
    data<-data[which(data$padj!=""),]
    P<-data$padj
    picname <-c(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),"-volcano","-padj-",FDR,"-FC-",FC, ".",opt$type,sep=""))
    type<-c("padj")
    thresh=FDR
  }
  log2FC<-data$log2FoldChange
  df <-data.frame(P,log2FC)
  tmp<-df[which(df$log2FC!="Inf"&df$log2FC!="-Inf"),]
  max=max(tmp$log2FC)
  min=min(tmp$log2FC)
  df$log2FC<-as.numeric(sub("-Inf",min,df$log2FC))
  df$log2FC<-as.numeric(sub("Inf",max,df$log2FC))
  ##pdf
  pdf(paste0(picname, ".pdf"), height=8, width=7)
  if ( !is.null( species ) ){
     pic_title = c(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),":",type,"<",thresh,"&& |log2FC|>",log2(FC)," for ",species))
  }else{
     pic_title = c(paste(condition,"_",as.character(group$case_name),"-vs-",as.character(group$control_name),":",type,"<",thresh,"&& |log2FC|>",log2(FC)))
  }
  fv(unlist(df$log2FC), unlist(df$P),  n1, n2, alpha=c(thresh),f1=c(log2(FC)),name=c(type), title= pic_title )
  dev.off()
  ##png
  png(paste0(picname, ".png"), height=2500, width=2200, res=300)
  fv(unlist(df$log2FC), unlist(df$P),  n1, n2, alpha=c(thresh),f1=c(log2(FC)),name=c(type), title=pic_title)
  dev.off()
  
}

##################################################################################
##                         main function module         ######
##################################################################################
#command line parameter definition#####
suppressPackageStartupMessages(library("optparse"))
option_list = list(
  make_option(c("-e", "--expression"), type="character", default=NULL,
              help="expression matrix file name", metavar="character"),
  make_option(c("-d", "--group"), type="character", default=NULL,
              help="[REQUIRED]comparison information file name", metavar="character"),
  make_option(c("-p", "--pvalue"), type="double",
              help="[REQUIRED]pvalue ratio threshold.Filtering can be performed using any one of (-p), (-f) at a time", metavar="double"),
  make_option(c("-f", "--fdr"), type="double",
              help="fdr ratio threshold.Filtering can be performed using any one of (-p), (-f) at a time", metavar="double"),
  make_option(c("-x", "--foldchange"),type="double", default=2.0,
              help="foldchange threshold [default %default]", metavar="double"),
  make_option( c("-c","--contrast"),type = "character",default = NULL,
               help = "[OPTIONAL]levels of a factor used to compare with for final differenetial results.The format is Factor:interesting_levle:reference_level."),
  make_option(c("-t", "--type"),type="character",default=NULL,
              help="[REQUIRED]transcript or gene type:mRNA,circRNA,lncRNA,miRNA,gene", metavar="character"),
  make_option(c("-m","--species"),type="character",default = "NULL",help="[Optional]if the analysis type is dualRNA,this would be host or guest"),
  make_option(c("-o", "--outputdir"), type="character", default="./",
              help="output directory for results", metavar="character"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##################################################################################
## get this R script file's own directory
args <- commandArgs(trailingOnly = F)
script.dir <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

#######import library###################################################################################################
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
########parameter check
if (is.null(opt$group)){
  print_help(opt_parser)
  stop("Sample groupping file must be supplied", call.=FALSE)
}else{
  group<-read.table(normalizePath(opt$group),header=T,sep="\t",check.names=FALSE,stringsAsFactors=FALSE,quote="")
}
if (! is.null(opt$pvalue)  &  ! is.null(opt$fdr)){
  print_help(opt_parser)
  stop("Filtering can be performed using any one of (-p), (-f) at a time.", call.=FALSE)
}
if (! is.null(opt$expression)){
  expression<-read.table(normalizePath(opt$expression),header = T,check.names=FALSE,quote="")
  colnames(expression)[1]<-paste0(opt$type,"_id")
  express_file_abpath<-normalizePath(opt$expression)
}else{
  cat("expression matrix file is not supplied\n")
}

if(file.exists(opt$outputdir)=="FALSE"){
  dir.create(opt$outputdir)
}
if ( !is.null(opt$species) ){
   species = opt$species
}

setwd(opt$outputdir)

###################
#if(file.exists(paste0(opt$type,"_diff_stat.xls"))==TRUE) {
#  suppressPackageStartupMessages(library(ggplot2))
#  file.remove(paste0(opt$type,"_diff_stat.xls"))
#}

#recontrust the groupping according to the contrast string
#get contrast elements,the contrast string format-factor:case_levele_in_factor:control_level_in_factor
contrasts = unlist( strsplit( opt$contrast,":",perl=T ) )
all_levels = unique(group[,contrasts[1]])
if ( contrasts[2] == "all" & contrasts[3] != "all" ){
    all_case_levels = all_levels[-which(all_levels==contrasts[3])] #without the control level
    all_comparisions = paste(contrasts[1],all_case_levels,contrasts[3],sep = ":")
}else if( contrasts[2] == "all" & contrasts[3] == "all" ){
  combine_of2 = combn(all_levels,2)
  all_comparisions = c( paste(contrasts[1],combine_of2[2,],combine_of2[1,],sep = ":"),paste(contrasts[1],combine_of2[1,],combine_of2[2,],sep = ":"))
}else if ( contrasts[2] != "all" & contrasts[3] == "all" ){
  all_control_levels = all_levels[-which(all_levels==contrasts[2])] #without the case level
    all_comparisions = paste(contrasts[1],contrasts[2],all_control_levels,sep = ":")
}else{
  all_comparisions = opt$contrast
}
print(all_comparisions)
#for each parsed the actual contrast or the comparision,reconstruct the groupping
mgroups = data.frame("case" = character(),"case_name" = character(),"control" = character(),"control_name" = character(),stringsAsFactors =F)
for ( comps in all_comparisions ){
    contrast_elements = unlist( strsplit( comps, ":", perl=T ) )
    control_sample = paste(group[which(group[,contrast_elements[1]]==contrast_elements[3]),"sample"], collapse = ",")
    case_sample = paste(group[which(group[,contrast_elements[1]]==contrast_elements[2]),"sample"],collapse = ",")
    mgroups[nrow(mgroups)+1,] = c( case_sample,contrast_elements[2],control_sample,contrast_elements[3])
}
print(mgroups)


###################diff_stat_barplot####################################################################################
print("barplot  for DEG statistic results is beginning:")
diff_stat_barplot(paste0(opt$type,"_diff_stat.xls"),opt$type,species)
print("barplot for DEG statistic results is OK:")
###################heatmap#############################
if (! is.null(opt$expression)) {
  registerDoParallel(cores=length(rownames(mgroups)))
  if(length(rownames(mgroups))>0)
  {   print("heatmap is beginning:")
    foreach(i=1:nrow(mgroups))%dopar%heatmap_plot(contrasts[1],expression,mgroups[i,],opt$foldchange,opt$pvalue,opt$fdr,species)
    print("heatmap is OK")
    doParallel::stopImplicitCluster()
  }
}
###################volcano##############################################################################################
registerDoParallel(cores=length(rownames(mgroups)))
print("volcano is beginning:")
foreach(i=1:nrow(mgroups))%dopar%volcano_plot(contrasts[1],mgroups[i,],opt$foldchange,opt$pvalue,opt$fdr,species)
print("volcano is OK")
doParallel::stopImplicitCluster()
