##' A universal enrichment analyzer for GO/KEGG using clusterProfile
##' @title Run_clusterProfile.r
##' @author YXF
##' @date 20180919
##' @mail yxfhenu@163.com
##' @md AnLau
##' @date 20181018
##' @mail anlau62@gmail.com


##加载R包
suppressMessages(library(clusterProfiler))
suppressMessages(library(docopt))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessag es(library(stringr))
suppressMessages(library(grid))
suppressMessages(library(oebio))
suppressMessages(library(Cairo))


########################################################################################################################
## configuration for docopt
doc <- paste0("
Usage:
	enricher.r -i <gene> -b <gobk> -c <category> -k <keggbk> -l <kegglevel>  [-m <minsize>] -p <prefix> [-o <outdir>] [-n <name>]
	enricher.r go -i <gene> -b <gobk>  -c <category>  [-m <minsize>] -p <prefix> [-o <outdir>] [-n <name>]
	enricher.r kegg -i <gene> -k <keggbk> -l <kegglevel> [-m <minsize>]  -p <prefix> [-o <outdir>] [-n <name>]
	enricher.r easyrich -i <gene> -s <background> -p <prefix> [-m <minsize>] [-o <outdir>] [-n <name>] 
	enricher.r goplot -br <goresult> -p <prefix> [-o <outdir>] 
	enricher.r keggplot -kr <keggresult> -p <prefix> -l <kegglevel> [-o <outdir>] 
	enricher.r easyplot -er <easyresult> -p <prefix> [-o <outdir>] 
	
Options:
	-i --gene  <gene> gene list ,eg:list1,list2.
	-m --minsize minimal size of genes annotated for testing term, terms smaller than this are excluded from the analysis. 
		Usually,using the default parameter:5 .If go annotation information is not complete enough,can selectly change to 2. [default: 5]
	-b --gobk <gobk> go annotation file,without header.eg:go.backgroud.txt.
	-c --category <category> go term category ,eg:category.txt.
	-k --keggbk <keggbk> kegg annotation file, without header. eg:kegg.backgroud.txt.
	-l --kegglevel <kegglevel> kegg level database, withouth header. eg:KEGGpathway_three_levels_v2.xls.
	-p --prefix <prefix> output prefix name,eg:prefix1,prefix2.
	-o --outdir <outdir> output directory [default: enrichment].
	-br --gorsult <goresult> go enrichment result file, eg:enrichment-GO-prefix.xls.
	-kr --keggrsult <keggresult> kegg enrichment result file, eg:enrich_KEGG-prefix.xls.
	-s --background <background> easyrich background file, ps:must contain 3 columns like keggbackground.
	-n --name <name> outputfile last column's header.[default: geneID]
	-er --easyresult <easyresult> easy enrichment result file, eg:enrichment-prefix.xls.

")

## docopt parsing
opt <- docopt(doc)
print(opt)
########################################################################################################################
##GO/KEGG条目长度处理
strLenLimit <- function(string,lenNum) {
  string <- as.character(string)
  if(str_length(string)>lenNum){string <- sub("[^ ]+$", "...", substr(string,1,(lenNum-3)))}
  return(string)	
}

#########################################################################################################################
##根据某一列反向提取dataframe
minimum <- function(dataframe, key, num) { 
  ascend_order = order(dataframe[key],decreasing=F)
  return(head(dataframe[ascend_order,], num))
}

##########################################################################################################################
##拆分特定列，并转换为数值类型
col_split <- function(col, spr, num, names){
  tmp = as.vector(col)
  colsplited = as.data.frame(matrix(unlist(strsplit(tmp, split = spr)),ncol=num, byrow = T))
  colnames(colsplited) = names
  for(i in names){
    colsplited[,i] = as.numeric(as.character(colsplited[,i]))
  }
  return(colsplited)
}

########################################################################################################################
## 解析背景文件
parse_bk <- function(bk){
  bk <- read.delim(bk, header = F,sep="\t")
  bk$V2 <- gsub(";", "|", gsub(",", "|", bk$V2))
  bk_parse <- separate_rows(bk, V2, V3, sep = "\\|", convert = TRUE)
  return(bk_parse)
}

#########################################################################################################################
##GO条形图绘制
go_bar_plot <- function(go_enrich,prefix,outdir){
  go_enrich <- go_enrich[order(go_enrich[,"pvalue"]),]
  
  top_bp <- head(go_enrich[which(go_enrich[,'category']=="biological_process"),],10)
  top_cc <- head(go_enrich[which(go_enrich[,'category']=="cellular_component"),],10)
  top_mf <- head(go_enrich[which(go_enrich[,'category']=="molecular_function"),],10)
  
  top_go_enrich <- rbind(top_bp,top_cc,top_mf)
  write.table(top_go_enrich, paste0(outdir,"/enrichment-GO-",prefix,"-top10.xls"),sep="\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
  top_go_enrich["Description"] <- apply(top_go_enrich["Description"], 1, strLenLimit,lenNum=50) 
  # top_go_enrich[,"Description"] <- factor(top_go_enrich[,"Description"], levels = top_go_enrich[,"Description"])
  
  p <- ggplot(data=top_go_enrich, aes(x=Description, y=-log(pvalue,10))) +
    geom_bar(fill="tomato4",stat="identity",position=position_dodge(0.7),width=0.5) +
    theme_bw()+
    theme(text=element_text(size=12),
          plot.title = element_text(hjust = 0.5, size=18, family = "ArialMT"),
          panel.grid =element_blank(),
          legend.position = "none",
          strip.background = element_rect(fill="white"))+
    facet_grid(category~.,scales="free_y")+
    labs(title = "GO Enrichment Top10", x="", y=expression('-log'[10]*' Pvalue'))+
    coord_flip()
  
  CairoPDF(paste0(outdir,"/enrichment-GO-",prefix,"-top10.pdf"),width = 12, height = 8)
  plot(p)
  dev.off()
  
  CairoPNG(paste0(outdir,"/enrichment-GO-",prefix,"-top10.png"),width = par('din')[1]*1.5, height = par('din')[2], units='in', dpi=600)
  plot(p)
  dev.off()
  
  file.remove("Rplots.pdf")
}

#########################################################################################################################
##KEGG_level1_top10气泡图绘制
kegg_level1_top10_plot <- function(kegg_enrich, kegglevel,prefix,outdir){
  kegglevel <- read.delim(kegglevel, header = F, sep="\t")
  colnames(kegglevel) <- c("ID","level1","level2","level3")
  
  oldlevel1 <- c("Cellular Processes","Environmental Information Processing",
                 "Genetic Information Processing","Human Diseases",
                 "Metabolism","Organismal Systems")
  newlevel1 <- c("CP","EIP","GIP","HD","Meta.","OS")
  
  for(i in 1:length(oldlevel1)){
    kegglevel$level1 <- sub(oldlevel1[i], newlevel1[i], kegglevel$level1)
  }
  kegg_enrich$ID <- sub("[a-zA-Z]+","ko",kegg_enrich$ID)
  kegg_enrich_level <- merge(kegg_enrich, kegglevel,by="ID")
  
  
  top_kegg_enrich_level = data.frame()
  for(level1 in newlevel1){
    min_level1 <- minimum(kegg_enrich_level[kegg_enrich_level$level1 == level1,],"pvalue",10)
    top_kegg_enrich_level <- rbind(top_kegg_enrich_level,min_level1)
  }
  write.table(top_kegg_enrich_level, paste0(outdir,"/enrichment-KEGG-",prefix,"-level1-top10.xls"),sep="\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
  xlab <- "Enrichment Score"
  size.lab <- "Number"
  p <- qplot(EnrichmentScore, Description, data=top_kegg_enrich_level, size=ListHits, 
             colour=pvalue, xlab=xlab, ylab="")+
			 ggtitle("KEGG Enrichment Level1 Top 10")+
			 theme_bw()+
			 theme(plot.title = element_text(hjust = 0.5))+
			scale_colour_gradient(low="grey80",high="darkred") +
			labs(size=size.lab) + labs(colour="pvalue") +
			facet_grid(level1~.,scales = "free_y")
  CairoPDF(paste0(outdir,"/enrichment-KEGG-",prefix,"-level1-top10.pdf"),width = 12, height = 8)
  plot(p)
  dev.off()
  
  CairoPNG(paste0(outdir,"/enrichment-KEGG-",prefix,"-level1-top10.png"),width = par('din')[1]*1.5, height = par('din')[2], units='in', dpi=600)
  plot(p)
  dev.off()
  
  file.remove("Rplots.pdf")
}

#########################################################################################################################
##easy富集绘制气泡图
easy_plot <- function(easy_enrich, prefix, outdir){
  easy_enrich <- as.data.frame(easy_enrich)
  top20 <- head(easy_enrich[order(easy_enrich[, "pvalue"]),], 20)
  if (nrow(top20) == 0) {
    print("top20 items = 0,program exit!")
    q()
  }
  write.table(top20[, c(1, 2, 8, 7, 10, 3)], paste0(outdir, "/", "enrichment-",prefix, ".top20.xls"), sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)
  # top20[, "ID"] <- sub("^path:", "", paste(top20[, 1], top20[, 2], sep = ": "))
  top20$EnrichmentScore <- sapply(top20$EnrichmentScore, function(x) eval(parse(text = x)))
  xlab <- "EnrichmentScore"
  title <- paste0(prefix, ": ", "Enrichment top 20")
  size.lab <- "Number"
  p = qplot(ListHits, ID, data = top20, size = ListHits,
            colour = pvalue, xlab = xlab, ylab = "") +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradientn(colours = rainbow(6)) +
    labs(size = size.lab) +
    labs(colour = "pvalue")
  ggsave(paste0(outdir, "/", "enrichment-",prefix, ".top20.pdf"), height = 8, width = 10, plot = p)
  ggsave(paste0(outdir, "/", "enrichment-",prefix, ".top20.png"), type = "cairo-png", height = 8, width = 10, plot = p)
  print(paste0(outdir, "/", "enrichment-",prefix,".top20.png(pdf) is OK"));
}

########################################################################################################################
##使用clusterProfiler包对非模式生物进行GO富集分析+绘图
enrich_go <- function(gene, gobk, category, minsize, prefix, outdir,name){
  print(paste0("GO for ", prefix, " is begin:"))
  go_term2gene <- data.frame(gobk$go_term, gobk$gene)
  go_term2name <- data.frame(gobk$go_term, gobk$name)

  gene_list <- read.delim(gene, header=F, sep = "\t")
  category <- read.delim(category, header = F, sep = "\t", quote = "", row.names = 1)
  
  go_enrich <- enricher(gene = gene_list$V1,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = as.numeric(minsize),
    maxGSSize = Inf,
    pAdjustMethod = "BH",
    TERM2GENE = go_term2gene,
    TERM2NAME = go_term2name)
  
  if(!is.null(go_enrich)){
    if (file.exists(outdir)== "FALSE") {dir.create(outdir)}
    
    go_enrich <- as.data.frame(go_enrich)
    go_enrich["category"] <- category[go_enrich[, 1],]
    GeneRatio = col_split(go_enrich$GeneRatio,"/",2,c("ListHits","ListTotal"))
    BgRatio = col_split(go_enrich$BgRatio,"/",2,c("PopHits","PopTotal"))
    
    go_enrich <- cbind(go_enrich, GeneRatio,BgRatio)
    go_enrich <- transform(go_enrich, EnrichmentScore=((ListHits/ListTotal)/(PopHits/PopTotal)),qvalue=NULL, GeneRatio=NULL, BgRatio=NULL,Count=NULL)
    
    tailcol <- c("pvalue","p.adjust","geneID")
    go_enrich <- go_enrich[c(setdiff(colnames(go_enrich),tailcol),tailcol)]
    names(go_enrich)[length(names(go_enrich))] <- name
      
    write.table(as.data.frame(go_enrich), paste0(outdir, "/", "enrichment-GO-",prefix, ".xls"), sep = "\t", quote = FALSE,row.names = F)
    go_bar_plot(go_enrich, prefix, outdir)
    
    print(paste0("GO for ", prefix, " is done!"))
  } else{
    print(paste0("GO for ", prefix, " is NULL!"))
  }
}

########################################################################################################################
##使用clusterProfiler包对非模式生物进行KEGG富集分析+绘图
enrich_kegg <- function(gene, keggbk, kegglevel, minsize, prefix, outdir,name) {
  print(paste0("KEGG for ", prefix, " is begin:"))
  kegg_term2gene <- data.frame(keggbk$ko_term, keggbk$gene)
  kegg_term2name <- data.frame(keggbk$ko_term, keggbk$name)
  
  gene_list <- read.delim(gene, header=F, sep = "\t")
  
  kegg_enrich <- enricher(gene = gene_list$V1,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = as.numeric(minsize),
    maxGSSize = Inf,
    TERM2GENE = kegg_term2gene,
    TERM2NAME = kegg_term2name)
  
  if (!is.null(kegg_enrich)){
    if (file.exists( outdir)== "FALSE") {dir.create(outdir)}
    kegg_enrich <- as.data.frame(kegg_enrich)
    GeneRatio = col_split(kegg_enrich$GeneRatio,"/",2,c("ListHits","ListTotal"))
    BgRatio = col_split(kegg_enrich$BgRatio,"/",2,c("PopHits","PopTotal"))
    
    kegg_enrich <- cbind(kegg_enrich, GeneRatio,BgRatio)
    kegg_enrich <- transform(kegg_enrich, EnrichmentScore=((ListHits/ListTotal)/(PopHits/PopTotal)),qvalue=NULL, GeneRatio=NULL, BgRatio=NULL,Count=NULL)
    
    tailcol <- c("pvalue","p.adjust","geneID")
    kegg_enrich <- kegg_enrich[c(setdiff(colnames(kegg_enrich),tailcol),tailcol)]
    names(kegg_enrich)[length(names(kegg_enrich))] <- name
    
    write.table(kegg_enrich, paste0(outdir, "/", "enrichment-KEGG-",prefix, ".xls"), sep = "\t",quote = FALSE, row.names = F)
    kegg_level1_top10_plot(kegg_enrich, kegglevel, prefix, outdir)
    print(paste0("KEGG for ", prefix, " is done!"))
  } else{
    print(paste0("KEGG for ", prefix, " is NULL!"))
  }
}

########################################################################################################################
##使用clusterProfiler包进行富集分析+绘制气泡图
enrich_easy <- function(gene,background,minsize, prefix, outdir,name){
  print(paste0("easyenrich for ", prefix, " is begin:"))
  easy_term2gene <- data.frame(background$term, background$gene)
  easy_term2name <- data.frame(background$term, background$name)
  
  gene_list <- read.delim(gene, header=F, sep = "\t")
  
  easy_enrich <- enricher(gene = gene_list$V1,
                          pvalueCutoff = 1,
                          qvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          minGSSize = as.numeric(minsize),
                          maxGSSize = Inf,
                          TERM2GENE = easy_term2gene,
                          TERM2NAME = easy_term2name)
  
  if (!is.null(easy_enrich)){
    if (file.exists( outdir)== "FALSE") {dir.create(outdir)}
    easy_enrich <- as.data.frame(easy_enrich)
    
    GeneRatio = col_split(easy_enrich$GeneRatio,"/",2,c("ListHits","ListTotal"))
    BgRatio = col_split(easy_enrich$BgRatio,"/",2,c("PopHits","PopTotal"))
    
    easy_enrich <- cbind(easy_enrich, GeneRatio,BgRatio)
    easy_enrich <- transform(easy_enrich, EnrichmentScore=((ListHits/ListTotal)/(PopHits/PopTotal)),qvalue=NULL, GeneRatio=NULL, BgRatio=NULL,Count=NULL)
    
    tailcol <- c("pvalue","p.adjust","geneID")
    easy_enrich <- easy_enrich[c(setdiff(colnames(easy_enrich),tailcol),tailcol)]
    names(easy_enrich)[length(names(easy_enrich))] <- name
    
    write.table(easy_enrich, paste0(outdir, "/", "enrichment-",prefix, ".xls"), sep = "\t",quote = FALSE, row.names = F)
    easy_plot(easy_enrich,prefix,outdir)
    print(paste0("enrichment for ", prefix, " is done!"))
  }else{
    print(paste0("enrichment for ", prefix, " is NULL!"))
  }
}

########################################################主函数模块########################################################
########GO富集分析+绘图
if (! is.null(opt$gobk) & ! is.null(opt$gene) & ! is.null(opt$category) & is.null(opt$keggrsult) & is.null(opt$gorsult)) {
    Sys.time()
    file_list <- c(unlist(strsplit(opt$gene, ",")))
    prefix <- c(unlist(strsplit(opt$prefix, ",")))
    gobk <- parse_bk(opt$gobk)
    colnames(gobk) <- c("gene","go_term","name")
    ##设置最大并行数
    if (length(file_list) < 10) {
        registerDoParallel(cores = length(file_list))}
    else {
        registerDoParallel(cores = 10)}
    if (file.exists(opt$outdir) == "FALSE") {
        dir.create(opt$outdir)}
    if (file.exists( paste0(opt$outdir,"/GO_enrichment"))== "FALSE") {
        dir.create(paste0(opt$outdir,"/GO_enrichment"))}
    ##执行GO富集分析函数
    print("enrichment for GO is beginning:")
    foreach(i = 1 : length(file_list)) %dopar% enrich_go(file_list[i], gobk, opt$category, opt$minsize, prefix[i], paste0(opt$outdir,"/GO_enrichment","/",prefix[i]),opt$name)
    print("enrichment for GO is done!")
    doParallel::stopImplicitCluster()
    Sys.time()
}
######KEGG富集分析+绘图
if (! is.null(opt$keggbk) & !is.null(opt$kegglevel) & ! is.null(opt$gene) & is.null(opt$keggrsult) & is.null(opt$gorsult)) {
    Sys.time()
    file_list <- c(unlist(strsplit(opt$gene, ",")))
    prefix <- c(unlist(strsplit(opt$prefix, ",")))
    keggbk <- parse_bk(opt$keggbk)
	  colnames(keggbk) <- c("gene","ko_term","name")
    ##设置最大并行数
    if (length(file_list) < 10) {
        registerDoParallel(cores = length(file_list))}
    else {
        registerDoParallel(cores = 10)}
    if (file.exists(opt$outdir) == "FALSE") {
        dir.create(opt$outdir)}
    if (file.exists( paste0(opt$outdir,"/KEGG_enrichment")) == "FALSE") {
        dir.create(paste0(opt$outdir,"/KEGG_enrichment"))}
    ##执行KEGG富集分析函数
    print("enrichment for KEGG is beginning:")
    foreach(i = 1 : length(file_list)) %dopar% enrich_kegg(file_list[i], keggbk, opt$kegglevel, opt$minsize, prefix[i], paste0(opt$outdir,"/KEGG_enrichment","/",prefix[i]),opt$name)
    print("enrichment for KEGG  is done!")
    doParallel::stopImplicitCluster()
    Sys.time()
}
######GO绘图
if (! is.null(opt$gorsult)){
  Sys.time()
  go_enrich <- read.delim(opt$gorsult)
  prefix <- opt$prefix
  outdir <- opt$outdir
  if (file.exists(outdir)== "FALSE") {dir.create(outdir)}
  print(paste0("GO plot for ",prefix, " is begining "))
  go_bar_plot(go_enrich, prefix, outdir)
  print(paste0("GO plot for ",prefix, " is ending "))
}
######KEGG绘图
if (!is.null(opt$keggrsult) & !is.null(opt$kegglevel)){
  Sys.time()
  kegg_enrich <- read.delim(opt$keggrsult)
  kegglevel <- opt$kegglevel
  prefix <- opt$prefix
  outdir <- opt$outdir
  if (file.exists( outdir)== "FALSE") {dir.create(outdir)}
  print(paste0("KEGG plot for ",prefix, " is begining "))
  kegg_level1_top10_plot(kegg_enrich, kegglevel, prefix, outdir)
  print(paste0("KEGG plot for ",prefix, " is ending "))
}
######easy富集分析
if (! is.null(opt$background) & ! is.null(opt$gene)) {
  Sys.time()
  file_list <- c(unlist(strsplit(opt$gene, ",")))
  prefix <- c(unlist(strsplit(opt$prefix, ",")))
  background <- parse_bk(opt$background)
  colnames(background) <- c("gene","term","name")
  ##设置最大并行数
  if (length(file_list) < 10) {
    registerDoParallel(cores = length(file_list))}
  else {
    registerDoParallel(cores = 10)}
  if (file.exists(opt$outdir) == "FALSE") {
    dir.create(opt$outdir)}
  ##执行KEGG富集分析函数
  print("easy enrichment is beginning:")
  foreach(i = 1 : length(file_list)) %dopar% enrich_easy(file_list[i], background, opt$minsize, prefix[i], opt$outdir,opt$name)
  print("easy enrichment is done!")
  doParallel::stopImplicitCluster()
  Sys.time()
}
######easy富集绘图
if (!is.null(opt$easyresult)){
  Sys.time()
  easy_enrich <- read.delim(opt$easyresult)
  prefix <- opt$prefix
  outdir <- opt$outdir
  if (file.exists( outdir)== "FALSE") {dir.create(outdir)}
  print(paste0("easy plot for ",prefix, " is begining "))
  easy_plot(easy_enrich, prefix, outdir)
  print(paste0("easy plot for ",prefix, " is ending "))
}
########################################################################################################################

