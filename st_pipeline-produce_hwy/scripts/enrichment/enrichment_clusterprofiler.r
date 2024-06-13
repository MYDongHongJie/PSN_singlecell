##' A universal enrichment analyzer for enrichment using clusterProfile
##' @title clusterProfile.r
##' @author YXF
##' @date 20180919
##' @mail yxfhenu@163.com
##' @md AnLau
##' @date 20181018
##' @date 20181122
##' @mail liuan_62@hotmail.com


##加载R包
suppressMessages(library(clusterProfiler))
suppressMessages(library(docopt))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(grid))
#suppressMessages(library(oebio))
#suppressMessages(library(oeRtools))
suppressMessages(library(Cairo))
#library("oeRtools")

########################################################################################################################
## configuration for docopt
doc <- paste0("
Usage:
	enricher.r -i <gene> -j <gobk> -c <category> -k <keggbk>   [-m <minsize>] -p <prefix> [-o <outdir>] [-n <name>]
	
Options:
	-i --gene  <gene> gene list ,eg:list1,list2.
	-m --minsize minimal size of genes annotated for testing term, terms smaller than this are excluded from the analysis. 
		Usually,using the default parameter:5 .If go annotation information is not complete enough,can selectly change to 2. [default: 5]
	-j --gobk <gobk> go annotation file,without header.eg:go.backgroud.txt.
	-c --category <category> go term category ,eg:category.txt.
	-k --keggbk <keggbk> kegg annotation file, without header. eg:kegg.backgroud.txt.
	-p --prefix <prefix> output prefix name,eg:prefix1,prefix2.
	-o --outdir <outdir> output directory [default: enrichment].
	-s --background <background> easyrich background file, ps:must contain 3 columns like keggbackground.
	-n --name <name> outputfile last column's header.[default: geneID]
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
  print("first step done")
  bk$V2 <- gsub(";", "|", gsub(",", "|", bk$V2))
  print("2rd step done")
  bk_parse <- separate_rows(bk, V2, V3, sep = "\\|", convert = TRUE)
  print("final step done")
  return(bk_parse)
}

#########################################################################################################################

########################################################################################################################

########################################################################################################################
##使用clusterProfiler包对非模式生物进行GO富集分析+绘图
enrich_go <- function(gene, gobk, category, minsize, prefix, outdir,name){
  print(paste0("GO for ", prefix, " is begin:"))
  go_term2gene <- data.frame(gobk$go_term, gobk$gene)
  go_term2name <- data.frame(gobk$go_term, gobk$name)

  gene_list <- read.delim(gene, header=F, sep = "\t")
  categoryfile <- read.delim(category, header = F, sep = "\t", quote = "", row.names = 1)
  #go_level2_file = paste0(dirname(category),"/GOterm_three_levels_v2.xls",sep="")
  #go_level2 <- read.delim(go_level2_file, header=F, sep = "\t")  #level2 info
  #names(go_level2) = c("id","Category","GO_classify2","Term")
  go_enrich <- enricher(gene = gene_list$V1,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = as.numeric(minsize),
    maxGSSize = Inf,
    pAdjustMethod = "BH",
    TERM2GENE = go_term2gene,
    TERM2NAME = go_term2name)
  
  if(!is.null(go_enrich)){
    if (file.exists(outdir)== "FALSE") {dir.create(outdir,recursive = TRUE)}
    
    go_enrich <- as.data.frame(go_enrich)
    go_enrich["Category"] <- categoryfile[go_enrich[, 1],]
    GeneRatio = col_split(go_enrich$GeneRatio,"/",2,c("ListHits","ListTotal"))
    BgRatio = col_split(go_enrich$BgRatio,"/",2,c("PopHits","PopTotal"))
    
    go_enrich <- cbind(go_enrich, GeneRatio,BgRatio)
    go_enrich <- transform(go_enrich, EnrichmentScore=((ListHits/ListTotal)/(PopHits/PopTotal)),qvalue=NULL, GeneRatio=NULL, BgRatio=NULL,Count=NULL)
    
    tailcol <- c("pvalue","p.adjust","geneID")    
    go_enrich <- go_enrich[c(setdiff(colnames(go_enrich),tailcol),tailcol)]
    go_enrich$geneID <- gsub("/",";",go_enrich$geneID)
    names(go_enrich)[length(names(go_enrich))] <- name
    names(go_enrich)[names(go_enrich)=="ID"]="id"
    names(go_enrich)[names(go_enrich)=="pvalue"]="p-value"
    names(go_enrich)[names(go_enrich)=="p.adjust"]="q-value"
    names(go_enrich)[names(go_enrich)=="EnrichmentScore"]="Enrichment_score"
    names(go_enrich)[names(go_enrich)=="Description"]="Term"
    #go_enrich_reslut <- merge(go_enrich,go_level2[,c(1,3)],by.x="id",by.y="id",all.x=T)
    #go_enrich <- go_enrich_reslut[c(1,2,3,12,4,5,6,7,9,10,8,11)]
    go_enrich <- go_enrich[c(1,2,3,4,5,6,7,9,10,8,11)]
    write.table(as.data.frame(go_enrich), paste0(outdir, "/", "enrichment-go-",prefix, ".xls"), sep = "\t", quote = FALSE,row.names = F)
    #go_bar_plot(go_enrich, prefix, outdir)
    
    print(paste0("GO for ", prefix, " is done!"))
  } else{
    print(paste0("GO for ", prefix, " is NULL!"))
  }
}

########################################################################################################################
##使用clusterProfiler包对非模式生物进行KEGG富集分析+绘图
enrich_kegg <- function(gene, keggbk, minsize, prefix, outdir,name,category) {
  print(paste0("KEGG for ", prefix, " is begin:"))

  kegg_term2gene <- data.frame(keggbk$ko_term, keggbk$gene)
  kegg_term2name <- data.frame(keggbk$ko_term, keggbk$name)
	print("333")
  TLA =unique(gsub("\\d","",as.character(head(keggbk$ko_term)))) #三字母缩写
  gene_list <- read.delim(gene, header=F, sep = "\t")
      print("222")
  kegg_level2_file = paste0(dirname(category),"/KEGGpathway_three_levels_v2.xls",sep="")
  print("111")
  kegg_level2 <- read.delim(kegg_level2_file, header=F, sep = "\t")  #level2 info
  names(kegg_level2) = c("id","Classification_level1","Classification_level2","Term")
  kegg_level2['id'] = gsub("ko",TLA,kegg_level2$id)
  print("000")
  kegg_enrich <- enricher(gene = gene_list$V1,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = as.numeric(minsize),
    maxGSSize = Inf,
    TERM2GENE = kegg_term2gene,
    TERM2NAME = kegg_term2name)

  if (!is.null(kegg_enrich)){
    if (file.exists( outdir)== "FALSE") {dir.create(outdir,recursive = TRUE)}
    kegg_enrich <- as.data.frame(kegg_enrich)
    GeneRatio = col_split(kegg_enrich$GeneRatio,"/",2,c("ListHits","ListTotal"))
    BgRatio = col_split(kegg_enrich$BgRatio,"/",2,c("PopHits","PopTotal"))

    kegg_enrich <- cbind(kegg_enrich, GeneRatio,BgRatio)
    kegg_enrich <- transform(kegg_enrich, EnrichmentScore=((ListHits/ListTotal)/(PopHits/PopTotal)),qvalue=NULL, GeneRatio=NULL, BgRatio=NULL,Count=NULL)
    
    tailcol <- c("pvalue","p.adjust","geneID")
    kegg_enrich <- kegg_enrich[c(setdiff(colnames(kegg_enrich),tailcol),tailcol)]
    kegg_enrich$geneID <- gsub("/",";",kegg_enrich$geneID)

    names(kegg_enrich)[length(names(kegg_enrich))] <- name
    names(kegg_enrich)[names(kegg_enrich)=="ID"]="id"
    names(kegg_enrich)[names(kegg_enrich)=="pvalue"]="p-value"
    names(kegg_enrich)[names(kegg_enrich)=="p.adjust"]="q-value"
    names(kegg_enrich)[names(kegg_enrich)=="EnrichmentScore"]="Enrichment_score"
    names(kegg_enrich)[names(kegg_enrich)=="Description"]="Term"
    kegg_enrich_reslut <- merge(kegg_enrich,kegg_level2[,c(1,2,3)],by.x="id",by.y="id",all.x=T)
    kegg_enrich <- kegg_enrich_reslut[c(1,2,11,12,3,4,5,6,8,9,7,10)]
    prefix_ <- gsub("-Total","",prefix)
    prefix_ <- gsub("-Up","",prefix_)
    prefix_ <- gsub("-Down","",prefix_)
    # kegg_enrich$hyperlink_only_excel <-  paste0( "=HYPERLINK(\"../../3.KEGG_map/",prefix_,"/",kegg_enrich[ ,"id"],".html\"", ",\"",kegg_enrich[ ,"id"] ,"\")")
    write.table(kegg_enrich, paste0(outdir, "/", "enrichment-kegg-",prefix, ".xls"), sep = "\t",quote = FALSE, row.names = F)
    print(paste0("KEGG for ", prefix, " is done!"))
  } else{
    print(paste0("KEGG for ", prefix, " is NULL!"))
  }
}

########################################################################################################################
##使用clusterProfiler包进行富集分析+绘制气泡图
########################################################主函数模块########################################################
########GO富集分析+绘图
if (! is.null(opt$gobk) & ! is.null(opt$gene) & ! is.null(opt$category) & is.null(opt$keggrsult) & is.null(opt$gorsult)) {
    Sys.time()
    file_list <- c(unlist(strsplit(opt$gene, ",")))
    prefix <- c(unlist(strsplit(opt$prefix, ",")))
    gobk <- parse_bk(opt$gobk)
	print('bkdone')
    colnames(gobk) <- c("gene","go_term","name")
    ##设置最大并行数
    if (length(file_list) < 10) {
        registerDoParallel(cores = length(file_list))}
    else {
        registerDoParallel(cores = 10)}
    if (file.exists(opt$outdir) == "FALSE") {
        dir.create(opt$outdir,recursive = TRUE)}

    ##执行GO富集分析函数
    print("enrichment for GO is beginning:")
    foreach(i = 1 : length(file_list)) %dopar% enrich_go(file_list[i], gobk, opt$category, opt$minsize, prefix[i], opt$outdir,opt$name)
    print("enrichment for GO is done!")
    doParallel::stopImplicitCluster()
    Sys.time()
}
######KEGG富集分析+绘图
if (! is.null(opt$keggbk) & ! is.null(opt$gene) & is.null(opt$keggrsult) & is.null(opt$gorsult)) {
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
        dir.create(opt$outdir,recursive = TRUE)}

    ##执行KEGG富集分析函数
    print("enrichment for KEGG is beginning:")
    #foreach(i = 1 : length(file_list)) %dopar% enrich_kegg(file_list[i], keggbk, opt$minsize, prefix[i], paste0(opt$outdir,"/2.KEGG_enrichment","/",prefix[i]),opt$name,opt$category)
    foreach(i = 1 : length(file_list)) %dopar% enrich_kegg(file_list[i], keggbk, opt$minsize, prefix[i], opt$outdir,opt$name,opt$category)
    print("enrichment for KEGG  is done!")
    doParallel::stopImplicitCluster()
    Sys.time()
}
######单独GO绘图
######easy富集分析
########################################################################################################################

