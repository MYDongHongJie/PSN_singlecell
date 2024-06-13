#!/usr/bin/env Rscript
library("optparse")

option_list = list(
	make_option(c("-d", "--diff"), type="character", default=NULL, help="diff result file", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
	make_option(c("-j", "--gobg"), type="character", default=NULL, help="go backgroud file", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript topGO.r -d diff-Group2-vs-Group1-Total.xls -m Total -j go.backgroud.xls -o outdir/");
opt = parse_args(opt_parser);
if (is.null(opt$diff) | is.null(opt$gobg) | is.null(opt$outpath) | is.null(opt$mark)){
	print_help(opt_parser)
	stop("--diff --gobg --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

library(topGO)
diff_genes<-read.delim(opt$diff, header=T, sep="\t", quote="")
geneID2GO<-readMappings(opt$gobg)
interesting_genes=factor(diff_genes[,1])
geneNames<- names(geneID2GO)
geneList <- factor(as.integer (geneNames %in% interesting_genes))
names(geneList)=geneNames

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata,classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 20)

anno_gene <-max(allRes[,4])
if (anno_gene>2){
	png(paste0(opt$outpath, "/topGO_MF_", opt$mark, ".png"), width=3000, height=3000, res=300)
	res<-try(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'def'))
	if(inherits(res, "try-error")){
	} else{
	dev.off()
	pdf(paste0(opt$outpath, "/topGO_MF_", opt$mark, ".pdf"))
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'def')
	dev.off()
	print(paste0(opt$outpath, "/topGO_MF_", opt$mark, ".pdf(png) is OK"))}
}

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata,classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 20)
anno_gene <-max(allRes[,4])
if (anno_gene>2){
	png(paste0(opt$outpath, "/topGO_BP_", opt$mark, ".png"), width=3000, height=3000, res=300)
	res<-try(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'def'))
	if(inherits(res, "try-error")){
	} else{
	dev.off()
	pdf(paste0(opt$outpath, "/topGO_BP_", opt$mark, ".pdf"))
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'def')
	dev.off()
	print(paste0(opt$outpath, "/topGO_BP_", opt$mark, ".pdf(png) is OK"))}
}

GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata,classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 20)

anno_gene <-max(allRes[,4])
if (anno_gene>2){
	png(paste0(opt$outpath, "/topGO_CC_", opt$mark, ".png"), width=3000, height=3000, res=300)
	res<-try(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'def'))
	if(inherits(res, "try-error")){
	} else{
	dev.off()
	pdf(paste0(opt$outpath, "/topGO_CC_", opt$mark, ".pdf"))
	showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'def')
	dev.off()
	print(paste0(opt$outpath, "/topGO_CC_", opt$mark, ".pdf(png) is OK"))}
}

