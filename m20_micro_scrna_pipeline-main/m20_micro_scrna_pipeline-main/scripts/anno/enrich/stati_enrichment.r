#!/usr/bin/env Rscript
# by The Coder, 20160812
# Packing script(add optparse), 20180111
library("optparse")

option_list = list(
	make_option(c("-j", "--goEnrichDir"), type="character", default="GO_enrichment", help="go enrichment directory. default: GO_enrichment", metavar="character"),
	make_option(c("-k", "--keggEnrichDir"), type="character", default="KEGG_enrichment", help="kegg enrichment directory. default: KEGG_enrichment", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript stati_enrichment.r -j GO_enrichment -k KEGG_enrichment");
opt = parse_args(opt_parser);
if(!file.exists(opt$goEnrichDir) | !file.exists(opt$keggEnrichDir)){
	print_help(opt_parser)
	stop(paste(opt$goEnrichDir, opt$keggEnrichDir, "is not exists!", sep=" "))
}
opt$goEnrichDir<-gsub("/$", "", opt$goEnrichDir)
opt$keggEnrichDir<-gsub("/$", "", opt$keggEnrichDir)

stati_go <- function(f) {
	name <- sub("^enrichment-go-", "", basename(f))
	name <- sub(".xls$", "", name)
	d <- read.delim(paste0(opt$goEnrichDir ,"/", f), header=T, sep="\t")
	pval0.05 <- nrow(d[which(d$pval < 0.05),])
	pval0.01 <- nrow(d[which(d$pval < 0.01),])
	padj0.05 <- nrow(d[which(d$padj < 0.05),])
	padj0.01 <- nrow(d[which(d$padj < 0.01),])
	return(c(name, nrow(d), pval0.05, pval0.01, padj0.05, padj0.01))
}

files <- list.files(path=opt$goEnrichDir, pattern="^enrichment-go.+.xls$", recursive=T)
stati <- matrix(unlist(lapply(files, stati_go)), nrow=6)
stati <- t(stati)
colnames(stati) <- c("groups", "Tested term", "pvalue<0.05", "pvalue<0.01", "padj<0.05", "padj<0.01")
write.table(stati, paste0(opt$goEnrichDir, "/enrichment_go.xls"), sep="\t", row.names=F, quote=F)


stati_kegg <- function(f) {
	name <- sub("^enrichment-kegg-", "", basename(f))
	name <- sub(".xls$", "", name)
	d <- read.delim(paste0(opt$keggEnrichDir, "/", f), header=T, sep="\t")
	pval0.05 <- nrow(d[which(d$pval < 0.05),])
	pval0.01 <- nrow(d[which(d$pval < 0.01),])
	padj0.05 <- nrow(d[which(d$padj < 0.05),])
	padj0.01 <- nrow(d[which(d$padj < 0.01),])
	return(c(name, nrow(d), pval0.05, pval0.01, padj0.05, padj0.01))
}

files <- list.files(path=opt$keggEnrichDir, pattern="^enrichment-kegg.+.xls$", recursive=T)
stati <- matrix(unlist(lapply(files,stati_kegg)), nrow=6)
stati <- t(stati)
colnames(stati) <- c("groups", "Tested term", "pvalue<0.05", "pvalue<0.01", "padj<0.05", "padj<0.01")
write.table(stati, paste0(opt$keggEnrichDir, "/enrichment_kegg.xls"), sep="\t", row.names=F, quote=F)

