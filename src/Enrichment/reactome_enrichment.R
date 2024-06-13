#!/usr/bin/env Rscript
# by The Coder, 20160726
# Packing script(add optparse), 20180111
library("optparse")

option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-j", "--gobg"), type="character", default=NULL, help="go backgroud file", metavar="character"),
	make_option(c("-c", "--gocategory"), type="character", default=NULL, help="go category file", metavar="character"),
	make_option(c("-k", "--keggbg"), type="character", default=NULL, help="kegg backgroud file", metavar="character"),
	make_option(c("-r", "--reactomebg"), type="character", default=NULL, help="Reactome backgroud file", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript enrichment.r -i a-vs-b-Down.xls,a-vs-b-Total.xls,a-vs-b-Up.xls -j go.backgroud.xls -c category.xls -k kegg.backgroud.xls -r reactome.backgroud.xls -o enrichment");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$gobg) | is.null(opt$gocategory) |is.null(opt$reactomebg) | is.null(opt$keggbg) | is.null(opt$outpath)){
        print_help(opt_parser)
        stop("--input --gobg --gocategory --keggbg --reactomebg --outpath must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

category <- read.delim(opt$gocategory, header=F, sep="\t", quote="", row.names=1)
go.bg <- read.delim(opt$gobg, header=F, sep="\t", quote="")
kegg.bg <- read.delim(opt$keggbg, header=F, sep="\t", quote="")
reactome.bg<- read.delim(opt$reactomebg, header=F, sep="\t", quote="")

map <- function(diff.bg) {
	id <- strsplit(as.character(diff.bg[2]), ",")[[1]]
	term <- strsplit(as.character(diff.bg[3]), "|", fix=T)[[1]]
	return(as.data.frame(cbind(diff.bg[1], id, term)))
}

enrichment <- function (diff, bg) {	
	bg.total <- unlist(lapply(as.vector(bg[[2]]), function(x) strsplit(x, ",") ))
	bg.total <- as.data.frame(table(bg.total))
	diff.bg <- bg[!is.na(match(bg[,1], diff[,1])),]
	if(nrow(diff.bg)==0 ){ return(NULL) }

	d <- apply(diff.bg, 1, map)
	d <- Reduce(rbind, d)
	enrich <- as.data.frame( table(d[,2]) )
	colnames(enrich)[c(1,2)] <- c("id", "ListHits")

	for(i in 1:nrow(enrich)) {
		enrich[i, "term"] <- d[which(d[,2]==enrich[i,1]), 3][1]
		enrich[i, "Gene"] <- paste(d[which(d[,2]==enrich[i,1]), 1], collapse="; ")
		enrich[i, "PopHits"] <- bg.total[which(bg.total[,1]==as.vector(enrich[i,1])), 2]
	}
	enrich["ListTotal"] <- nrow(diff.bg)
	enrich["PopTotal"] <- nrow(bg)
	enrich["pval"] <- phyper(enrich[,"ListHits"], enrich[,"PopHits"], 
		enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)

	enrich["padj"] <- p.adjust(enrich[,"pval"], method="fdr")
	enrich["Enrichment_score"] <- 
		(enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
	enrich <- enrich[order(enrich["pval"]), ]
	return(enrich[c(1,3,2,6,5,7,8,9,10,4)])
}

infile<-strsplit(opt$input,split=",")[[1]]
for(i in 1:length(infile)) {
	f <- infile[i]
	diff <- read.delim(f, header=T, sep="\t", quote="")

	d <- enrichment(diff, go.bg)
	d["category"] <- category[as.vector(d[,1]), ]
	d <- d[, c(1, 2, ncol(d), 3:(ncol(d)-1))]

	if(!is.null(d)) { 
		write.table(d, paste0(opt$outpath, "/enrichment-go-", basename(f)), sep="\t", row.names=F, quote=F) 
	}

	d <- enrichment(diff, kegg.bg)
	if(!is.null(d)) { 
		write.table(d, paste0(opt$outpath, "/enrichment-kegg-", basename(f)), sep="\t", row.names=F, quote=F) 
	}
	d <- enrichment(diff, reactome.bg)
	if(!is.null(d)) { 
		write.table(d, paste0(opt$outpath, "/enrichment-Reactomebg-", basename(f)), sep="\t", row.names=F, quote=F) 
	}
}
