#!/usr/bin/env Rscript
library("optparse")

option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top20_KEGG.r -i enrichment-kegg-Group1-vs-Group2-Down.txt -m Down -o outdir/");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$outpath) | is.null(opt$mark)){
	print_help(opt_parser)
	stop("--input --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

library(ggplot2);

enrich <- read.delim(opt$input, head=TRUE, sep="\t", quote="")
enrich <- enrich[which(enrich["ListHits"]>2), ]

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-kegg-", "", basename(opt$input)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
	groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
	groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
	groupname <- gsub("-Total$", "(Total)", groupname)
}

top20 <- head(enrich[order(enrich[,"pval"]),], 20)
#stopifnot(nrow(top20)>0)
if(nrow(top20)==0){
    print("top20 KEGG terms is equal to zero, program exit!")
    q()
}
write.table(top20[,c(1,2,7,9,10)], paste0(opt$outpath, "/KEGG.top.", opt$mark, ".xls"), sep="\t",
	quote=FALSE, col.names=TRUE, row.names=FALSE)
top20[, "term"] <- sub("^path:", "", paste(top20[,1], top20[,2], sep=": "))

xlab <- "Enrichment Score"
title <- paste0(groupname, ": ", "KEGG Enrichment top 20")
size.lab <- "Number"
p=qplot(Enrichment_score, term, data=top20, size=ListHits, 
	colour=pval, xlab=xlab, ylab="")+ggtitle(title)+
	theme_bw()+
	theme(plot.title = element_text(hjust = 0.5))+
	scale_colour_gradientn(colours=rainbow(6)) +
	labs(size=size.lab) + labs(colour="pvalue")+guides(size = guide_legend(order = 1))
ggsave(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".pdf"), height=8, width=10, plot=p)
ggsave(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png"), type="cairo-png", height=8, width=10, plot=p)
print(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png(pdf) is OK"));
