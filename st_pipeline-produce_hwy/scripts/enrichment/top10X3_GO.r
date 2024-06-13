#!/usr/bin/env Rscript
# by The Coder, 20160527
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top10X3_GO.r -i enrichment-go-Group1-vs-Group2-Down.txt -m Down -o outdir/");
opt = parse_args(opt_parser);
if(is.null(opt$input) | is.null(opt$outpath) | is.null(opt$mark)){
	print_help(opt_parser)
	stop("--input --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$input)){
	print("input is not exists")
	q()
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-go-", "", basename(opt$input)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Total$", "(Total)", groupname)
}

limit70 <- function(s,df) {
	k <- as.character(s)
	if(str_length(s)>55){
	  go  <- subset(df, df['Term'] == s)
	  go <- go$id[1]
	  term <- sub("[^ ]+$", "...", substr(k,1,55))
	  k <- paste0(go,":",term)
	}
	return(k)	
}

top10 <- function(i) { return(i[order(head(i, 10)["p-value"]),]) }

d <- read.delim(opt$input, sep="\t", header=T, quote="",check.names=F)
# d <- d[which(d[,"ListHits"]>2),]
d <- d[order(d[,"p-value"]),]

#stopifnot(nrow(d)>0)
if(nrow(d)==0){
    print("d items = 0, program exit!")
    q()
}

dp <- rbind(top10(d[d["Category"]=="biological_process", ]),
	top10(d[d["Category"]=="cellular_component", ]),
	top10(d[d["Category"]=="molecular_function", ]))
write.table(dp[,c(1,2,3,4,8,10,11)], paste0(opt$outpath, "/GO.top.", opt$mark, ".xls"), sep="\t", quote=FALSE,
	col.names=TRUE, row.names=FALSE)

d1<-dp[which(dp[,3]=="biological_process"),]
#d1_s<-d1[order(d1[,4],decreasing=F),]
d2<-dp[which(dp[,3]=="cellular_component"),]
#d2_s<-d2[order(d2[,4],decreasing=F),]
d3<-dp[which(dp[,3]=="molecular_function"),]
#d3_s<-d3[order(d3[,4],decreasing=F),]
d_l<-rbind(d1,d2,d3)
d_l[which(dp$Category=="biological_process"), "color"] <- "#4DAF4A"
d_l[which(dp$Category=="cellular_component"), "color"] <- "#377EB8"
d_l[which(dp$Category=="molecular_function"), "color"] <- "#E41A1C"

#d_l["Term"] <- apply(d_l["Term"], 1, limit70)
d_l$Term <- factor(d_l$Term, levels=d_l$Term)
colnames(d_l)[8] = "pvalue"
print(d_l$pvalue)
if( max(d_l$pvalue) == "0" ) {
	min_y = 0
	max_y = 50
	p=ggplot(data=d_l, aes(x=Term, y=-log(pvalue,10), width=0.6, fill=Category,space=0.6,cex.main=3,cex.lab=2)) +
	    geom_bar(stat="identity",position=position_dodge(0.7),width=0.5) + ylim(min_y,max_y) +
	    labs(x="", y=expression('-log'[10]*' p-value'), title = paste0(groupname,": ","Top 30 GO Term")) +
	    theme_bw() + scale_fill_manual(values=unique(d_l$color))  +
	    theme(axis.text.x=element_text(angle = 45, hjust=1, size=14,color=d_l$color))+
	    theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
	    theme(legend.text=element_text(size=10))+
	    theme(text=element_text(size=14)) +
	    theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
	    theme(plot.margin=unit(c(2,2,2,18), "lines"))+
	    theme(panel.grid =element_blank()) 
} else {
	p=ggplot(data=d_l, aes(x=Term, y=-log(pvalue,10), width=0.6, fill=Category,space=0.6,cex.main=3,cex.lab=2)) +
	    geom_bar(stat="identity",position=position_dodge(0.7),width=0.5)  +
	    labs(x="", y=expression('-log'[10]*' p-value'), title = paste0(groupname,": ","Top 30 GO Term")) +
	    theme_bw() + scale_fill_manual(values=unique(d_l$color))  +
	    theme(axis.text.x=element_text(angle = 45, hjust=1, size=14,color=d_l$color))+
	    theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
	    theme(legend.text=element_text(size=10))+
	    theme(text=element_text(size=14)) +
	    theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
	    theme(plot.margin=unit(c(2,2,2,18), "lines"))+
	    theme(panel.grid =element_blank()) 
}

ggsave(paste0(opt$outpath, "/GO.top.", opt$mark, ".pdf"), height=10, width=18, plot=p)
ggsave(paste0(opt$outpath, "/GO.top.", opt$mark, ".png"), type="cairo-png", height=10, width=18, plot=p)
print(paste0(opt$outpath, "/GO.top.", opt$mark, ".png(pdf) is OK"));
