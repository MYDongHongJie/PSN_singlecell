#!/usr/local/bin/r
library("optparse")

option_list = list(
	make_option(c("-b", "--bpGOterm"), type="character", default=NULL, help="input file", metavar="character"),
	make_option(c("-m", "--mfGOterm"), type="character", default=NULL, help="input file", metavar="character"),
	make_option(c("-c", "--ccGOterm"), type="character", default=NULL, help="input file", metavar="character"),
	make_option(c("-l",  "--level2"), type="character", default=NULL, help="input file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$bpGOterm) | is.null(opt$mfGOterm) | is.null(opt$ccGOterm) | is.null(opt$level2)){
	print_help(opt_parser)
	stop("--bpGOterm --mfGOterm --ccGOterm --level2 must be supplied", call.=FALSE)
}

library("GO.db")

xx.bp <- as.list(GOBPANCESTOR)
bp.id<-names(xx.bp)
xx.cc <- as.list(GOCCANCESTOR)
cc.id<-names(xx.cc)
xx.mf <- as.list(GOMFANCESTOR)
mf.id<-names(xx.mf)

endnode<-read.csv(opt$level2, header=F)
endnode<-as.matrix(endnode)

d.BP<-as.matrix(read.delim(opt$bpGOterm,header=T,sep="\t"))
d.MF<-as.matrix(read.delim(opt$mfGOterm,header=T,sep="\t"))
d.CC<-as.matrix(read.delim(opt$ccGOterm,header=T,sep="\t"))

BP<-d.BP[,1]
MF<-d.MF[,1]
CC<-d.CC[,1]

pos.bp<-match(BP,bp.id)
pos.cc<-match(CC,cc.id)
pos.mf<-match(MF,mf.id)

ancestor.bp<-xx.bp[pos.bp]
ancestor.cc<-xx.cc[pos.cc]
ancestor.mf<-xx.mf[pos.mf]

####childrens ancestor
matrix.bp<-matrix(c("term","ancestor"),ncol=2)
for(i in 1:length(ancestor.bp)){
	d<-match(ancestor.bp[[i]],endnode[,1])
	pos<-d[!is.na(d)]
	term.bp<-names(ancestor.bp[i])
	bp.ancestor<-paste(endnode[pos,1],collapse=",",sep="")
	matrix.bp1<-c(term.bp,bp.ancestor)
	matrix.bp<-rbind(matrix.bp,matrix.bp1)
	}
bpPrefix<-opt$bpGOterm
bpPrefix<-gsub(".\\w*$","",bpPrefix,perl=TRUE)
write.table(matrix.bp,paste(bpPrefix,".ancestor.xls",sep=""),col.names=F,row.names=F,sep="\t",quote=F)

matrix.cc<-matrix(c("term","ancestor"),ncol=2)
for(i in 1:length(ancestor.cc)){
	d<-match(ancestor.cc[[i]],endnode[,1])
	pos<-d[!is.na(d)]
	term.cc<-names(ancestor.cc[i])
	cc.ancestor<-paste(endnode[pos,1],collapse=",",sep="")
	matrix.cc1<-c(term.cc,cc.ancestor)
	matrix.cc<-rbind(matrix.cc,matrix.cc1)
	}
ccPrefix<-opt$ccGOterm
ccPrefix<-gsub(".\\w*$","",ccPrefix,perl=TRUE)
write.table(matrix.cc,paste(ccPrefix,".ancestor.xls",sep=""),col.names=F,row.names=F,sep="\t",quote=F)

matrix.mf<-matrix(c("term","ancestor"),ncol=2)
for(i in 1:length(ancestor.mf)){
	d<-match(ancestor.mf[[i]],endnode[,1])
	pos<-d[!is.na(d)]
	term.mf<-names(ancestor.mf[i])
	mf.ancestor<-paste(endnode[pos,1],collapse=",",sep="")
	matrix.mf1<-c(term.mf,mf.ancestor)
	matrix.mf<-rbind(matrix.mf,matrix.mf1)
	}
mfPrefix<-opt$mfGOterm
mfPrefix<-gsub(".\\w*$","",mfPrefix,perl=TRUE)
write.table(matrix.mf,paste(mfPrefix,".ancestor.xls",sep=""),col.names=F,row.names=F,sep="\t",quote=F)

