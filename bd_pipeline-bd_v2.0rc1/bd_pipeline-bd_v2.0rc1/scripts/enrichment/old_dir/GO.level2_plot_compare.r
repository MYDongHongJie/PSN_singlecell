#!/usr/bin/env Rscript
library("optparse")

option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="col3,col4 names in input file, eg:Up,Down", metavar="character"),
	make_option(c("-n", "--group"), type="character", default=NULL, help="group name, eg: Group2-vs-Group1", metavar="character"),
	make_option(c("-p", "--prefix"), type="character", default=NULL, help="outfile prefix", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript GO.level2_plot_compare.r -i up_down_GO.level2.stat.xls -n Group2-vs-Group1 -m Up,Down -p up_down_GO.level2.stat -o outdir/");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$outpath) | is.null(opt$mark) | is.null(opt$prefix) | is.null(opt$group)){
	print_help(opt_parser)
	stop("--input --outpath --mark --group --prefix must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

d<-read.table(opt$input, header=F, sep="\t", quote="", check.names=F, comment.char = "#", colClasses = "character")
line12<-read.table(opt$input, header=F, sep="\t", quote="", check.names=F, comment.char = "", colClasses = "character", nrows= 2)

mark1 <- strsplit(as.character(opt$mark), ",")[[1]][1]
mark2 <- strsplit(as.character(opt$mark), ",")[[1]][2]

numbers1 <- as.numeric(line12[2,3])
numbers2 <- as.numeric(line12[2,4])

three_numbers1<-floor(numbers1*0.1) 
two_numbers1<-floor(numbers1*0.01)
one_numbers1<-floor(numbers1*0.001)

three_numbers2<-floor(numbers2*0.1)
two_numbers2<-floor(numbers2*0.01)
one_numbers2<-floor(numbers2*0.001)

d[which(d[,3]<0.1),3]<-0.1
d[which(d[,4]<0.1),4]<-0.1
d1<-d[which(d[,1]=="biological process"),]
d1_s<-d1[order(d1[,2],decreasing=F),]
d2<-d[which(d[,1]=="cellular component"),]
d2_s<-d2[order(d2[,2],decreasing=F),]
d3<-d[which(d[,1]=="molecular function"),]
d3_s<-d3[order(d3[,2],decreasing=F),]
d_l<-rbind(d1_s,d2_s,d3_s)

d_l$V3<-as.numeric(d_l$V3)
d_l$V4<-as.numeric(d_l$V4)
d_l$V3<-log10(d_l$V3+0.001)+1
d_l$V4<-log10(d_l$V4+0.001)+1
dd<-d_l[,c(2,3,4)]
names(dd)=c("term",mark1,mark2)

if (numbers1<100){
a=119.7
b=119.7
c=119.7
d=120.5
} else if ((numbers1<1000)&(numbers1>=100)){
a=119.7
b=119.7
c=120.5
d=121.4
} else if ((numbers1<10000)&(numbers1>=1000)){
a=119.7
b=120.5
c=121.4
d=122.4
} else {
a=120.5
b=121.4
c=122.4
d=123.1
}

rownames(dd)<-dd$term
data<-dd[,c(2,3)]
s<-t(data)
s<-as.matrix(s)
if(grepl("up,down", tolower(opt$mark), perl=T)){
	colors<-c("#FF6A6A", "#00FF7F")
} else {
	colors<-c("#377EB8", "#E41A1C")
}

#png(paste(opt$outpath, "/", opt$prefix, ".png", sep=""),height=1000,width=2000)
png(paste(opt$outpath, "/", opt$prefix, ".png", sep=""),height=4180,width=8360,res=300)
par(mar=c(30,10,8,8),mgp=c(4.5,1,0))
bp<-barplot(s,width=0.8,space=c(0,0.2),beside=TRUE,border=NA,ylim=c(0,3),las=3,col=colors, args.legend=list(bty="n",cex=2),horiz = FALSE,legend.text=rownames(s),main=paste0("Gene Ontology Classification","(",opt$group,")"),ylab="percentage of genes",axes=F,cex.names=1.5,cex.axis=2,cex.main=3,cex.lab=2)
#bp<-barplot(s,beside=TRUE,border=NA,ylim=c(0,3),las=3,col=colors, args.legend=list(bty="n",cex=2),horiz = FALSE,legend.text=rownames(s),main=paste0("Gene Ontology Classification","(",opt$group,")"),ylab="percentage of genes",axes=F,cex.names=1.5,cex.axis=2,cex.main=3,cex.lab=2)

axis(2,at=c(0:3),las=1,labels=c(0.1,1,10,100),cex.axis=2)
lines(c(1.1,1.1),c(-1.15,-2.58),xpd=T,col="#4DAF4A",lwd=2.5)
lines(c(1.1,39.7),c(-2.58,-2.58),xpd=T,col="#4DAF4A",lwd=2.5)
lines(c(39.7,39.7),c(-2.58,-1.4),xpd=T,col="#4DAF4A",lwd=2.5)
text(29,-2.72,"Biological Process",pos=2,xpd=T,cex=2.4,col="#4DAF4A")

lines(c(41.5,41.5),c(-0.3,-2.58),xpd=T,col="#377EB8",lwd=2.5)
lines(c(41.5,75.1),c(-2.58,-2.58),xpd=T,col="#377EB8",lwd=2.5)
lines(c(75.1,75.1),c(-2.58,-0.65),xpd=T,col="#377EB8",lwd=2.5)
text(68,-2.72,"Cellular Component",pos=2,xpd=T,cex=2.4,col="#377EB8")

lines(c(76.9,76.9),c(-1.3,-2.58),xpd=T,col="#E41A1C",lwd=2.5)
lines(c(76.9,112),c(-2.58,-2.58),xpd=T,col="#E41A1C",lwd=2.5)
lines(c(112,112),c(-2.58,-1.1),xpd=T,col="#E41A1C",lwd=2.5)
text(103,-2.72,"Molecular Function",pos=2,xpd=T,cex=2.4,col="#E41A1C")

text(a,-0.2,one_numbers2,pos=2,xpd=T,cex=2,col=colors[2])
text(b,0.8,two_numbers2,pos=2,xpd=T,cex=2,col=colors[2])
text(c,1.8,three_numbers2,pos=2,xpd=T,cex=2,col=colors[2])
text(d,2.8,numbers2,pos=2,xpd=T,cex=2,col=colors[2])

axis(4,at=c(0:3),las=1, labels=c(one_numbers1,two_numbers1,three_numbers1,numbers1),cex.axis=2,col.axis=colors[1])
mtext("gene number",side=4,line=6,cex=2)
dev.off()
print(paste(opt$outpath, "/", opt$prefix, ".png is OK", sep=""))

pdf(paste(opt$outpath, "/", opt$prefix, ".pdf", sep=""),height=14,width=28)
par(mar=c(30,10,8,8),mgp=c(4.5,1,0))
bp<-barplot(s,width=0.8,space=c(0,0.2),beside=TRUE,border=NA,ylim=c(0,3),las=3,col=colors, args.legend=list(bty="n",cex=2),horiz = FALSE,legend.text=rownames(s),main=paste0("Gene Ontology Classification","(",opt$group,")"),ylab="percentage of genes",axes=F,cex.names=1.5,cex.axis=2,cex.main=3,cex.lab=2)
axis(2,at=c(0:3),las=1,labels=c(0.1,1,10,100),cex.axis=2)
lines(c(1.1,1.1),c(-1.15,-2.58),xpd=T,col="#4DAF4A",lwd=2.5)
lines(c(1.1,39.7),c(-2.58,-2.58),xpd=T,col="#4DAF4A",lwd=2.5)
lines(c(39.7,39.7),c(-2.58,-1.4),xpd=T,col="#4DAF4A",lwd=2.5)
text(29,-2.72,"Biological Process",pos=2,xpd=T,cex=2.4,col="#4DAF4A")

lines(c(41.5,41.5),c(-0.3,-2.58),xpd=T,col="#377EB8",lwd=2.5)
lines(c(41.5,75.1),c(-2.58,-2.58),xpd=T,col="#377EB8",lwd=2.5)
lines(c(75.1,75.1),c(-2.58,-0.65),xpd=T,col="#377EB8",lwd=2.5)
text(68,-2.72,"Cellular Component",pos=2,xpd=T,cex=2.4,col="#377EB8")

lines(c(76.9,76.9),c(-1.3,-2.58),xpd=T,col="#E41A1C",lwd=2.5)
lines(c(76.9,112),c(-2.58,-2.58),xpd=T,col="#E41A1C",lwd=2.5)
lines(c(112,112),c(-2.58,-1.1),xpd=T,col="#E41A1C",lwd=2.5)
text(103,-2.72,"Molecular Function",pos=2,xpd=T,cex=2.4,col="#E41A1C")

text(a,-0.2,one_numbers2,pos=2,xpd=T,cex=2,col=colors[2])
text(b,0.8,two_numbers2,pos=2,xpd=T,cex=2,col=colors[2])
text(c,1.8,three_numbers2,pos=2,xpd=T,cex=2,col=colors[2])
text(d,2.8,numbers2,pos=2,xpd=T,cex=2,col=colors[2])

axis(4,at=c(0:3),las=1, labels=c(one_numbers1,two_numbers1,three_numbers1,numbers1),cex.axis=2,col.axis=colors[1])
mtext("gene number",side=4,line=6,cex=2)
dev.off()
print(paste(opt$outpath, "/", opt$prefix, ".pdf is OK", sep=""))
