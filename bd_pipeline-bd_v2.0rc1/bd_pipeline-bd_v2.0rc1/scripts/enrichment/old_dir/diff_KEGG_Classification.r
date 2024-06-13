#!/usr/bin/env Rscript
library("optparse")

option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="mark is shown in the picture title, eg:Group2-vs-Group1(Up)", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, help="output file prefix", metavar="character")
);
opt_parser = OptionParser(option_list=option_list, epilogue = "Rscript diff_KEGG_Classification.r -i diff-KEGG_Classification.xls -m Group2-vs-Group1(Up) -o out/KEGG_Classification");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$output) | is.null(opt$mark)){
	print_help(opt_parser)
	stop("--input --output --mark  must be supplied", call.=FALSE)
}

library(ggplot2)
data <- read.delim(opt$input, sep="\t", header=T, quote="")

data$Classification_level2 <- factor(data$Classification_level2, levels=data$Classification_level2)
mylabel=paste(data$Classification_level1, data$Classification_level2, sep="--")
p=ggplot(data=data, aes(x=Classification_level2, y=percentage, fill=Classification_level1)) +coord_flip()+
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.7) +
  geom_text(aes(label=gene_number),hjust=-0.2, vjust=0.4,size=2.5) + 
  #scale_x_discrete(breaks=data$Classification_level2, labels=mylabel)+
  labs(x="", y="Percent of Genes(%)" ,title=paste0(opt$mark, ": ","KEGG Pathway Classification")) +
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=9, family = "ArialMT")) +
  theme_bw() + theme(panel.grid =element_blank()) +
  theme(legend.title=element_blank())+
  theme(legend.key.height=unit(0.5,"cm"),legend.text=element_text(size=8))+
  theme(axis.text.y=element_text(size=8,color="black")) + 
  theme(axis.text.x=element_text(hjust=1, size=8,color="black"))
  #legend.position="none")

ggsave(paste0(opt$output, ".pdf"), height=6, width=11, plot=p)
ggsave(paste0(opt$output, ".png"), type="cairo-png", height=6, width=11, plot=p)
print(paste0(opt$output, ".png(pdf) is OK"));
