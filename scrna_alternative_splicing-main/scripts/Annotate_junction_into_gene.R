
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

#######################
# 导入自定义函数
#######################

source("scripts/MARVEL_functions.R")

#######################
# 命令行参数
#######################

option_list <- list(
  make_option("--junction", type="character", default=NULL, help="样本文件"),
  make_option("--gtf", type="character", default=NULL, help="gtf文件"),
  make_option("--configfile", type="character", default="config/config.yaml", help="config/config.yaml文件，我们直接从文件中读取参数"),
  make_option("--out", type="character", default="features.anno.tsv", help="输出已经注释好的feature文件")
)

args <- parse_args(OptionParser(option_list=option_list))

junctionfile <- args$junction
gtffile <- args$gtf
outfile <- args$out
configfile <- args$configfile

config <- yaml.load_file(configfile)
Species <- config$report$Species
genename_or_geneid <- config$params$splicing$prepare$genename_or_geneid

#######################
message("读取 geneInfo.tab 文件，获得每个特征的 gene_type")
#######################

if (tolower(Species) == "human"){
  M20FFPE_refgenome <- config$database$Homo$m20$ref_genome
} else if (tolower(Species) == "mouse"){
  M20FFPE_refgenome <- config$database$Mus$m20$ref_genome
} else {
  message("您没有正确设置物种，目前只支持 Human 和 Mouse")
  quit(save = "no", status = 1, runLast = TRUE)
}

geneInfo_file <- paste0(M20FFPE_refgenome, "/STAR_index/geneInfo.tab")
geneInfo <- read.table(geneInfo_file, skip = 1, header = F, sep = "\t", col.names = c("ensembl_id", "gene_name", "gene_type") )



#######################
message("根据gtf文件注释junction对应的基因")
#######################

gtf <- as.data.frame(fread(gtffile), sep="\t", header=FALSE, stringsAsFactors=FALSE)
if (dim(gtf)[2] != 9){
  message("GTF文件没有正确读取，请检查")
  quit(save = "no", status = 1, runLast = TRUE)
}
## MARVEL内部函数调用的时候以以下名字调用，所以需要修改掉默认的名字
names(gtf) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")

## 解析gtf文件
gtf <- gtf[which(gtf$V3=="exon"), ]
if (genename_or_geneid == "gene_name"){
    # Retrieve gene names
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("gene_name", x, value=TRUE))
    . <- gsub("gene_name", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)

    gtf$gene_short_name <- .   
}

if (genename_or_geneid == "gene_id"){
    # Retrieve gene names        
    . <- strsplit(gtf$V9, split=";")
    . <- sapply(., function(x) grep("gene_id", x, value=TRUE))
    . <- gsub("gene_id", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)

    gtf$gene_short_name <- .
}

gtf$V4 <- gtf$V4 - 1
gtf$V5 <- gtf$V5 + 1

logical <- grepl("chr", gtf$V1[1], fixed=TRUE)

# Collapse start position
# Keep unique entries
gtf.small <- gtf[, c("V1", "V5", "gene_short_name")]

if(logical==FALSE) {
    gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V5, sep="")
} else if(logical==TRUE) {
    gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V5, sep="")
}
    
gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
gtf.small <- unique(gtf.small)

# Collapse
gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
. <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)

# Save as new object
gtf.small.collapsed.start <- gtf.small.collapsed
names(gtf.small.collapsed.start) <- paste(names(gtf.small.collapsed.start), ".start", sep="")

# Collapse end position
# Keep unique entries
gtf.small <- gtf[, c("V1", "V4", "gene_short_name")]

if(logical==FALSE) {
    gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V4, sep="")
} else if(logical==TRUE) {
    gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V4, sep="")
}

gtf.small <- gtf.small[,c("chr.pos", "gene_short_name")]
gtf.small <- unique(gtf.small)

# Collapse
gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
. <- by(gtf.small[,"gene_short_name"], gtf.small[,"chr.pos"], function(x) {paste(x, collapse="|")})
gtf.small.collapsed <- data.frame("chr"=as.character(names(.)), "gene_short_name"=as.character(.), stringsAsFactors=FALSE)

# Save as new object
gtf.small.collapsed.end <- gtf.small.collapsed
names(gtf.small.collapsed.end) <- paste(names(gtf.small.collapsed.end), ".end", sep="")

#######################
message("读取 junction 文件")
#######################

junction_data <- read.table(junctionfile, header = F, sep = "\t")
junction_data <- junction_data[,1:3]
names(junction_data) <- c("chr", "start", "end")
junction_data$coord.intron <- paste(junction_data$chr, junction_data$start, junction_data$end, sep=":")

junction_data$chr.start <- paste(junction_data$chr, ":", junction_data$start, sep="")
junction_data$chr.end <- paste(junction_data$chr, ":", junction_data$end, sep="")

# 将解析的 gtf 内容和 junction 合并，仅保留成功注释的 junction
junction_data <- join(junction_data, gtf.small.collapsed.start, by="chr.start", type="left")
junction_data <- join(junction_data, gtf.small.collapsed.end, by="chr.end", type="left")
junction_data <- na.omit(junction_data)

junction_data



