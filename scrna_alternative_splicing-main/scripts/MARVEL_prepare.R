# .libPaths("/public/dev_scRNA/zhangchunyuan/tools/R4.2.0_MARVEL/lib/R/library/")
suppressMessages(library(MARVEL))

# Load adjunct packages for selected MARVEL features
# General data processing, plotting
suppressMessages(library(ggnewscale))
suppressMessages(library(ggrepel))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(stringr))
suppressMessages(library(textclean))

# Gene ontology analysis
suppressMessages(library(AnnotationDbi))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))

# ad hoc gene candidate gene analysis
suppressMessages(library(gtools))

# Visualising splice junction location
suppressMessages(library(GenomicRanges))
suppressMessages(library(IRanges))
suppressMessages(library(S4Vectors))
suppressMessages(library(wiggleplotr))

# Load adjunct packages for this tutorial
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(patchwork))
suppressMessages(library(purrr))
suppressMessages(library(scCustomize))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(yaml))
# 从外部文件中导入函数
source("scripts/MARVEL_functions.R")

################################
################################

option_list <- list(
  make_option("--samplefile", type="character", default="config/samples.csv", help="样本文件"),
  make_option("--gtf", type="character", default=NULL, help="gtf文件"),
  make_option("--coordinate", type="character", default=NULL, help="降维聚类文件"),
  make_option("--h5seurat", type="character", default="result/Count_QC/filtered.h5seurat", help="原项目中h5seurat对象，从中读取meta.data表格"),
  make_option("--cluster", type="character", default=NULL, help="聚类结果文件"),
  make_option("--SJ_datapath", type="character", default="result/Splicing/STARSolo", help="可变剪切矩阵结果文件夹"),
  make_option("--Gene_datapath", type="character", default="", help="基因矩阵结果路径"),
  make_option("--scale_factor", type="numeric", default=1e6, help="可变剪切结果文件夹"),
  make_option("--out_MARVEL_data", type="character", default="result/Splicing/MARVEL.RData", help="输出创建好的marvel数据"),
  make_option("--configfile", type="character", default="config/config.yaml", help="config/config.yaml文件，我们直接从文件中读取参数")
)

args <- parse_args(OptionParser(option_list=option_list))

#############
samplefile <- args$samplefile
gtffile <- args$gtf
coordinatefile <- args$coordinate
clusterfile <- args$cluster
h5seurat <- args$h5seurat
scale_factor <- args$scale_factor
Gene_datapath <- args$Gene_datapath
SJ_datapath <- args$SJ_datapath
out_MARVEL_data <- args$out_MARVEL_data
configfile <- args$configfile
config <- yaml.load_file(configfile)

show_protein_coding_only <- config$params$splicing$plot$show_protein_coding_only
######

genename_or_geneid <- config$params$splicing$prepare$genename_or_geneid
use_routine_genearray <- config$params$splicing$prepare$use_routine_genearray

if(genename_or_geneid == "gene_name"){genecolumn = 2}
if(genename_or_geneid == "gene_id"){genecolumn = 1}

Species <- config$report$Species

check_paramaters(genename_or_geneid)
check_paramaters(use_routine_genearray)
check_paramaters(Species)
check_paramaters(show_protein_coding_only)

#######################
# 读取 geneInfo.tab 文件，获得每个特征的 gene_type
#######################

if (tolower(Species) == "human"){
  M20FFPE_refgenome <- config$database$Homo$m20$ref_genome
} else if (tolower(Species) == "mouse"){
  M20FFPE_refgenome <- config$database$Mus$m20$ref_genome
} else {
  message("您没有正确设置物种，目前只支持 Human 和 Mouse")
  quit(save = "no", status = 1, runLast = TRUE)
}

# geneInfo_file <- paste0(M20FFPE_refgenome, "/STAR_index/geneInfo.tab")
# geneInfo <- read.table(geneInfo_file, skip = 1, header = F, sep = "\t", col.names = c("ensembl_id", "gene_name", "gene_type") )
# protein_coding_genes <- geneInfo[geneInfo$gene_type == "protein_coding",]$gene_name


# samplefile <- "config/samples.csv"
# gtffile <- "/data/database/reference/M20FFPE_refgenome/refdata-GRCh38-2020-A-2.0.0/genome.gtf"
# coordinatefile <- "result/Clustering/umap_Dimension_Reduction/umap_Dimension_Reduction_coordination.csv"
# clusterfile <- " result/Clustering/clustering_results.csv"
# scale_factor <- 1e6
# datapath <- "result/Splicing"
# out_MARVEL_data <- "result/Splicing/MARVEL.RData"

################################
################################

sample_meta <- read.table(samplefile, sep=",", header = T)
samplelist <- sample_meta$sampleid
samplelist.paste <- paste0(samplelist, collapse = ", ")
message(paste0("本次分析共有样本 ",length(samplelist)," 个, 分别是：", samplelist.paste))

##################

message("从result/Count_QC/filtered.h5seurat文件中读取meta.data表格")
Seurat_metadata <- LoadH5Seurat(h5seurat)
raw_metadata <- Seurat_metadata@meta.data
raw_metadata$cell.id <- paste(raw_metadata$sampleid, str_split_fixed(raw_metadata$rawbc, pattern = "-", n=2)[,1], sep="_" )
remove(Seurat_metadata)
message("读取meta.data完成")

##################

message("开始读取gtf文件")
gtf <- as.data.frame(fread(gtffile), sep="\t", header=FALSE, stringsAsFactors=FALSE)
if (dim(gtf)[2] != 9){
  message("GTF文件没有正确读取，请检查")
  quit(save = "no", status = 1, runLast = TRUE)
}
## MARVEL内部函数调用的时候以以下名字调用，所以需要修改掉默认的名字
names(gtf) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")

# if(FALSE){
#   gtf.gene <- gtf[gtf$V3 == "gene",]
#   gtf.gene$n <- unlist(lapply(gtf.gene$V9, function(x) {length(unlist(str_split(x, pattern = ";")))}))

#   # 从gtf文件中将 gene_id 和 gene_name 列提取出来，用于后续使用
#   if(all(grepl("gene_id", gtf.gene$V9), grepl("gene_name", gtf.gene$V9))){
    
#     gtf.gene$gene_name <- apply(gtf.gene, 1, function(x) { grep( x =str_split_fixed(x["V9"], pattern = ";", n=as.numeric(x["n"])), pattern = "gene_name", value=T )} ) 
#     gtf.gene$gene_name <- gsub(" ","",gsub("\"","", gsub("gene_name", "", gtf.gene$gene_name)))
#     gtf.gene$gene_id <- apply(gtf.gene, 1, function(x) { grep( x =str_split_fixed(x["V9"], pattern = ";", n=as.numeric(x["n"])), pattern = "gene_id", value=T )} ) 
#     gtf.gene$gene_id <- gsub(" ","",gsub("\"","", gsub("gene_id", "", gtf.gene$gene_id)))
    
#   } else {
#     message("您的gtf文件中缺少 gene_name 或者缺少 gene_id, 请检查")
#     quit(save = "no", status = 1, runLast = TRUE)
#   }
#   ID_table <- gtf.gene[,c("gene_id", "gene_name")]
#   remove(gtf.gene)
#   if( any(c(grepl(x = ID_table$gene_id, pattern = "_"), grepl(x = ID_table$gene_name, pattern = "_"))) ){
#     message("#######################################")
#     message("#######################################")
#     message("")
#     message("Warning: 您的gtf文件中gene_id或者gene_name中有下划线 “_”, 这会导致注释错误")
#     message("")
#     message("#######################################")
#     message("#######################################")
#   }
# }

message("读取gtf文件完成")

##################

message("开始读取聚类坐标文件")
# 不同平台、不同降维方法使用的列名不同，但是内容是相同的，因此我们直接给每列重命名，避免因名字不同造成的报错
coord.table <- read.table(coordinatefile, sep = ",", header = TRUE)
names(coord.table) <- c("Barcode", "x", "y") 
cluster.table <- read.table(clusterfile, sep = ",", header = TRUE)
names(cluster.table) <- c("Barcode", "sampleid", "clusters", "group")
df.coord <- join(coord.table, cluster.table, by="Barcode",  type="left")

df.coord$sample_number <- data.frame(str_split_fixed(df.coord$Barcode, pattern = "-", n = 2))[,2]
df.coord$real_Barcode <- data.frame(str_split_fixed(df.coord$Barcode, pattern = "-", n = 2))[,1]

df.coord$cell.id <- paste(df.coord$sampleid, df.coord$real_Barcode, sep = "_")
df.coord <- subset(df.coord, select = c(cell.id, x, y, sampleid, clusters, group ))
message("读取聚类坐标文件完成")


#############

get_genearray_datapath <- function(){
  if(use_routine_genearray == TRUE){
    data_path <- paste0("result/Splicing/STARSolo/", sample, "/", sample, "_Solo.out/Gene_upstream/filtered/")
  }else if(use_routine_genearray == FALSE){
    if (dir.exists("result/cellranger")){
      data_path <- paste0("result/Splicing/STARSolo/", sample, "/" , sample, "_Solo.out/Gene/filtered/")
    } else if (dir.exists("result/1.STARsolo")) {
      data_path <- paste0("result/Splicing/STARSolo/", sample, "/" , sample, "_Solo.out/GeneFull/filtered/")
    } else {
      message("请检查您的文件路径")
      quit(save = "no", status = 1, runLast = TRUE)
    }
  }else{
    message("请检查config文件中 use_routine_genearray 参数是否正确设置")
  }
  return(data_path)
}


message("开始读取基因矩阵原始数据")
seurat_list <- list()
barcode_list <- list()  # 将样本的barcode保存起来，SJ矩阵利用barcode提取细胞
for (sample in samplelist) {
  data_path <- get_genearray_datapath()
  # message(paste0("从 ", data_path, " 读取基因矩阵"))
  message(paste0("正在处理样本 ",sample))
  tenX_martix <- Read10X(data_path, gene.column = genecolumn)
  barcode_list[[sample]] <- colnames(tenX_martix)
  # message("将矩阵转化为 Seurat 对象")
  seurat_obj <- CreateSeuratObject(tenX_martix)
  seurat_obj$orig.ident <- sample
  seurat_list[[sample]] <- seurat_obj
  remove(tenX_martix)
  remove(seurat_obj)
}
message("样本读取完成，开始合并样本")
seurat_merged_filtered <- Merge_Seurat_List(seurat_list, add.cell.ids = samplelist, merge.data = TRUE)
# 仅保留蛋白编码的基因
# genes_in_seurat <- seurat_merged_filtered@assays$RNA@counts@Dimnames[[1]]
# keep_genes <- intersect(genes_in_seurat, protein_coding_genes)
# seurat_merged_filtered <- subset(seurat_merged_filtered, features =  keep_genes)  # 以基因为单位提取seurat对象
df.gene.count <- as(as.sparse(seurat_merged_filtered@assays$RNA@counts),"dgTMatrix") 
df.gene.count.pheno <- data.frame("cell.id" = row.names(seurat_merged_filtered@meta.data))
df.gene.count.feature <- data.frame("gene_short_name" = unlist(seurat_merged_filtered@assays$RNA@data@Dimnames[1]))
remove(seurat_list)
message("读取基因矩阵原始数据完成")

#################
# 记录被ReadX改变的基因名
#################

featurefile <- paste0(data_path, "/features.tsv.gz")
ID_table <- read.table(gzfile(featurefile), header = F, sep="\t", col.names = c("gene_id", "gene_real_name", "V3"))
ID_table$gene_name <- seurat_merged_filtered@assays$RNA@counts@Dimnames[[1]]
ID_table$V3 <- NULL

################

message("开始读取基因矩阵，并进行未经log转换的标准化")
# message("我们使用gene_id而非gene_name作为差异分析的基本单位，避免一个gene_name对应多个gene_id时可能引入的错误")
seurat_list <- list()
for (sample in samplelist) {
  data_path <- get_genearray_datapath()
  # message(paste0("从 ", data_path, " 读取基因矩阵"))
  message(paste0("正在处理样本 ", sample))
  tenX_martix <- Read10X(data_path, gene.column = genecolumn)
  # message("将矩阵转化为 Seurat 对象")
  seurat_obj <- CreateSeuratObject(tenX_martix)
  # NormalizeData函数共3个标准化的方法，但是只有RC是不经过LOG的
  seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "RC", scale.factor = scale_factor )
  seurat_obj$orig.ident <- sample
  seurat_list[[sample]] <- seurat_obj
  remove(tenX_martix)
  remove(seurat_obj)
}
message("样本读取完成，开始合并样本")
seurat_merged_normalization <- Merge_Seurat_List(seurat_list, add.cell.ids = samplelist, merge.data = TRUE)
# 仅保留蛋白编码的基因
# genes_in_seurat <- seurat_merged_normalization@assays$RNA@counts@Dimnames[[1]]
# keep_genes <- intersect(genes_in_seurat, protein_coding_genes)
# seurat_merged_normalization <- subset(seurat_merged_normalization, features =  keep_genes)   # 以基因为单位提取seurat对象
df.gene.norm <- as(as.sparse(seurat_merged_normalization@assays$RNA@counts), "dgTMatrix")
df.gene.norm.pheno <- data.frame("cell.id" = row.names(seurat_merged_normalization@meta.data), 
                                 "donor.id" = seurat_merged_normalization@meta.data$orig.ident)
df.gene.norm.feature <- data.frame("gene_short_name" = unlist(seurat_merged_normalization@assays$RNA@data@Dimnames[1]))
remove(seurat_list)
message("读取基因矩阵数据并标准化完成")


#################

message("开始读取SJ数据")

seurat_list <- list()
for (sample in samplelist) {
  data_path <- paste0(SJ_datapath, "/", sample, "/", sample,"_Solo.out/SJ/raw/")
  # message(paste0("从 ",data_path, " 读取SJ矩阵数据"))
  tenX_martix <- Read10X(data_path, gene.column = 1)

  # message("读取SJ数据完成，开始注释匹配外显子的junction")
  keep_junctions <- extract_juntions_by_annotate(gtf=gtf,
                                junctions=rownames(tenX_martix),
                                genename_or_geneid=genename_or_geneid
                                )
                                
  if ( length(keep_junctions) ==  0 ){
    message(paste0("样本 ",sample, " 没有 junction 注释到外显子上"))
    quit(save = "no", status = 1, runLast = TRUE)
  }

  message(paste0("样本 ",sample, " 原始SJ矩阵共有 ", length(rownames(tenX_martix))," 个junction，其中注释到外显子的 junction 共 ", length(keep_junctions), " 个，占检测到junction总量的 ", length(keep_junctions)/length(rownames(tenX_martix)) ))
  tenX_martix <- tenX_martix[keep_junctions,]

  barcode_from_gene_martix <- barcode_list[[sample]]
  message(paste0("样本 ", sample," 原始SJ矩阵共有 ", length(colnames(tenX_martix)), " 个细胞，其中高质量的细胞共有 ", length(barcode_from_gene_martix), " 个，占检测到细胞总量的 ",  length(barcode_from_gene_martix)/length(colnames(tenX_martix)) ))
  tenX_martix <- tenX_martix[,barcode_from_gene_martix]
  # junction_size <- dim(tenX_martix)
  # message(paste0("提取 junction 完成，共有junction ", junction_size[1], " 个，共有细胞 " , junction_size[2],"个，将矩阵转化为 Seurat 对象"))
  seurat_obj <- CreateSeuratObject(tenX_martix)  # SJ 矩阵只能使用第一列
  seurat_obj$orig.ident <- sample
  seurat_list[[sample]] <- seurat_obj
  remove(tenX_martix)
  remove(seurat_obj)
}

message("Seurat 对象创建完毕，开始合并样本")
seurat_merged_SJ <- Merge_Seurat_List(seurat_list, add.cell.ids = samplelist, merge.data = TRUE)
message("样本合并完成，开始依据基因-细胞矩阵提取SJ矩阵中高质量细胞")

seurat_merged_SJ <- subset(seurat_merged_SJ, cells=df.coord$cell.id)    # 以细胞为单位提取seurat对象，SJ矩阵没有对应的 gene_name
message("高质量提取细胞完成，将 Seurat 对象转为 marvel 所需数据格式")

df.sj.count <- as(as.sparse(seurat_merged_SJ@assays$RNA@counts), "dgTMatrix")
df.sj.count.pheno <- data.frame("cell.id" = row.names(seurat_merged_SJ@meta.data))
df.sj.count.feature <- data.frame("coord.intron" = unlist(seurat_merged_SJ@assays$RNA@data@Dimnames[1]))
remove(seurat_list)
message("读取SJ数据完成")


##################

message("创建MARVEL对象")
marvel <- CreateMarvelObject.10x(gene.norm.matrix=df.gene.norm,
                                 gene.norm.pheno=df.gene.norm.pheno,  # marvel$sample.metadata
                                 gene.norm.feature=df.gene.norm.feature,
                                 gene.count.matrix=df.gene.count,
                                 gene.count.pheno=df.gene.count.pheno,
                                 gene.count.feature=df.gene.count.feature,
                                 sj.count.matrix=df.sj.count,
                                 sj.count.pheno=df.sj.count.pheno,
                                 sj.count.feature=df.sj.count.feature,
                                 pca=df.coord,
                                 gtf=gtf
)
message("MARVEL对象创建成功")

##################

message("调整sample.metadata表格")
cellnumber.1 <- length(marvel$sample.metadata$cell.id)
marvel$sample.metadata <- plyr::join(marvel$sample.metadata, marvel$pca, type="inner") # 当有多列内容相同时，不使用by，by会导致相同的内容重复出现
marvel$sample.metadata <- plyr::join(marvel$sample.metadata, raw_metadata, type="inner")
cellnumber.2 <- length(marvel$sample.metadata$cell.id)
message(paste0("sample.metadata 包含细胞数目为 ", cellnumber.1, " 与降维位置文件以 inner 合并后的细胞数目为 ", cellnumber.2))
marvel$sample.metadata$clusters <- as.character(marvel$sample.metadata$clusters)
message("调整sample.metadata表格完成")

##################

message("调整gene.metadata表格")
marvel$gene.metadata <- cbind(marvel$gene.metadata, ID_table)
if ( !all(marvel$gene.metadata$gene_name == marvel$gene.metadata$gene_short_name )){
  message("您的基因-细胞矩阵和junction-细胞矩阵对应不上，请核查")
  quit(save = "no", status = 1, runLast = TRUE)
}
message("调整gene.metadata表格完成")

##################

message("开始进行数据分析前准备")
### 基因注释
marvel <- AnnotateGenes.10x.ZCY(MarvelObject=marvel)
# ### 给SJ注释
marvel <- AnnotateSJ.10x.ZCY(MarvelObject=marvel)
### Only splice junctions whose start and end mapped to same gene are retained.
marvel <- ValidateSJ.10x(MarvelObject=marvel)
### Subset CDS genes
if(show_protein_coding_only){
  marvel <- FilterGenes.10x(MarvelObject=marvel, gene.type="protein_coding")
}
### ensures that our data is ready for further downstream analysis including differential expression analysis, functional annotation, and candidate gene analysis.
marvel <- CheckAlignment.10x(MarvelObject=marvel)
message("数据分析前准备完成")

message("保存已经处理好的MARVEL对象")
save(marvel, file=out_MARVEL_data)
message("保存完成")
