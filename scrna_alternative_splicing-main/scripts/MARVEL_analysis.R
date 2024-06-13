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
suppressMessages(library(patchwork))
suppressMessages(library(purrr))
suppressMessages(library(scCustomize))

suppressMessages(library(optparse))

suppressMessages(library(yaml))
suppressMessages(library(cowplot))

##########################################
# 从外部文件中导入函数
source("scripts/MARVEL_functions.R")

############################

option_list <- list(
  make_option("--marvelproject", type="character", default="result/Splicing/MARVEL.RData", help="MARVEL文件"),
  make_option("--diff_type", type="character", default=NULL, help="差异比较的类型，可以是样本之间，也可以是cluster之间，等。目前只能在样本之间计算"),
  make_option("--diff_control", type="character", default=NULL, help="差异比较的control"),
  make_option("--diff_case", type="character", default=NULL, help="差异比较的case"),
  make_option("--outpath", type="character", default="result/Splicing", help="输出文件所在位置" ),
  make_option("--configfile", type="character", default="config/config.yaml", help="config/config.yaml文件，我们直接从文件中读取参数")
)

args <- parse_args(OptionParser(option_list=option_list))

marvelproject <- args$marvelproject
diff_type <- args$diff_type
diff_control <- args$diff_control
diff_case <- args$diff_case
outpath <- args$outpath
configfile <- args$configfile
config <- yaml.load_file(configfile)

########

filter_hard <- config$params$splicing$analysis$filter_hard
filter_gene_in_percent_cell <- config$params$splicing$analysis$filter_gene_in_percent_cell
filter_sj_in_percent_cell <- config$params$splicing$analysis$filter_sj_in_percent_cell
iterations <- config$params$splicing$analysis$iterations
min_gene_expression <- config$params$splicing$analysis$min_gene_expression
genename_or_geneid <- config$params$splicing$prepare$genename_or_geneid

check_paramaters(marvelproject)
check_paramaters(diff_type)
check_paramaters(diff_control)
check_paramaters(diff_case)
check_paramaters(configfile)
check_paramaters(filter_hard)
check_paramaters(filter_gene_in_percent_cell)
check_paramaters(filter_sj_in_percent_cell)
check_paramaters(iterations)
check_paramaters(min_gene_expression)
check_paramaters(genename_or_geneid)

#############################
# 创建文件夹
#############################



diff_group <- paste0(diff_type, "_", black2underline(diff_case), "-vs-", black2underline(diff_control) )
create_folder(paste0(outpath,"/", diff_group ))

#############################
#############################

message("加载MARVEL文件")
load(marvelproject)
message("加载完成")

##########################################

message("开始做可变剪切分析")

##########################################

if ( filter_hard == FALSE ){
  message("以样本为单位统计基因在细胞中的表达情况，并绘制表达密度分布图")
  # 返回数据保存在 marvel$pct.cells.expr$Gene$Data
  marvel <- PlotPctExprCells.Genes.10x.ZCY(MarvelObject=marvel,
                                      diff_type=diff_type,
                                      g1=diff_control,
                                      g2=diff_case,
                                      min.pct.cells=5  # 密度图绘制之前不对数据做任何筛选
                                      )
  result_gene_exp <- plot_density(marvel$pct.cells.expr$Gene$Data, diff_control, diff_case, "Gene")
  max_height_location.gene <- result_gene_exp$number  # 该值会传递给函数 CompareValues.SJ.10x()，用于差异可变剪切基因的筛选
  message(paste0("we observe majority of genes to be expressed in ~", max_height_location.gene, "% of cells"))
  plot.gene <- result_gene_exp$figure
  table(marvel$pct.cells.expr$Gene$Data$cell.group)
  pdfname <- paste0(outpath,"/", diff_group, "/diff_SJ_exp/distribution_gene_expression.pdf")
  ggsave(plot = plot.gene, width = 6, height = 5, filename = pdfname )
  pngname <- paste0(outpath,"/", diff_group, "/diff_SJ_exp/distribution_gene_expression.png")
  ggsave(plot = plot.gene, width = 6, height = 5, filename = pngname )
  message("基因表达情况统计完成")


  message("统计 SJ 表达情况" )
  marvel <- PlotPctExprCells.SJ.10x.ZCY(MarvelObject=marvel,
                                    diff_type=diff_type,
                                    g1=diff_control,
                                    g2=diff_case,
                                    min.pct.cells.genes=5,  # 在百分之5以上细胞中表达的 gene 才保留, 仅用于当前画图，不影响后续差异可变剪切的筛选
                                    min.pct.cells.sj=5,     # 在百分之5以上细胞中表达的 junction 才保留, 仅用于当前画图，不影响后续差异可变剪切的筛选
                                    downsample=FALSE,  # 若要按照比例抽取 SJ 矩阵中的细胞，则先抽取，在统计表达 SJ 细胞百分比
                                    downsample.pct.sj=max_height_location.gene  # 只有 downsample = TRUE 时，该参数才有意义
                                    )
  # min.pct.cells.genes=5 和 min.pct.cells.sj=5 的初始值不能设置太小，若设置为0，则细胞表达密度图的最高值可能位于1或者2附近，但是如果基因仅在1%或者2%的细胞中表达，该基因就没有分析意义了

  result_SJ_exp <- plot_density(marvel$pct.cells.expr$SJ$Data, diff_control, diff_case, "SJ")
  max_height_location.SJ <- result_SJ_exp$number  # 该值会传递给函数 CompareValues.SJ.10x()，用于差异可变剪切基因的筛选
  message(paste0("we observe majority of splice junctions to be expressed in ~", max_height_location.SJ , "% of cells."))
  plot.SJ <- result_SJ_exp$figure
  table(marvel$pct.cells.expr$SJ$Data$cell.group)

  pdfname <- paste0(outpath,"/", diff_group, "/diff_SJ_exp/distribution_SJ_expression.pdf")
  ggsave(plot = plot.SJ, width = 6, height = 5, filename = pdfname )
  pngname <- paste0(outpath,"/", diff_group , "/diff_SJ_exp/distribution_SJ_expression.png")
  ggsave(plot = plot.SJ, width = 6, height = 5, filename = pngname )
  message("SJ 表达情况统计完成" )
} else if ( filter_hard == TRUE & !is.null(filter_gene_in_percent_cell) & !is.null(filter_sj_in_percent_cell) ){
  max_height_location.gene <- filter_gene_in_percent_cell
  max_height_location.SJ <- filter_sj_in_percent_cell
  message(paste0("筛选至少在百分之",max_height_location.gene, "的细胞中表达的基因，筛选至少在百分之", max_height_location.SJ, "的细胞中表达的junction"))
}else{
  message("Error: 您没有正确设置参数 filter_hard、 filter_gene_in_percent_cell、 filter_sj_in_percent_cell")
  quit(save = "no", status = 1, runLast = TRUE)
}

message("找可变剪切，并计算表达量，求PSI值")

marvel <- CompareValues.SJ.10x.ZCY(MarvelObject=marvel,
                              g1=diff_control,
                              g2=diff_case,
                              diff_type=diff_type,
                              min.pct.cells.genes=max_height_location.gene,   # 至少在百分之多少的细胞中表达的gene才保留
                              min.pct.cells.sj=max_height_location.SJ,        # 至少在百分之多少的细胞中表达的junction才保留
                              min.gene.norm=min_gene_expression,              # 该参数从外部传参进来，用于按照筛选样本总体平均水平表达量筛选基因，将过于低表达的基因筛掉。常规单细胞不筛选，我们这里也设置为0
                              n.iterations=iterations,
                              downsample=FALSE,  # 抽取细胞，让两个分组间的细胞数量相等，默认不抽取，也没有设置参数接口
                              show.progress=TRUE
                              )

message("计算差异可变剪切对应基因的表达量")
marvel <- CompareValues.Genes.10x(MarvelObject=marvel,
                                  method="wilcox",
                                  log2.transform=TRUE,
                                  show.progress=TRUE
                                  )

############################
############################

marvelfile <- paste0(outpath,"/", diff_group, "/", "MARVEL.RData")
message(paste0("保存完成差异分析的MARVEL对象到: ",  marvelfile))
save(marvel, file=marvelfile)
message("保存MARVEL对象完成")

totaltable <- marvel$DE$SJ$Table

totaltable$gene_id <- marvel$gene.metadata[match(totaltable$gene_short_name, marvel$gene.metadata$gene_short_name),]$gene_id
totaltable$gene_name <- marvel$gene.metadata[match(totaltable$gene_short_name, marvel$gene.metadata$gene_short_name),]$gene_name

if (all(totaltable$gene_name == totaltable$gene_short_name)){

  totaltable$gene_short_name <- NULL
  totaltable <- dplyr::rename(totaltable, "log2FoldChange.sj(psi)" = "log2fc", 
                                        "p-value.sj" = "pval", 
                                        "log2FoldChange.gene" = "log2fc.gene.norm", 
                                        "p-value.gene" = "pval.gene.norm", 
                                        "q-value.gene" = "pval.adj.gene.norm"
                                        )
  
  totaltable <- dplyr::select(totaltable, "coord.intron","gene_name", "gene_id", everything())  

  tablename <- paste0(outpath, "/", diff_group, "/diff_SJ_exp/", "splicing_junction_differect_expression.total.csv")
  message(paste0("保存差异分析总表到：", tablename))
  write.table(totaltable, file = tablename, quote = F, col.names = T, row.names = F, sep=",")
  message("保存差异分析总表完成")

} else {
  message("基因名对应不上，请检查")
  quit(save = "no", status = 1, runLast = TRUE)
}


message("")
message("")
message("!!! 以下内容为运行过程中产生的 warning 信息 !!!")
message("")
message("")

warnings()

