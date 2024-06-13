# .libPaths("/public/dev_scRNA/zhangchunyuan/tools/R4.2.0_MARVEL/lib/R/library/")
options(warn = 0)  # 当有warning时R会执行非0退出，这会导致snakemake将正确输出的结果删除，设置该参数让warning保持正常退出即可

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
suppressMessages(library(hash))
suppressMessages(library(optparse))
suppressMessages(library(yaml))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))

# 从外部文件中导入函数
source("scripts/MARVEL_functions.R")

##########################################

option_list <- list(
  make_option("--marvelproject", type="character", help="MARVEL文件"),
  make_option("--diff_type", type="character", help="差异比较的类型，可以是样本之间，
                                              也可以是cluster之间，或者细胞类型之间，只要是在前期已经完成的分类即可进行分析"),
  make_option("--diff_control", type="character", help="差异比较的control"),
  make_option("--diff_case", type="character", help="差异比较的case"),
  make_option("--outpath", type="character", default="result/Splicing", help="输出文件所在位置" ),
  make_option("--analysis_type", type="character", help="用于指定分析类型，
                              当前的分析类型有：
                              volcano：绘制火山图
                              diff_SJ_exp：统计各类型差异junction，并绘制占比饼图
                              scatter： 绘制差异junction在每个细胞中表达量的散点图
                              tabulate： 绘制junction差异表达点图
                              structure： 绘制差异junction在现有基因注释下的结构图，用以提示差异junction的位置
                              sashimi：绘制sashimi图，用以展示可变剪切
                              "), # 候选列表：volcano, diff_SJ_exp, scatter, tabulate, structure, outdiffgenelist
  make_option("--configfile", type="character", default="config/config.yaml", help="config/config.yaml文件，我们直接从文件中读取参数"),
  # 中间不允许有空白行
  # 不能解析亚参数
  # 以下参数仅当 --analysis_type sashimi 时起作用
  make_option("--gff", type="character", default="/data/database/reference/M20FFPE_refgenome/refdata-GRCh38-2020-A-2.0.0/genome.gff3", help="绘制sashimi时需要的注释文件。【注意】不是gtf" ),
  make_option("--intron_s", type="integer", default=1, help='绘制sashimi图时，内含子缩减倍数，--intron_s 5 意味着内含子的长度为 100/5 = 20，本流程默认值为 10 ' ),
  make_option("--exon_s", type="integer", default=1, help='绘制sashimi图时，外显子缩减倍数，--exon_s 5 意味着外显子的长度为 100/5 = 20，本流程默认值为 1 '),
  make_option("--min_count", type="integer", default=0, help="绘制sashimi图时，最少reads数量，少于该值时不绘制，默认值为 0"),
  make_option("--figheight", type="integer", default=7, help="图片高度，默认值为 7"),
  make_option("--figwidth", type="integer", default=8, help="图片宽度，默认值为 8")
)

args <- parse_args(OptionParser(option_list = option_list))

# marvel对象
marvelproject <- args$marvelproject
# 差异分组类型，如样本间
diff_type <- args$diff_type
# 指定差异分析的样本
diff_control <- args$diff_control
diff_case <- args$diff_case
# 指定输出文件位置
outpath <- args$outpath
# 指定分析类型
analysis_type <- args$analysis_type
analysis_type <- gsub(" ", "", analysis_type)
analysis_type <- unlist(strsplit(analysis_type, ","))

# 指定config文件
###################
configfile <- args$configfile
config <- yaml.load_file(configfile)
###################

###################
# 当绘制sashimi图的时候需要调用亚参数
###################
if ('sashimi' %in% analysis_type){

  gfffile <- args$gff
  exon_s <- args$exon_s
  intron_s <- args$intron_s
  min_count <- args$min_count
  figheight <- args$figheight
  figwidth <- args$figwidth

  check_paramaters(gfffile)
  check_paramaters(exon_s)
  check_paramaters(intron_s)
  check_paramaters(min_count)
  check_paramaters(figheight)
  check_paramaters(figwidth)

}

warningfile <- config$warningfile
warningfile <- file(warningfile, open="a")

check_paramaters(marvelproject)
check_paramaters(diff_type)
check_paramaters(diff_control)
check_paramaters(diff_case)
check_paramaters(outpath)
check_paramaters(analysis_type)
check_paramaters(configfile)
check_paramaters(warningfile)


diff_group <- paste0(diff_type, "_", black2underline(diff_case), "-vs-", black2underline(diff_control))

if (length(analysis_type) == 0){
  message("Error: 您没有指定分析类型")
  quit(save = "no", status = 1, runLast = TRUE)
}

#############################
# 互斥函数判断
#############################

min_gene_expression <- config$params$splicing$analysis$min_gene_expression
min_gene_expression_name <- paste0("min.log2EP-", min_gene_expression)

pvalue_sj <- config$params$splicing$plot$pvalue_sj

# 将 delta_sj 和 log2fc_sj 设置为互斥参数，若设置错误，则执行非0退出
delta_sj <- config$params$splicing$plot$delta_sj
log2fc_sj <- config$params$splicing$plot$log2fc_sj

if( is.null(delta_sj) & !is.null(log2fc_sj) ){
  delat_or_log2fc <- paste0("log2fc.sj-",log2fc_sj)
} else if( !is.null(delta_sj) & is.null(log2fc_sj) ){
  delat_or_log2fc <- paste0("delta.sj-",delta_sj)
} else {
  message("Error: delta_sj 和 log2fc_sj 二选一")
  quit(save = "no", status = 1, runLast = TRUE)
}

pvalue_gene <- config$params$splicing$plot$pvalue_gene
qvalue_gene <- config$params$splicing$plot$qvalue_gene
log2fc_gene <- config$params$splicing$plot$log2fc_gene

# 将pvalue_gene 和 qvalue_gene 设置为互斥参数
if (is.null(pvalue_gene) & !is.null(qvalue_gene)) {
  pval_or_padj <- paste0("qvalue.gene-", qvalue_gene)
  diff_long_name <- paste0("pvalue.sj-", pvalue_sj, "-", delat_or_log2fc, "-", pval_or_padj, "-log2fc.gene-", log2fc_gene, "-", min_gene_expression_name )
} else if(!is.null(pvalue_gene) & is.null(qvalue_gene)){
  pval_or_padj <- paste0("pvalue.gene-", pvalue_gene)
  diff_long_name <- paste0("pvalue.sj-", pvalue_sj, "-", delat_or_log2fc, "-", pval_or_padj, "-log2fc.gene-", log2fc_gene, "-", min_gene_expression_name )
} else {
  message("Warning: 您没有设置 gene 的筛选标准")
}

psi_or_expr <- config$params$splicing$plot$psi_or_expr
dimension_type <- config$params$splicing$prepare$dimension_type
diff_exp_junction_only <- config$params$splicing$plot$diff_exp_junction_only
rescale_introns <- config$params$splicing$plot$rescale_introns
show_protein_coding_only <- config$params$splicing$plot$show_protein_coding_only
genename_or_geneid <- config$params$splicing$prepare$genename_or_geneid
top_DE_junctions <- config$params$splicing$plot$top_DE_junctions

check_paramaters(min_gene_expression)
check_paramaters(delta_sj)
check_paramaters(log2fc_sj)
check_paramaters(pvalue_gene)
check_paramaters(qvalue_gene)
check_paramaters(log2fc_gene)
check_paramaters(psi_or_expr)
check_paramaters(dimension_type)
check_paramaters(diff_exp_junction_only)
check_paramaters(rescale_introns)
check_paramaters(show_protein_coding_only)
check_paramaters(genename_or_geneid)
check_paramaters(top_DE_junctions)

#############################
# 加载marvel文件
#############################

message("加载MARVEL文件")
load(marvelproject)
message("加载完成")

#############################
# 从config文件中读取参数：是按照自定义的细胞表达过滤标准，还是按照最大细胞表达百分比作为细胞过滤标准
#############################

if (config$params$splicing$analysis$filter_hard == TRUE){
  filter_gene_in_percent_cell <- config$params$splicing$analysis$filter_gene_in_percent_cell
  filter_sj_in_percent_cell <- config$params$splicing$analysis$filter_sj_in_percent_cell
  if( is.null(filter_gene_in_percent_cell) | is.null(filter_sj_in_percent_cell) ){
    message("Error: 您需要指定基因和junction的细胞最低表达百分比： filter_gene_in_percent_cell 和 filter_sj_in_percent_cell")
    quit(save = "no", status = 1, runLast = TRUE)
  }
}else{
  result_gene_exp <- plot_density(marvel$pct.cells.expr$Gene$Data, diff_control, diff_case, "Gene")
  max_height_location.gene <- result_gene_exp$number
  filter_gene_in_percent_cell <- max_height_location.gene
  result_SJ_exp <- plot_density(marvel$pct.cells.expr$SJ$Data, diff_control, diff_case, "SJ")
  max_height_location.SJ <- result_SJ_exp$number
  filter_sj_in_percent_cell <- max_height_location.SJ
}

#########################
# 绘制火山图，火山图的绘制并不涉及到基因的上下调关系
#########################

if ("volcano" %in% analysis_type){
  message("开始绘制火山图")
  path <- paste0(outpath, "/", diff_group, "/volcano")
  create_folder(path)
  volcan_long_name <- paste0("pvalue.sj-", pvalue_sj, "-", delat_or_log2fc, "-", min_gene_expression_name)

  # write_table(data=marvel$DE$SJ$Table, outtablename="test.csv")

  if(!is.null(log2fc_sj)){
    message("绘制差异火山图 - log2fc")
    marvel <- PlotDEValues.SJ.10x.ZCY(MarvelObject=marvel,
                                diff_type=diff_type,
                                g1=diff_control,
                                g2=diff_case,
                                min.pct.cells.gene=filter_gene_in_percent_cell,
                                min.pct.cells.sj=filter_sj_in_percent_cell,
                                pval=pvalue_sj,
                                log2fc=log2fc_sj,
                                delta=NULL,
                                min.gene.norm=min_gene_expression,
                                anno=FALSE
    )
  }

  if(is.null(log2fc_sj) & !is.null(delta_sj)){
    message("绘制差异火山图 - delta")
      marvel <- PlotDEValues.SJ.10x.ZCY(MarvelObject=marvel,
                                diff_type=diff_type,
                                g1=diff_control,
                                g2=diff_case,
                                min.pct.cells.gene=filter_gene_in_percent_cell,
                                min.pct.cells.sj=filter_sj_in_percent_cell,
                                pval=pvalue_sj,
                                delta=delta_sj,
                                min.gene.norm=min_gene_expression,
                                anno=FALSE
    )
  }

  pngname <- paste0(path, "/", "diff_splicing_volcano." , volcan_long_name, ".png")
  ggsave(plot = marvel$DE$SJ$VolcanoPlot$SJ$Plot, filename = pngname, width = 6, height = 5)
  pdfname <- paste0(path, "/", "diff_splicing_volcano." , volcan_long_name, ".pdf")
  ggsave(plot = marvel$DE$SJ$VolcanoPlot$SJ$Plot, filename = pdfname, width = 6, height = 5)
  remove(path)
  message("差异火山图绘制完成")
}

#########################
# 根据config.yaml文件中设置的显著性标准，筛选显著差异的可变剪切，并统计可变剪切的种类
#########################


IsoSwitch <- function(){
  if(!is.null(config$params$splicing$plot$log2fc_sj)){
    message("统计可变剪切的种类 - log2fc")
    marvel <<- IsoSwitch.10x.ZCY(MarvelObject=marvel,
                            pval.sj=pvalue_sj,
                            log2fc.sj=log2fc_sj,
                            pval.gene=pvalue_gene,
                            log2fc.gene=log2fc_gene,
                            min.gene.norm=min_gene_expression,
                            min.pct.cells.gene=filter_gene_in_percent_cell,
                            min.pct.cells.sj=filter_sj_in_percent_cell
                            )
  }

  if(is.null(log2fc_sj) & !is.null(delta_sj)){
    message("统计可变剪切的种类 - delta")
    marvel <<- IsoSwitch.10x.ZCY(MarvelObject=marvel,
                            pval.sj=pvalue_sj,
                            delta.sj=delta_sj,
                            pval.gene=pvalue_gene,
                            log2fc.gene=log2fc_gene,
                            min.gene.norm=min_gene_expression,
                            min.pct.cells.gene=filter_gene_in_percent_cell,
                            min.pct.cells.sj=filter_sj_in_percent_cell
                            )
  }
}

#########################
# 统计每种类型
#########################


if ("diff_SJ_exp" %in% analysis_type){

  message("开始统计差异junction的类型")

  path <- paste0(outpath, "/", diff_group, "/diff_SJ_exp")
  create_folder(path)

  IsoSwitch()

  DEtable <- marvel$SJ.Gene.Cor$Data
  DEtable <- subset(DEtable, select = -c(cor))

  DEtable$gene_id <- marvel$gene.metadata[match(DEtable$gene_short_name, marvel$gene.metadata$gene_short_name),]$gene_id
  DEtable$gene_name <- marvel$gene.metadata[match(DEtable$gene_short_name, marvel$gene.metadata$gene_short_name),]$gene_name
  DEtable$gene_short_name <- NULL

  DEtable <- dplyr::rename(DEtable, "log2FoldChange.sj(psi)" = "log2fc", 
                                          "p-value.sj" = "pval", 
                                          "log2FoldChange.gene" = "log2fc.gene.norm", 
                                          "p-value.gene" = "pval.gene.norm", 
                                          "q-value.gene" = "pval.adj.gene.norm", 
                                          "classification" = "cor.complete",
                                          "Regulation.sj" = "sig.sj",
                                          "Regulation.gene" = "sig.gene"
                                          )

  DEtable <- dplyr::select(DEtable, "coord.intron","gene_name", "gene_id", everything())  

  pngname <- paste0(path, "/", "diff_SJ_exp.", diff_long_name, ".png")
  if (!file.exists(pngname)){
    ggsave(plot = marvel$SJ.Gene.Cor$Proportion$Plot, filename = pngname, width = 6, height = 5)
  }
  pdfname <- paste0(path, "/", "diff_SJ_exp.", diff_long_name, ".pdf")
  if (!file.exists(pdfname)){
    ggsave(plot = marvel$SJ.Gene.Cor$Proportion$Plot, filename = pdfname, width = 6, height = 5)
  }
  tablename = paste0(path, "/", "diff_SJ_exp.", diff_long_name, ".csv")
  if (!file.exists(tablename)){
    write.table(DEtable, file = tablename, quote = F, col.names = T, row.names = F, sep=",")
  }
  
  remove(path)
  remove(DEtable)
  message("统计完成")
}


#######################################
# 以细胞为单位绘制基因和junction表达散点图
#######################################
    

if ( "scatter" %in% analysis_type ) {

  scatter_plot_funtion <- function(MarvelObject, g1, g2, min_gene_expression, dimension_type, diff_type, gene){
    MarvelObject <- MarvelObject
    g1 <- g1
    g2 <- g2
    min_gene_expression <- min_gene_expression
    dimension_type <- dimension_type
    diff_type <- diff_type
    gene <- gene

    ######

    group_plot <- MarvelObject$adhocPlot$PCA$CellGroup  # 需要先在外部执行完毕 PlotValues.PCA.CellGroup.10x.ZCY
    MarvelObject <- PlotValues.PCA.Gene.10x.ZCY(MarvelObject=MarvelObject,
                                                g1=g1,
                                                g2=g2,
                                                gene_short_name=gene,
                                                dimension_type=dimension_type,
                                                diff_type=diff_type,
                                                # color.gradient=c("grey","cyan","green","yellow","red"),
                                                log2.transform=TRUE
                                                )
    gene_plot <- MarvelObject$adhocPlot$PCA$Gene # 需要先在外部执行完毕 PlotValues.PCA.Gene.10x.ZCY

    scatter_junction_plot <- list()
    junctions <- MarvelObject$SJ.Gene.Cor$Data[(MarvelObject$SJ.Gene.Cor$Data$gene_short_name == gene),]$coord.intron # 需要先在外面执行完 IsoSwitch.10x.ZCY
    
    if (psi_or_expr == "psi"){ 
      log2.transform = FALSE
      color.gradient = c("grey","cyan","green","yellow","red")
    }
    if (psi_or_expr == "expr"){
      log2.transform = TRUE
      color.gradient = c("grey90","blue","red")
    }

    for(junction in junctions){
      MarvelObject <- PlotValues.PCA.PSI.10x.ZCY(MarvelObject=MarvelObject,
                                      g1=g1,
                                      g2=g2,
                                      coord.intron=junction,
                                      psi_or_expr=psi_or_expr,
                                      # min.gene.count=min_gene_expression, 表达量过低的细胞不显示
                                      log2.transform=log2.transform,
                                      color.gradient=color.gradient,
                                      dimension_type=dimension_type,
                                      diff_type=diff_type
                                      )
      scatter_junction_plot[[junction]] <- MarvelObject$adhocPlot$PCA$PSI
    }

    junction_figures <- ggarrange(plotlist = scatter_junction_plot, 
      nrow = 1, 
      common.legend = TRUE, 
      legend="right",
      align = "v"
      )

    figures <- ggarrange(group_plot, gene_plot, junction_figures , 
      labels = c("A", "B", "C"),  nrow = 1, align = "v",
      widths = c(1, 1, length(scatter_junction_plot)),
      heights = c(1,1,1)
      ) 
    result <- list()
    result[["figures"]] <- figures
    result[["Count_DE_junction"]] <- length(scatter_junction_plot)
    return(result)
  }


  gave_scatter_figures <- function(gene, figures, Count_DE_junction, path){
    gene <- gene
    figures <- figures
    Count_DE_junction <- Count_DE_junction
    path <- path

    ###

    create_folder(path)
    figure_width <- 10 + 5*Count_DE_junction
    pdfname <- paste0(path,"/", gene, ".pdf")
    ggsave(plot = figures, filename = pdfname, width = figure_width, height = 5.5, limitsize = FALSE, bg = "white")
    pngname <- paste0(path,"/", gene, ".png")
    ggsave(plot = figures, filename = pngname, width = figure_width, height = 5.5, limitsize = FALSE, bg = "white")
  }

  ###################

  message("开始绘制散点图，将表达信息映射到降维结果上")

  IsoSwitch()
  create_folder(paste0(outpath,"/", diff_group, "/scatter"))
  scatter_path <- paste0(outpath, "/", diff_group, "/scatter/", diff_long_name)
  create_folder(scatter_path)

  which_dimension_type(dimension_type)

  message("绘制分组散点图")
  marvel <- PlotValues.PCA.CellGroup.10x.ZCY(MarvelObject=marvel,
                                        diff_type=diff_type,
                                        g1=diff_control,
                                        g2=diff_case,
                                        dimension_type=dimension_type,                                        
                                        legendtitle <- ""
                                        )
  message("绘制分组散点图完成")


  # 返回hash形式存储的候选基因
  candidate_genes <- get_candidate_genes(MarvelObject=marvel, 
                                          top_DE_junctions=top_DE_junctions, 
                                          delta=delta_sj
                                          )

  for (DE_type in keys(candidate_genes)) {
    for (gene in candidate_genes[[DE_type]]) {

      message(paste0("正在绘制类型 ", DE_type, " 的基因 ", gene, " 表达量散点图"))
      DE_type_path <- paste0(scatter_path, "/", DE_type)
      create_folder(DE_type_path)
      figures_result <- scatter_plot_funtion(MarvelObject=marvel,
                                          g1=diff_control,
                                          g2=diff_case,
                                          min_gene_expression=min_gene_expression,
                                          dimension_type=dimension_type,
                                          diff_type=diff_type,
                                          gene=gene
                                          )
      figures <- figures_result$figures
      Count_DE_junction <- figures_result$Count_DE_junction
      gave_scatter_figures(gene=gene, figures=figures, Count_DE_junction=Count_DE_junction, path=DE_type_path)
      
    }
  }
  donefile <- paste0(scatter_path, "/scatter.done" )
  file.create(donefile)
}


#######################################
# 绘制基因表达量的面板图
#######################################


if ("tabulate" %in% analysis_type){

  IsoSwitch()
  
  message("开始绘制基因和junction的表达量面板图")

  path <- paste0(outpath, "/", diff_group, "/tabulate/", diff_long_name)
  create_folder(path)

  # 返回候选基因
  candidate_genes <- get_candidate_genes(MarvelObject=marvel, 
                                          top_DE_junctions=top_DE_junctions, 
                                          delta=delta_sj
                                          )

  gene_SJ_exp_table <- data.frame()

  for (DE_type in keys(candidate_genes)) {
    for (gene in candidate_genes[[DE_type]]) {

      message(paste0("正在绘制类型 ", DE_type, " 的基因 ", gene, " 表达量散点图"))
      create_folder(paste0(path, "/", DE_type))

      marvel <- adhocGene.TabulateExpression.Gene.10x.ZCY(MarvelObject=marvel,
                                                    gene_short_name=gene,
                                                    diff_type=diff_type,
                                                    g1=diff_control,
                                                    g2=diff_case,
                                                    min.pct.cells=filter_gene_in_percent_cell
                                                    )

      gene_exp_plot <- marvel$adhocGene$Expression$Gene$Plot

      # 从config文件中读取参数：是否仅绘制差异表达（目标）的junction，如果是FALSE，则绘制原始数据中包含的所有junction，如果为真，则仅绘制差异表达的junction
      # junction按照在细胞中表达的百分比排序，给junction重命名

      marvel <- adhocGene.TabulateExpression.PSI.10x.ZCY(MarvelObject=marvel,
                                                    min.pct.cells=filter_sj_in_percent_cell,
                                                    diff_exp_junction_only=diff_exp_junction_only
                                                    )
      junction_exp_plot <- marvel$adhocGene$Expression$PSI$Plot

      tabulate_figure <- plot_grid(gene_exp_plot, junction_exp_plot, nrow = 1, rel_widths = c(1, 2), labels = c("A", "B"), align = "h")

      pngname <- paste0(path, "/", DE_type, "/", gene, ".png")
      ggsave(plot = tabulate_figure, filename = pngname, width = 10, height = 5)
      pdfname <- paste0(path, "/", DE_type, "/", gene, ".pdf")
      ggsave(plot = tabulate_figure, filename = pdfname, width = 10, height = 5)

      table_gene_exp <- marvel$adhocGene$Expression$Gene$Table
      table_SJ_exp <- marvel$adhocGene$Expression$PSI$Table
      names(table_gene_exp)[which(colnames(table_gene_exp) == "mean.expr")] <- "mean.expr.gene" 
      names(table_gene_exp)[which(colnames(table_gene_exp) == "pct.cells.expr")] <- "pct.cells.expr.gene"

      . <- plyr::join(table_gene_exp, table_SJ_exp, by="cell.group", type = "full")
      .$classification <- DE_type

      gene_SJ_exp_table <- rbind(gene_SJ_exp_table, .)
      remove(.)

    }
  }
  gene_SJ_exp_table$figure.column <- NULL
  gene_SJ_exp_table <- dplyr::select(gene_SJ_exp_table, "coord.intron","gene_name", "gene_id", "cell.group", everything())
  tablename <- paste0(path, "/", "gene_SJ_expression.csv")
  write.table(x = gene_SJ_exp_table, file = tablename, col.names = T, row.names = F, quote = F, sep = ",")
  remove(path)
  message("基因和junction表达量面板图绘制完成")

}


#######################################
# 以基因为单位输出组合图，包含基因表达量比较、junction表达量比较、、基因结构示意图
#######################################

if("structure" %in% analysis_type){

  IsoSwitch()
  message("开始绘制基因和junction结构图")

  path <- paste0(outpath, "/", diff_group, "/structure/", diff_long_name)
  create_folder(path)

  structure_junction_table <- data.frame()
  structure_transtript_table <- data.frame()

  # 返回候选基因，后续基因以哈希的形式存储，基因类型为键，基因列表为值
  candidate_genes <- get_candidate_genes(MarvelObject=marvel, 
                                          top_DE_junctions=top_DE_junctions, 
                                          delta=delta_sj)

  for (DE_type in keys(candidate_genes)) {
    for (gene in candidate_genes[[DE_type]]) {

      create_folder(paste0(path, "/", DE_type))
      message(paste0("正在绘制类型 ", DE_type, " 的基因 ", gene, " 结构图"))
      
      marvel <- adhocGene.PlotSJPosition.10x.ZCY(MarvelObject=marvel,
                                    gene_short_name=gene,
                                    diff_exp_junction_only=diff_exp_junction_only,
                                    rescale_introns=rescale_introns,
                                    show.protein.coding.only=show_protein_coding_only
                                    )

      grange.exon.list <- marvel$adhocGene$SJPosition$exonfile 
      height <- 2 + length(grange.exon.list) * 0.4

      structure_figure <- marvel$adhocGene$SJPosition$Plot.2
      pngname <- paste0(path, "/", DE_type, "/", gene, ".png")
      ggsave(plot = structure_figure, filename = pngname, width = 10, height = height)
      pdfname <- paste0(path, "/", DE_type, "/", gene, ".pdf")
      ggsave(plot = structure_figure, filename = pdfname, width = 10, height = height)

      temp.structure_junction_table <- marvel$adhocGene$SJPosition$junction_table
      # temp.structure_junction_table$gene_name <- gene
      structure_junction_table <- rbind(temp.structure_junction_table, structure_junction_table)
      remove(temp.structure_junction_table)

      temp.structure_transtript_table <- marvel$adhocGene$SJPosition$structure_table
      # temp.structure_transtript_table$gene_name <- gene
      structure_transtript_table <- rbind(temp.structure_transtript_table, structure_transtript_table)
      remove(temp.structure_transtript_table)
    }
  }
  message("所有结构图绘制完成，输出结构表")

  structure_junction_table.name <- paste0(path, "/structure_junctions.csv")
  message(paste0("保存junction结构表 ", structure_junction_table.name ))
  write.table(x=structure_junction_table, file = structure_junction_table.name, col.names = T, row.names = F, quote = F, sep = ",")
  structure_transtript_table.name <- paste0(path, "/structure_transtripts.csv")
  message(paste0("保存转录本结构表 ", structure_transtript_table.name))
  write.table(x=structure_transtript_table, file = structure_transtript_table.name, col.names = T, row.names = F, quote = F, sep = ",")
  remove(path)
  message("结构图绘制完成")

}

#######################################
# 输出差异基因列表
#######################################

if ("sashimi" %in% analysis_type){
  IsoSwitch()
  # message("输出差异表达的基因列表")
  path <- paste0(outpath, "/", diff_group, "/sashimi")
  create_folder(path)

  # 数据集里面cluster以数字的形式存储，输出时需要将Cluster字符添加进去
  if (diff_type == "clusters"){
    diff_case <- paste0('Cluster', diff_case)
    diff_control <- paste0('Cluster', diff_control)
  }

  # 获得差异基因
  candidate_genes <- get_candidate_genes(MarvelObject=marvel, 
                                          top_DE_junctions=top_DE_junctions, 
                                          delta=delta_sj)

  # 处理gff文件，从gff文件中获取基因的基本信息
  # gff <- as.data.frame(fread(gfffile), sep="\t", header=FALSE, stringsAsFactors=FALSE)  #fread读取gff文件会丢失信息
  gff <- read.table(gfffile, header = F, sep = "\t")

  if (dim(gff)[2] != 9){
    message("GFF文件没有正确读取，请检查")
    quit(save = "no", status = 1, runLast = TRUE)
  }

  names(gff) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
  gff <- subset(gff, gff$V3 == 'gene')

  parsing_gtf <- function(annfile, ... ){
    annfile <- annfile
    #########
    for (subinfor in list(...)) {
      . <- strsplit(gff$V9, split=";")
      . <- sapply(., function(x) grep(subinfor, x, value=TRUE))
      . <- gsub(subinfor, "", .)
      . <- gsub(" ", "", .)
      . <- gsub("\"", "", .)
      . <- gsub("=", "", .)
      annfile[,subinfor] <- .
    }
    return(annfile)
  }

  gff <- parsing_gtf(gff, "gene_name")


  for (DE_type in keys(candidate_genes)) {

    subpath <- paste0(path, "/", DE_type)
    create_folder(subpath)



    for (gene in candidate_genes[[DE_type]]) {

      subgff <- subset(gff, gff$gene_name == gene)

      if(dim(subgff)[1] == 1){

        gene_name <- subgff$gene_name
        
        if(gene_name == gene){
          chr <- subgff$V1
          start <- subgff$V4
          end <- subgff$V5
          strand <- subgff$V7

          #######
          case_bam <- paste0(outpath, '/' , diff_group, '/sashimi/' , diff_group, '.g1.bam')
          control_bam <- paste0(outpath, '/' , diff_group, '/sashimi/' , diff_group, '.g2.bam')
          sashimi_location <- paste0( chr, ':', strand, ':', start, ':', end, ':', gfffile)



          com_setting <- paste0('rmats2sashimiplot ',
                        ' --b1 ', case_bam, 
                        ' --b2 ', control_bam,
                        ' -c ',  sashimi_location, 
                        ' --l1 ', diff_case,
                        ' --l2 ', diff_control,
                        ' -o ' , paste0(subpath, '/', gene),
                        ' --event-type SE ',
                        ' --intron_s ', as.character(intron_s),
                        ' --exon_s ', as.character(exon_s),
                        ' --no-text-background ',
                        ' --min-counts ', as.character(min_count),
                        ' --fig-height ', as.character(figheight), 
                        ' --fig-width ', as.character(figwidth)
                        )

          message(com_setting)
          system(com_setting)

          ########

          end = end + 1
          sashimi_event <- paste0( chr, '_', start, '_', end, '_', strand )
          index_path <- paste0( outpath, '/', diff_group, '/sashimi/', DE_type, '/', gene, '/Sashimi_index' )
          setting_file <- paste0( index_path, '/sashimi_plot_settings.txt' )
          outputdir <- paste0( outpath, '/', diff_group, '/sashimi/', DE_type )

          com_plotting <- paste0('sashimi_plot ',
                        ' --plot-event ', sashimi_event, ' ', index_path, ' ', setting_file,
                        ' --output-dir ' , outputdir,
                        ' --no-posteriors ',
                        ' --plot-title ', gene,
                        ' --plot-label ', gene
          )

          message(com_plotting)
          system(com_plotting)
          unlink(paste0(outputdir, "/", gene))

        } else {
          message(paste0('基因 ', gene, '有问题，请查证'))
        }
      } else {
        message(paste0('基因 ', gene, ' 在gff文件中不唯一'))
      }

    }
  }
  donefile <- paste0(path, "/sashimi.done" )
  message(paste0('生成done文件：',donefile))
  file.create(donefile)
}


#######################################
# 
#######################################



message("")
message("")
message("!!! 以下内容为运行过程中产生的 warning 信息 !!!")
message("")
message("")

warnings()

close(warningfile)

