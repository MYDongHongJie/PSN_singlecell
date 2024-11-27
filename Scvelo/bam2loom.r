#!/usr/bin/env Rscript
# Title     : bam2loom.r
# Objective : creat loom file
# Created by: donghongjie
# Created on: 2024/5/14
library("optparse")
library(jsonlite)

option_list <- list(
	make_option(c("-i", "--input"), help="cellranger path "),
	#make_option(c("-i", "--idents"), help=" the colnames in sce@meta.data "),
	make_option(c("-s", "--species"), help="species"),
	make_option(c("-g", "--gtf"), help="if is.null species ,please provide gtf path")
)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

folder_path = opt$input
species = opt$species
# 列出文件夹中的所有文件和文件夹
contents <- list.files(folder_path, full.names = TRUE, recursive = FALSE)

# 定义一个函数来提取子文件夹名称和父文件夹名称
extract_subfolder_and_parent <- function(path) {
  # 使用正则表达式或字符串操作函数提取子文件夹名称和父文件夹名称
  result <- gsub("(.*/)([^/]+/outs/possorted_genome_bam.bam).*", "\\1\\2", path)
  return(result)
}

# 应用函数提取子文件夹名称和父文件夹名称
subfolders_and_parents <- sapply(contents, extract_subfolder_and_parent)

# 去重
subfolders_and_parents <- unique(subfolders_and_parents)

velocy='/PERSONALBIO/work/singlecell/s04/Test/donghongjie/Miniconda/envs/scvelo/bin/velocyto'
json_dic<-read_json("/PERSONALBIO/work/singlecell/s00/software/script/1.source/ref.json")

if (!is.null(species)){
	dir =  json_dic[[species]]['dir']
	gtffile = file.path(dir,'genes/genes.gtf')
	if (!file.exists(gtffile)){
		print('该物种的genes.gtf路径不对或者注释文件没有解压')
	}
}else {
	 gtffile = opt$gtf
}


# 循环处理每个样本
for (i in subfolders_and_parents) {
  sample_name <- tail(strsplit(i, '/')[[1]], 1)
  
  # 构建SBATCH脚本内容
  sbatch_script <- c(
    "#!/bin/bash",
    paste("#SBATCH -J", sample_name),
    paste("#SBATCH -e", paste0(sample_name, ".err")),
    paste("#SBATCH -o", paste0(sample_name, ".out")),
    "#SBATCH -p Batch3,Batch2",
		"#SBATCH --mem 60G",
    paste(velocy, "run10x", i, gtffile)
  )
  
  # 将内容写入到文件中
  writeLines(sbatch_script, con = paste0(sample_name, ".sh"))
	system(glue::glue('sbatch {sample_name}.sh'))
}

