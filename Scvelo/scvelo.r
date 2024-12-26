#!/usr/bin/env Rscript
# Title     : scvelo.r
# Objective : scvelo
# Created by: donghongjie
# Created on: 2024/5/14
library("optparse")

option_list <- list(
	make_option(c("-i", "--rds"), help="input rds "),
	#make_option(c("-i", "--idents"), help=" the colnames in sce@meta.data "),
	make_option(c("-l", "--loom"), help="the path of loom file"),
	make_option(c("-g", "--groupby"), help="The feature of metadata,split by ,"),
	make_option(c("-r", "--reduction"), help="TSNE or UMAP",default="umap"),
	make_option(c("--subset"), help="subsetcell,eg:celltype %in% c('CD4','CD8')"),
	make_option(c("-o", "--output"),help='output path',default="./")
)

opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#wkdir = getwd()
temp_dir = file.path(opt$output,'scvelo_temp')
if (!file.exists(temp_dir)){
	dir.create(temp_dir,recursive = TRUE)
}
if (!file.exists(file.path(temp_dir,'data.h5Seurat'))){
	library(SeuratDisk)
	library(Seurat)
	data_ob = readRDS(opt$rds)
	if ("SCT" %in% names(data_ob@assays)){
		for (redution in names(data_ob@reductions)){
			temp_sct = data_ob[[redution]]
			data_ob[[redution]] =NULL
			data_ob[[redution]] <- temp_sct
			data_ob[[redution]]@assay.used <- "RNA"
		}
	}
	data_ob[["RNA2"]] <- as(object = data_ob[["RNA"]], Class = "Assay")
	data_ob[["RNA"]] = data_ob[["RNA2"]]
	data_ob[["RNA2"]] =NULL
	for (group in unlist(strsplit(opt$groupby, ","))){
			data_ob@meta.data[,group] <- as.character(data_ob@meta.data[,group])
	}
	for (character_meta in c("sample","group")){
		if (character_meta %in% names(data_ob@meta.data)){
			data_ob@meta.data[,character_meta] <- as.character(data_ob@meta.data[,character_meta])
		}
	}
	if ( !is.null(opt$subset) ){
      df = data_ob@meta.data
			print(opt$subset)
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      data_ob = data_ob[, rownames(desired_cells)]
  }
	SaveH5Seurat(data_ob, filename = file.path(temp_dir,'data.h5Seurat'),overwrite=FALSE)
	Convert(file.path(temp_dir,'data.h5Seurat'), dest = "h5ad")

}


pyscvelo = '/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/Scvelo/scvelocity.py'

#cmd = glue::glue("python {pyscvelo} --input scvelo_temp/data.h5Seurat --loom_dir {opt$loom} --metadata scvelo_temp/meta_temp.xls --output {opt$output}  --groupby {opt$groupby} --basis {opt$reduction}")
cmd = glue::glue("source /PERSONALBIO/work/singlecell/s04/Test/donghongjie/Miniconda/bin/activate scvelo &&  python {pyscvelo} --input {temp_dir}/data.h5ad --loom_dir {opt$loom}  --output {opt$output}  --groupby {opt$groupby} --basis {opt$reduction}")
system(cmd)
if (file.exists(file.path(opt$output,'Scvelocity结果说明.docx'))){
	print('删除中间文件')
	system(glue::glue('rm -rf  {opt$loom}/scvelo_temp'))
}