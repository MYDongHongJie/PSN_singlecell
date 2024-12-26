# .libPaths(c("/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/sc/lib/R/library","/PERSONALBIO/work/singlecell/s00/software/R.library/library"))
library(Seurat)
#rds <- "/PERSONALBIO/work/singlecell/s02/single_obj_doing/M02/S202212006/03.label/03.sub/Epi/sub.rds"
suppressPackageStartupMessages( library("argparse") )
parser <- ArgumentParser(description = "scenic  parameters",
                        usage = "%(single_obg)s [global options]" )
parser$add_argument("--rds", type="character", help="input rds file")
parser$add_argument("--outdir", type="character", help="output directory")
parser$add_argument("--subset", type = "character", default = NULL,
             help = "The conditional expression to subset cells used for subtyping.[default: %(default)s]")
parser$add_argument("--species", type = "character", default = "mmu",
             help = "The species for the project .[default: %(default)s]")
parser$add_argument("--mem", type = "character", default = "50G",
             help = "Script runtime memory.[default: %(default)s]")						 
opt = parser$parse_args()


rds = opt$rds
outdir = opt$outdir
species = opt$species
if ( !file.exists(outdir) ){
    dir.create( outdir, recursive = TRUE )
}

single_ob <- readRDS(rds)
if ( !is.null(opt$subset) ){
      df = single_ob@meta.data
      desired_cells= subset(df, eval( parse(text=opt$subset)))
      single_ob = single_ob[, rownames(desired_cells)]
}

#注意矩阵一定要转置，不然会报错
if(length(rownames(single_ob@meta.data))>20000){
    single_ob <- subset(single_ob,downsample=5000)
}
dir.create("temp")
if (single_ob@version >=5){
	write.csv(t(as.matrix(single_ob@assays$RNA$counts)),file = "temp/matrix.csv",quote=F)
}else{
	write.csv(t(as.matrix(single_ob@assays$RNA@counts)),file = "temp/matrix.csv",quote=F)
}

write.table(single_ob@meta.data,'temp/metadata.txt',sep='\t',quote=F)


  
# 构建SBATCH脚本内容
sbatch_script <- c(
	"#!/bin/bash",
	paste("#SBATCH -J pyscenic"),
	paste("#SBATCH -e pyscenic.err"),
	paste("#SBATCH -o pyscenic.out"),
	"#SBATCH -p Batch3,Batch2",
	paste("#SBATCH --mem", opt$mem),
	paste("source /PERSONALBIO/work/singlecell/s02/software/miniconda3/bin/activate CellID && python /PERSONALBIO/work/singlecell/s02/software/script/pyScenic/PyScenic.py --matrix temp/matrix.csv --meta temp/metadata.txt --outdir", outdir, "--species", species)  
)

# 构建sbatch脚本文件名,加一个随机数
#随机生成一个四位数
random_number <- sample(1:9999, 1)
sample_name <- paste0("scenic",random_number)

# 将sbatch脚本内容写入文件
writeLines(sbatch_script, con = paste0(sample_name, ".sh"))


# 将内容写入到文件中
writeLines(sbatch_script, con = paste0(sample_name, ".sh"))
system(glue::glue('sbatch {sample_name}.sh'))
