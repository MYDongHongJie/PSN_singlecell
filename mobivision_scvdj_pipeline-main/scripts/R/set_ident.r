if (!requireNamespace('argparse', quietly = TRUE)) {
  stop("Please install argparse to parse the command line paramters!")
}
suppressPackageStartupMessages( library("argparse") )
suppressPackageStartupMessages( library("magrittr") )
suppressPackageStartupMessages( library(Seurat) )
suppressPackageStartupMessages( library(SeuratDisk) )

parser = ArgumentParser(description = "single cell sequencing data manipulating toolsets.",
                        usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input_dir", type = "character",
             help = "The output directory of create step.")

opt = parser$parse_args()


prefix = paste0(opt$input_dir,"seurat")
filename = glue::glue("{prefix}.h5seurat")
obj <- LoadH5Seurat(filename)
Idents(obj) <- 1
SaveH5Seurat(obj,filename,overwrite = TRUE, verbose = FALSE)
