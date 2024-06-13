CreateX <-function (data.dir, assay.meta, source = NULL, assay, cell.delim = "-", 
    aggr = FALSE, add.meta = FALSE, outfile = NULL, use.fragment = FALSE, 
    gtf = NULL, chrom.size = NULL, ...) 
{
    if (!is.null(assay.meta)) {
        assay_metadata <- read.csv(normalizePath(assay.meta), 
            sep = ",", header = T, colClasses = c("character"))
        rownames(assay_metadata) <- assay_metadata[["sampleid"]]
    }
    if (add.meta == TRUE) {
        additional_metap <- file.path(assay_metadata[["addition.anno"]])
    }
    else {
        additional_metap <- NULL
    }
    if (!is.null(gtf)) {
        gtf <- rtracklayer::readGFF(gtf)
        if (!"gene_biotype" %in% colnames(gtf)) {
            gtf <- gtf %>% dplyr::rename(gene_biotype = "gene_type")
        }
        annotations.gr <- GenomicRanges::makeGRangesFromDataFrame(gtf, 
            keep.extra.columns = TRUE)
        if (!is.null(chrom.size)) {
            chrom.size = read.delim(chrom.size, sep = "\t", header = F, 
                row.names = 1)
            chrom.size = chrom.size[annotations.gr@seqinfo@seqnames, 
                1]
            annotations.gr@seqinfo@seqlengths = chrom.size
        }
    }
    else {
        annotations.gr <- NULL
    }
    data.dir <- normalizePath(sub("\\/$", "", data.dir, perl = T))
    mtx_exists <- length(Sys.glob(file.path(data.dir,"*", "outs","per_sample_outs" ,"*","count","sample_filtered_*", "barcodes*"))) > 0
    if (mtx_exists) {
        if (source == "mtx") {
            matrix_path <- unlist(lapply(assay_metadata$sampleid, 
              function(x) unique(dirname(Sys.glob(file.path(data.dir, 
                x, "outs","per_sample_outs" ,x,"count","sample_filtered_*", "barcodes*"))))))
            all_fragments <- NULL
            readf <- ReadMTX    
    }else {
        stop("NO filtered matrix FOUND! Check your input directory!")  
    }
    matrix_path <- unique(matrix_path)
    options(future.globals.maxSize = Inf)
    future::plan("multicore", workers = length(matrix_path))
    seu_list <- future.apply::future_lapply(seq_along(matrix_path), 
        function(idx) {
            seu_ob = readf(matrix_path[idx], assay = assay, 
              add.meta = additional_metap[idx], fragmentx = all_fragments[idx], 
              granges = annotations.gr,...)
            seu_ob[["rawbc"]] = paste(gsub("-\\d+$", "", 
              colnames(seu_ob)), idx, sep = "-")
            seu_ob <- Seurat::RenameCells(seu_ob, new.names = paste(idx, 
              colnames(seu_ob), sep = "-"))
            if (tolower(assay) == "spatial") 
              names(seu_ob@images) = assay_metadata$sampleid[idx]
            seu_ob
        }, future.seed = 2020)
    if (length(seu_list) > 1) {
        seurat_ob <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
    }
    else {
        seurat_ob <- seu_list[[1]]
    }
    seurat_ob[["orig.ident"]] = seurat_ob[["rawbc"]]
    seurat_ob <- UpdataCellMeta(seurat_ob, metadata = assay_metadata, 
        cell.delim = cell.delim)
    }
    if (!is.null(outfile)) {
        SaveX(seurat_ob, output = outfile, outformat = "h5seurat")
    }
    else {
        return(seurat_ob)
    }
}
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("OESingleCell"))
suppressPackageStartupMessages( library("argparse") )
suppressPackageStartupMessages( library("magrittr") )
parser = ArgumentParser(description = "single cell sequencing data manipulating toolsets.",
                        usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input", type = "character", default = NULL,
     help = "cellranger path ")
parser$add_argument("-f", "--informat", type = "character", default = NULL,
             help = "The format of data object, the possible choices can be:h5seurat,(seurat)rds,(sce)rds, loom.[default: %(default)s]")
parser$add_argument("-o", "--output", type = "character", default = "./",
             help = "the output directory of results."  )
parser$add_argument("-d", "--outformat", type = "character", default = "h5seurat",
             help = "the output format of data object, possible choices:h5seurat,seurat,anndata,sce,CellDataSet(monocle2)[default: %(default)s]")
parser$add_argument("--update", default= "TRUE", type="character",
             help="whether update the data in the object on disk using the newly produced results. Set this to FALSE when you use this script for subclustering! Only availiable for h5seurat input.[default TRUE]")
parser$add_argument("-s", "--source", type = "character", default = "mtx",
     help = "The source of data, possible choices:h5, mtx, xsv(including BD,tsv and csv), Visium, slideseq.[default: %(default)s]")
parser$add_argument("-m", "--metadata", type = "character",
     help = "the sample metadata which must include sample id in this assay design.")
parser$add_argument("-x", "--metrics", type = "character", default = "percent.mito",
     help = "the additional QC metrics list to calculate in advance.[default %(default)s]")
parser$add_argument("--gcolumn", type = "integer", default = 2,
     help = "Specify which column of genes.tsv or features.tsv to use for feature names.[default %(default)s]")
parser$add_argument("--aggr", type = "character", default = "FALSE",
     help = "[OPTIONAL] wether the results derived from cellranger aggr.[default %(default)s]")
parser$add_argument("--feature.meta", type = "character", default = NULL,
     help = "[OPTIONAL]the reference genome annotaion file from 10X.")
parser$add_argument("--chrom.size", type = "character", default = NULL,
     help = "[OPTIONAL]length of each chromosome of the genome file from 10X.")
parser$add_argument("--cell.meta", type = "character", default = "FALSE",
     help = "[OPTIONAL] logical indication wether there is cell annotation file for, 'singlecell.csv' for in ATAC-seq for example.[default: FALSE]")
parser$add_argument("--assay", type = "character", default = NULL,
             help = "the main assay in data object to use. When it comes to multimodal assay, this is the assay used to initialize the object, all the other assays will merged into it.") 
parser$add_argument("--fragment", type = "character", default = "FALSE",
     help = "[OPTIONAL] use the fragment annotion file 'fragments.tsv.gz' in outs directory from cellranger-atac count results for each sample, Only valid for scATAC-seq assay.[default FALSE]")
parser$add_argument("--transpose", type = "character", default = "FALSE",
     help = "[OPTIONAL]wether to transpose the count matrix when comes to the cell-feature matrix.[default %(default)s]")
parser$add_argument("--gset", type = "character", default = NULL,
     help = "[OPTIONAL]the gmt formated customized genes set file for QC metrics' calculation.")
parser$add_argument("--prefix", type = "character", default = "seurat",
             help = "the prefix of output file without file extension.[default %(default)s]")
opt = parser$parse_args()
data_ob <- CreateX(data.dir = opt$input,
                                  assay.meta = opt$metadata,
                                  source = opt$source, assay = opt$assay,
                                  aggr = as.logical(opt$aggr),
                                  add.meta = as.logical(opt$cell.meta),
                                  use.fragment = as.logical(opt$fragment),
                                  gene.column = opt$gcolumn,
                                  gtf = opt$feature.meta,
                                  chrom.size = opt$chrom.size)
    
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir,recursive = T)
    }
}
output_dir = normalizePath(output_dir ) 

if ( !is.null(opt$gset) ){
  # gmt = GSEABase::getGmt(con=opt$gset)
  gmt_list = unlist(strsplit(opt$gset,",",perl = T))
  #gset_list = GSEABase::geneIds(GSEABase::getGmt(con=opt$gset))
  gset_list <- lapply(gmt_list, function(gmtfile){
                      gset_list <- GSEABase::geneIds(GSEABase::getGmt(con=gmtfile))  
                      return(gset_list)
  })
  for (i in c(1:length(gset_list))){
      data_ob = OESingleCell::CalculateMetrics(data_ob, metrics = NULL, gset = gset_list[[i]])
  }
} else if (!is.null(opt$metrics)) {
  # metrics = c("percent.mito", "CC.difference"),
  metrics = unlist( strsplit(opt$metrics, ","))
  data_ob = OESingleCell::CalculateMetrics(data_ob, metrics = metrics)
}

if ( !is.null(opt$feature.meta) ){
  data_ob = OESingleCell::addFeatureAnno( data_ob,
                                          assay = opt$assay,
                                          gtf = opt$feature.meta)
}

# Now the h5seurat format is used as the standard data object in disk or better interpre
if ( tolower(opt$outformat) == "h5seurat" ){
  prefix = file.path(output_dir, opt$prefix)
  SeuratDisk::SaveH5Seurat(data_ob, filename = glue::glue("{prefix}.h5seurat") ,
                            overwrite = TRUE, verbose = FALSE)
}else{
  OESingleCell::SaveX(data_ob, output = output_dir,
                      outformat = opt$outformat, prefix = opt$prefix, update = FALSE )
}


