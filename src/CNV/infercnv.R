#!/usr/bin/env Rscript
# Title     : infercnv.R
# Objective : infer the copu number alteration(CNA) from the tumor and normal single cell sequencing data.
# Created by: hanmin
# Created on: 2019/8/29

rm(list=ls())
suppressWarnings({
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("Matrix") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("ggplot2") )
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages( library("OESingleCell"))
    suppressPackageStartupMessages( library("infercnv"))
    suppressPackageStartupMessages(library("dplyr"))
})


#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i"), type = "character", default = "TRUE",
        help = "the processed seurat object saved as R object in RDS format."),
    make_option( c("--informat", "-f" ), type = "character", default = "seurat",
        help = "The indication of type of input expression matrix, the possible type can be:
                            seurat: the seurat object from the clustering results."),
    make_option( c("--gene_order", "-g"), type = "character",
        help = "[OPTIONAL]data file containing the positions of each gene along each chromosome in the genome.
                The gene_order_file, contains chromosome, start, and stop position for each gene, tab-delimited without header:
                        chr   start  stop
                DDX11L1 chr1  11869  14412
                WASH7P  chr1  14363  29806
                FAM138A chr1  34554  36081
                OR4F5   chr1  69091  70008
                OR4F29  chr1 367640 368634
                "),
    make_option( c("--malignant", "-t"), type = "character",
        help = "the cell type name to be assumed as cancer cells"),
    make_option( c("--celltype", "-l"), type = "character",
        help = "the cell type annotation column name to use in seurat metadata."),
    make_option( c("--mode", "-u"), type = "character", default = "samples",
        help = "Grouping level for image filtering or HMM predictions. Options can be:samples,subclusters,cells.
                samples: fastest, but subclusters is ideal.
                subclusters: detect the subclusters  of tumor"),
    make_option( c("--assay", "-e"), type = "character", default = "RNA",
                 help = "[OPTIONAL] The array result to use in case of multimodal analysis."),
    make_option( c("--clusting2use", "-m"), type = "character", default = "ward.D2",
        help = "the hierarchical clustering methods of cells to use.
                Valid choices are: 'ward.D', 'ward.D2', 'single', 'complete',
                'average', 'mcquitty', 'median', 'centroid'."),
    make_option( c("--refgroup", "-r"), type = "character",
        help = "a comma seperated list containing the classifications of the reference (normal) cells to use for infering cnv.
                If the normal reference cell is supplied from customized source, all the reference cell type will be named as 'normal'."),
    make_option( c("--refexp"), type = "character",
        help = "[OPTIONAL]the reference (normal) cells expression count matrix for inferring cnv."),
    make_option( c("--gtf", "-z"), type = "character",
        help = "the exact gtf file for the alignment for this project."),
    # make_option( c("--chr2exclude", "-x"), type = "character", default = "chrM",
    #     help = "list of chromosomes in the reference genome annotations that should be excluded from analysis."),
    make_option( c("--cutoff", "-c"), type = "double", default = 0.1,
        help = "min average read counts per gene among reference cells.
                cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics."),
    make_option( c("--doHMM"),type="logical", default = FALSE,
        help="wether to using HMM when predicting the CNV for each cell. If no CNV gene prediction needed, Set it False to save time." ),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--pval","-p"),type="double", default = 0.05,
        help="max p-value for defining a significant tumor subcluster." ),
    make_option( c("--down","-d"),type="double",
        help="down sample percentage, down normal cell" ),
    make_option( c("--ncores", "-j" ), type="integer", default = 10,
        help="the number of CPUs used to improve the performace.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir)
    }
}
output_dir = normalizePath(output_dir )


if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){
        # the input is a seurat object produced by previous analysis
        seurat_ob = readRDSMC( opt$input, cores = 10)
        print(table(seurat_ob@meta.data[,opt$celltype]))
        ## down sample
        if ( !is.null(opt$down) ){
            if ( !is.null( opt$refgroup ) ){
                refcells = unlist(strsplit(opt$refgroup, ",", perl = T))
                sampled_cellmeta = seurat_ob@meta.data[which(seurat_ob@meta.data[,opt$celltype] %in% refcells),] %>% tibble::rownames_to_column() %>%
                                         dplyr::group_by( .dots= "clusters" ) %>%
                                         dplyr::sample_frac( size = (1-opt$down),replace = F) %>%
                                         tibble::column_to_rownames()
                subsetcell=setdiff(colnames(seurat_ob),rownames(sampled_cellmeta))
                seurat_ob = seurat_ob[, subsetcell]
            }else if ( !is.null( opt$malignant ) ){
                malignant = unlist(strsplit(opt$malignant, ",", perl = T))
                sampled_cellmeta = seurat_ob@meta.data[which(! seurat_ob@meta.data[,opt$celltype] %in% malignant),] %>% tibble::rownames_to_column() %>%
                                         dplyr::group_by( .dots= "clusters" ) %>%
                                         dplyr::sample_frac( size = (1-opt$down),replace = F) %>%
                                         tibble::column_to_rownames()
                subsetcell=setdiff(colnames(seurat_ob),rownames(sampled_cellmeta))
                seurat_ob = seurat_ob[, subsetcell]
            }
        print("after down sample")
        print(table(seurat_ob@meta.data[,opt$celltype]))
        }
        # if the input seurat object version is less than 3, upgrade it to version 3
        if ( seurat_ob@version < 3 ){
        seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }

        # change the default assay for reduction if necessary
        if ( !is.null( opt$assay) ){
            DefaultAssay(seurat_ob) = opt$assay
        }else{
            DefaultAssay(seurat_ob) = "RNA"
        }
    }
}

# prepare the gene annotation file.
if ( !is.null( opt$gene_order) ){
    gene_order_f = opt$gene_order
}else if ( !is.null( opt$gtf ) ){
    gtf = plyranges::read_gff(opt$gtf)
    gene.chr = gtf %>% plyranges::filter(type == "gene" & gene_name %in% rownames(seurat_ob)) %>%
        as.data.frame() %>%
        dplyr::select(gene_name, seqnames, start, end) %>%
        dplyr::distinct(gene_name, .keep_all=T) %>%
        dplyr::mutate( seqnames = paste0("chr", seqnames))
}else{
    stop("NO gene coordination annotation file is not available!")
}

# if ( !is.null(opt$chr2exclude) ){
#     chr2exclude = unlist(strsplit(opt$chr2exclude, ",", perl =T))
# }else{
#     chr2exclude = "chrM"
# }

# split the whole seurat to equal parts with almost same normal cells and tumor cells for pallalization.
# for the tumor cells, equally seperate them to even parts.
# for the normal cells, we tend to substract all the normal cells out and downsample to the desired number
# for high efficiency.
# then combine the same normal cells with different tumor cells parts.

# make the celltype annotation file for CNV run
#a description of the cells, indicating the cell type classifications.
# The annotations_file, containing the cell name and the cell type classification, tab-delimited without header.
#             V1                   V2
# 1 MGH54_P2_C12 Microglia/Macrophage
# 2 MGH36_P6_F03 Microglia/Macrophage
# 3 MGH53_P4_H08 Microglia/Macrophage
# 4 MGH53_P2_E09 Microglia/Macrophage
# 5 MGH36_P5_E12 Oligodendrocytes (non-malignant)
# 6 MGH54_P2_H07 Oligodendrocytes (non-malignant)
# ...
# 179  93_P9_H03 malignant
# 180 93_P10_D04 malignant
# 181  93_P8_G09 malignant
# 182 93_P10_B10 malignant
# 183  93_P9_C07 malignant
# 184  93_P8_A12 malignant
cellanno = FetchData(seurat_ob, vars = opt$celltype ) %>% tibble::rownames_to_column(var = "cellbarcode")
if ( !is.null( opt$refgroup ) ){
    refcells = unlist(strsplit(opt$refgroup, ",", perl = T))
    count_mat = GetAssayData(seurat_ob, "counts")
}else if ( !is.null( opt$malignant ) ){
    malignant = unlist(strsplit(opt$malignant, ",", perl = T))
    refcells = setdiff( unique(seurat_ob@meta.data[, opt$celltype]), malignant)
    count_mat = GetAssayData(seurat_ob, "counts")
}else{
    print("NO reference normal cells or malignant cells are specified!
           The internal customized normal reference data will be used!")
    refexp = readRDSMC(opt$refexp, cores = 10)
    refcell_anno = data.frame(cellbarcode = colnames(refexp), celltype = "normal")
    cellanno = rbind(cellanno, refcell_anno )
    com.genes = intersect( rownames(seurat_ob), rownames(refexp) )
    count_mat = cbind(GetAssayData(seurat_ob, "counts")[com.genes,], refexp[com.genes,])
}

# =================================================================================
# main workflow
# =================================================================================
tempdir = tempdir()
cnv_celltyping = file.path(tempdir, "cnv_celltype_group.xls")
write.table(cellanno, cnv_celltyping, sep = "\t", col.names = F,row.names = F, quote = F)
gene_order_f = file.path( tempdir, "gene_order_file.xls" )
write.table(gene.chr, gene_order_f, col.names =F, row.names =F, sep = "\t", quote =F )

infercnv_obj = CreateInfercnvObject(raw_counts_matrix= count_mat,
                                    annotations_file= cnv_celltyping,
                                    delim="\t",
                                    gene_order_file= gene_order_f,
                                    ref_group_names=refcells)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= opt$cutoff,
                             analysis_mode= opt$mode, # detect the subclusters  of tumor
                             tumor_subcluster_pval=opt$pval,
                             hclust_method = opt$clusting2use,
                             out_dir= output_dir,
                             num_threads=opt$ncores,
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             no_plot = T,
                             no_prelim_plot = F,
                             HMM=opt$doHMM)

# visualize the cnv value for each cell in heatmap for optimal cluster numbers
pdf( file.path(output_dir, "raw_heatmap.pdf"), width = 18, height = 12 )
ComplexHeatmap::Heatmap(
                     t(as.matrix(infercnv_obj@expr.data)),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_row_names =F,
                     show_column_names = F,
                     name="CNV level",
                     use_raster=TRUE,
                     raster_quality=4 )
dev.off()

