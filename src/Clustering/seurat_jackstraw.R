#!/usr/bin/env Rscript
## Perform the Seurat JackStraw analysis, which can be useful
## for deciding on the number of PCA components to include in the
## downstrem analysis

rm(list=ls())
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("Matrix") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("scatterplot3d"))
suppressPackageStartupMessages( library("oebio"))
        
#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
                 help = "[Required]The input exprssion matrix in several possible format."),
    make_option( c("--informat", "-f" ), type = "character", default = "tenx",
                 help = "[Required]The indication of type of input expression matrix, the possible type can be:
                        tenx:the directory of cellranger count/aggr results with sampleid as its subdirectory.
                        raw_dense: the raw gene expression matrix file, it can be very large.
                        raw_sparse: the raw gene expression matrix in sparse format.
                        seurat: the seurat object from the clustering results."),
    make_option(c("--numreplicate","-r"), type="integer", default=100,
                help="Number of replicates"),
    make_option(c("--PCs", "-c"), type = "integer", default = 20,
                help="Total Number of PCs to compute and store."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory.", metavar="character"),
    make_option( c("--ident2use", "-t" ), type = "character", default = NULL,
        help = "[OPTIONAL]The column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
        help = "[OPTIONAL]The subset of cluster ids used for subtyping."),
    make_option( c("--metadata", "-m" ), type="character", default = NULL,
        help="[Required]the sample metadata which must include sample id in this assay design.", metavar="character")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
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


if ( is.null(opt$ident2use) ){
    print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
    ident2use = "clusters"
}else{
    ident2use = opt$ident2use
}

#parse the assay design metadata
if ( !file.exists( opt$metadata) ){
    print_help(opt_parser)
    stop( "Warning:the metadata of this single cell assay is NOT AVAILABLE!")
}else{
    assay_metadata = read.csv(opt$metadata,sep=",",header =T )
}
# ####################################################################
#read in the 10X data in different format
# ###################################################################
##======1.1 Setup the Seurat Object (Required)=========================
if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){# the input is a seurat object which may contain more than one sample
        seurat_ob = readRDS( opt$input )

        if ( is.null(seurat_ob@meta.data$clusters) ){
            seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
        }else{
            seurat_ob = SetAllIdent( seurat_ob, id = "clusters")
        }

        #get the subset of cells used for visualization
        if ( !is.null(opt$which_cells)){
            cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
            seurat_ob = SubsetData(seurat_ob,cells.use = FastWhichCells(seurat_ob,
                                    group.by = ident2use, subset.value = cluster_list ),subset.raw = T)
        }
    }else{
        if ( opt$informat == "tenx" ){
            # Initialize the Seurat object with the raw (non-normalized data).
            # Keep all genes expressed in >= 3 cells (~0.1% of the data).
            #Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X.
            #A vector or named vector can be given in order to load several data directories of several samples.
            #If a named vector is given, the cell barcode names will be prefixed with the name.
            tenx_path = sub("\\/$","",opt$input,perl=T)
            tenx_path = normalizePath(tenx_path)
            #the gene-barcode matrix result directory for all samples from 10X
            outs_filtered_path = "outs/filtered_gene_bc_matrices"
            matrix_path = apply( assay_metadata,1,
                            function(samplex) file.path(tenx_path,samplex["sampleid"],outs_filtered_path,samplex["specie"]))
            countMatrixSparse <- Read10X(matrix_path)
        }
        #read data from the count matrix file, row as genes and columns as cells
        #the input matrix can be stored as sparse format or text format
        if ( opt$informat == "raw_dense" ){
            mycountmatrix = read.table(opt$input, header =T, row.name = 1)
            #notice that this may take to much memory
            #transform the density count matrix to the sparse form
            countMatrixSparse = Matrix(as.matrix(mycountmatrix), sparse=T)
            #remove the original matrix to reduce memory usage
            rm(mycountmatrix)
        }
        if ( opt$informat == "raw_sparse" ){
            countMatrixSparse = readMM(opt$input)
            rm(mycountmatrix)
        }

        #add in the metadata of all cells
        #preserve the sample group infomation into the seurat object,the order of the sample
        #is the same as the order in the metadata list
        cellnames = colnames(countMatrixSparse)
        if ( length(grep("\\d_[A-Z]*",cellnames,perl=T) ) < length(cellnames) ){
            #the barcodes of cells in first sample don't prefix with sample index number
            #so we add it for convenince
            firstsample = cellnames[grep("^[A-Z]",cellnames,perl=T)]
            cellnames[1:length(firstsample)] = gsub("^","1_",firstsample,perl=T)
            colnames(countMatrixSparse) = cellnames
        }
        #the barcodes of cells in first sample don't prefix with sample index number
        #so we add it for convenince
        sampleidx =  gsub("_.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
        #integrate the metadata from the assay design
        cell_meta = vector( )
        for ( colidx in colnames(assay_metadata) ){
            cell_meta= cbind(cell_meta, as.vector(assay_metadata[sampleidx, colidx]))
        }
        colnames(cell_meta) = colnames(assay_metadata)
        rownames(cell_meta) = cellnames
        cell_meta = as.data.frame(cell_meta)
        cell_meta$orig.ident = cell_meta$sampleid
        #the minimium cell number one gene is detected
        if ( is.null(opt$min.cell) ){
            min.cell_N = 0
        }else{
            min_cell4gene = opt$min.cell
            min.cell_N = round(min_cell4gene * ncol(countMatrixSparse))
        }
        #construct the seurat object using the meta data above
        seurat_ob <- CreateSeuratObject( raw.data = countMatrixSparse,min.cells=min.cell_N,
                                        names.field = 1, meta.data = cell_meta,
                                        names.delim = "_" )
        seurat_ob@meta.data$orig.ident = as.factor(cell_meta$orig.ident)
        rm(countMatrixSparse)
    }
}

#=====1.3 Normalizing the data========
seurat_ob <- NormalizeData(object = seurat_ob, normalization.method = "LogNormalize",scale.factor = 10000)

#========1.4 find the highly variable genes
seurat_ob = FindVariableGenes(object = seurat_ob, mean.function = ExpMean, dispersion.function = LogVMR,do.plot = F,
                              x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
#  1.5 Seurat Clustering
#=====1.5 Scaling the data and removing unwanted sources of variation=====
mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(x = seurat_ob@data), value = T)
percent.mito <- Matrix::colSums(seurat_ob@raw.data[mito.genes, ])/Matrix::colSums(seurat_ob@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
seurat_ob <- AddMetaData(object = seurat_ob,
                         metadata = percent.mito, col.name = "percent.mito")
seurat_ob <- ScaleData(object = seurat_ob,vars.to.regress = c("nUMI", "percent.mito")) #takes some time


#the genes in seurat_ob@var.genes are used as input, but can be defined using pc.genes. 
#We have typically found that running dimensionality reduction on highly variable genes can improve performance.
seurat_ob <- RunPCA(object = seurat_ob,pc.genes = seurat_ob@var.genes,pcs.compute=opt$PCs, do.print = F)
# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
seurat_ob <- ProjectPCA(object = seurat_ob, do.print = FALSE)

pca_outdir = file.path(output_dir,"PCA_Dimension_Reduction")
if ( !dir.exists(pca_outdir) ){
  dir.create(pca_outdir)
}

PCHeatmap(object = seurat_ob, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
ggsave(file.path(pca_outdir,"seurat_One_PC_PCHeatmap.png"),dpi=1000)

PCHeatmap(object = seurat_ob, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
                            label.columns = FALSE, use.full = FALSE)
ggsave(file.path(pca_outdir,"seurat_many_PC_PCHeatmap.png"),dpi=1000)

#======1.7 Determine statistically significant principal components=====
# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
#It takes long time
nPCs <- min(dim(seurat_ob@dr$pca@cell.embeddings)[2],30)
seurat_ob<- JackStraw(seurat_ob, num.replicate=opt$numreplicate,
               num.pc = nPCs)
JackStrawPlot(object = seurat_ob, PCs = 1:nPCs)
ggsave(file.path(pca_outdir,"JackStrawPlot.pdf"))

PCElbowPlot(object = seurat_ob, num.pc = nPCs)
ggsave(file.path(pca_outdir,"OCElbowPlot.pdf"))
