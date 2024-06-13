#!/usr/bin/env Rscript
rm(list=ls())
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("Matrix") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("scatterplot3d"))
suppressPackageStartupMessages( library("destiny"))
suppressPackageStartupMessages( library("OESingleCell"))

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
        help = "The input exprssion matrix in several possible format."),
    make_option( c("--informat", "-f" ), type = "character", default = "tenx",
        help = "The indication of type of input expression matrix, the possible type can be:
                                tenx:the directory of cellranger count/aggr results with sampleid as its subdirectory.
                                raw_dense: the raw gene expression matrix file, it can be very large.
                                raw_sparse: the raw gene expression matrix in sparse format.
                                seurat: the seurat object from the clustering results."),
    make_option( c("--components", "-t"), type="integer", default=20,
        help="the appropriate number of statistically significant components to use for clustering,
                         which derived from the JackStraw result."),
    make_option( c("--min.cell","-x" ),type="double", default = 0,
        help="the minimium proportion of cell number one gene detected", metavar = "minimium proportion"),
    make_option( c("--batch", "-b" ), type = "logical", default = F,
        help = "Wether to remove batch effect using combat from sva package."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--metadata", "-m" ), type="character", default = NULL,
        help="the sample metadata which must include sample id in this assay design."),
    make_option( c("--cluteringuse", "-c" ), type="character", default = "snn",
        help="the clustering methods used to produce cell clusters.The option can be Kmeans,snn" ),
    make_option( c("--from"), type = "character", default = NULL,
        help="the original cluster id used in calibrate the clustering results."),
    make_option( c("--to"), type = "character", default = NULL,
        help="the adjusted cluster id used in calibrate the clustering results."),
    make_option( c("--resolution","-r"), type = "double", default = 1,
        help = "vaule used to set the resolution of cluster distiguish,
                         use a value above(below)1.0 if you want to obtain a larger(smaller) number of communities."),
    make_option( c("--outformat", "-v"), type = "character", default = "RDS",
        help = "the seurat object saved as R object in RDS format in case of reanalysis."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
        help = "The subset of cluster ids used for subtyping.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$resolution) ){
    print("The clustering resolution value will be set to 0.8 as default.")
    resolution = 0.8
}else{
    resolution = opt$resolution
}

#check the components to use
if ( is.null(opt$components) ){
    print("NO components number is available, the default 10 will be used")
    components_num = 20
}else{
    components_num  = opt$components
}

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

#parse the assay design metadata
if ( !file.exists( opt$metadata) ){
    print_help(opt_parser)
    stop( "Warning:the metadata of this single cell assay is NOT AVAILABLE!")
}else{
    assay_metadata = read.csv(opt$metadata,sep=",",header =T )
}

output_dir = normalizePath(output_dir )
# ####################################################################
#read in the 10X data in different format
# ####################################################################
#======1.1 Setup the Seurat Object (Required)=========================
if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){
        # the input is a seurat object produced by previous analysis
        seurat_ob = readRDS( opt$input )

        if ( !is.null(opt$from) & !is.null(opt$to) ){
            from_ident = unlist( strsplit( opt$from, ",", perl =T ) )
            to_ident = unlist( strsplit(opt$to), ",", perl = T )
            seurat_ob@ident = plyr::mapvalues(x = seurat_ob@ident, from = from_ident, to = to_ident )
        }

        if ( !is.null(opt$which_cells)){
            cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
            seurat_ob = SubsetData(seurat_ob,cells.use = WhichCells(seurat_ob, ident = cluster_list))
        }
    }else{
        if ( opt$informat == "tenx" ){
            # Initialize the Seurat object with the raw (non-normalized data).
            # Keep all genes expressed in >= 3 cells (~0.1% of the data).
            #Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X.
            #A vector or named vector can be given in order to load several data directories of several samples.
            #If a named vector is given, the cell barcode names will be prefixed with the name.
            #the gene-barcode matrix result directory for all samples from 10X
            outs_filtered_path = "outs/filtered_gene_bc_matrices"
            tenx_path = sub("\\/$","",opt$input,perl=T)
            tenx_path = normalizePath(tenx_path)
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

#batch effect remove using Combat from sva package
if ( opt$batch == T){
    suppressPackageStartupMessages(library(sva))
    m = as.matrix(seurat_ob@data)
    com = ComBat(dat=m, batch=seurat_ob@meta.data$batchid,
    prior.plots=FALSE, par.prior=TRUE)
    seurat_ob@data = Matrix(com)
    rm(m)
    rm(com)
}

origin_metadata = seurat_ob@meta.data
#=========dimension reduction and clustering
#the genes in seurat_ob@var.genes are used as input, but can be defined using pc.genes.
#We have typically found that running dimensionality reduction on highly variable genes can improve performance.
seurat_ob <- RunPCA(object = seurat_ob,pc.genes = seurat_ob@var.genes,do.print = F)
# node.scores <- AssessNodes(seurat_ob)
# node.scores = node.scores[order(node.scores$oobe, decreasing = TRUE), ]
##################################################
## ********* USER DEFINED SECTION ***************
##################################################
max_allowed_oobe = 0.10   # any tree branches exceeding this OOBE are collapsed.
## must re-play this section for each merge event, since only one merge per 'play'
# node_to_merge = node.scores$node[node.scores$oobe > max_allowed_oobe]
# if (length(node_to_merge) > 0) {
#   # have at least one node to merge
#   message("merging high OOBE clusters")
#   # do merge for top set:
#   seurat2obj = MergeNode(object=seurat_ob, node.use=node_to_merge[1])
#   ####################################################################################
#   # now redo the few steps we performed just above - in this new post-merge clustering
#
#   # examine the new tree post-merge
#   seurat_ob <- BuildClusterTree(seurat_ob, do.reorder = TRUE, reorder.numeric = TRUE)
#   TSNEPlot(object = seurat_ob, do.label = TRUE) # tSNE plot as a reference to the cluster ids

#   # reexamine classification error in the context of the updated clusters:
#   node.scores <- AssessNodes(seurat_ob)
#   node.scores = node.scores[order(node.scores$oobe, decreasing = TRUE), ]
#   print(head(node.scores))
# } else {
#   message("no clusters with high OOBE to merge.")
# }

#after v2.3.4,seurat use the destiny package as the default diffusion map provider
#so make sure the destiny R package is installed before running this script
message("Beginning Diffusion Map Dimension Reduction")
diffusion_outdir = file.path(output_dir,"diffusion_Dimension_Reduction")
if ( !dir.exists(diffusion_outdir) ){
    dir.create(diffusion_outdir)
}

#dimension reduction using diffusion method
seurat_ob <- RunDiffusion(seurat_ob,genes.use = seurat_ob@var.genes)
message("Beginning to cluster the cell using Diffusion Map reduction results")
seurat_ob <- FindClusters(object = seurat_ob,
                            reduction.type = "pca",
                            dims.use = 1:components_num,
                            resolution = resolution,  #orginal 0.6,Further subdivisions 0.8
                            print.output = 0,
                            save.SNN = TRUE)
#keep the clutering results using the specified resolution to seurat object metadata
diffusion_cluster_result_colname2d = paste(reductmethod,"2D","res",resolution,sep = ".")
# from.id = levels(seurat_ob@ident)
# to.id = as.numeric(from.id)+1
# seurat_ob@ident = plyr::mapvalues(x=seurat_ob@ident,from=from.id,to=to.id)
seurat_ob <- StashIdent(object = seurat_ob, save.name = diffusion_cluster_result_colname2d)
cell_count_diffusion = table(seurat_ob@ident)
diffusion_cell_count_labels = paste(paste(names(cell_count_diffusion),cell_count_diffusion,sep="-")," cells")

# note that you can set do.label=T to help label individual clusters
pp = DMPlot(object = seurat_ob, dim.1 = 1, dim.2 = 2, do.label = T,group.by=diffusion_cluster_result_colname2d,do.return =T)
pp + scale_color_discrete(breaks=levels(seurat_ob@ident),labels = diffusion_cell_count_labels)
ggsave(file.path(diffusion_outdir,paste0("diffusion_groupby_cluster","_resolution",resolution,"_plot.pdf",collapse="")))
ggsave(file.path(diffusion_outdir,paste0("diffusion_groupby_cluster","_resolution",resolution,"_plot.png")),dpi=1000)
if ( opt$groupby == T ){
    DMPlot(object = seurat_ob, dim.1 = 1, dim.2 = 2, do.label = T,group.by="orig.ident")
    ggsave(file.path(diffusion_outdir,paste0("diffusion_groupby_design","_resolution",resolution,"_plot.pdf",collapse="")))
    ggsave(file.path(diffusion_outdir,paste0("diffusion_groupby_design","_resolution",resolution,"_plot.png",collapse="")),dpi=1000)
}
diffusion2d <- as.data.frame(seurat_ob@dr$dm@cell.embeddings)
diffusion2d_data <- merge(diffusion2d, origin_metadata, by=0)

rownames(diffusion2d_data) <- diffusion2d_data$Row.names
diffusion2d_data$Row.names <- NULL
diffusion2d_data$barcode <- row.names(diffusion2d_data)
diffusion2d_data$cluster = seurat_ob@meta.data[,diffusion_cluster_result_colname2d]
## save the annotated diffusion data frame
write.table(diffusion2d_data,
            file.path(diffusion_outdir,paste0("diffusion.2D.resolution",resolution,".output.xls",collapse="")),
            sep="\t", quote=FALSE, row.names=FALSE)

#3D diffusion
seurat_ob <- RunDiffusion(seurat_ob,genes.use = seurat_ob@var.genes, max.dim = 3)
diffusion_cluster_result_colname3d = paste(reductmethod,"3D","res",resolution,sep = ".")
#keep the clutering results using the specified resolution to seurat object metadata
# from.id = levels(seurat_ob@ident)
# to.id = as.numeric(from.id)+1
# seurat_ob@ident = plyr::mapvalues(x=seurat_ob@ident,from=from.id,to=to.id)
seurat_ob <- StashIdent(object = seurat_ob, save.name = diffusion_cluster_result_colname3d)
cell_count_diffusion3d = table(seurat_ob@ident)
diffusion_cell_count_labels3d = paste(paste(names(cell_count_diffusion3d),cell_count_diffusion3d,sep="-")," cells")

DCs = data.frame(DC1=seurat_ob@dr$dm@cell.embeddings[,1],
DC2=seurat_ob@dr$dm@cell.embeddings[,2],
DC3=seurat_ob@dr$dm@cell.embeddings[,3],
cluster = seurat_ob@ident)
cluster_num = length(levels(seurat_ob@ident))
col_from_seurat = hcl(h = seq(15, 375, length = cluster_num+1 ), l = 65, c = 100)[1:cluster_num]
all_col = factor(DCs$cluster,level=levels(DCs$cluster),labels = col_from_seurat)
pdf(file.path(diffusion_outdir,
paste0("3D_diffusion_groupby_cluster","_resolution",resolution,"_plot.pdf")),
width=20,height = 20)
plot3d = with(DCs,scatterplot3d(DC1, DC2, DC3, highlight.3d=F,
mar = c(3,3,1,10)+0.1, cex.axis = 1,cex.symbols = 1.5,
color = all_col, col.grid="lightblue", grid =T,scale.y = 0.8,
box = T,angle = 50, main=NULL, pch=20))
legend("topright",pch = 19,legend =tsne_cell_count_labels3d,col = col_from_seurat,
ncol=1,inset = c(-0.25,0.2),cex = 1.3,xpd = T,horiz = F,bty = "n")
dev.off()
#legend(plot3d$xyz.convert(40,40,60),pch = 16,labels=diffusion_cell_count_labels3d,
#		legend = levels(diffusion$cluster),col = as.character(all_col),ncol = 2)
#save the output
diffusion3d <- as.data.frame(seurat_ob@dr$dm@cell.embeddings)
diffusion3d_data <- merge(diffusion3d, origin_metadata, by=0)

rownames(diffusion3d_data) <- diffusion3d_data$Row.names
diffusion3d_data$Row.names <- NULL
diffusion3d_data$barcode <- row.names(diffusion3d_data)
diffusion3d_data$cluster = seurat_ob@meta.data[,diffusion_cluster_result_colname3d]
write.table(diffusion3d_data,
file.path(diffusion_outdir,paste0("diffusion_resolution",resolution,"_output.xls",collapse="")),
sep="\t", quote=FALSE, row.names=FALSE)
dev.off()
message("Diffusion Map based Dimension Reduction and Clutering finished")

metadata = cbind(rownames(seurat_ob@meta.data), seurat_ob@meta.data)
colnames(metadata) = c("barcode",colnames(seurat_ob@meta.data))
write.table(metadata, file.path(output_dir,paste0("clustering_results.resolution",resolution,".xls", collapse="")),
sep="\t",row.names = F,col.names=T, quote = F)
if (!is.null(opt$outformat) & opt$outformat == "RDS"){
saveRDSMC(seurat_ob,file.path(dirname(output_dir),paste0("seurat_object.clustering_resolution",resolution,".rds", collapse="")))
