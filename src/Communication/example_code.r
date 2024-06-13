#!/usr/bin/env Rscript
# this is the script used to infer the cell-cell communication
# using the high expressed genes or differential expressed genes
# in each cell. The key point is to find the Ligand-Receptor pairs
# in the selected genes.

suppressPackageStartupMessages( library("iTALK") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("Seurat") )

#=command line parameters setting=============================
optioption_list = list(
        make_option( c("--input", "-i" ), type = "character",
             help = "The input exprssion matrix in several possible format."),
        make_option( c("--informat", "-f" ), type = "character", default = "seurat",
             help = "The indication of type of input expression matrix, the possible type can be:
                    sce: the SingleCellExperiment object with celltype annotation.,
                    seurat: the seurat object with celltype annotation."),
        make_option( c("--customized_LR", "-l"), type = "character", default = NULL,
                      help = "the additional customized Ligand-Receptor pairs annotation to use
                        for the inference. If NULl, the default database in the package will be used."),
        make_option( c("--topn", "-n" ), type = "integer", default = 100,
                     help = "the top N highly expressed candidate genes for Ligand-Receptor pairs."),
        make_option( c("--celltype2use", "-c"), type="character", default="celltype",
            help="the cell type annotation column in the cell annotation meta data."),
        make_option( c("--npairs2plot", "-m"), type = "integer", default = 50,
                     help = "the number of Ligand-Receptor pairs to plot on the circos plot."),
        make_option( c("--diffexp", "-e" ), type = "character", default = NULL,
                     help = "[OPTIONAL]The differential expressed genes list 
                            user specified as the canididate Ligand-Receptor pairs."),
        make_option( c("--data2use", "-d"), type="character", default="raw",
            help="the expression data type used as input, options can be raw or normalized."),
        make_option( c("--genetype", "-g"), type="character", default="hvg",
            help="the genes used as candidate Ligand-Receptor pairs, options can be: 
                  highly variable genes(hvg) or all genes(all)."),
        make_option( c("--output","-o"),type="character", default = "./CellCommunication",
            help="the output directory of Clustering results." )
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

npairs = opt$npairs2plot
#=================================================================================
# read in the expression matrix 
#=================================================================================
if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){
        # the input is a seurat object produced by previous analysis
        seurat_ob = readRDS( opt$input )
        # if the input seurat object version is less than 3, upgrade it to version 3
        seurat_ob = UpdateSeuratObject(seurat_ob)
    }else if ( opt$informat == "sce" ){
      # the input is a SingleCellExperiment object
        sc_object = readRDS( opt$input )
        seurat_ob = as.SingleCellExperiment(sc_object)
    }
    # read the data, it should be a cell(row)-gene(column) expression matrix in data.frame format
    # the expression matrix should integrate the cell type annotation stored in the matrix as 
    # cell_type.
    if ( !is.null( opt$diffexp) ){
      genes_in = read.table( opt$genes, sep = "\t", header = T)
      candidate_genes = genes_in$gene
    }else{
      if ( opt$genetype == "hvg" ){ # high variable genes prefered
          candidate_genes = VariableFeatures(seurat_ob)
      }else{ # all genes
          candidate_genes = rownames(seurat_ob)
      }
    }
  
    if ( opt$data2use == "raw" ){
        data = as.data.frame(t(as.matrix(GetAssayData(seurat_ob, slot = "counts")[candidate_genes,])))
    }else{
        data = as.data.frame(t(as.matrix(GetAssayData(seurat_ob, slot = "data")[candidate_genes,])))
    }
    data$cell_type = tryCatch({as.factor(FetchData(seurat_ob, vars = opt$celltype2use)[,1]) },
          error = function(err){
            warnings( "NO cell type annotation Found in the cell metadata slot,
                      the clusters annotation will be used instead.")
            seurat_ob = StashIdent(seurat_ob, save.name = "clusters" )
            clusters_ann = as.factor(FetchData(seurat_ob, vars = "clusters" )[,1])
            return( clusters_ann)
          })
}

# find the candidata genes for Ligand-Receptor pairs from the 
# the highly expressed genes or differential expressed genes.
highly_exprs_genes<-rawParse(data,top_genes=opt$topn,stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(data$cell_type))
res<-NULL
for(comm_type in comm_list){
    # this function loads the highly expressed genes or differential expressed genes
    # as a dataframe. Significant interactions are found through mapping these
    # genes to our ligand-receptor database.
    res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
    #plot by ligand category overall network plot
    pdf(file.path(output_dir,paste0("cell-cell communication netwotk of genes in ", comm_typem,".pdf")))
    NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    dev.off()
    png(file.path(output_dir,paste0("cell-cell communication netwotk of genes in ", comm_type,".png")))
    NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    dev.off()
    #top N ligand-receptor pairs
    pdf(file.path(output_dir,paste0("cell-cell communication circos plot of genes in ", comm_typem,".pdf")))
    LRPlot(res_cat[1:20,],datatype='mean count',
           cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],
           link.arr.width=res_cat$cell_to_mean_exprs[1:20])
    title(comm_type)
    dev.off()
    png(file.path(output_dir,paste0("cell-cell communication circos plot of genes in ", comm_typem,".png")))
    LRPlot(res_cat[1:npairs,],datatype='mean count',
           cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:npairs],
           link.arr.width=res_cat$cell_to_mean_exprs[1:npairs])
    title(comm_type)
    dev.off()
    res<-rbind(res,res_cat)
}
write.table(res, file = file.path(output_dir,"all Ligand-Receptor cell communication results.xls"), sep = "\t", col.names = 1)
# plot for all Ligand-Receptor pairs
pdf(file.path(output_dir,"cell-cell communication circos plot of genes in all types.pdf"))
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
dev.off()
png(file.path(output_dir,"cell-cell communication circos plot of genes in all types.png"))
LRPlot(res[1:npairs,],datatype='mean count',
       cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:npairs],
       link.arr.width=res$cell_to_mean_exprs[1:npairs])
dev.off()