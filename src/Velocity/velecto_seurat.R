#!/usr/bin/env Rscript
.libPaths("/home/lipeng/R/x86_64-conda_cos6-linux-gnu-library/3.6")

#.libPaths("/data/software/conda_envs/seurat3.1/lib/R/library")
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
library("optparse")

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
        help = "The input exprssion matrix in seurat object."),
    make_option( c("--informat", "-f" ), type = "character", default = "seurat",
        help = "The indication of type of input expression matrix, the possible type can be:
               seurat or sce."),
    make_option( c("--loom", "-l" ), type = "character", 
        help = "The comma seperated loom file list from run velocyto."),
    make_option( c("--output","-o"),type = "character", default = "./",
        help = "the output directory of Clustering results." ),
    make_option( c("--metadata", "-m" ), type = "character", default = NULL,
        help = "[optional]the sample metadata which must include sample id in this assay design."),
    make_option( c("--pointsize", "-s"), type = "double", default = 0.5,
        help = "[optional]the point size in the plot."),
    make_option( c("--reduct1"),type = "character",default = "pca",
        help = "[optional] the (primary) reduction methods used, here is to calculate cell distance."),
    make_option( c("--reduct.vis","-r"),type = "character",default = "tsne",
        help = "[optional]the results of computed reduction methods used as input for secondary reduction.
             choice can be ica,cca,pca,rpca,mnn,tsne,umap etc."),
    make_option( c("--ident2use", "-q" ), type = "character", default = NULL,
        help = "[optional]the column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--assay" ), type = "character", default = "RNA",
        help = "[OPTIONAL]the assay used to calulation in case of multimodal data."),
    make_option( c("--groupby", "-g"),type = "character",default = "clusters",
        help = "[optional] which cloumn for plot color in meta.data ,like clusters,groups et"),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
        help = "[optional]the subset of cluster ids used for subtyping."),
    make_option(c("--downsample", "-e"),type = "character", default = "30000",
        help = "the downsample number of cells "
      ),
    make_option( c("--replot" ), type = "logical", default = FALSE,
        help = "[optional]replot.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
suppressWarnings({
  #suppressPackageStartupMessages(library("optparse"))
  suppressPackageStartupMessages(library("velocyto.R"))
  suppressPackageStartupMessages(library("Seurat"))
  suppressPackageStartupMessages(library("SeuratWrappers"))
  suppressPackageStartupMessages(library("OESingleCell"))
})

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
        seurat_ob = readRDS( opt$input )
        # if the input seurat object version is less than 3, upgrade it to version 3
        if ( seurat_ob@version < 3 ){ seurat_ob = UpdateSeuratObject(seurat_ob) }
    }else if ( opt$informat == "sce" ){
      # the input is a SingleCellExperiment object
        sc_object = readRDS( opt$input )
        seurat_ob = as.SingleCellExperiment(sc_object)
    }
  
    # change the default assay for reduction if necessary
    if ( !is.null( opt$assay) ){
        DefaultAssay(seurat_ob) = opt$assay
    }else{
        DefaultAssay(seurat_ob) = "RNA"
    }
}

if (is.null(opt$downsample)) {
  downsample <- 30000
} else {
  downsample <- opt$downsample
}
if ( !opt$reduct1 %in% names(seurat_ob@reductions) ) stop(paste0("Primary reduction ",opt$reduct1," is not present in the imput file." ))
if ( !opt$reduct.vis %in% names(seurat_ob@reductions) ) stop(paste0("Secondary reduction ",opt$reduct.use," is not present in the imput file." ))

#update the metedata in the seurat_ob@meta.data with new additional sample metadata
if ( !is.null(opt$metadata) ){
    additional_metadata = read.csv(opt$metadata,sep=",",header =T )
    rownames(additional_metadata) = additional_metadata$sampleid
    cellnames = Cells(seurat_ob)
    eampleidx =  gsub("_|-[ATGC]{16,}.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
    #integrate the additional metadata from the assay design
    additional_cell_meta = vector( )
    for ( colidx in colnames(additional_metadata) ){
        additional_cell_meta = cbind(additional_cell_meta, as.vector(additional_metadata[sampleidx, colidx]))
    }
    colnames(additional_cell_meta) = colnames(additional_metadata)
    rownames(additional_cell_meta) = cellnames
    additional_cell_meta = as.data.frame(additional_cell_meta)
    seurat_ob = AddMetaData( seurat_ob, metadata = additional_cell_meta)
}

#get the subset of cells used for visualization if necessay
if ( !is.null(opt$which_cells)){
    if ( is.null(opt$ident2use ) ){
        print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
        ident2use = "clusters"
    }else{
        ident2use = opt$ident2use
    }
    cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
    #cells used as interested cells
    seurat_ob = SubsetData(seurat_ob,cells = 
                       OldWhichCells( seurat_ob, subset.name= ident2use, accept.value = cluster_list))
#    seurat_ob = subset(seurat_ob,cells = 
#                       OldWhichCells( seurat_ob, subset.name= ident2use, accept.value = cluster_list))
}

if (ncol(seurat_ob) > 50000) {
# library(sampling) need install
ratio <- as.numeric(opt$downsample) / ncol(seurat_ob)
metadata_temp <- as.data.frame(seurat_ob@meta.data)
# strata(metadata_temp,stratanames="clusters",ratio,description=FALSE)
cells_sample <- c()
for (i in unique(seurat_ob$clusters)) {
  cells_temp <- rownames(metadata_temp)[which(metadata_temp$clusters == i)]
  cells_temp_sample <- sample(cells_temp, ceiling(length(cells_temp) * ratio), replace = FALSE, prob = NULL)
  cells_sample <- append(cells_sample, cells_temp_sample)
}
seurat_ob <- subset(seurat_ob, cells = cells_sample)
}

if (opt$replot == F) {
    # This is generated from the Velocyto python command line tool.
    # You need a loom file before you can proceed
    # ldat <- read.loom.matrices( opt$loom )
    ldat_seurat_list = lapply( unlist(strsplit(opt$loom, ",")), function(loomx){
        ldat = ReadVelocity(loomx)
        ldat_seurat = as.Seurat(ldat)
        # This is a little bit of foo-magic that needs to be adjusted on a per-sample
        # basis depending on the cell names and how you ran the pipeline. Each cell
        # stored in the loom object and seurat have an ID, make sure these are the same.
        # What this step does is essentially this:
        # > head(Cells(ldat_seurat))
        # [1] "possorted_genome_bam_XL2S3:AAAGATGCATACTACGx"
        # [2] "possorted_genome_bam_XL2S3:ACCTTTATCTTTAGTCx"
        # [3] "possorted_genome_bam_XL2S3:AAGGAGCCACGCATCGx"
    
        cellids = gsub(":", "_",Cells(ldat_seurat))
        cellids = gsub("x$", "", cellids) # trim the suffix "x" at the end of the cell names
        # Now the names in emat and nmat will match up to the cell names used in my seurat object
        # > head(Cells(ldat_seurat))
        # [1] "possorted_genome_bam_XL2S3_AAAGATGCATACTACG"
        # [2] "possorted_genome_bam_XL2S3_ACCTTTATCTTTAGTC"
        # [3] "possorted_genome_bam_XL2S3_AAGGAGCCACGCATCG"
        ldat_seurat = RenameCells(ldat_seurat, new.names = cellids)
    })
    if (length(ldat_seurat_list) >1 ) {
        ldat_seurat = Seurat:::merge.Seurat( ldat_seurat_list[[1]],
                                    y = ldat_seurat_list[2:length(ldat_seurat_list)],
                                    merge.data = T)
    }else {
        ldat_seurat = ldat_seurat_list[[1]]
    }
    ldat_seurat = SubsetData(ldat_seurat, cells = Cells(seurat_ob))
    seurat_ob[["spliced"]] = ldat_seurat[["spliced"]]
    seurat_ob[["unspliced"]] = ldat_seurat[["unspliced"]]
    
    # Gather the spliced and unspliced estimates
    # emat <- ldat$spliced
    # nmat <- ldat$unspliced
    
    seurat_ob <- RunVelocity(seurat_ob, deltaT = 1, kCells = 25, fit.quantile = 0.02,reduction = opt$reduct1)
    saveRDSMC(seurat_ob, file.path(output_dir, "seurat_Velocity.rds"))
}

# Estimate the cell-cell distances 
# cell.dist <- as.dist(1-armaCor(t( emb )))

# Main velocity estimation
# rvel.cd <- gene.relative.velocity.estimates(emat,nmat, deltaT=2,
#                                             kCells=10, cell.dist=cell.dist,
#                                             fit.quantile=fit.quantile, n.cores=24)

# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
# gg <- TSNEPlot(object)
# ggplot_build(gg)$data
# colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
# names(colors) <- rownames(emb)
# ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seurat_ob)))
 seurat_ob = SetIdent(seurat_ob,value= opt$groupby )
 ident.colors <- CustomCol2(1:length(unique(Idents(seurat_ob))))
 #names(ident.colors) <- levels(seurat_ob)
 names(ident.colors) <- as.character(unique(Idents(seurat_ob)))
 cell.colors <- ident.colors[Idents(seurat_ob)]
 names(cell.colors) <- colnames(seurat_ob)
 leg_name  = levels(Idents(seurat_ob))
 leg_color = ident.colors
# take embedding from the Seurat data object
# NOTE: This assumes you have a seurat data object loaded
# into R memory prior to using this script. STOP and rerun seurat
# pipeline if you do not have this loaded. In my case, my seurat object is simply myData
emb = Embeddings( seurat_ob, reduction = opt$reduct.vis)
arrow.scale=(range(emb[,1])[2]-range(emb[,1])[1])/40*4 # grid.n=40. Here the max arrow.scale is twice the length of one grid cell
pdf(paste0(output_dir,"/velocity_plot_by_",opt$groupby,".pdf"), width = 10, height = 10)
print(opt$pointsize)
print(arrow.scale)

par(mai=c(1,1,1,2))

show.velocity.on.embedding.cor(emb, 
                                vel = Tool(object = seurat_ob, slot = "RunVelocity"),
                                n = 200, 
                                scale = "sqrt", min.arrow.size = 0,
                                cell.colors = ac(x = cell.colors, alpha = 0.8),
                                cex = opt$pointsize , arrow.scale = arrow.scale, show.grid.flow = TRUE,
                                min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                                do.par = F,n.cores=24, cell.border.alpha = 0)
								
legend(max(seurat_ob@reductions[[opt$reduct.vis]]@cell.embeddings[,1])+0.9,1.5, xpd=TRUE,legend=leg_name, col=leg_color, pch=19 ) 

#legend("topright", inset=c(-0.2,0), ,legend=leg_name, col=leg_color, pch=15 ) 
dev.off()
if(!file.exists(file.path(output_dir, "RNA Velocity分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/RNA Velocity分析说明.docx",
  file.path(output_dir, "RNA Velocity分析说明.docx"))
}

