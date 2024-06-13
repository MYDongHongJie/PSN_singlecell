#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("monocle"))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))

#=command line parameters setting=============================
option_list = list(
    make_option( c( "--INEXPRESS","-i"), type = "character",
        help = "the input expression matrix to be analyzed, the format
                          can be specified by the option -f INFORMAT."
    ),
    make_option( c("--INFORMAT", "-f" ), type = "character",default = NULL,
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
    make_option( c("--perplexity", "-p"), type="integer", default=30,
        help="The value of the perplexity used for tSNE"),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of Clustering results." ),
    make_option( c("--metadata", "-m" ), type="character", default = NULL,
        help="the sample metadata which must include sample id in this assay design."),
    make_option( c("--from"), type = "character", default = NULL,
        help="the original cluster id used in calibrate the clustering results."),
    make_option( c("--to"), type = "character", default = NULL,
        help="the adjusted cluster id used in calibrate the clustering results."),
    make_option( c("--resolution","-r"), type = "double", default = 1,
        help = "vaule used to set the resolution of cluster distiguish,
                     use a value above(below)1.0 if you want to obtain a larger(smaller) number of communities."),
    make_option( c("--CORES", "-j" ), type = "integer", default = 20,
        help = "the core number used to run this script." ),
    make_option( c("--colorby", "-C" ), type = "character", default = "clusters,sampleid,group,State",
        help = "[Otional]visualize cells' metadata by coloring cells in different color according to cell grouping list."),
    make_option( c("--groupby", "-c" ), type = "character", default = "orig.ident",
        help = "[Otional]the column name in cell grouping metadata used with which_cell parameter for subseting cell groups."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
        help = "The subset of cluster ids used for subtyping."),
    make_option( c("--pointsize", "-s" ), type = "double",  default = 1,
        help = "[Otional]the point size in the plot."),
    make_option( c("--design", "-d" ), type = "character",  default = "clusters",
        help = "[Otional]The group design to find ordering genes."),
    make_option(c("--downsample", "-e"),
        type = "character", default = "30000",
        help = "the downsample number of cells "
      ),
    make_option( c("--assay"), type = "character", default = "RNA",
                 help = "[OPTIONAL] The array result to use in case of multimodal analysis."),
    make_option(c("--topn", "-n"),
        type = "integer", default = NULL,
        help = "the subset of ordering gene "
      ),
    make_option( c("--rds" ), type = "character",  default = NULL,
        help = "[Otional]The pseudotime_result.rds input to replot.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$downsample)) {
  downsample <- 30000
} else {
  downsample <- opt$downsample
}


if ( is.null(opt$output)){
    print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
    output_dir = getwd()
}else{
    output_dir = opt$output
    if ( !file.exists(output_dir) ){
        dir.create(output_dir)
    }
}

#check the components to use
if ( is.null(opt$components) ){
    print("NO components number is available, the default will be used")
}

if ( is.null( opt$groupby ) ){
    print("The groupping information is not specified. The sampleid will be used for plot")
    groupfactor = "orig.ident"
}else{
    groupfactor = opt$groupby
}

#check the value of perplexity
if ( is.null(opt$perplexity) ){
    print("NO perplexity number is specified, the default will be used!")
}

if ( is.null( opt$design ) ){
    print("NO group design is provided! The clusters will be used!")
    design = "clusters"
}else{
    design = opt$design
}

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}


gbm.PlotAbundances <- function (x = seurat, prop.by = "res.1", group.by = "sampleid", 
    split.by = NULL, method = "barplot", ncol = NULL, cols = CustomCol(1:length(unique(x[, 
        prop.by])))) 
{
    md <- x@phenoData@data
    field4emeta = colnames(md)
    if (is.null(prop.by) | !prop.by %in% field4emeta) {
        stop("prop.by should be one of the clustering results name in meta.data.")
    }
    if (is.null(group.by) | !group.by %in% field4emeta) {
        stop("One of sampleid and clusters field name should be provided for groupping!")
    }
    cluster_ids <- md[, prop.by]
    sample_ids <- md[, group.by]
    counts <- table(cluster_ids, sample_ids)
    df <- melt(t(round(t(counts)/colSums(counts) * 100, 2)), 
        varnames = c(prop.by, group.by), value.name = "freq")
    df[, group.by] = factor(df[, group.by])
    df[, prop.by] = factor(df[, prop.by])
    if (!is.null(split.by)) {
        xm = unique(md[, c(group.by, split.by)])
        m <- match(df[, group.by], xm[, group.by])
        df[, split.by] = xm[m, setdiff(names(xm), names(df))]
        df[, split.by] = factor(df[, split.by])
    }
    p = switch(method, barplot = ggplot(df) + geom_bar(aes_string(y = "freq", 
        x = group.by, fill = prop.by), position = "stack", stat = "identity") + 
        scale_fill_manual(values = cols) + scale_y_continuous(expand = c(0, 
        0), labels = seq(0, 100, 25)) + labs(x = NULL, y = "Propotion [%]") + 
        theme_bw() + theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
            color = NA), panel.border = element_blank(), axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 90, 
            hjust = 1, vjust = 0.5)), boxplot = ggplot(df) + 
        theme_bw() + geom_boxplot(aes_string(y = "freq", x = group.by, 
        color = group.by, fill = group.by), position = position_dodge(), 
        alpha = 0.25, outlier.color = NA) + scale_fill_manual(values = cols) + 
        geom_point(position = position_jitter(width = 0.25), 
            aes_string(x = group.by, y = "freq", color = group.by)) + 
        guides(fill = FALSE) + labs(x = NULL, y = "Proportion [%]") + 
        theme(strip.background = element_rect(fill = NA, color = NA), 
            panel.spacing = unit(0.2, "cm"), strip.text = element_text(face = "bold"), 
            panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
            axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
            axis.text.x = element_text(angle = 90, hjust = 1, 
                vjust = 0.5)), pie = ggplot(df, aes(x = 0, y = freq)) + 
        geom_bar(aes_string(fill = prop.by), width = 1, stat = "identity") + 
        theme_void() + scale_fill_manual(values = cols) + theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
        coord_polar(theta = "y", start = 0) + facet_wrap(eval(expr(~!!ensym(group.by))), 
        ncol = ncol) + theme(strip.background = element_rect(fill = NA, 
        color = NA), strip.text = element_text(face = "bold")))
    if (!is.null(split.by)) {
        p = p + facet_grid(eval(expr(~!!ensym(split.by))), scales = "free_x", 
            space = "free")
    }
    return(p)
}


# #####################################################################
#
#  2. Load and normalize data
#  https://www.bioconductor.org/packages/devel/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf
#
# ####################################################################
#whatever format the input file has, all the input expression matrix
#should be transformed to the seurat object finally.
if ( !is.null(opt$INEXPRESS) ){
    # the input is a seurat object which may contain more than one sample
    if ( opt$INFORMAT == "seurat" ){
        seurat_ob = readRDS( opt$INEXPRESS )
        if ( seurat_ob@version < 3 ){
            seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }
        
        if ( !is.null( opt$assay) ){
            DefaultAssay(seurat_ob) = opt$assay
        }else{
            DefaultAssay(seurat_ob) = "RNA"
        }
 
        if ( is.null(seurat_ob@meta.data$clusters) ){
            seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
        }else{
            seurat_ob = SetIdent( seurat_ob, value = "clusters")
        }

        if ( !is.null(opt$from) & !is.null(opt$to) ){
            from_ident = unlist( strsplit( opt$from, ",", perl =T ) )
            to_ident = unlist( strsplit(opt$to), ",", perl = T )
            seurat_ob@ident = plyr::mapvalues(x = seurat_ob@ident, from = from+ident, to = to_ident )
        }

        #subset the seurat object using the specified cells of clusters from the
        #command line if available
        if ( !is.null(opt$which_cells)){
            cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
            seurat_ob = SubsetData(seurat_ob,cells = OldWhichCells(seurat_ob,
                                        subset.name = groupfactor, accept.value = cluster_list ))
        }

        #check the components to use
        if ( is.null(opt$components) ){
            if ( is.null(Misc(seurat_ob, "optimal_pc")) ){
                print( "NO previous components calculation is AVAILABLE, 20 will be used as default.")
                components_num = 20
            }else{
                components_num = Misc(seurat_ob, "optimal_pc")
            }
        }else{
            components_num  = opt$components
        }

        #check the value of perplexity
        if ( is.null(opt$perplexity) ){
            if ( is.null(Misc(seurat_ob, "perplexity")) ){
                print( "NO previous perplexity is AVAILABLE, 30 will be used as default.")
                perplexity = 30
            }else{
                perplexity = Misc(seurat_ob, "perplexity")
            }
        }else{
            perplexity = opt$perplexity
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
          print(dim(seurat_ob))
    }else{
        if ( opt$informat == "tenx" ){
            # Initialize the Seurat object with the raw (non-normalized data).
            # Keep all genes expressed in >= 3 cells (~0.1% of the data).
            #Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X.
            #A vector or named vector can be given in order to load several data directories of several samples.
            #If a named vector is given, the cell barcode names will be prefixed with the name.
            #the gene-barcode matrix result directory for all samples from 10X
            if ( !file.exists( opt$metadata) ){
               stop( "Warning:the metadata of this single cell assay is NOT AVAILABLE!")
            }else{
               assay_metadata = read.table(opt$metadata,sep="\t",header =T, row.names = F,check.names = )
            }
            outs_filtered_path = "outs/filtered_gene_bc_matrices"
            tenx_path = sub("\\/$","",opt$INEXPRESS,perl=T)
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
        seurat_ob <- NormalizeData(object = seurat_ob, normalization.method = "LogNormalize", scale.factor = 10000)

        #========1.4 find the highly variable genes
        seurat_ob = FindVariableGenes(object = seurat_ob, mean.function = ExpMean,
        dispersion.function = LogVMR,do.plot = F,
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
        if ( opt$BATCH == T){
            suppressPackageStartupMessages(library(sva))
            m = as.matrix(seurat_ob@data)
            com = ComBat(dat=m, batch=seurat_ob@meta.data$batchid,
            prior.plots=FALSE, par.prior=TRUE)
            seurat_ob@data = Matrix(com)
            rm(m)
            rm(com)
        }

        seurat_ob <- RunPCA(object = seurat_ob,pc.genes = seurat_ob@var.genes,do.print = F)

        seurat_ob <- FindClusters(object = seurat_ob, reduction.type = "pca",
        dims.use = 1:components_num, resolution = resolution,
        print.output = 0, save.SNN = TRUE,force.recalc=T)
        pca_cluster_result_colname = paste("pca","2D","res",resolution,sep = ".")
        #as default the seurat will store the clustering results in the seurat_object@ident
        #to keep the clutering results using the specified resolution to seurat_object@metadata for reuse
        seurat_ob <- StashIdent(object = seurat_ob, save.name = pca_cluster_result_colname)
    }
}

# if (ncol(seurat_ob@data)<100) {
    # seurat_ob <- RunTSNE(object = seurat_ob, dims.use = 1:components_num, perplexity=10,
                        # check_duplicates = F, do.fast = T, dim.embed = 2)
# } else {
    # seurat_ob <- RunTSNE(object = seurat_ob, dims.use = 1:components_num, perplexity=perplexity,
                        # check_duplicates = F, do.fast = T, dim.embed = 2)
# }
if (is.null(opt$rds)){
    #transform the seurat object to the monocle CellDataSet
    # gbm_cds = importCDS( seurat_ob, import_all = F )
    gbm_cds = as.CellDataSet( seurat_ob )
    #Estimate size factors and dispersions
    gbm_cds = estimateSizeFactors(gbm_cds)
    gbm_cds = estimateDispersions(gbm_cds)
    #Filtering low-quality cells
    gbm_cds = detectGenes(gbm_cds,min_expr = 1)
    #keep only the genes expressed in at least 10 cells of the data set.
    expressed_genes = row.names(subset(fData(gbm_cds),num_cells_expressed>10))
 
    pData(gbm_cds)$Total_mRNA = Matrix::colSums(exprs(gbm_cds))
    # gbm_cds = clusterCells(gbm_cds)
    #Find ordering genes
    clustering_DEGs = differentialGeneTest(gbm_cds,fullModelFormulaStr = paste0("~", design),cores = opt$CORES)
    featureData(gbm_cds)@data[rownames(clustering_DEGs),"pval"]=clustering_DEGs$pval
    featureData(gbm_cds)@data[rownames(clustering_DEGs),"qval"]=clustering_DEGs$qval
    ordering_genes <- row.names (subset(clustering_DEGs, qval < 0.01))
    if(!is.null(opt$topn)){
        topn=opt$topn
        data <- featureData(gbm_cds)@data
        if(opt$assay=="RNA"){
            ordering_genes <- subset(data, qval < 0.01) %>% tibble::rownames_to_column(var="gene" ) %>% dplyr::arrange(desc(vst.variance.standardized))  %>% dplyr::top_n(as.numeric(topn),vst.variance.standardized)
        }else if(opt$assay=="SCT"){
            ordering_genes <- subset(data, qval < 0.01) %>% tibble::rownames_to_column(var="gene" ) %>% dplyr::arrange(desc(sct.variance))  %>% dplyr::top_n(as.numeric(topn),sct.variance)  
        } 
        write.table(ordering_genes,file.path(output_dir,'ordering_genes.xls'),sep='\t',row.names=F,quote=F)   
        ordering_genes <- ordering_genes$gene
    }
    # ordering.genes = row.names(fData(gbm_cds))[1:100]
    gbm_cds = setOrderingFilter(gbm_cds,ordering_genes = ordering_genes)
    #the dependent package irlba should not be less than 2.3.2, or the following step will fail
    #if the size of cells is big.
    gbm_cds = reduceDimension(gbm_cds,max_components = 2,verbose = T,check_duplicates = F)
    gbm_cds = orderCells(gbm_cds,reverse = F)

    #Remove the state with 0 cell
    for (i in which(table(gbm_cds$State)==0)){
    gbm_cds$State=as.numeric(gbm_cds$State)
    gbm_cds$State[gbm_cds$State>i]=gbm_cds$State[gbm_cds$State>i]-1
    gbm_cds$State=as.factor(gbm_cds$State)
    }

    saveRDSMC(gbm_cds,file.path(output_dir,"pseudotime_results.rds"))

    # Pseudotime on tsne plot
    seurat_ob@meta.data$Pseudotime = ""
    seurat_ob@meta.data[rownames(gbm_cds@phenoData),]$Pseudotime = gbm_cds$Pseudotime
    seurat_ob@meta.data$Pseudotime = as.numeric(seurat_ob@meta.data$Pseudotime)
    # p = FeaturePlot(seurat_ob, features= "Pseudotime", reduction = "tsne", pt.size= opt$pointsize) + scale_colour_viridis_c(option = "inferno")
    # ggsave(file.path(output_dir,paste0("Pseudotime_on_tsne_plot.pdf",collapse = "")), plot = p)
    # ggsave(file.path(output_dir,paste0("Pseudotime_on_tsne_plot.png",collapse = "")), plot = p, dpi=1000)
    Pseudotime_State= gbm_cds$State
    seurat_ob <- AddMetaData(object = seurat_ob,metadata = Pseudotime_State, col.name = "Pseudotime_State")
    if ( !is.null(opt$INEXPRESS) ){
        if ( !is.null(opt$which_cells)){
            saveRDS(seurat_ob,file.path(output_dir,paste0("subseted_",groupfactor,"_",paste0(cluster_list,collapse="-"),"_seurat_withpseudotime.rds")))                                               
        } else{
            saveRDS(seurat_ob,file.path(output_dir,"seurat_withpseudotime.rds"))
        }
    }
} else  gbm_cds = readRDS(opt$rds)
# cell trajectory plot
# p = plot_cell_trajectory(gbm_cds, color_by = "Pseudotime", cell_size = opt$pointsize, show_branch_points = F) + scale_colour_viridis_c(option = "inferno")
# ggsave(file.path(output_dir,"cell_trajectory_color_by_Pseudotime.pdf"), plot = p )
# ggsave(file.path(output_dir,"cell_trajectory_color_by_Pseudotime.png"), plot = p, dpi=1000)

nlevel = length(unique(gbm_cds@phenoData@data[,"State"]))
pp = plot_cell_trajectory(gbm_cds, color_by = "State", cell_size = opt$pointsize, show_branch_points = F) + scale_colour_manual( values = CustomCol2(1:nlevel)) + guides(colour = guide_legend(override.aes = list(size=2)))
ggsave(file.path(output_dir,"cell_trajectory_color_by_State.pdf"), plot = pp )
ggsave(file.path(output_dir,"cell_trajectory_color_by_State.png"), plot = pp, dpi=1000)

color_by = unlist(strsplit(opt$colorby,",",perl=T))

for ( color_group in color_by ){
    if(color_group =="clusters"){
        nlevels = sort(unique(gbm_cds@phenoData@data[,"clusters"]))
        user_color_pal = CustomCol2(nlevels)
        nlevel = length(unique(gbm_cds@phenoData@data[,"clusters"]))
        }else if(paste0(color_group,"_col") %in% colnames(gbm_cds@phenoData@data)){
            groupby_col = paste0(color_group,"_col")
            nlevel = length(unique(gbm_cds@phenoData@data[,color_group]))
            nlevel_list = sort(as.character(unique(gbm_cds@phenoData@data[,color_group])))
            tmp_df <- unique(gbm_cds@phenoData@data[c(color_group, groupby_col)])
            new_celltype_pal <- as.vector(tmp_df[,groupby_col])
            names(new_celltype_pal) <-  as.vector(tmp_df[,color_group])
            new_celltype_pal = as.list(new_celltype_pal)
            user_color_pal = new_celltype_pal[nlevel_list]
    }else {
        nlevel = length(unique(gbm_cds@phenoData@data[,color_group]))
        user_color_pal = CustomCol2(1:nlevel)
    }

    # 获取总字符数
    total_char = sum(nchar(as.character(unique(gbm_cds@phenoData@data[,color_group]))))

    if (total_char < 40){
        p = plot_cell_trajectory(gbm_cds, color_by = color_group, cell_size = opt$pointsize, show_branch_points = F) + scale_colour_manual( values = user_color_pal) + guides(colour = guide_legend(override.aes = list(size=2),nrow=as.numeric((ifelse(nlevel >4,"2","1")))))
        ggsave(file.path(output_dir,paste0("cell_trajectory_color_by_",color_group,".pdf",collapse = "")), plot = p)
        ggsave(file.path(output_dir,paste0("cell_trajectory_color_by_",color_group,".png",collapse = "")), plot = p, dpi=1000)
        # facet_wrap
        p_facet = plot_cell_trajectory(gbm_cds, color_by = color_group, cell_size = opt$pointsize, show_branch_points = F ) +
              facet_wrap(eval(expr(~!!ensym(color_group))),scales = "free_x", ncol = as.numeric((ifelse(nlevel >20,"4","2")))) +scale_colour_manual( values = user_color_pal) + 
              guides(colour = guide_legend(override.aes = list(size=2),nrow=as.numeric((ifelse(nlevel >4,"2","1")))))
        ggsave(file.path(output_dir,paste0("cell_trajectory_group_by_",color_group,".pdf",collapse = "")), plot = p_facet,height= ceiling(nlevel/2)*4)
        try(ggsave(file.path(output_dir,paste0("cell_trajectory_group_by_",color_group,".png",collapse = "")), plot = p_facet, dpi=1000,height= ceiling(nlevel/2)*4))
    } else {
        # 获取最大字符数
        nmax = max(nchar(as.character(unique(gbm_cds@phenoData@data[,color_group]))))
        # 计算nrow数
        nrows = ceiling(nlevel * nmax / 40)

        p = plot_cell_trajectory(gbm_cds, color_by = color_group, cell_size = opt$pointsize, show_branch_points = F) + scale_colour_manual( values = user_color_pal) + guides(colour = guide_legend(override.aes = list(size=2),nrow=nrows))
        ggsave(file.path(output_dir,paste0("cell_trajectory_color_by_",color_group,".pdf",collapse = "")), plot = p)
        ggsave(file.path(output_dir,paste0("cell_trajectory_color_by_",color_group,".png",collapse = "")), plot = p, dpi=1000)
        # facet_wrap
        p_facet = plot_cell_trajectory(gbm_cds, color_by = color_group, cell_size = opt$pointsize, show_branch_points = F ) +
              facet_wrap(eval(expr(~!!ensym(color_group))),scales = "free_x", ncol = as.numeric((ifelse(nlevel >20,"4","2")))) +scale_colour_manual( values = user_color_pal) + 
              guides(colour = guide_legend(override.aes = list(size=2),nrow=nrows))
        ggsave(file.path(output_dir,paste0("cell_trajectory_group_by_",color_group,".pdf",collapse = "")), plot = p_facet,height= ceiling(nlevel/2)*as.numeric((ifelse(nlevel >20,"2","4"))),width=as.numeric((ifelse(nlevel >20,"12","7"))),limitsize = FALSE)
        try(ggsave(file.path(output_dir,paste0("cell_trajectory_group_by_",color_group,".png",collapse = "")), plot = p_facet, dpi=1000,height= ceiling(nlevel/2)*as.numeric((ifelse(nlevel >20,"2","4"))),width=as.numeric((ifelse(nlevel >20,"12","7"))),limitsize = FALSE))
    }
}


# count summary
DATA <- as.data.frame( gbm_cds@phenoData@data[,c("sampleid", "group","State")] ) %>%
    group_by( .dots= c( "State","sampleid","group")) %>%
    dplyr::summarize(cell_number = n())
write.table(as.data.frame(DATA),file.path(output_dir,"Pseudotime_State_count_info.xls"),sep="\t",col.names=T, row.names =F)

# Pie plot
pie_by_state = gbm.PlotAbundances(gbm_cds, prop.by = "sampleid", group.by = "State", method = "pie", ncol = 4, cols = CustomCol2(1:length(table(gbm_cds@phenoData@data[,"sampleid"]))))
ggsave(file.path(output_dir,"groupby_Pseudotime_State_pie_plot.pdf"), plot = pie_by_state)
ggsave(file.path(output_dir,"groupby_Pseudotime_State_pie_plot.png"), plot = pie_by_state, dpi=1000)

pie_by_sampleid = gbm.PlotAbundances(gbm_cds, prop.by = "State", group.by = "sampleid", method = "pie", ncol = 4, cols = CustomCol2(1:length(table(gbm_cds@phenoData@data[,"State"]))))
ggsave(file.path(output_dir,"groupby_sampleid_pie_plot.pdf"),plot = pie_by_sampleid)
ggsave(file.path(output_dir,"groupby_sampleid_pie_plot.png"),plot = pie_by_sampleid, dpi=1000)

if(!file.exists(file.path(output_dir, "拟时序分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/拟时序分析说明.docx",
  file.path(output_dir, "拟时序分析说明.docx"))
}

