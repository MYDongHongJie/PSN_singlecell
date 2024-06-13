#!/usr/bin/env Rscript
#Because of the variablity of single cell RNA-seq expression ,there are many variables could
# result in the hetergenesity.
#This is the script used to visualize the differet QC metrics of the gene expression matrix.
# The possible QC metrics can be:
# nGene: the number of genes for each cells
# nUMI : the number of UMIs for each cells
# percent.mito： the percentage of mitochondrial derived transcript of each cell
# percent.ribo: the percentage of ribosome transcript of each cell
# cell cycle: the cell cycle inference
# batch effect: the assay batch derived variance
# etc.
# By visualizing we could decide which variable should be preprocessed before further analysis.

#========import packages=====================================
rm(list=ls())
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("optparse"))

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-d" ), type = "character",
                 help = "The input exprssion matrix in several possible format."),
    make_option( c("--informat", "-f" ), type = "character", default = "tenx",
                 help = "The indication of type of input expression matrix, the possible type can be:
                        tenx:the directory of cellranger count/aggr results with sampleid as its subdirectory.
                        sce: the singlecellexperiment format of single cell expression matrix.
                        seurat: the seurat object from the clustering results."),
    make_option( c("--metadata", "-m" ), type="character", default = NULL,
        help="the sample metadata which must include sample id in this assay design.", metavar="character"),
    make_option( c("--colorby", "-c" ), type = "character", default = "orig.ident",
        help = "[Otional]visualize cells' metadata by coloring cells in different color according to cell grouping levels."),
    make_option( c("--prefix", "-p" ), type = "character", default = "beforeQC",
        help = "[Optional]the prefix of the output file."),
    make_option( c("--metrics", "-v" ), type = "character", default = "nGene,nUMI,percent.mito,percent.ribo,cellcycle",
        help = "The QC metrics list of source of hetergenesity to visualize.[default:nGene,nUMI,percent.mito,percent.ribo,cellcycle]"),
    make_option( c("--perplexity", "-P"), type="integer", default=30,
        help="[Optional]The value of the perplexity used for tSNE"),
    make_option( c("--components", "-n"), type="integer", default=20,
        help="[Optional]the appropriate number of statistically significant components to use for clustering,
                 which derived from the JackStraw result."),
    make_option( c("--ptsize" ), type="character", default = 1,
        help="[OPTIONAL]setting the dot point size on the violin plot. Set this to NA/0 when you prefer no points."),
    make_option( c("--pointsize", "-s"), type = "double", default = 0.5,
        help = "[OPTIONAL]the point size in the plot."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of QC results.", metavar="outputdir")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null( opt$prefix ) ){
    prefix = "beforeQC"
}else{
    prefix = opt$prefix
}

if ( is.null( opt$colorby ) ){
    print("The groupping information is not specified. The sampleid will be used for plot")
    groupfactor = "orig.ident"
}else{
    groupfactor = opt$colorby
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

output_dir = normalizePath(output_dir )

if ( is.null( opt$metrics ) ){
    print( "NO QC metrics list AVAILABLE! The default nGene,nUMI,percent.mito,percent.ribo,cellcycle will be choosed!")
    QC_metrics = unlist( strsplit( "nGene,nUMI,percent.mito,percent.ribo,cellcycle", ",", perl = T) )
}else{
    QC_metrics = unlist( strsplit( opt$metrics, ",", perl = T) )
}

#wether to plot the dot on the voilin plot and the size of the dot to use
if ( is.null(opt$ptsize) ){
    pt.size = 1
}else if ( as.numeric(opt$ptsize) == 0 | opt$pt.size == "NA"){
    pt.size = 0
}else{
    pt.size = as.numeric(opt$ptsize)
}

######################################################################
#read in the 10X data
# ####################################################################
#======1.1 Setup the Seurat Object (Required)=========================
if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){# the input is a seurat object which may contain more than one sample
        seurat_ob = readRDS( opt$input )
        # In case of filtered seurat object, there seurat@data are the filtered matrix, in which case the cells
        # in seurat@data should be used.
        seurat_ob@raw.data = seurat_ob@raw.data[,colnames(seurat_ob@data)]

        #check the components to use
        if ( is.null(opt$components) ){
            if ( is.null(seurat_ob@calc.params$RunPCA$pcs.compute) ){
                print( "NO previous components calculation is AVAILABLE, 20 will be used as default.")
                components_num = 20
            }else{
                components_num = seurat_ob@calc.params$RunPCA$pcs.compute
            }
        }else{
            components_num  = opt$components
        }

        #check the value of perplexity
        if ( is.null(opt$perplexity) ){
            if ( is.null(seurat_ob@calc.params$BuildSNN$k.param) ){
                print( "NO previous perplexity is AVAILABLE, 30 will be used as default.")
                perplexity = 30
            }else{
                perplexity = seurat_ob@calc.params$BuildSNN$k.param
            }
        }else{
            perplexity = opt$perplexity
        }

        seurat_ob <- NormalizeData(object = seurat_ob, normalization.method = "LogNormalize",scale.factor = 10000)
        seurat_ob = FindVariableGenes(object = seurat_ob, mean.function = ExpMean,
        dispersion.function = LogVMR,do.plot = F,
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
        seurat_ob <- ScaleData(object = seurat_ob,display.progress = F,num.cores = 10,do.par= T) #takes some time

        seurat_ob <- RunPCA(object = seurat_ob,do.print = F)
        if (ncol(seurat_ob@data)<100) {
            seurat_ob <- RunTSNE(object = seurat_ob, dims.use = 1:components_num, perplexity= perplexity,
                                check_duplicates = F, do.fast = T, dim.embed = 2)
        } else {
            seurat_ob <- RunTSNE(object = seurat_ob, dims.use = 1:components_num, perplexity=perplexity,
                                check_duplicates = F, do.fast = T, dim.embed = 2)
        }
        # The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
        # For non-UMI data, nUMI represents the sum of the non-normalized values within a cell.
        # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
        # ,since this represents non-transformed and non-log-normalized counts.
        # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
        if ( "nGene" %in% QC_metrics ){
            vln4nGene = VlnPlot(object = seurat_ob, features.plot = "nGene", group.by = groupfactor, do.return = T,
                        nCol = 1, x.lab.rot = T,point.size.use = pt.size)
            ggsave(file.path(output_dir,paste0(prefix,"_total_genes_per_cell.png", collapse = "")),dpi=1000, plot = vln4nGene)
            ggsave(file.path(output_dir,paste0(prefix,"_total_genes_per_cell.pdf", collapse = "")), plot = vln4nGene)
            pdf(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_PCA_plot.pdf", collapse = "")), width = 8)
            FeaturePlot(seurat_ob, features.plot="nGene", reduction.use = "pca",
                            cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="nGene", reduction.use = "pca",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            pdf(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_tSNE_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot=c("nGene"), reduction.use = "tsne",
                            cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_tSNE_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot=c("nGene"), reduction.use = "tsne",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
        }

        if ( "nUMI" %in% QC_metrics ){
            vln4nUMI = VlnPlot(object = seurat_ob, features.plot = "nUMI" ,group.by = groupfactor, nCol = 1,
                                do.return = T, x.lab.rot = T,point.size.use = pt.size)
            ggsave(file.path(output_dir,paste0(prefix,"_total_UMIs_per_cells.png")),dpi=1000, plot = vln4nUMI)
            ggsave(file.path(output_dir,paste0(prefix,"_total_UMIs_per_cells.pdf")), plot = vln4nUMI)
            pdf(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="nUMI", reduction.use = "pca",
                                cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="nUMI", reduction.use = "pca",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()

            pdf(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_tSNE_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="nUMI", reduction.use = "tsne",
                                cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"nUMI4each_cell_on_tSNE_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="nUMI", reduction.use = "tsne",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
        }

        if ( "batch" %in% QC_metrics ){
            #to visulaize the metadata of each cell on the PCA plot
            pca4batch = PCAPlot(object = seurat_ob, dim.1 = 1, dim.2 = 2,
                                pt.size = opt$pointsize, group.by = "batchid", do.return = T)
            ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_PCA_plot.pdf",collapse="")))
            ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_PCA_plot.png",collapse="")),
                                dpi=1000, plot = pca4batch)
            tsne4batch = TSNEPlot(object = seurat_ob, do.label = F, pt.size = opt$pointsize, group.by = "batchid", do.return = T )
            ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_tSNE_plot.pdf", collapse = "")), plot = tsne4batch)
            ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_tSNE_plot.png", collapse = "")),
                                dpi=1000, plot = tsne4batch)
        }

        # We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
        if ( "percent.mito" %in% QC_metrics ){
            mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(x = seurat_ob@data), value = T,perl=T)
            percent.mito <- Matrix::colSums(seurat_ob@raw.data[mito.genes, ])/Matrix::colSums(seurat_ob@raw.data)
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
            vln4mito = VlnPlot(object = seurat_ob, features.plot = "percent.mito",group.by = groupfactor, nCol = 1,do.return = T,
                        x.lab.rot = T,point.size.use = pt.size)
            ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_violin_plot.png")),
                        dpi=600, plot = vln4mito )
            ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_violin_plot.pdf" )), plot = vln4mito)

            pdf(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="percent.mito", reduction.use = "pca",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="percent.mito", reduction.use = "pca",
            cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            pdf(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot=c('percent.mito'), reduction.use = "tsne",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_tSNE_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot=c('percent.mito'), reduction.use = "tsne",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
        }

        # We calculate the percentage of ribosome genes here and store it in percent.ribo using AddMetaData.
        if ( "percent.ribo" %in% QC_metrics ){
            #Calculate percent ribosomal genes.
            ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", rownames(seurat_ob@data), value = TRUE,ignore.case = T)
            percent.ribo <- Matrix::colSums(seurat_ob@raw.data[ribo.genes, ])/Matrix::colSums(seurat_ob@raw.data)
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.ribo, col.name = "percent.ribo")
            vln4ribo = VlnPlot(object = seurat_ob, features.plot = "percent.ribo",group.by = groupfactor, nCol = 1,
                                do.return = T, x.lab.rot = T,point.size.use = pt.size)
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_violin_plot.png")),
                                dpi=600, plot = vln4ribo )
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_violin_plot.pdf" )), plot = vln4ribo)

            pdf(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="percent.ribo", reduction.use = "pca",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="percent.ribo", reduction.use = "pca",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            pdf(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot='percent.ribo', reduction.use = "tsne",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_tSNE_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot='percent.ribo', reduction.use = "tsne",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
        }

        #wether to regress out the cell cycle effect
        if ( "cellcycle" %in% QC_metrics ){
            s.genes = R.utils::capitalize(tolower(s.genes))
            g2m.genes = R.utils::capitalize(tolower(g2m.genes))
            seurat_ob <- CellCycleScoring(object = seurat_ob, s.genes = s.genes, g2m.genes = g2m.genes,set.ident = TRUE)
            seurat_ob <- RunPCA(object = seurat_ob, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
            pca4cellcycle = PCAPlot(object = seurat_ob, dim.1 = 1, dim.2 = 2, pt.size = opt$pointsize, do.return = T)
            ggsave(file.path(output_dir,paste0(prefix, "_visualize_cell_cycle4each_cell_PCA_plot.pdf",collapse="")),
                    plot = pca4cellcycle)
            ggsave(file.path(output_dir,paste0(prefix, "_visualize_cell_cycle4each_cell_PCA_plot.png",collapse="")),
                    dpi=1000, plot = pca4cellcycle)
            # seurat_ob <- RunTSNE(object = seurat_ob, dims.use = 1:components_num, perplexity=perplexity,force.recalc = T,
            #                     check_duplicates = F, do.fast = T, dim.embed = 2,genes.use =c(s.genes, g2m.genes))
            # tsne4cellcycle = TSNEPlot(object = seurat_ob, do.label = F, pt.size = opt$pointsize, do.return = T )
            # ggsave(file.path(output_dir,paste0(prefix, "_visualize_cell_cycle4each_cell_tSNE_plot.pdf",collapse="")),
            #         plot = tsne4cellcycle)
            # ggsave(file.path(output_dir,paste0(prefix, "_visualize_cell_cycle4each_cell_tSNE_plot.png",collapse="")),
            #         dpi=1000, plot = tsne4batch)
        }

        if ( "top50" %in% QC_metrics ){
            top50_gene = apply( seurat_ob@raw.data, 2, function(x) sum(x[order(x, decreasing = T)][1:50])/sum(x))
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = top50_gene, col.name = "top50")
            vln4top50 = VlnPlot(object = seurat_ob, features.plot = "top50",group.by = groupfactor, nCol = 1,
                                do.return = T, x.lab.rot = T,point.size.use = pt.size)
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_violin_plot.png")),
                                dpi=600, plot = vln4top50 )
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_violin_plot.pdf" )),
                                plot = vln4top50)
            pdf(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="top50", reduction.use = "pca",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="top50", reduction.use = "pca",
                        cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            pdf(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="top50", reduction.use = "tsne",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_tSNE_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="top50", reduction.use = "tsne",
                            cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
        }

        if ( "UMI.per.gene" %in% QC_metrics ){
            tmp = seurat_ob@meta.data$nUMI/seurat_ob@meta.data$nGene
            names(tmp) = row.names(seurat_ob@meta.data)
            seurat_ob = AddMetaData(seurat_ob, metadata = tmp, "UMI.per.Gene")
            vln4UPG = VlnPlot(object = seurat_ob, features.plot = "UMI.per.Gene",group.by = groupfactor,
                        nCol = 1, x.lab.rot = T,point.size.use = pt.size, do.return = T)
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_violin_plot.png")),
            dpi=600, plot = vln4top50 )
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_violin_plot.pdf" )),
            plot = vln4top50)
            pdf(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="UMI.per.Gene", reduction.use = "pca",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="UMI.per.Gene", reduction.use = "pca",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()

            pdf(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.pdf", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="UMI.per.Gene", reduction.use = "tsne",
                                    cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
            png(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.png", collapse = "")))
            FeaturePlot(seurat_ob, features.plot="UMI.per.Gene", reduction.use = "tsne",
                            cols.use = c("grey","red"),pt.size=opt$pointsize,no.legend = F)
            dev.off()
        }
    }
}
# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
# par(mfrow = c(1, 2))
# nUMI2mito = GenePlot(object = seurat_ob, gene1 = "nUMI", gene2 = "percent.mito",col.use = "red")
# ggsave(file.path(output_dir,"nUMI_vs_percent.mito_gene_feature_plot.pdf"),plot = nUMI2mito)
# nUMI2nGene = GenePlot(object = seurat_ob, gene1 = "nUMI", gene2 = "nGene", col.use="green")
# ggsave(file.path(output_dir,"gene_feature_plot.png"), plot= nUMI2nGene,dpi=600)


#save the filtered data into disk by sample seperately if there are more than one sample
#for this project
# for ( idx in 1:nrow(assay_metadata) ){
#     sample_id = assay_metadata[idx,"sampleid"]
#     message(paste0("Save QC results for Sample ",sample_id))
#     sample_matrix = file.path(output_dir,sample_id,"outs","filtered_gene_bc_matrices",assay_metadata[idx,"specie"])
#     if (!dir.exists(sample_matrix)){
#         dir.create(sample_matrix,recursive = T)
#     }
#     #subset the seurat object for the specified specie
#     # sample_seurat = SubsetData(seurat_ob, cells.use = seurat_ob@cell.names[seurat_ob@meta.data$orig.ident==sample_id])
#     sample_seurat = SubsetData(seurat_ob, cells.use = FastWhichCells(seurat_ob, group.by = "sampleid", subset.value = sample_id ),subset.raw = T)
#
#     #visualize the data after filtering
#     VlnPlot(object = sample_seurat, features.plot = "nGene", nCol = 1,group.by="orig.ident", point.size.use = pt.size)
#     ggsave(file.path(output_dir,sample_id, "afterQC_total_genes_per_cell.png"),dpi=1000)
#     ggsave(file.path(output_dir,sample_id, "afterQC_total_genes_per_cell.pdf"))
#
#     VlnPlot(object = sample_seurat, features.plot = "nUMI" , nCol = 1,group.by="orig.ident", point.size.use = pt.size)
#     ggsave(file.path(output_dir, sample_id, "afterQC_total_UMIs_per_cells.png"),dpi=1000)
#     ggsave(file.path(output_dir, sample_id, "afterQC_total_UMIs_per_cells.pdf"))
#
#     VlnPlot(object = sample_seurat, features.plot = "percent.mito", nCol = 1,group.by="orig.ident", point.size.use = pt.size)
#     ggsave(file.path(output_dir,sample_id, "afterQC_total_mitochondroin_per_cells.png"), dpi=600)
#     ggsave(file.path(output_dir,sample_id, "afterQC_total_mitochondroin_per_cells.pdf"))
# }
