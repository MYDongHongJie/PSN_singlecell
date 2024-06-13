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
suppressWarnings({
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("optparse"))
    suppressPackageStartupMessages(library("OESingleCell"))
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages(library("ggplot2"))
})

#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
                 help = "The input exprssion matrix in several possible format."),
    make_option( c("--informat", "-f" ), type = "character", default = "tenx",
                 help = "The indication of type of input expression matrix, the possible type can be:
                        seurat: the seurat object from the clustering results."),
    make_option( c("--metadata", "-m" ), type="character", default = NULL,
        help="[Optional]the additional sample metadata in case of assay design changes.", metavar="character"),
    make_option( c("--colorby", "-c" ), type = "character", default = "sampleid",
        help = "[Otional]visualize cells' metadata by coloring cells in different color according to cell grouping levels."),
    make_option( c("--prefix", "-d" ), type = "character", default = "beforeQC",
        help = "[Optional]the prefix of the output file."),
    make_option( c("--metrics", "-v" ), type = "character", default = "nGene,nUMI,percent.mito,percent.ribo,batch,cellcycle",
        help = "The QC metrics list of source of hetergenesity to visualize.[default:nGene,nUMI,percent.mito,percent.ribo,cellcycle]"),
    make_option( c("--perplexity", "-P"), type="integer", default=30,
        help="[Optional]The value of the perplexity used for tSNE"),
    make_option( c("--components", "-n"), type="integer", default=20,
        help="[Optional]the appropriate number of statistically significant components to use for clustering,
                 which derived from the JackStraw result."),
    make_option( c("--ptsize" ), type="double", default = 1,
        help="[OPTIONAL]setting the dot point size on the violin plot. Set this to NA/0 when you prefer no points."),
    make_option( c("--alpha" ), type="double", default = 0.6,
        help="[OPTIONAL]setting the dot opacity on the violin plot."),
    make_option( c("--assay" ), type="character", default = "RNA",
        help="[OPTIONAL]the assay used to calulation in case of multimodal data."),
    make_option( c("--ncpus", "-p" ), type = "integer", default = 10,
        help = "[Optional]the number of CPUs used to parallize this job."),
    make_option( c("--pointsize", "-s"), type = "double", default = 0.5,
        help = "[OPTIONAL]the point size in the plot."),
    make_option( c("--colshema"),type="character", default = "red,grey",
        help="the color schema for groupping coloring."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory of QC results.", metavar="outputdir"),
    make_option( c("--mito", "-M" ), type = "character", default = NULL,
                  help = "mitochondrial gene list.") );
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
    groupfactor = "sampleid"
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
}else if ( as.numeric(opt$ptsize) == 0 | opt$ptsize == "NA"){
    pt.size = 0
}else{
    pt.size = as.numeric(opt$ptsize)
}

assay2use = opt$assay
alpha2use = opt$alpha
######################################################################
#read in the 10X data
# ####################################################################
#======1.1 Setup the Seurat Object (Required)=========================
if ( !is.null(opt$informat) & !is.null(opt$input) ){
    if ( opt$informat == "seurat" ){# the input is a seurat object which may contain more than one sample
        seurat_ob = readRDSMC( opt$input, cores = availableCores() )
        if ( seurat_ob@version < 3 ){
            seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
        }

        #set the default assay as the specified
        DefaultAssay(seurat_ob) = assay2use

        #update the metedata in the singlecell_ob@meta.data with new additional sample metadata
        if ( !is.null(opt$metadata) ){
            additional_metadata = read.csv(opt$metadata,sep=",",header =T )
            rownames(additional_metadata) = additional_metadata$sampleid
            cellnames = Cells(seurat_ob)
            sampleidx =  gsub("(-|_)[ATGC]{16,}.*","",cellnames,perl=T) #the index order is the same as the row index of the assay metadata
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

        # setting the maxumium chunck mermory usage much bigger in case of big data
        #options(future.globals.maxSize= Inf)
        #plan("multiprocess", workers = opt$ncpus) # parallization start from here
        #if ( !length(VariableFeatures(seurat_ob))){ # if calculated, ignore it.
        #    seurat_ob <- NormalizeData(object = seurat_ob,
        #                            normalization.method = "LogNormalize", scale.factor = 10000)
        #    seurat_ob = FindVariableFeatures(object = seurat_ob,
        #                            mean.function = "FastExpMean",
        #                            dispersion.function = "FastLogVMR",do.plot = F )
        #    seurat_ob <- ScaleData(object = seurat_ob, verbose = F) #takes some time
        #}

        ##check the components to use
        #suppressWarnings({
        #    if ( is.null(opt$components) ){
        #        if ( is.null(Misc(seurat_ob, "components_num")) ){
        #            print( "NO previous components calculation is AVAILABLE, 20 will be used as default.")
        #            Misc(seurat_ob, "components_num") = 20
        #        }else{
        #            components_num = Misc(seurat_ob, "components_num")
        #        }
        #    }else{
        #         Misc(seurat_ob, "components_num") = opt$components
        #        components_num = opt$components
        #    }
        #})

        ##check the value of perplexity
        #suppressWarnings({
        #    if ( is.null(opt$perplexity) ){
        #        if ( is.null(Misc(seurat_ob, "perplexity")) ){
        #            print( "NO previous perplexity is AVAILABLE, 30 will be used as default.")
        #            Misc(seurat_ob, "perplexity") = 30
        #        }else{
        #            perplexity = Misc(seurat_ob, "perplexity")
        #        }
        #    }else{
        #        perplexity = opt$perplexity
        #        Misc(seurat_ob, "perplexity") = opt$perplexity
        #    }
        #})

        #if (!"pca" %in% names(Key(seurat_ob))){ # if calcuated, ignore it
        #    seurat_ob <- RunPCA(object = seurat_ob, features = VariableFeatures(seurat_ob),verbose = F)
        #}
        #if ( !"tsne" %in% names(Key(seurat_ob)) ){# if calcuated, ignore it
        #    seurat_ob <- RunTSNE(object = seurat_ob, dims = 1:components_num,
        #                        reduction = "pca",
        #                        # assay = assay2use,
        #                        features = VariableFeatures(seurat_ob) ,
        #                        perplexity=perplexity,force.recalc = T, check_duplicates = F,
        #                        do.fast = T, dim.embed = 2)
        #}
        # The number of genes and UMIs are automatically calculated for every object by Seurat.
        # For non-UMI data, nUMI represents the sum of the non-normalized values within a cell.
        # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
        # ,since this represents non-transformed and non-log-normalized counts.
        # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
        if ( "nGene" %in% QC_metrics ){
            nfeatures = paste0( "nFeature_", assay2use)
            vln4nGene = VlnPlot(object = seurat_ob, features = "nFeature_RNA",x.lab.rot = T, alpha = alpha2use,
                            group.by = groupfactor,ncol = 1, pt.size = pt.size) + ggtitle("nGene")+labs(xlab = "") +
                            theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_total_genes4each_cell_on_violin_plot.png", collapse = "")),
                    dpi=1000, plot = vln4nGene)
            ggsave(file.path(output_dir,paste0(prefix,"_total_genes4each_cell_on_violin_plot.pdf", collapse = "")),
                    plot = vln4nGene)
            # pdf(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_PCA_plot.pdf", collapse = "")), width = 8)
            # gg_feature_pca = FeaturePlot(seurat_ob, features= nfeatures,
            #                                 reduction = "pca", cols = c("grey","red"),
            #                                 pt.size=opt$pointsize) + ggtitle("nGene")
            # ggsave(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_PCA_plot.pdf", collapse = "")), plot = gg_feature_pca)
            # ggsave(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_PCA_plot.png", collapse = "")),
            #         dpi = 1000, plot = gg_feature_pca)
            # gg_feature_tsne = FeaturePlot(seurat_ob, features = nfeatures,
            #                             reduction = "tsne", cols = c("grey","red"),
            #                             pt.size=opt$pointsize) + ggtitle("nGene")
            # ggsave(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_tSNE_plot.pdf", collapse = "")),
            # plot = gg_feature_tsne)
            # ggsave(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_tSNE_plot.png", collapse = "")),
            # dpi = 1000, plot = gg_feature_tsne)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features = nfeatures, reduction = "pca",
            #             cols = c("grey","red"),pt.size=opt$pointsize) + ggtitle("nGene")
            # dev.off()
            # pdf(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_tSNE_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features = nfeatures, reduction = "tsne",
            #                 cols = c("grey","red"),pt.size=opt$pointsize) + ggtitle("nGene")
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_nGene4each_cell_on_tSNE_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features = nfeatures, reduction = "tsne",
            #             cols = c("grey","red"),pt.size=opt$pointsize) + ggtitle("nGene")
            # dev.off()
        }

        if ( "nUMI" %in% QC_metrics ){
            nCount = paste0( "nCount_", assay2use)
            vln4nUMI = VlnPlot(object = seurat_ob, features = nCount ,alpha = alpha2use,group.by = groupfactor, ncol = 1,
                                x.lab.rot = T,pt = pt.size) + ggtitle("nUMI")+labs(xlab = "") +
                                theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_total_UMIs4each_cell_on_violin_plot.png")),
                    dpi=1000, plot = vln4nUMI)
            ggsave(file.path(output_dir,paste0(prefix,"_total_UMIs4each_cell_on_violin_plot.pdf")), plot = vln4nUMI)
            # gg_nCount_pca = FeaturePlot(seurat_ob, features= nCount, reduction = "pca", cols = c("grey","red"),
            #             pt.size=opt$pointsize) + ggtitle("nUMI")
            # ggsave(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.pdf", collapse = "")), plot = gg_nCount_pca)
            # ggsave(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.png", collapse = "")),
            #         dpi = 1000, plot = gg_nCount_pca)
            # gg_nCount_tsne = FeaturePlot(seurat_ob, features= nCount, reduction = "tsne", cols = c("grey","red"),
            #             pt.size=opt$pointsize) + ggtitle("nUMI")
            # ggsave(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_tSNE_plot.pdf", collapse = "")), plot = gg_nCount_tsne)
            # ggsave(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.png", collapse = "")),
            #         dpi = 1000, plot = gg_nCount_tsne)
            # pdf(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features= nCount, reduction = "pca", cols = c("grey","red"),
            #             pt.size=opt$pointsize) + ggtitle("nUMI")
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features= nCount, reduction = "pca", cols = c("grey","red"),
            #             pt.size=opt$pointsize) + ggtitle("nUMI")
            # dev.off()

            # pdf(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_tSNE_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features= nCount, reduction = "tsne", cols = c("grey","red"),
            #             pt.size=opt$pointsize) + ggtitle("nUMI")
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_nUMI4each_cell_on_tSNE_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features= nCount, reduction = "tsne", cols = c("grey","red"),
            #             pt.size=opt$pointsize) + ggtitle("nUMI")
            # dev.off()
        }

        # visualize to check wether there is batch effect among the samples.
        if ( "batch" %in% QC_metrics ){
            #to visulaize the metadata of each cell on the PCA plot
            # pca4batch = DimPlot(object = seurat_ob, dims=c(1,2), reduction = "pca",
            #                     label = F, pt.size = opt$pointsize, group.by = "batchid") +
            #                    theme( plot.title = element_text(hjust = 0.5))
            # ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_PCA_plot.pdf",collapse="")))
            # ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_PCA_plot.png",collapse="")),
            #                     dpi=1000, plot = pca4batch)
            tsne4batch = DimPlot(object = seurat_ob, dims=c(1,2), reduction = "tsne",
                                label = F,pt.size = opt$pointsize, group.by = "batchid")+
                                theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_tSNE_plot.pdf", collapse = "")), plot = tsne4batch)
            ggsave(file.path(output_dir,paste0(prefix,"_visualize_batch_effect_on_tSNE_plot.png", collapse = "")),
                                dpi=1000, plot = tsne4batch)

            # if ( !is.null( seurat_ob@dr$cca ) ){ #NO CCA analysis before
            #     p3 <- MetageneBicorPlot(seurat_ob, grouping.var = "batchid", dims.eval = 1:30, display.progress = FALSE)
            #     ggsave(file.path(output_dir,"metageneBicorPlot_for_CCA_batch_effect_remove.png"),dpi=1000,plot = p3)
            #     ggsave(file.path(output_dir,"metageneBicorPlot_for_CCA_batch_effect_remove.pdf"), plot = p3)
            # }
        }

        # We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
        if ( "percent.mito" %in% QC_metrics ){
		    if ( !is.null(opt$mito) ){
				mito.genes <- read.delim(opt$mito)
				raw.counts = GetAssayData(seurat_ob, slot = "counts")
				percent.mito <- Matrix::colSums(raw.counts[mito.genes[,1], ])/Matrix::colSums(raw.counts)
			}else{
				mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(seurat_ob), value = T,perl=T)

				raw.counts = GetAssayData(seurat_ob, slot = "counts")
				percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)
			}
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
            vln4mito = VlnPlot(object = seurat_ob, features = "percent.mito",
                        group.by = groupfactor, ncol = 1,alpha = alpha2use,
                        x.lab.rot = T,pt.size = pt.size)+labs(xlab= "") +
                        theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_violin_plot.png")),
                        dpi=600, plot = vln4mito )
            ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_violin_plot.pdf" )), plot = vln4mito)

            # gg_mito_pca = FeaturePlot(seurat_ob, features="percent.mito",
            #                             reduction = "pca", cols = c("grey","red"),
            #                             pt.size=opt$pointsize)
            # ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")),
            #         plot = gg_mito_pca)
            # ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")),
            #         dpi = 1000, plot = gg_mito_pca)
            # gg_mito_tsne = FeaturePlot(seurat_ob, features=c('percent.mito'),
            #                             reduction = "tsne", cols = c("grey","red"),
            #                             pt.size=opt$pointsize)
            # ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")),
            #             plot = gg_mito_tsne)
            # ggsave(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")),
            #         dpi = 1000, plot = gg_mito_tsne)
            # pdf(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features="percent.mito", reduction = "pca", cols = c("grey","red"),
            #             pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features="percent.mito", reduction = "pca", cols = c("grey","red"),
            #             pt.size=opt$pointsize)
            # dev.off()
            # pdf(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features=c('percent.mito'), reduction = "tsne", cols = c("grey","red"),
            #             pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_mitochondroin_transcript_ratio_in_each_cell_tSNE_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features=c('percent.mito'), reduction = "tsne", cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
        }

		if ( "log10GenesPerUMI" %in% QC_metrics ){
			seurat_ob@meta.data$log10GenesPerUMI <- log10(seurat_ob@meta.data$nFeature_RNA)/log10(seurat_ob@meta.data$nCount_RNA)  
            vln4log10 = VlnPlot(object = seurat_ob, features = "log10GenesPerUMI",
                        group.by = groupfactor, ncol = 1,alpha = alpha2use,
                        x.lab.rot = T,pt.size = pt.size)+labs(xlab= "") +
                        theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_log10GenesPerUMI_transcript_ratio_in_each_cell_violin_plot.png")),
                        dpi=600, plot = vln4log10 )
            ggsave(file.path(output_dir,paste0(prefix,"_log10GenesPerUMI_transcript_ratio_in_each_cell_violin_plot.pdf" )), plot = vln4log10)
		}

        # We calculate the percentage of ribosome genes here and store it in percent.ribo using AddMetaData.
        if ( "percent.ribo" %in% QC_metrics ){
            #Calculate percent ribosomal genes.
            ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", rownames(seurat_ob), value = TRUE,ignore.case = T)
            percent.ribo <- PercentageFeatureSet(seurat_ob, features = ribo.genes )
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.ribo, col.name = "percent.ribo")
            vln4ribo = VlnPlot(object = seurat_ob, features = "percent.ribo",
                                group.by = groupfactor, ncol = 1,alpha = alpha2use,
                                x.lab.rot = T,pt.size = pt.size)+labs(xlab = "") +
                                theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_violin_plot.png")),
                                dpi=600, plot = vln4ribo )
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_violin_plot.pdf" )), plot = vln4ribo)

            gg_ribo_pca = FeaturePlot(seurat_ob, features="percent.ribo", reduction = "pca",
                                    cols = c("grey","red"),pt.size=opt$pointsize)
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")),
                    plot = gg_ribo_pca)
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")),
                     dpi = 1000, plot = gg_ribo_pca)
            gg_ribo_tsne = FeaturePlot(seurat_ob, features='percent.ribo', reduction = "tsne",
                                    cols = c("grey","red"),pt.size=opt$pointsize)
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")),
                    plot = gg_ribo_tsne)
            ggsave(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_tSNE_plot.png", collapse = "")),
                    dpi = 1000, plot = gg_ribo_tsne)
            # pdf(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features="percent.ribo", reduction = "pca",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features="percent.ribo", reduction = "pca",
            #             cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # pdf(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features='percent.ribo', reduction = "tsne",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_ribosome_transcript_ratio_in_each_cell_tSNE_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features='percent.ribo', reduction = "tsne",
            #             cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
        }

        if ( "top50" %in% QC_metrics ){
            top50_gene = apply( seurat_ob@raw.data, 2, function(x) sum(x[order(x, decreasing = T)][1:50])/sum(x))
            seurat_ob <- AddMetaData(object = seurat_ob, metadata = top50_gene, col.name = "top50")
            vln4top50 = VlnPlot(object = seurat_ob, features = "top50",
                                    group.by = groupfactor, ncol = 1,alpha = alpha2use,
                                 x.lab.rot = T,pt.size = pt.size)+labs(xlab = "") +
                                    theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_violin_plot.png")),
                                dpi=600, plot = vln4top50 )
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_violin_plot.pdf" )),
                                plot = vln4top50)
            gg_top50_pca = FeaturePlot(seurat_ob, features="top50", reduction = "pca",
                                    cols = c("grey","red"),pt.size=opt$pointsize)
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_PCA_plot.png",
                    dpi = 1000,collapse = "")))
            FeaturePlot(seurat_ob, features="top50", reduction = "tsne",
                                    cols = c("grey","red"),pt.size=opt$pointsize)
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            ggsave(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_tSNE_plot.png",
                        dpi = 1000, collapse = "")))
            # pdf(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_PCA_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features="top50", reduction = "pca",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features="top50", reduction = "pca",
            #             cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # pdf(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_tSNE_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features="top50", reduction = "tsne",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_transcript4top50_gene_ratio_in_each_cell_tSNE_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features="top50", reduction = "tsne",
            #                 cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
        }

        if ( "UMI.per.gene" %in% QC_metrics ){
            tmp = seurat_ob@meta.data$nCount_RNA/seurat_ob@meta.data$nFeature_RNA
            names(tmp) = row.names(seurat_ob@meta.data)
            seurat_ob = AddMetaData(seurat_ob, metadata = tmp, "UMI.per.Gene")
            vln4UPG = VlnPlot(object = seurat_ob, features = "UMI.per.Gene",
                        group.by = groupfactor,alpha = alpha2use,
                        ncol = 1, x.lab.rot = T,pt.size = pt.size)+labs(xlab = "") +
                        theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_violin_plot.png")),
            dpi=600, plot = vln4top50 )
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_violin_plot.pdf" )),
            plot = vln4top50)
            gg_umi4gene_pca = FeaturePlot(seurat_ob, features="UMI.per.Gene", reduction = "pca",
                                    cols = c("grey","red"),pt.size=opt$pointsize)
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.pdf", collapse = "")),
                    plot = gg_umi4gene_pca)
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.png", collapse = "")),
                    dpi = 1000, plot = gg_umi4gene_pca)
            gg_umi4gene_tsne = FeaturePlot(seurat_ob, features="UMI.per.Gene", reduction = "tsne",
                                    cols = c("grey","red"),pt.size=opt$pointsize)
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.pdf", collapse = "")),
                    plot = gg_umi4gene_tsne)
            ggsave(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.png", collapse = "")),
                        dpi = 1000, plot = gg_umi4gene_tsne)
            # pdf(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features="UMI.per.Gene", reduction = "pca",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features="UMI.per.Gene", reduction = "pca",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            #
            # pdf(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.pdf", collapse = "")))
            # FeaturePlot(seurat_ob, features="UMI.per.Gene", reduction = "tsne",
            #                         cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
            # png(file.path(output_dir,paste0(prefix,"_mean_UMIs4each_gene_in_each_cell_PCA_plot.png", collapse = "")))
            # FeaturePlot(seurat_ob, features="UMI.per.Gene", reduction = "tsne",
            #                 cols = c("grey","red"),pt.size=opt$pointsize)
            # dev.off()
        }

        #visualize to check wether regress out the cell cycle effect
        if ( "cellcycle" %in% QC_metrics ){
            print("calulating cell cycle.")
            genes.inuse = rownames(GetAssayData(seurat_ob, slot="counts"))
            s.genes = CaseMatch(search = cc.genes$s.genes, match = genes.inuse)
            g2m.genes = CaseMatch(search = cc.genes$g2m.genes, match = genes.inuse)
            seurat_ob <- CellCycleScoring(object = seurat_ob, s.features = s.genes,
                                            g2m.features = g2m.genes,set.ident = T)
            # seurat_obx <- RunPCA(object = seurat_obx, features = c(s.genes, g2m.genes),
            # assay = assay2use, verbose = F)
            # pca4cellcycle = DimPlot(object = seurat_obx, label = F, pt.size = opt$pointsize)+
            #                 theme( plot.title = element_text(hjust = 0.5))
            # ggsave(file.path(output_dir,
            #             paste0(prefix, "_visualize_cell_cycle4each_cell_PCA_plot.pdf",collapse="")),
            #             plot = pca4cellcycle)
            # ggsave(file.path(output_dir,
            #         paste0(prefix, "_visualize_cell_cycle4each_cell_PCA_plot.png",collapse="")),
            #         dpi = 1000, plot = pca4cellcycle)
            #
            # seurat_obx <- RunTSNE(object = seurat_obx, dims = 1:components_num,
            # reduction = "pca", assay = assay2use,features =c(s.genes, g2m.genes),
            # perplexity=perplexity,force.recalc = T, check_duplicates = F,
            # do.fast = T, dim.embed = 2)
            # tsne4cellcycle = DimPlot(object = seurat_obx, label = F, pt.size = opt$pointsize)+
            #                     theme( plot.title = element_text(hjust = 0.5))
            # ggsave(file.path(output_dir,paste0(prefix, "_visualize_cell_cycle4each_cell_tSNE_plot.pdf",collapse="")),
            # plot = tsne4cellcycle)
            # ggsave(file.path(output_dir,paste0(prefix, "_visualize_cell_cycle4each_cell_tSNE_plot.png",collapse="")),
            #         dpi = 1000, plot = tsne4cellcycle)
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
