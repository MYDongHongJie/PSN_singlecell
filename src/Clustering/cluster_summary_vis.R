#!/usr/bin/env Rscript
# This script is mainly used to visualize the summary or annotations of cells derived
# from the previous clustering and QC in different plot types, including: barplot,
# pie plot, boxplot and feature plot

rm(list=ls())
#========import packages=====================================
suppressWarnings({
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("docopt"))
    suppressPackageStartupMessages(library("OESingleCell"))
    suppressPackageStartupMessages(library("future"))
    suppressPackageStartupMessages(library("dplyr"))
    suppressPackageStartupMessages(library("ggplot2"))
})

# command line parameters settings using docopt
'usage:cluster_summary_vis.R
    cluster_summary_vis.R elbow [options]
    cluster_summary_vis.R jackstraw [options]
    cluster_summary_vis.R dimvis [options]
    cluster_summary_vis.R summarize [options]
    cluster_summary_vis.R visfeature [options]

options:
    global_options:
        -i <input>, --input <input>                   # The input exprssion matrix in several possible format.
                                                         Only Seurat object is supported currently.
        -o <outdir>, --outdir <outdir>                # the output directory of results.[default: .]
        -d <reduct>, --reduct <reduct>                # the dimension reduction result to use in this run.[default: tsne]
        -t <idents>, --idents <idents>                # [OPTIONAL]The column name in cell metadata used as identity
                                                         of each cell combined with --cell2use.[default NULL]
        -u <cell2use>, --cell2use <cell2use>          # [OPTIONAL]The subset of group/cluster ids used for visualizing.[default:NULL]
        -b <splitby>, --splitby <splitby>             # [OPTIONAL]visualize cells in seperate plot split by this groupping variable, only 2D supported.[default: NULL]
        --assay <assay>                               # [OPTIONAL]the assay to use in case of multimodal data.[default: RNA]
        --slot <slot>                                 # [OPTIONAL]Slot to pull data from, should be one of counts, data, or scale.data.[default: data]
        --from <from>                                 # [OPTIONAL]the original cluster id used in calibrate the clustering results.[default: NULL]
        --ncores <ncores>                             # [OPTIONAL]the number of CPUs used to improve the performace.[default: 10]
        --to <to>                                     # [OPTIONAL]the adjusted cluster id used in calibrate the clustering results.[default: NULL]
    subcommand:
        elbow                                         # run elbow to evaluate the dimension reduction result for the optimal dimensions
        jackstraw                                     # run jackstraw to evaluate the dimension reduction result for the significant dimensions
            options:
                --components <components>             # [OPTIONAL]the dimensions used for visualization or calibration.[default: 20]
                --overwrite <overwrite>               # [OPTIONAL]overwrite the input data object.
        visfeature                                    # visualize the cell annotation on the specified reduction feature plot
            options:
                --features <features>                 # the feature list seperated by comma.
                --feature1 <feature1>                 # [OPTIONAL] First feature to plot. Typically feature expression but can also be metrics,
                                                         PC scores, etc. - anything that can be retreived with FetchData
                --feature2 <feature2>                 # [OPTIONAL] Second feature to plot.
                --feature3 <feature3>                 # [OPTIONAL] third feature to plot.
        dimvis                                         # visualize the dimension reduction and clustering results.
            options:
                -r <resolution>, --resolution <resolution>  # [OPTIONAL]resolution used in cell clustering.
                -c <groupby>, --groupby <groupby>             # [OPTIONAL]visualize cell metadata by coloring cells
                                                                 in different color according to cell clustering metadata.[default: clusters]
                -s <ptsize>, --ptsize <ptsize>                # [OPTIONAL]the point size in the plot.[default: 0.1]
                -k <dims>, --dims <dims>                      # [OPTIONAL]the 2D/3D of the plot.Currently 3D plot is only support in plotly style.[default: 2]
                --dosummary <summary>                         # wether to plot the summary statistics of cells clustering.[default: T]
        summarize
            options:
                --plot <plottype>                       # [OPTIONAL]the plot type used to summarize the cell clustering results.Options can be:
                                                                  barplot, pie, boxplot.[default: barplot]
                --propby <propby>                       # [OPTIONAL]calculate the proportion of cell in the specified groups.[default: group]
                --groups <groupby>                      # [OPTIONAL]visualize cell metadata by coloring cells
                                                            in different color according to cell clustering metadata.[default: clusters]
' -> doc
opt = docopt(doc)

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}

#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$outdir) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    }else{
        output_dir = opt$outdir
        dir.create(output_dir, recursive = T)
    }
}

output_dir = normalizePath(output_dir)
######################################################################
if ( !is.null(opt$input) ){
    seurat_ob = readRDSMC( opt$input, cores = availableCores() )
    if ( seurat_ob@version < 3){
        #make sure the seurat object work with the latest seurat package
        seurat_ob = UpdateSeuratObject(seurat_ob)
    }
    DefaultAssay(seurat_ob) = opt$assay # setting the default assay in this run

    # get or set the clusters actually in use
    if ( is.null(seurat_ob@meta.data$clusters) ){
        seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
    }else{
        seurat_ob = SetIdent( seurat_ob, value = "clusters")
    }

    # extract the preserved resolution setting
    if ( is.null(opt$resolution) ){
        findcluster_record = Command(seurat_ob, command = "FindClusters")
        resolution = findcluster_record$resolution
    }else{
        resolution = opt$resolution
    }

    if ( is.null(opt$idents) ){
        print("NO cell identity column name AVAILABLE! The clusters column will be used as default.")
        ident2use = "clusters"
    }else{
        ident2use = opt$idents
    }

    if ( !is.null( opt$splitby) & opt$splitby != "NULL" ){
        facetby = unlist( strsplit( opt$splitby, ",", perl =T ) )
    }else{
        facetby = NULL
    }

    seurat_ob = StashIdent(seurat_ob, save.name = "clusters")
    cell_count_tsne = table(Idents(seurat_ob))
    tsne_cell_count_labels = paste(paste(names(cell_count_tsne),cell_count_tsne,sep="-")," cells")

    if ( !is.null(opt$from) & !is.null(opt$to) ){
        from_ident = unlist( strsplit( opt$from, ",", perl =T ) )
        to_ident = unlist( strsplit(opt$to, ",", perl = T ) )
        Idents(seurat_ob) = plyr::mapvalues(x = Idents(seurat_ob),
                                            from = from_ident, to = to_ident )
    }

    #get the subset of cells used for visualization
    if ( !is.null(opt$which_cells)){
        cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
        seurat_ob = SubsetData(seurat_ob,cells =
                           OldWhichCells( object,
                                        subset.name= ident2use,
                                        accept.value = cluster_list))
    }

    # =================================================================================
    # run each subcommand from here
    # =================================================================================
    # to visulaize the metadata of each cell on the tsne plot
    if ( opt$elbow ){
        ggelb = ElbowPlot2(seurat_ob)
        ggsave(file.path(output_dir,"Elbow_plot.pdf"),plot=ggelb)
        ggsave(file.path(output_dir,"Elbow_plot.png"),plot=ggelb, dpi = 600)
        message("the number of optimal pc used is: ",seurat_ob@misc$optimal_pc)
        quit()
    }
    if ( opt$jackstraw ){
        # setting the cores for parallization
        suppressPackageStartupMessages(library("future"))
        options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
        plan("multicore", workers = min(availableCores(), opt$ncores)) # parallization using specified CPUs start from here
        components2use = min( opt$components, ncol(Embeddings(seurat_ob[[opt$reduct]])))
        seurat_ob <- JackStraw(seurat_ob , num.replicate = 100,
                               reduction = opt$reduct, dims = components2use)
        seurat_ob <- ScoreJackStraw(seurat_ob, dims = components2use, reduction = opt$reduct)
        ggjack = JackStrawPlot(object = seurat_ob, dims = 1:components2use)
        ggsave(file.path(output_dir,"jackstraw_plot.pdf"),plot=ggjack)
        ggsave(file.path(output_dir,"jackstraw_plot.png"),plot=ggjack, dpi = 600)
        if ( opt$overwrite ){
            saveRDSMC(seurat_ob,opt$input )
        }else{
            saveRDSMC( seurat_ob, file.path( output_dir, basename(opt$input)) )
        }
       quit()
    }

    if ( opt$dimvis ){
		if (dim(seurat_ob)[2] < 500){
			pointsize=1.5
		}else{
			pointsize=as.numeric(opt$ptsize)
		}
        if ( is.null( opt$groupby ) ){
            print("The groupping information is not specified. The sampleid will be used for plot")
            cellfeature2vis = "clusters"
        }else{
            cellfeature2vis = unlist(strsplit(opt$groupby, ",", perl = T))
        }
        for (groupfactor in cellfeature2vis ){
            output_dir = file.path(output_dir, paste0( "visualize_cluster_by_", groupfactor, collapse = ""))
            if ( !file.exists(output_dir) ){
                dir.create(output_dir, recursive = T)
            }
            color_counter = 0
            if ( !is.null(facetby) ){
               for ( facetbyx in facetby ){
                   if ( !is.null(facetbyx) ){
                       nrow = ceiling(length(unique(seurat_ob@meta.data[,facetbyx]))/2)
                   }
                   nlevel = length(unique(seurat_ob@meta.data[,groupfactor]))
                   color_counter = nlevel
                   seurat_ob = SetIdent( seurat_ob, value = groupfactor)
                   groupvis_split = DimPlot(object = seurat_ob,
                                    dims = c(1,2),reduction = opt$reduct,
                                    pt.size = pointsize, ncol = 2,
                                    group.by = groupfactor, split.by = facetbyx)+
                       theme( plot.title = element_text(hjust = 0.5))
                   groupvis_split = groupvis_split + scale_colour_manual( values = CustomCol2(1:nlevel))
                   ggsave(file.path(output_dir,paste0("splitby-",facetbyx,"_resolution",resolution,"_split_plot.pdf",collapse="")),
                   limitsize = FALSE, plot = groupvis_split, height = 6*nrow, width = 14)
                   ggsave(file.path(output_dir,paste0("splitby-",facetbyx,"_resolution",resolution,"_split_plot.png",collapse="")),
                   limitsize = FALSE, plot = groupvis_split, width = 10, height = 4*nrow )


                   nfacet = length(unique(seurat_ob@meta.data[,facetbyx]))
                   groupvis_all = DimPlot(object = seurat_ob, dims = c(1,2),reduction = opt$reduct,
                   pt.size = pointsize, group.by = facetbyx)+
                       theme( plot.title = element_text(hjust = 0.5)) +
                       scale_colour_manual( values = CustomCol2(1:nfacet))
                   ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution",resolution,"_contrast_plot.pdf",collapse="")),
                   limitsize = FALSE, plot = groupvis_all, height = 7, width = 7)
                   ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution",resolution,"_contrast_plot.png",collapse="")),
                   limitsize = FALSE, plot = groupvis_all, width = 7, height = 7 )
                   color_counter = color_counter + nfacet

               }
            }else{
                if ( as.numeric(opt$dims) == 2 ){
                    groupvis = DimPlot(object = seurat_ob, dims = c(1,2),reduction = opt$reduct,
                    pt.size = pointsize, group.by = groupfactor)+
                        theme( plot.title = element_text(hjust = 0.5)) +
                        #scale_colour_manual( values = CustomCol2(nlevel:nlevel*2))
                        scale_colour_manual( values = scale_colour_manual( values = CustomCol2(1:nlevel))) ###keep the same color whether split.by or not
                    ggsave(file.path(output_dir,paste0("groupby-",groupfactor,"_resolution",resolution,"_contrast_plot.pdf",collapse="")),
                    limitsize = FALSE, plot = groupvis, height = 10, width = 12)
                    ggsave(file.path(output_dir,paste0("groupby-",groupfactor,"_resolution",resolution,"_contrast_plot.png",collapse="")),
                    limitsize = FALSE, plot = groupvis, width = 15, height = 10 )
                }else{
                    suppressPackageStartupMessages(library("SeuratPlotly"))
                    groupvis = DimPlotly3d(seurat_ob, reduction= opt$reduct,
                                        pt_size = pointsize,
                                        grouping_var = groupfactor,
                                        plot_grid = F)
                    htmlwidgets::saveWidget(groupvis, file = file.path(output_dir,paste0("groupby-",groupfactor,"_resolution",resolution,"_contrast_3D_plot.html",collapse="")))
                }
            }
        }

        if ( opt$dosummary ){
            seurat_ob = SetIdent( seurat_ob, value = opt$groupby )
            DATA <- as.data.frame( seurat_ob@meta.data[,c("sampleid", opt$groupby)] ) %>%
                group_by( .dots= c("sampleid", opt$groupby)) %>%
                dplyr::summarize(cell_number = n()) %>%
                mutate(freq = (cell_number / sum(cell_number)) * 100)
            write.table(as.data.frame(DATA),
                    file.path(output_dir,file="clust_cond_freq_info.xls"),
                    sep="\t",col.names=T, row.names =F)
            # visulaize the summary statistics of cells clusters in each groupping factor
            if ( !is.null(facetby) ){
               for ( facetbyx in facetby ){
                    if ( length(unique(seurat_ob@meta.data[[facetbyx]]))  != length(levels(seurat_ob@meta.data[[facetbyx]])) ) {
                        seurat_ob@meta.data[[facetbyx]] <- factor( seurat_ob@meta.data[[facetbyx]],levels=unique(sort(seurat_ob@meta.data[[ facetbyx]] )) )
                    }
                    clust_sum_all = PlotAbundances(seurat_ob, prop.by = opt$groupby , group.by = facetbyx, method = "barplot",
                    cols= CustomCol2(1:length(unique(Idents(seurat_ob))))
                    )
                    ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution-",resolution,"_summary_plot.pdf",collapse="")),plot=clust_sum_all)
                    ggsave(file.path(output_dir,paste0("groupby-",facetbyx,"_resolution-",resolution,"_summary_plot.png",collapse="")),dpi=1000, plot = clust_sum_all)
                }
            }
            seurat_ob = SetIdent(seurat_ob, value = "sampleid")
            clust_sum_all2 = PlotAbundances(seurat_ob, prop.by = "sampleid", group.by = opt$groupby, method = "barplot",
            cols= CustomCol2(1:length(unique(Idents(seurat_ob))))
            )
            ggsave(file.path(output_dir,paste0("groupby-",opt$groupby,"_resolution-",resolution,"_summary_plot.pdf",collapse="")),plot=clust_sum_all2)
            ggsave(file.path(output_dir,paste0("groupby-",opt$groupby,"_resolution-",resolution,"_summary_plot.png",collapse="")),dpi=1000, plot = clust_sum_all2)
        }
        quit()
    }

        # summarize the cell type composition in each sample and group
    if ( opt$summarize ){
        DATA <- as.data.frame( seurat_ob@meta.data[,c(opt$groups, opt$propby)] ) %>%
                            group_by( .dots= c(opt$groups, opt$propby)) %>%
                            dplyr::summarize(cell_number = n()) %>%
                            mutate(freq = (cell_number / sum(cell_number)) * 100)
        write.table(as.data.frame(DATA),
                    file.path(output_dir,file="clust_cond_freq_info.xls"),
                    sep="\t",col.names=T, row.names =F)
        clust_sum_all = PlotAbundances(seurat_ob, prop.by = opt$propby, group.by = opt$groups,
                                split.by = facetby, method = opt$plot, ncol = ifelse(opt$plot=="pie",4,1),
                                cols= CustomCol2(1:length(unique(seurat_ob@meta.data[,opt$propby])))
                                )
        ggsave(file.path(output_dir,paste0("groupby-",opt$groups,"_resolution-",resolution,"summary_plot.pdf",collapse="")),plot=clust_sum_all)
        ggsave(file.path(output_dir,paste0("groupby-",opt$groups,"_resolution-",resolution,"summary_plot.png",collapse="")),dpi=1000, plot = clust_sum_all)
        quit()
    }

    # visulaize the feature of cell on the dimension reduction plot, the feature can be anything pulled by FetechData
    if ( opt$visfeature ){
        DefaultAssay(seurat_ob) = opt$assay
        feature_list = unlist(strsplit(opt$features, ",") )
        for ( featurex in feature_list ){
            ggfeature = FeaturePlot(seurat_ob, features = featurex,cols = c("grey","red"),
                                    split.by = facetby, reduction= reduct,
                                    ncol = 1, pt.size = as.numeric(opt$ptsize)) +
                                    theme( plot.title = element_text(hjust = 0.5))
            ggsave(file.path(output_dir, paste0(featurex, "_visulazation.pdf")),plot=ggfeature)
            ggsave(file.path(output_dir, paste0(featurex, "_visulazation.png")),dpi = 600, plot=ggfeature)
        }

        if ( !is.null(opt$feature1) & !is.null(opt$feature2) ){
            ggscater = FeatureScatter(seurat_ob, feature1 = opt$feature1, slot = opt$slot,
                                    feature2 = opt$feature2, group.by = opt$groupby )
            ggsave(file.path(output_dir, paste0(feature1,"_VS_",feature2, "_visulazation.pdf")),plot=ggscater)
            ggsave(file.path(output_dir, paste0(feature1,"_VS_",feature2, "_visulazation.png")),dpi = 600, plot=ggscater)
        }
        quit()
    }
}
