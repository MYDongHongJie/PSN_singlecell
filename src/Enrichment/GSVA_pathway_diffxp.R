#=================================================================================
# library loading
#=================================================================================
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))

CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
  return(my_palette[n])
}
option_list = list(
    make_option(c("--input", "-i"), type = "character",
            help = "the GSVA enrichment score matrix."),
    make_option(c("--seurat", "-v"), type = "character",
            help = "the seurat object of gene expression matrix."),
    make_option( c( "--metadata", "-m"), type = "character",
            help = "[Optional]the metadata of each sample/cell in the input file."),
    make_option( c("--contrast", "-c"),type = "character",default = NULL,
            help = "[Required]levels of a factor used to compare with for final differenetial results.
                    The format is Factor:interesting_level:reference_level."),
    make_option( c( "--fdr", "-F" ), type = "double", default = NULL,
            help = "[Optional]the adjust P value threshold."),
    make_option( c( "--pval", "-p" ), type = "double", default = 0.05,
            help = "[Optional]the P value threshold."),
    make_option( c( "--topn", "-n" ), type = "integer", default = 10,
            help = "the maximium number of ranked pathway terms on the top for each cluster "),
    make_option( c( "--matrix", "-x" ), type = "character", default = NULL,
            help = "the GSVA diffexp result matrix."),
    #make_option( c("--groupby", "-g" ), type = "character",default = "clusters",
    #        help = "metadata column, a factor used to compare "),
    make_option( c("--term", "-t" ), type = "character",default = NULL,
            help = "attention term info, can be list or dataframe, but first column must be termid and column id is geneset. "),
    make_option( c("--display_p", "-d" ), type = "character",default = TRUE,
            help = "Whether to show a significant P value"),
    make_option( c("--outdir", "-o" ), type = "character",
            help = "the output directory of the analysis results."),
    make_option( c("--cell_heatmap"), type = "character",default = FALSE,
            help = "cell heatmap."),
    make_option( c("--WHICH_GROUP", "-q" ), type = "character",
            help = "[Optional]select the groupping column in metadata used to do subsettig of cell if necessary."),                                                                          
    make_option( c("--WHICH_CELLS", "-u" ), type = "character",
            help = "[Optional] the level id list in selected cell groupping used to celltyping.For some cell
                   clusters with high hertergensity this may be useful for sub-celltyping
                   combined with the option -l/--LEVEL with single mode.
                   If not specified with cell clusters's ID, all cells will be used.")

);
opt_parser = OptionParser(option_list=option_list, usage = "Rscripit GSVA_pathway_diffxp.R -i gsva_results -v RDS -m metadata.xls -c group:variable1:variable2 -p 0.05 -n 10 -o GSVA/");
opts = parse_args(opt_parser);
#=================================================================================
# check the command line parameters
#=================================================================================

if ( is.null(opts$outdir) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opts$outdir) ){
        output_dir = opts$outdir
    }else{
        output_dir = opts$outdir
        dir.create(output_dir)
    }
}
output_dir = normalizePath(output_dir )

if ( is.null( opts$fdr ) ){
    FDR_cutoff = NULL
    if ( is.null(opts$pval) ){
        stop( "Warning:NO FDR or P value is specified!")
    }else{
        pval_cutoff = opts$pval
    }
}else{
    FDR_cutoff = opts$fdr
}

if ( is.null(opts$matrix) ){
    suppressPackageStartupMessages(library(GSEABase))
    suppressPackageStartupMessages(library(Seurat))

    if ( !is.null(opts$input )){
        #gsva_scores = read.table(opts$input, sep ="\t", header = T, row.names = 1,check.names = F, quote='"')
		gsva_scores = fread(opts$input, sep ="\t", header = T,check.names = F, quote='"') %>% as.data.frame()
        row.names(gsva_scores) = gsva_scores$geneset
        gsva_scores = gsva_scores[,-1]
    }else{
        stop("NO GSVA results is AVAILABLE!")
    }

    if ( !is.null(opts$seurat )){
        seurat_ob = readRDS( opts$seurat )
        # if the input seurat object version is less than 3, upgrade it to version 3
        if ( seurat_ob@version < 3 ){
            seurat_ob = UpdateSeuratObject(seurat_ob) # make sure the seurat object match with the latest seurat package
        }
        if ( !is.null(opts$WHICH_CELLS)){
            cluster_list = unlist(strsplit( opts$WHICH_CELLS,",",perl = T))
            seurat_ob = SubsetData(seurat_ob,  subset.name = opts$WHICH_GROUP, accept.value = cluster_list)
            gsva_scores = gsva_scores[,intersect(colnames(gsva_scores),colnames(seurat_ob))]
        }
        gsva_scores = gsva_scores[,intersect(colnames(gsva_scores),colnames(seurat_ob))]
        seurat_ob = subset(seurat_ob, cells=intersect(colnames(gsva_scores),colnames(seurat_ob)) )
    }else{
        stop("NO seurat object is AVAILABLE!")
    }

    # update the metedata in the seurat_ob@meta.data with new additional sample metadata
    if ( !is.null(opts$metadata) ){
        additional_metadata = read.csv(opts$metadata,sep=",",header =T )
        rownames(additional_metadata) = additional_metadata$sampleid
        cellnames = Cells(seurat_ob)
        sampleidx =  gsub("(_|-)[ATGC]{16,}.*","",cellnames,perl=T) # the index order is the same as the row index of the assay metadata
        # integrate the additional metadata from the assay design
        additional_cell_meta = vector( )
        for ( colidx in colnames(additional_metadata) ){
            additional_cell_meta = cbind(additional_cell_meta, as.vector(additional_metadata[sampleidx, colidx]))
        }
        colnames(additional_cell_meta) = colnames(additional_metadata)
        rownames(additional_cell_meta) = cellnames
        additional_cell_meta = as.data.frame(additional_cell_meta)
        seurat_ob = AddMetaData( seurat_ob, metadata = additional_cell_meta)
    }

    assay_metadata = seurat_ob@meta.data

    #=================================================================================
    # parse the contrast for this assay
    #=================================================================================
    # if available,the contrast string from the user must obey the right format:
    # the_interested_factor:the_interested_levels_in_this_factor:the_reference_level_in_this_factor

    if ( is.null(opts$contrast ) ){ #no contrast is provided
        stop("NO contrast string is AVAILABLE!")
    }else{
        contrast = opts$contrast
    }

    contrasts = unlist( strsplit(contrast,":",perl = T) )
    all_levels = as.vector(unique(assay_metadata[,contrasts[1]]))
    groupby = contrasts[1]
    if ( contrasts[2] == "all" & contrasts[3] != "all" ){
        case_levels = all_levels[-which(all_levels==contrasts[3])] #delete the reference level
        all_comparisions = paste(contrasts[1],case_levels,contrasts[3],sep = ":")
    }else if ( contrasts[2] != "all" & contrasts[3] == "all" ){
        ref_levels = all_levels[-which(all_levels==contrasts[2])] #delete the interested level
        all_comparisions = paste(contrasts[1],contrasts[2],ref_levels,sep = ":")
    }else if ( contrasts[2] == "all" & contrasts[3] == "all" ){
        all_comparisions = lapply(all_levels,
                    function(x) paste(contrasts[1],x,paste0(all_levels[-which(all_levels==x)],collapse = ","),sep = ":"))
        all_comparisions = unlist(all_comparisions)
    }else{
        all_comparisions = contrast
    }

    #=================================================================================
    # add the GSVA results to seurat object
    #=================================================================================
    seurat_ob[['GSVA']] = CreateAssayObject(counts = gsva_scores)
    seurat_ob = ScaleData(seurat_ob, assay = "GSVA")
    #saveRDS(seurat_ob,"seurat_GSVA.rds")


    #=================================================================================
    # differential Gene set analysis using limma function definition
    #=================================================================================
    if ( contrasts[2] == "all" & contrasts[3] == "all" ){
        # find the significant differential pathway for each Interested_group against all other Interested_group
        gsva_results = c()
        for ( contrastx in all_comparisions ){
            contrastsx = unlist(strsplit(contrastx,':',perl =T))
            DEG_gsva_tmp = FindMarkers( seurat_ob, 
                                        ident.1 = as.character(unlist( strsplit(contrastsx[2],",",perl = T) )), 
                                        ident.2 = as.character(unlist( strsplit(contrastsx[3],",",perl = T) )), 
                                        test.use = "limma", 
                                        min.pct = 0, 
                                        group.by = contrastsx[1], 
                                        logfc.threshold = -Inf, 
                                        assay = "GSVA", 
                                        slot = "data" )
            DEG_gsva_tmp$cluster = as.character(contrastsx[2])
            DEG_gsva = DEG_gsva_tmp %>% rownames_to_column(var = "geneset")
            gsva_results = rbind(gsva_results, DEG_gsva) 
        }
        gsva_results = gsva_results %>% dplyr::select(geneset, logFC, AveExpr, t, p_val, adj.P.Val, cluster)
        colnames(gsva_results) = c("geneset", "logFC", "avgExp", "t", "pval", "FDR", "Interested_group")
        write.table(gsva_results, file.path(output_dir,"GSVA_results.xls"), quote=F, sep="\t", col.names=T, row.names=F)

        # visualize topn the expressed pathway using heatmap
        
        # 如果seurat里的metadata列是因子类型，则使用seurat里的顺序，如果不是，则转换成因子型
        if (is.factor(seurat_ob@meta.data[[contrasts[1]]])) {
            gsva_results$Interested_group = factor(gsva_results$Interested_group,levels = levels(seurat_ob@meta.data[[contrasts[1]]]))
        } else {
            gsva_results$Interested_group = as.factor(gsva_results$Interested_group)
        }

        if ( is.null(FDR_cutoff) ){
            plot_term = gsva_results %>% subset( pval < pval_cutoff&t > 0 ) %>% 
                group_by( Interested_group ) %>% top_n(opts$topn,t) %>% arrange(Interested_group)
                
        }else{
            plot_term = gsva_results %>% subset( FDR < FDR_cutoff&t > 0 ) %>% 
                group_by( Interested_group ) %>% top_n(opts$topn,t) %>%  arrange(Interested_group)
        }
        write.table(plot_term, file.path(output_dir,paste0("GSVA_top",opts$topn,"_results.xls")), 
                    quote=F, sep="\t", col.names=T, row.names=F)
        #准备差异显著的P值加*号文件 
        if(opts$display_p){
            if ( is.null(FDR_cutoff) ){
                plot_term_p = gsva_results %>% subset(geneset %in% unique(plot_term$geneset)) %>% subset( pval < pval_cutoff ) %>%  group_by( Interested_group )  %>% arrange(Interested_group)
            }else{
                plot_term_p = gsva_results %>% subset(geneset %in% unique(plot_term$geneset)) %>% subset( FDR < FDR_cutoff ) %>%  group_by( Interested_group ) %>%  arrange(Interested_group)
            }

            plot_data_p = spread(as.data.frame(plot_term_p[,c("geneset","Interested_group","t")]),Interested_group,t) %>% column_to_rownames(var = "geneset")
            plot_data_p = plot_data_p[as.vector(unique(plot_term$geneset)),]
            display_numbers = ifelse(is.na(plot_data_p) , '' , '*')
        }else{
            display_numbers = FALSE
        }
        
        plot_data = spread(as.data.frame(gsva_results[,c("geneset","Interested_group","t")]),Interested_group,t) %>% 
                        column_to_rownames(var = "geneset")
        plot_data = plot_data[as.vector(unique(plot_term$geneset)),]
        p = pheatmap( plot_data,
                    scale = "row",
                    color = colorRampPalette(c("#406AA8", "white","#D91216"))(100),
                    show_colnames = T,
                    cluster_cols = F,
                    cluster_rows = F,
                    border_color = "white",
                    fontsize_row = 8, cellheight = 12, cellwidth = 12, display_numbers = display_numbers)
        ggsave(file.path(output_dir, paste0("top", opts$topn, "_gsva_term.pdf")), plot=p, width=9, height=7*opts$topn*0.1*length(unique(plot_term$Interested_group))/3)
        ggsave(file.path(output_dir, paste0("top", opts$topn, "_gsva_term.png")), plot=p, width=9, height=7*opts$topn*0.1*length(unique(plot_term$Interested_group))/3, dpi=1000)
        write.table(as.data.frame(plot_data) %>% rownames_to_column(var = "geneset"),
                    file.path(output_dir,paste0("top", opts$topn, "_term_t_value.xls")),quote=F,col.names=T,row.names=F,sep="\t")
        if( as.logical(opts$cell_heatmap)){
            cell_order = rownames(assay_metadata[order(assay_metadata[,groupby]), ])
            color_map = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
            names=unique(seurat_ob@meta.data[,groupby])
            names(color_map) <- names
            color_use = list()
            color_use[[groupby]]=color_map
            groupby_anno = assay_metadata[order(assay_metadata[,groupby]), groupby, drop = FALSE]
            # print(head(cell_plotdata))
            if ( !is.null(opts$term )){
                plot_term <- read.delim(opts$term,sep="\t",header=T)
                colnames(plot_term)[1]="geneset"
                print('有输入的通路，使用输入通路绘图')
                print(plot_term)
                plot_term$geneset = intersect(plot_term$geneset,rownames(seurat_ob[['GSVA']]@scale.data))
                plot_term$Interested_group='x'
                print(plot_term)
            }
            cell_plotdata = seurat_ob[['GSVA']]@scale.data[as.vector(unique(plot_term$geneset)),cell_order]
            bks <- unique(c(seq(-2.5,0, length=50),  seq(0,2.5, length=50)))
            # palette <- colorRampPalette(c("#163E69","#009BC3", "white", "#B7181F","#8A0F1E"))(n=300)
            # palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(300)
            p = pheatmap(cell_plotdata,
                    # color=palette,
                    # cex=1,
                    # border=F,
                    # angle_col=45,
                    # # treeheight_row=36, treeheight_col=36,
                    # lwd=1.1,
                    # scale=ifelse(dim(plotdata)[2]==2,"none","row"),
                    breaks=bks,
                    show_colnames=FALSE,
                    cluster_col=FALSE,
                    cluster_row=FALSE,
                    annotation_colors = color_use,
                    # gaps_row=gaps_row,
                    annotation_col = groupby_anno)
                    # ,
                    # annotation_row = groupby2)
            ggsave(file.path(output_dir, paste0("top", opts$topn, "_gsva_term_cellheatmap.pdf")), plot=p, width=(9+max(nchar(rownames(cell_plotdata)))*0.1), height=7*opts$topn*0.1*length(unique(plot_term$Interested_group))/3)
            # ggsave(file.path(output_dir, paste0("top", opts$topn, "_gsva_term_cellheatmap.png")), plot=p, width=(9+max(nchar(rownames(cell_plotdata)))*0.1), height=7*opts$topn*0.1*length(unique(plot_term$Interested_group))/3, dpi=1000)
            print("Convert pdf to png...")
            plotname=paste0("top", opts$topn, "_gsva_term_cellheatmap")
            setwd(output_dir)
            system(paste0("module purge && module load OESingleCell/3.0.d && /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 800 -trim  ",plotname,".pdf  -quality 500  -flatten   ",plotname,".png "))
            write.table(as.data.frame(cell_plotdata) %>% rownames_to_column(var = "geneset"),
                        file.path(output_dir,paste0("top", opts$topn, "_term_cell_score.xls")),quote=F,col.names=T,row.names=F,sep="\t")
        }
    }else{
        for ( contrastx in all_comparisions ){
            contrastsx = unlist(strsplit(contrastx,':',perl =T))
            pathway_diffexp = FindMarkers( seurat_ob, 
                                        ident.1 = as.character(unlist( strsplit(contrastsx[2],",",perl = T) )), 
                                        ident.2 = as.character(unlist( strsplit(contrastsx[3],",",perl = T) )), 
                                        test.use = "limma", 
                                        min.pct = 0, 
                                        group.by = contrastsx[1], 
                                        logfc.threshold = -Inf, 
                                        assay = "GSVA", 
                                        slot = "data" )
            pathway_diffexp = pathway_diffexp %>% rownames_to_column(var = "geneset") %>% 
                            dplyr::select(geneset, logFC, AveExpr, t, p_val, adj.P.Val)
            colnames(pathway_diffexp) = c("geneset", "logFC", "avgExp", "t", "pval", "FDR")
            write.table(pathway_diffexp,file.path(output_dir,
                        paste0("diffexp_genesets_GSVA_score4",contrastsx[1],"_",contrastsx[2],
                        "-vs-",contrastsx[3],".xls")), quote = F,sep = '\t',row.names = F)
            
            # visualize all the expressed pathway using divergent barplot with text
            # labels on both directions of the y axis
            pathway_diffexp$just = ifelse( pathway_diffexp$t<0,0,1)
            pathway_diffexp$is.just = pathway_diffexp$just==1 
            if ( is.null(FDR_cutoff)){
                pathway_diffexp$is.sig = pathway_diffexp$pval<pval_cutoff
            }else{
                pathway_diffexp$is.sig = pathway_diffexp$FDR<FDR_cutoff
            }
            # pick equal number of pathways in both up regulated pathway
            # and down regulated pathway to visualize in case of too many
            # pathways
            pathway2vis = pathway_diffexp %>% filter( pval < pval_cutoff ) %>% group_by( is.just ) %>% arrange(abs(t)) %>% top_n(opts$topn,abs(t))

            pp = ggplot(pathway2vis, aes(reorder(geneset, t), t)) +
                    geom_col(aes(fill=is.just)) +
                    scale_fill_manual(values=c( "#6CC570","#2A5078")) +
                    coord_flip() +
                    labs(x="Pathway", y="t value of GSVA" ) +
                    theme_minimal() +
                    geom_text( aes(x= geneset, y=0, label = geneset), hjust = pathway2vis$just, size = 3.5 )+
                    theme(axis.text.y=element_blank()) +
                    theme(panel.grid =element_blank())+ylim(c(-max(abs(pathway2vis$t))-5,max(abs(pathway2vis$t))+5))+
                    labs( fill = paste0("t value > 0") )
            #theme(panel.grid =element_blank())+ylim(c(min(pathway2vis$t)-5,max(pathway2vis$t)+5))+
            max.nchar = max(apply(pathway2vis["geneset"],1,nchar))
            if(grepl("GO", pathway2vis$geneset[1])){
                ggsave(file.path(output_dir,
                    paste0("diffexp_genesets_GSVA_score4",contrastsx[1],"_",contrastsx[2],"-vs-",contrastsx[3],"_barplot.pdf")),
                    width = max.nchar/6, height = 8, plot = pp)
                ggsave(file.path(output_dir,
                    paste0("diffexp_genesets_GSVA_score4",contrastsx[1],"_",contrastsx[2],"-vs-",contrastsx[3],"_barplot.png")),
                    width = max.nchar/6, height = 8, plot = pp) }else{
                        ggsave(file.path(output_dir,
                            paste0("diffexp_genesets_GSVA_score4",contrastsx[1],"_",contrastsx[2],"-vs-",contrastsx[3],"_barplot.pdf")),
                            width = max.nchar/5, height = 8, plot = pp)
                        ggsave(file.path(output_dir,
                            paste0("diffexp_genesets_GSVA_score4",contrastsx[1],"_",contrastsx[2],"-vs-",contrastsx[3],"_barplot.png")),
                            width = max.nchar/5, height = 8, plot = pp) 
                    }
            
        }
    }

    if(!file.exists(file.path(output_dir, "GSVA分析说明.docx"))){
    file.copy("/public/scRNA_works/Documents/GSVA分析说明.docx",
    file.path(output_dir, "GSVA分析说明.docx"))
    }

}else{
    gsva_results <- read.delim(opts$matrix,sep="\t",header=T)
    if ( !is.null(opts$term )){
        plot_term <- read.delim(opts$term,sep="\t",header=T)
    }else{
        stop("NO attention termid  is AVAILABLE!")
    }
    if( "Interested_group" %in% colnames(gsva_results) ) {
        #if ( !is.null(opts$groupby )){
        #    contrasts = opts$groupby
        #}else{
        #    contrasts = "clusters"
        #}
        #if ( contrasts == "clusters" ){
        #    gsva_results$Interested_group = factor(gsva_results$Interested_group, 
        #                                    levels = sort( unique(as.numeric( gsva_results$Interested_group))))
        #}else{
        #    gsva_results$Interested_group = as.factor(gsva_results$Interested_group)
        #}
        if ("FALSE" %in% grepl("^\\d+\\.\\d+$|^\\d+$", unique(gsva_results$Interested_group))){
            print("The Interested_group is character type")
            gsva_results$Interested_group = as.factor(gsva_results$Interested_group)
        }else{
            print("The Interested_group is numeric type")
            gsva_results$Interested_group = factor(gsva_results$Interested_group, 
                                levels = sort( unique(as.numeric( gsva_results$Interested_group))))
        }
        ##准备差异显著的P值加*号文件        
        if(opts$display_p){
            if ( is.null(FDR_cutoff) ){
                plot_term_p = gsva_results %>% subset(geneset %in% unique(plot_term$geneset)) %>% subset( pval < pval_cutoff ) %>% group_by( Interested_group ) %>% arrange(Interested_group)
            }else{
                plot_term_p = gsva_results %>% subset(geneset %in% unique(plot_term$geneset)) %>% subset( FDR < FDR_cutoff ) %>%  group_by( Interested_group ) %>%  arrange(Interested_group)
            }
            
            plot_data_p = spread(as.data.frame(plot_term_p[,c("geneset","Interested_group","t")]),Interested_group,t) %>% column_to_rownames(var = "geneset")
            plot_data_p = plot_data_p[as.vector(unique(plot_term$geneset)),]
            display_numbers = ifelse(is.na(plot_data_p) , ' ' , '*')
        }else{
            display_numbers = FALSE
        }
        
        plot_data = spread(as.data.frame(gsva_results[,c("geneset","Interested_group","t")]),Interested_group,t) %>%
                        column_to_rownames(var = "geneset") 
        cell_length=ifelse(length(unique(plot_term$geneset)) > 20, 12, 20 )
        height=ifelse(length(unique(plot_term$geneset)) > 60, 7*log2(length(unique(plot_term$geneset)))*0.4, 7*log2(length(unique(plot_term$geneset)))*0.3 )
        plot_data = plot_data[as.vector(unique(plot_term$geneset)),]

        p = pheatmap( plot_data,
                    scale = "row",
                    color = colorRampPalette(c("#406AA8", "white","#D91216"))(100),
                    show_colnames = T,
                    cluster_cols = F,
                    cluster_rows = F,
                    border_color = "white",
                    fontsize_row = 8, cellheight = cell_length , cellwidth = cell_length, display_numbers = display_numbers )

        ggsave(file.path(output_dir, "sig_gsva_geneset.pdf"), width=8, height=height , plot=p,limitsize = FALSE)
        ggsave(file.path(output_dir, "sig_gsva_geneset.png"), width=8, height=height,  plot=p, dpi=1000,limitsize = FALSE)

        write.table(as.data.frame(plot_data) %>% rownames_to_column(var = "geneset"),
        file.path(output_dir,paste0("sig_term_t_value.xls")),quote=F,col.names=T,row.names=F,sep="\t")
    }else{
        pathway_diffexp <- gsva_results
        pathway_diffexp <- pathway_diffexp %>% filter( geneset %in% as.vector(unique(plot_term$geneset)) )
        # visualize all the expressed pathway using divergent barplot with text
        # labels on both directions of the y axis
        pathway_diffexp$just = ifelse( pathway_diffexp$t<0,0,1)
        pathway_diffexp$is.just = pathway_diffexp$just==1 
        # pick equal number of pathways in both up regulated pathway
        # and down regulated pathway to visualize in case of too many
        # pathways
        pathway2vis = pathway_diffexp %>% group_by( is.just ) %>% arrange(abs(t))
        pp = ggplot(pathway2vis, aes(reorder(geneset, t), t)) +
            geom_col(aes(fill=is.just)) +
            scale_fill_manual(values=c( "#6CC570","#2A5078")) +
            coord_flip() +
            labs(x="Pathway", y="t value of GSVA" ) +
            theme_minimal() +
            geom_text( aes(x= geneset, y=0, label = geneset), hjust = pathway2vis$just, size = 3.5 )+
            theme(axis.text.y=element_blank()) +
            theme(panel.grid =element_blank())+ylim(c(-max(abs(pathway2vis$t))-5,max(abs(pathway2vis$t))+5))+
            labs( fill = paste0("t value > 0") )
            #theme(panel.grid =element_blank())+ylim(c(min(pathway2vis$t)-5,max(pathway2vis$t)+5))+
        max.nchar = max(apply(pathway2vis["geneset"],1,nchar))
        contrasts_vs = gsub(".xls","",basename(opts$matrix))
        if(grepl("GO", pathway2vis$geneset[1])){
            ggsave(file.path(output_dir,
                paste0("sig_",contrasts_vs,"_barplot.pdf")),
                width = max.nchar/6, height = 8, plot = pp)
            ggsave(file.path(output_dir,
                paste0("sig_",contrasts_vs,"_barplot.png")),
                width = max.nchar/6, height = 8, plot = pp)
        }else{
            ggsave(file.path(output_dir,
                paste0("sig_",contrasts_vs,"_barplot.pdf")),
                width = max.nchar/5, height = 8, plot = pp)
            ggsave(file.path(output_dir,
                paste0("sig_",contrasts_vs,"_barplot.png")),
                width = max.nchar/5, height = 8, plot = pp) 
        }
    }
}
