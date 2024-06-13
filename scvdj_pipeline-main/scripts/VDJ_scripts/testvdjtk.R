# module load immunoinformatics/1.0.0
suppressPackageStartupMessages( library("immunarch") )
suppressPackageStartupMessages( library("ggplot2") )
suppressPackageStartupMessages( library("optparse"))

option_list = list(
    make_option( c( "--inputdir","-i"), type = "character",default = "result/cellranger",
        help = "the input dir."),
    make_option( c("--output","-o"),type="character", default = "result/VDJ_aggr",
        help="the output directory of Clustering results." ),
    make_option( c("--vdj","-t"),type="character", default = "TCR",
        help="type of vdj info." ),
    make_option( c("--metadata","-m"),type="character", default = "config/samples.csv",
        help="the sample metadata which must include sample id in this assay design.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###########################参数读取########
output_dir = opt$output
type = opt$vdj
#output_dir = file.path(rawdir,type)
metadata = opt$metadata
input_path = normalizePath(opt$inputdir)


if ( type == "TCR" ){
    type_path = "vdj_t"
}else if( type == "BCR"){
    type_path = "vdj_b"
}else{
    stop("NO right vdj type found, please check vdj type whether TCR or BCR")
}


assay_metadata = data.table::fread(metadata,header =T )
assay_metadata = assay_metadata %>% rename(Sample = sampleid ) %>% as.data.frame()
rownames(assay_metadata) = assay_metadata$Sample

#input_path = normalizePath(file.path(dir,"Cellranger"))

immdata = list() #an empty immdata list
# vdj_profile = Sys.glob(file.path(input_path,assay_metadata$Sample,"outs","filtered_contig_annotations.csv"))
vdj_profile = Sys.glob(file.path(input_path,assay_metadata$Sample,"outs","per_sample_outs",assay_metadata$Sample,type_path,"filtered_contig_annotations.csv"))

immdata$data = repLoad(vdj_profile, .format = "10x")


### summary

#opt$groupby="Sample,group"
if ( !is.null(opt$groupby) ){
    groupby = unlist(strsplit( opt$groupby, "," ))
}

if ( length(assay_metadata$Sample) > 1  ){
		summ_outdir = file.path( output_dir, "summary" )
		if ( !file.exists( summ_outdir ) ){
				dir.create( summ_outdir )
		}
    names(immdata$data) = assay_metadata$Sample
    immdata$meta = assay_metadata
    exp_cnt = repExplore(immdata$data, .method = "count")
    # for (group in groupby){
    #     p2 = vis(exp_cnt, .by = group, .meta = immdata$meta) + theme(panel.grid =element_blank()) + labs(color = group)
    #     ggsave( file.path(summ_outdir, paste0("clonotypes_count_rank_aboundance_groupby_",group, ".pdf")), plot = p2)
    #     ggsave( file.path(summ_outdir, paste0("clonotypes_count_rank_aboundance_groupby_",group, ".png")), dpi = 600, plot = p2)
    # }
    groupby="Sample"
    p2 = vis(exp_cnt, .by = groupby, .meta = immdata$meta) + theme(panel.grid =element_blank()) + labs(color = groupby)
    ggsave( file.path(summ_outdir, paste0("clonotypes_count_rank_aboundance_groupby_",groupby, ".pdf")), plot = p2)
    ggsave( file.path(summ_outdir, paste0("clonotypes_count_rank_aboundance_groupby_",groupby, ".png")), dpi = 600, plot = p2)
    write.table( exp_cnt, file.path(summ_outdir, "clonotypes_count_rank_aboundance.xls"), sep = "\t", col.names = T, row.names = T)

    # exp_vol = repExplore(immdata$data, .method = "volume")

    #     p3 = vis(exp_vol, .by = "Sample", .meta = immdata$meta ) + theme(panel.grid =element_blank()) 
    #     ggsave( file.path(summ_outdir, paste0("clonotypes_count_barplot_groupby_",group,".pdf")) , plot = p3)
    #     ggsave( file.path(summ_outdir, paste0("clonotypes_count_barplot_groupby_",group,".png")), dpi = 600, plot = p3)

    # quit()

    #output directory setting
    dividx_outdir = file.path( output_dir, "Diversity" )
    if ( !file.exists( dividx_outdir ) ){
        dir.create( dividx_outdir )
    }

###### chao1
        chao1_idx = repDiversity(.data = immdata$data , .method = 'chao1')
        # gg_chao1 = vis(chao1_idx, .by = groupby, .meta = immdata$meta) + 
        #     theme(panel.grid =element_blank(), axis.text.x = element_text(size = 10, angle = 45, hjust = 1) ) 
        groupby="Sample"
        gg_chao1 = vis(chao1_idx, .meta = immdata$meta) + 
            theme(panel.grid =element_blank(), axis.text.x = element_text(size = 10, angle = 90, hjust = 1) ) 


        ggsave( file.path(dividx_outdir, paste0("diversity.chao1_index_groupby_",groupby,".pdf")), plot = gg_chao1)
        ggsave( file.path(dividx_outdir, paste0("diversity.chao1_index_groupby_",groupby,".png")),
                dpi = 1000, plot = gg_chao1)
        chao1_idx = as.data.frame( chao1_idx )
        chao1_idx$sample = rownames( chao1_idx )
        write.table( chao1_idx,
                    file.path(dividx_outdir, "diversity.chao1_index.xls"),
                    sep = "\t", col.names = T, row.names = T)

        # Hill numbers
        hill_idx = repDiversity(.data = immdata$data, .method = 'hill', .max.q = 6,.min.q = 1, .do.norm = NA, .laplace = 0)
        gg_hill = vis(hill_idx, .by = groupby, .meta = immdata$meta)+ theme(panel.grid =element_blank()) + labs(color= groupby)
        ggsave( file.path(dividx_outdir, paste0("diversity.hill_index_groupby_",groupby,".pdf")), plot = gg_hill)
        ggsave( file.path(dividx_outdir, paste0("diversity.hill_index_groupby_",groupby,".png")),
                dpi = 600, plot = gg_hill)
        hill_idx = as.data.frame( hill_idx )
        write.table( hill_idx,
                    file.path(dividx_outdir, "diversity.hill_index.xls")
                    ,sep = "\t", col.names = T, row.names = F)


    dxx="d50"
    ratio_xx = as.numeric(gsub( "^d", "", dxx ))
    dxx_idx = repDiversity(.data = immdata$data, .method = dxx, .perc = ratio_xx)
    gg_dxx = vis( dxx_idx , .meta = immdata$meta)+ theme(panel.grid =element_blank(),axis.text.x = element_text(size = 10, angle = 90, hjust = 1) )
    ggsave( file.path(dividx_outdir, paste0("diversity.", dxx,"_groupby_", groupby, ".pdf")), plot = gg_dxx)
    ggsave( file.path(dividx_outdir, paste0("diversity.", dxx,"_groupby_", groupby, ".png")), dpi = 600, plot = gg_dxx)


### 绘制 3D barplot #### 

    suppressPackageStartupMessages(library(epade))
    suppressPackageStartupMessages(library(viridis))

    #output_dir = getwd()
    geneusage_outdir = file.path( output_dir, "Gene_usage" )
    if ( !file.exists( geneusage_outdir ) ){
        dir.create( geneusage_outdir )
    }

    for ( samplex in names(immdata$data) ){
        dat_long = immdata$data[[samplex]] %>%
                    dplyr::select( V.name, J.name) %>%
                    table() %>% as.data.frame.table() 


        data_wide = immdata$data[[samplex]] %>% 
                    dplyr::select( V.name, J.name) %>% 
                    table() %>% as.data.frame.table() %>% 
                    data.table::dcast(V.name~J.name, value.var = "Freq") %>%
                    tibble::column_to_rownames(var = "V.name") %>% 
                    as.matrix()
        # data_wide = immdata$data[[samplex]] %>%
        #         dplyr::select( V.name, J.name) %>%
        #         data.table::dcast(V.name~J.name, value.var = "J.name") %>%
        #         tibble::column_to_rownames(var = "V.name") %>%
        #         as.matrix()

        write.table( data_wide %>% as.data.frame() %>% tibble::rownames_to_column(var = "V.name"),
            file.path(geneusage_outdir, paste0( samplex,"_VJ_pairing_Usage4sample_counts.xls")),
            sep = "\t", col.names =T, row.names =F, quote =F )



        # 3D barplot of V_J gene usage with count number mapped to height of bars.
        pdf( file.path(geneusage_outdir, paste0( samplex,"_VJ_pairing_Usage4sample",  "_3D_barplot.pdf")))
        par(las=2, cex = 0.3)
        bar3d.ade(t(data_wide), wall=3,ylab = "Clonotypes",zw=0.3, xw = 0.8,
                lcol = F, xlab = "V genes", zlab = "J genes",
                col = viridis::viridis(dim(t(data_wide))[1]), alpha = 0.6,bgbox=F)
        dev.off()
    }

}else{
    ### 绘制 3D barplot #### 
    suppressPackageStartupMessages(library(epade))
    suppressPackageStartupMessages(library(viridis))

    #output_dir = getwd()
    geneusage_outdir = file.path( output_dir, "Gene_usage" )
    if ( !file.exists( geneusage_outdir ) ){
        dir.create( geneusage_outdir )
    }

    dat_long = immdata$data %>%
            dplyr::select( V.name, J.name) %>%
            table() %>% as.data.frame.table() 

    data_wide = immdata$data %>% 
                dplyr::select( V.name, J.name) %>% 
                table() %>% as.data.frame.table() %>% 
                data.table::dcast(V.name~J.name, value.var = "Freq") %>%
                tibble::column_to_rownames(var = "V.name") %>% 
                as.matrix()

    samplex=assay_metadata$Sample
    write.table( data_wide %>% as.data.frame() %>% tibble::rownames_to_column(var = "V.name"),
        file.path(geneusage_outdir, paste0( samplex,"_VJ_pairing_Usage4sample_counts.xls")),
        sep = "\t", col.names =T, row.names =F, quote =F )

    # 3D barplot of V_J gene usage with count number mapped to height of bars.
    pdf( file.path(geneusage_outdir, paste0( samplex,"_VJ_pairing_Usage4sample",  "_3D_barplot.pdf")))
    par(las=2, cex = 0.3)
    bar3d.ade(t(data_wide), wall=3,ylab = "Clonotypes",zw=0.3, xw = 0.8,
            lcol = F, xlab = "V genes", zlab = "J genes",
            col = viridis::viridis(dim(t(data_wide))[1]), alpha = 0.6,bgbox=F)
    dev.off()

}
