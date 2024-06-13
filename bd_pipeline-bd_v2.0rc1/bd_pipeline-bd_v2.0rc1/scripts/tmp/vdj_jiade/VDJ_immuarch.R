
# module load immunoinformatics/1.0.0
suppressPackageStartupMessages( library("immunarch") )
suppressPackageStartupMessages( library("ggplot2") )
library("optparse")

option_list = list(
    make_option( c( "--inputdir","-i"), type = "character",default = "result/VDJ_aggr",
        help = "the input dir."),
    make_option( c("--output","-o"),type="character", default = "result/VDJ_aggr",
        help="the output directory of Clustering results." ),
    make_option( c("--vdj","-t"),type="character", default = "BCR",
        help="type of vdj info." ),
    make_option( c("--metadata","-m"),type="character", default = "config/samples.csv",
        help="the sample metadata which must include sample id in this assay design.")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#1.get info 
assay_metadata = read.delim(opt$metadata,sep=",",header=T)
assay_metadata = assay_metadata %>% rename(Sample = sampleid ) %>% as.data.frame()
rownames(assay_metadata) = assay_metadata$Sample

input_path = normalizePath(opt$inputdir)

immdata = list() #an empty immdata list
vdj_profile = Sys.glob(file.path(input_path,opt$vdj,paste0("/Clonotypes/",assay_metadata$Sample,"_",opt$vdj,".xls")))

immdata$data = repLoad(vdj_profile,)

names(immdata$data) = assay_metadata$Sample
immdata$meta = assay_metadata

groupby="Sample"

### summary
output_dir = file.path(opt$output)
    #1.summary
    summ_outdir = file.path( output_dir, "summary" )
    if ( file.exists( summ_outdir ) ){
        summ_outdir = summ_outdir
    }else{
        dir.create( summ_outdir )
    }

    if ( !is.null(groupby) ){
        groupby = unlist(strsplit( groupby, "," ))
    }
    exp_cnt = repExplore(immdata$data, .method = "count")
    p2 = vis(exp_cnt, .by = groupby, .meta = immdata$meta) + theme(panel.grid =element_blank()) + labs(color = groupby)
    ggsave( file.path(summ_outdir, paste0("clonotypes_count_rank_aboundance_groupby_",groupby, ".pdf")), plot = p2)
    ggsave( file.path(summ_outdir, paste0("clonotypes_count_rank_aboundance_groupby_",groupby, ".png")), dpi = 600, plot = p2)
    write.table( exp_cnt, file.path(summ_outdir, "clonotypes_count_rank_aboundance.xls"), sep = "\t", col.names = T, row.names = T)

    #output directory setting
    dividx_outdir = file.path( output_dir, "Diversity" )
    if ( file.exists( dividx_outdir ) ){
        dividx_outdir = dividx_outdir
    }else{
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

