library(dplyr)
library(optparse)
#library(stringr)

option_list = list(
    make_option( c("--input", "-i"), type = "character", default = NULL,
                 help = "the list file of marker genes"),
    make_option( c("--output","-o"),type="character", default = "./",
                help="the output directory of results.", metavar="character"),
    make_option( c("--topn", "-n"), type="integer", default = 25,
                 help = "the number of top markers for each cluster to visualizse."),
    make_option( c("--topby", "-c"), type = "character", default = "gene_diff",
                 help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table."),
    make_option( c("--split", "-s"), type = "logical", default = FALSE,
                 help="the column used to pick top n marker gene to visulize.The
                 option can be one of the column in the input marker genes table.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#######
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir)
    }
}

topby = opt$topby
##################
global_DEGs = read.delim(opt$input,sep="\t")
if ("avg_logFC" %in% colnames(global_DEGs)){
    colnames(global_DEGs) = gsub("avg_logFC","avg_log2FC",colnames(global_DEGs))
    print("为兼容脚本，表格中的 avg_logFC 修改为 avg_log2FC。")
    if ( opt$topby == "avg_logFC" ){ topby = "avg_log2FC" }
}


topn_markers  = global_DEGs %>% dplyr::group_by(cluster) %>%
      dplyr::arrange(p_val,desc(avg_log2FC),desc(gene_diff)) %>%
      dplyr::top_n(opt$topn,.data[[topby]])

write.table(topn_markers,file = file.path(output_dir,paste0("top", opt$topn, "_markers_for_each_cluster.xls", collapse = "")),
              col.names =T,row.names = F,sep = "\t",quote=F)

if (opt$split == TRUE ){
    for (i in unique(topn_markers$cluster)){
        topn_markers_sub = subset(topn_markers, cluster == i)
        print(paste0("cluster ",i," 已筛出，文件名中非字字符将被替换为下划线。"))
        cluster_name = gsub("\\W+","_",i)
        write.table(topn_markers_sub, file = file.path(output_dir,paste0("top", opt$topn, "_markers_for_cluster",cluster_name,".xls", collapse = "")),
              col.names =T,row.names = F,sep = "\t",quote=F)
    }
}

