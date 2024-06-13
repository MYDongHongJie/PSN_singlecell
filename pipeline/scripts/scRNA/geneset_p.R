rm(list=ls())
## package
suppressWarnings({
	suppressPackageStartupMessages( library("Seurat") )
	suppressPackageStartupMessages( library("optparse") )
	suppressPackageStartupMessages( library("dplyr"))
	suppressPackageStartupMessages( library("OESingleCell"))
	suppressPackageStartupMessages( library("ggplot2") )
	suppressPackageStartupMessages( library("tibble") )
})

opt=list()

### rds 获取 #### 
opt$RDS="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-20625-m/SA20210518/new_celltype.no_unknown_sort_seurat_0.6.rds"

seurat_ob = readRDSMC(opt$RDS)
subset <- SubsetData(seurat_ob,subset.name="new_celltype", accept.value=c("Olfactory HBC","mOSN","SUS","BG","olfactory ensheathing glia","Vascular smooth muscle cells","Fibroblasts stromal cell","Pericytes","Respiratory ciliated cells","B cell","T cells","Myeloid cells") )
all <- seurat_ob
seurat_ob <- subset

seurat_ob@meta.data$new_celltype= factor(seurat_ob@meta.data$new_celltype, levels=c("Olfactory HBC","mOSN","SUS","BG","olfactory ensheathing glia","Vascular smooth muscle cells","Fibroblasts stromal cell","Pericytes","Respiratory ciliated cells","B cell","T cells","Myeloid cells") )

################ gene list 匹配
root_dir="./boxplot"
dir.create(root_dir)
opt=list()
topn_markers = data.frame()
opt$extraGene="/public/scRNA/works/liuhongyan/Project/scRNA/SA/HT2020-20625-m/SA20210809/genelist.xls"
extra_gene = read.table(opt$extraGene, sep="\t", header = T)
    if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
    formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
    match = CaseMatch(search = as.vector(formated_extra_gene$GENE),match = rownames(seurat_ob))
    filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
    if(length(filtered_gene)!=0){
        filtered_gene = as.data.frame(filtered_gene)
        colnames(filtered_gene) = "Gene"
        write.table(filtered_gene,file.path(root_dir,"genes_not_matched.xls"),quote = F,row.names=F)
        print("There are some mismatched gene symbol, Please check genes_not_matched.xls for the genename.")}
    formated_extra_gene = formated_extra_gene %>% filter(GENE %in% names(match)) %>% mutate(gene = match,folder_suffix = cluster) %>% select(cluster, gene,folder_suffix)
    topn_markers = rbind(topn_markers, formated_extra_gene)


groupby="new_celltype"
path4vis = file.path(root_dir,"geneset_visualization")
if ( file.exists( path4vis ) ){
	output_dir = path4vis
}else{
	output_dir = path4vis
	dir.create(output_dir)
}
seurat_ob = SetIdent( seurat_ob, value = groupby )
topn_markers2vis=list()
for ( clusterx in unique(topn_markers$folder_suffix) ){
	topn_markers2vis[[clusterx]] = subset(topn_markers,folder_suffix == clusterx)$gene
}
seurat_ob = AddModuleScore(seurat_ob,features=topn_markers2vis,name=names(topn_markers2vis))
colnames(seurat_ob@meta.data)[(dim(seurat_ob[[]])[2]-length(topn_markers2vis)+1):dim(seurat_ob[[]])[2]] = names(topn_markers2vis)

output_dir=root_dir
metadata = seurat_ob@meta.data
metadata= metadata[,c("clusters","group","sampleid",groupby,names(topn_markers2vis))]

data <- metadata %>% tibble::rownames_to_column("barcode")

write.table(data,file.path(output_dir,"AddModuleScore_genesetscore_data.xls"),sep="\t",row.names=F,quote=F)


groupby="group"
# p = ggplot(data,aes_string(x=groupby,y="Aging",fill=groupby))+
# 		stat_boxplot(geom = "errorbar",width = 0.5,lwd=0.5) + 
# 		geom_boxplot()+
# 		labs(x="",y="Gene Set Score" ) +
# 		theme(panel.grid.major =element_blank(),
# 				panel.grid.minor = element_blank(),
# 				panel.background = element_blank(),
# 				legend.position = "none",
# 				axis.line = element_line(color = "black"))

# ggsave(file.path(output_dir, "Aging_geneset0.png"), dpi=1000  )

library("ggplot2")
library("ggpubr")

# p <- ggboxplot(data, x = groupby, y = "Aging",
#           color = groupby, palette = "jco",
#           #add = "jitter",
#           # facet.by = "new_celltype", 
# 		  short.panel.labs = FALSE) +
# 		  facet_wrap(~new_celltype, ncol = 12, scales = "free_y")  + 
# 			labs(x="",y="Gene Set Score", title="Aging" ) +
#   			theme_bw() +
#   			theme( plot.title = element_text(hjust = 0.5, size=14),
# 			  axis.text.x = element_blank(),
# 			  axis.ticks.x = element_blank(), 
#          panel.grid = element_blank(),
#          legend.title=element_blank()) +
# 		 monocle:::monocle_theme_opts() 

# p  + stat_compare_means(aes(group = groupby), label = "p.signif")
# ggsave(file.path(output_dir, "Aging_geneset2.png"), dpi=1000, width=30, height=6 )






# p <- ggboxplot(data, x = groupby, y = "Aging",
#           color = groupby, palette = "jco",
#           #add = "jitter",
#           facet.by = "new_celltype", 
# 		  short.panel.labs = TRUE) +
# 		 stat_compare_means(method = "t.test", label =  "p.signif", label.x = 1.5) + 
# 			labs(x="",y="Gene Set Score", title="Aging" ) +
#   			theme_bw() +
#   			theme( plot.title = element_text(hjust = 0.5, size=14),
# 			  axis.text.x = element_blank(),
# 			  axis.ticks.x = element_blank(), 
#          panel.grid = element_blank(),
#          legend.title=element_blank())
# 		#   +
# 		#  monocle:::monocle_theme_opts() 

# # p  + stat_compare_means(aes(group = groupby), label = "p.signif")
# ggsave(file.path(output_dir, "Aging_geneset3.png"), dpi=1000, width=12, height=7 )




# CustomCol2 <- function(n){
#   my_palette=c(
#     "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
#     "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
#     "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
#     "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
#     "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17")
#   return(my_palette[n])
# }

# p <- ggboxplot(data, x = "new_celltype", y = "Aging",
#           color = groupby, palette = "jco",
#           #add = "jitter",
#           #facet.by = "new_celltype", 
# 		  short.panel.labs = TRUE) +
# 		 stat_compare_means(aes_string(group = groupby), method = "t.test", label =  "p.signif", label.x = 1.5) + 
# 			labs(x="",y="Gene Set Score", title="Aging" ) +
#   			theme_bw() +
#   			theme( plot.title = element_text(hjust = 0.5, size=20),
# 			#   axis.text.x = element_blank(),
# 			#   axis.ticks.x = element_blank(), 
#          panel.grid = element_blank(),
#          legend.title=element_blank())  +
# 		 scale_colour_manual(values=CustomCol2(1:length(unique(data[,groupby]))))
# 		#   +
# 		#  monocle:::monocle_theme_opts() 

# # p  + stat_compare_means(aes(group = groupby), label = "p.signif")
# ggsave(file.path(output_dir, "Aging_geneset4.png"), dpi=1000, width=20, height=6 )

for ( term in colnames(data)[6:10] ) {
	p <- ggboxplot(data, x = "new_celltype", y = term ,
          color = groupby, palette = "jco",
          #add = "jitter",
          #facet.by = "new_celltype", 
		  short.panel.labs = TRUE) +
		 stat_compare_means(aes_string(group = groupby), method = "t.test", label =  "p.signif", label.x = 1.5) + 
			labs(x="",y="Gene Set Score", title= term ) +
  			theme_bw() +
  			theme( plot.title = element_text(hjust = 0.5, size=20),
			#   axis.text.x = element_blank(),
			#   axis.ticks.x = element_blank(), 
         panel.grid = element_blank(),
         legend.title=element_blank())  +
		 scale_colour_manual(values=CustomCol2(1:length(unique(data[,groupby]))))
		#   +
		#  monocle:::monocle_theme_opts() 

# p  + stat_compare_means(aes(group = groupby), label = "p.signif")
	ggsave(file.path(output_dir, paste0(term,"_geneset.png") ), dpi=1000, width=20, height=8 )
	ggsave(file.path(output_dir, paste0(term,"_geneset.pdf") ), dpi=1000, width=20, height=8 )
}


> ggboxplot
function (data, x, y, combine = FALSE, merge = FALSE, color = "black",
    fill = "white", palette = NULL, title = NULL, xlab = NULL,
    ylab = NULL, bxp.errorbar = FALSE, bxp.errorbar.width = 0.4,
    facet.by = NULL, panel.labs = NULL, short.panel.labs = TRUE,
    linetype = "solid", size = NULL, width = 0.7, notch = FALSE,
    select = NULL, remove = NULL, order = NULL, add = "none",
    add.params = list(), error.plot = "pointrange", label = NULL,
    font.label = list(size = 11, color = "black"), label.select = NULL,
    repel = FALSE, label.rectangle = FALSE, ggtheme = theme_pubr(),
    ...)
{
    .opts <- list(combine = combine, merge = merge, color = color,
        fill = fill, palette = palette, title = title, xlab = xlab,
        ylab = ylab, bxp.errorbar = bxp.errorbar, bxp.errorbar.width = bxp.errorbar.width,
        facet.by = facet.by, panel.labs = panel.labs, short.panel.labs = short.panel.labs,
        linetype = linetype, size = size, width = width, notch = notch,
        select = select, remove = remove, order = order, add = add,
        add.params = add.params, error.plot = error.plot, label = label,
        font.label = font.label, label.select = label.select,
        repel = repel, label.rectangle = label.rectangle, ggtheme = ggtheme,
        ...)
    if (!missing(data))
        .opts$data <- data
    if (!missing(x))
        .opts$x <- x
    if (!missing(y))
        .opts$y <- y
    .user.opts <- as.list(match.call(expand.dots = TRUE))
    .user.opts[[1]] <- NULL
    for (opt.name in names(.opts)) {
        if (is.null(.user.opts[[opt.name]]))
            .opts[[opt.name]] <- NULL
    }
    .opts$fun <- ggboxplot_core
    if (missing(ggtheme) & (!is.null(facet.by) | combine))
        .opts$ggtheme <- theme_pubr(border = TRUE)
    p <- do.call(.plotter, .opts)
    if (.is_list(p) & length(p) == 1)
        p <- p[[1]]
    return(p)
}
<bytecode: 0x55d0eeb6f8d0>
<environment: namespace:ggpubr>