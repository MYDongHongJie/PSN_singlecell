casematch = function (search, match=rownames(seurat_ob),possibleOrigins=c(9606,10090),possibleTargets){
    search.match <- sapply(X = search, FUN = function(s) {
        s=sub(" ","",s)
        res= grep(pattern = paste0("^", s, "$"), x = match,
            ignore.case = TRUE, perl = TRUE, value = TRUE)
        if (length(res) != 0) {
            return(res)
        } else {
            res= try(autotranslate(genes=s, targetGenes=match, possibleOrigins = possibleOrigins , returnAllPossible = FALSE, db = homologene::homologeneData)[,2],TRUE)
            if (class(res)!="try-error") {
                return(res)
            } else {
                return(character())
            }
        }
    })
    return(unlist(x = search.match))
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

library(Seurat)
library(ggplot2)
seurat_ob<-readRDS("/public/scRNA/works/guokaiqi/HT2021-12657-16-human/houxu-20220706/rds/data_ob_v3.rds")
extra_gene=read.delim("/public/scRNA/works/guokaiqi/HT2021-12657-16-human/houxu-20220720/dotplot.txt")
groupby="new_celltype2"
output_dir="/public/scRNA/works/guokaiqi/project/scRNA/HT2021-12657-16-human/houxu-20220720/plots"
    if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
    formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
    match = casematch(search = as.vector(formated_extra_gene$GENE),match = rownames(seurat_ob))
    # filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
    # if(length(filtered_gene)!=0){
    #     filtered_gene = as.data.frame(filtered_gene)
    #     colnames(filtered_gene) = "Gene"
    #     write.table(filtered_gene,file.path(root_dir,"genes_not_matched.xls"),quote = F,row.names=F)
    #     print("There are some mismatched gene symbol, Please check genes_not_matched.xls for the genename.")}
    topn_markers = formated_extra_gene %>% dplyr::filter(GENE %in% names(match)) %>% dplyr::mutate(gene = match,folder_suffix = cluster) %>% dplyr::select(cluster, gene,folder_suffix)


topn_markers2vis = unique(as.vector(topn_markers$gene))
# cols <- c(
#   "Endothelial cells" = CustomCol2(1),
#   "Fibroblasts" = CustomCol2(2),
#   "Smooth muscle cells" = CustomCol2(3),
#   "Myeloid cells" = CustomCol2(4),
#   "Neutrophils" = CustomCol2(5),
#   "Mast cells" = CustomCol2(6),
#   "B cells" = CustomCol2(7),
#   "T_NK cells" = CustomCol2(8),
#   "Oligodendrocytes" = CustomCol2(9)
# )
# seurat_ob = SetIdent( seurat_ob, value = groupby )
# dot_plot = DotPlot(object = seurat_ob, features = topn_markers2vis,cols=cols ) + RotatedAxis() + scale_fill_manual(values = cols) + coord_flip()
# dot_plot = dot_plot + guides(color = guide_colorbar(order = 1, title = "Average Expression"))

# scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

# plot_data = dot_plot$data
# plot <- ggplot(data = plot_data, mapping = aes_string(y = "features.plot", x = "id")) + 
#         geom_point(aes(size = pct.exp, color= id, fill = id)) +
#         scale_size_area(max_size = 10) + scale_colour_manual(values = cols) + theme_cowplot() + 
#         # scale.func(range = c(0, 6)) + #,limits = c(scale.min, scale.max)
#         #   theme_minimal(base_family = "serif", base_size = 0) + ,vjust = 0.85,hjust = 0.75
#         theme(axis.text.x = element_text(angle = 90,hjust=1) ) + guides(size = guide_legend(title = "Percent Expressed")) + 
#         labs(y = "Features", x = "Identity") 



# ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.pdf", collapse = "_")),height=(0.3*length(topn_markers2vis)+2.5+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10))
# ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.png", collapse = "_")),dpi = 1000 ,limitsize = F ,height=(0.3*length(topn_markers2vis)+2.5+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10))


library(viridis)
groupby="new_celltype2"
topn_markers2vis=markers$genes
            data_ob = SetIdent( data_ob, value = groupby )
            levels(data_ob)=rev(levels(data_ob))
            dot_plot = DotPlot(object = data_ob, features = markers$genes ) + RotatedAxis() + scale_colour_viridis_c() 
            dot_plot + guides(color = guide_colorbar(order = 1, title = "Average Expression"))+ theme(axis.text.x = element_text(angle = 90,vjust=0.5) )
            ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.pdf", collapse = "_")),bg="white",height=4,width=(0.2*length(topn_markers2vis)+2.5+max(nchar(names(table(data_ob@meta.data[,groupby]))))/10))
            ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.png", collapse = "_")),dpi = 1000 ,bg="white",height=4,limitsize = F ,width=(0.2*length(topn_markers2vis)+2.5+max(nchar(names(table(data_ob@meta.data[,groupby]))))/10))
#### vlnplot
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
            p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
               xlab("") + ylab("") + ggtitle(feature) +
               coord_flip() +
               theme(legend.position = "none",
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
            #    axis.ticks.x = element_line(),
                axis.text.x = element_text(size = 8,angle = 0, hjust = 0.5),
               axis.title.x = element_text(size = 12, angle = 0 ),
               plot.margin = plot.margin )
            return(p)
            }

            StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
            plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
            plot_list[[1]]  <- plot_list[[1]] + 
            theme(axis.text.y=element_text(size = 12), axis.ticks.y = element_line() )
            p <- patchwork::wrap_plots(plotlist = plot_list, ncol = length(plot_list) )
            return(p)
            }
flip <- StackedVlnPlot(data_ob, topn_markers2vis, pt.size=0, cols=rev(CustomCol2(1:9)), group.by = groupby) + coord_flip()     
            pdf(file.path(output_dir,paste0("top", "marker_gene_flip_vlnplot_plot",".pdf", collapse = "_")),
                width = length(topn_markers2vis)*1 + max(nchar(as.character(data_ob@meta.data[,groupby])))/2 , height = 8)
            print(flip)
            dev.off()
            png(file.path(output_dir,paste0("top","marker_gene_flip_vlnplot_plot",".png",collapse = "_")),
                width = length(topn_markers2vis)*1+max(nchar(as.character(data_ob@meta.data[,groupby])))/2, height = 8 , res = 96, units = "in")
            print(flip)
            dev.off()

### circle
library(dplyr)
library(tibble)
library(ggrepel)
ditto = c("#E69F00","#56B4E9","#009E73", "#F0E442", "#0072B2" ,"#D55E00" ,"#CC79A7", "#666666" ,"#AD7700", "#1C91D4", 
          "#007756", "#D5C711", "#005685", "#A04700" ,"#B14380" ,"#4D4D4D","#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71" )
output_dir="/public/scRNA/works/guokaiqi/project/scRNA/HT2021-12657-16-human/houxu-20220720/CD8_T_cells"
seurat_ob<-ReadX(input="/public/scRNA/works/guokaiqi/project/scRNA/HT2021-12657-16-human/houxu-20220624/CD4_T_new_var_features2465/seurat.h5seurat",informat="h5seurat")
seurat_ob<-ReadX(input="/public/scRNA/works/guokaiqi/project/scRNA/HT2021-12657-16-human/houxu-20220708/CD8_T_cells/new_clusters.h5seurat",informat="h5seurat")

propby="new_clusters"
propby="Four_celltype_sub"
propby="clusters"
groupby="group"
groupby="sampleid"

    data = as.data.frame(seurat_ob@meta.data[,c(groupby,propby)]) %>% 
                    group_by(.dots = c(groupby,propby)) %>% 
                    dplyr::summarize(cell_number = n(), .groups = "drop_last") %>% 
                    mutate(freq = round((cell_number / sum(cell_number)) * 100,2)) 
     write.table(data,file.path(output_dir,paste0(groupby,'_grouped_pieplot_data.xls')), sep = '\t',col.names=T, row.names =F)

    pp = ggplot(data, aes(x = 0.6, y = freq , fill = get(propby))) + 
           geom_bar(stat = "identity" , position = "stack",color = "white" , width = 0.1 ) + 
           labs(x = "", y = "", title = "")+coord_polar(theta = "y",start = 1, direction = 1, clip = "on")+ 
           theme_bw()+scale_fill_manual(values=ditto[1:length(unique(seurat_ob@meta.data[[propby]]))])+
           geom_text(stat="identity",color = "black", size=3, position=position_stack(vjust = 0.6),
           aes(y = 0, x =0.2, label = ""))+
           geom_text(stat="identity",color = "black", size=2, check_overlap = TRUE, position=position_stack(vjust = 0.6), 
           aes(y = freq, x =0.7, label = scales::percent(freq/100)))+ #geom_text_repel(color= "grey20",size= 3,point.padding = NA) +
           facet_wrap(~get(groupby),ncol = 2)+theme_void()+theme(legend.title = element_blank())

           ggsave(file.path(output_dir,paste0('Grouped_',groupby,'_Column_pieplot.png')),bg="white",pp,width=7,height=10)
           ggsave(file.path(output_dir,paste0('Grouped_',groupby,'_Column_pieplot.pdf')),bg="white",pp,width=7,height=10)

### barplot

    DATA = as.data.frame(seurat_ob@meta.data[,c(propby,groupby)]) %>% 
                    group_by(.dots = c(propby,groupby)) %>% 
                    dplyr::summarize(cell_number = n(), .groups = "drop_last") %>% 
                    mutate(freq = (cell_number / sum(cell_number)) * 100)

    write.table(DATA,file.path(output_dir,paste0(groupby,'_grouped_barplot_data.xls')), sep = '\t',col.names=T, row.names =F)

lev=length(unique(seurat_ob@meta.data[[groupby]]))
    p = ggplot(DATA , aes(x = get(propby) , weight = freq ,fill = get(groupby) ))+ 
                coord_flip()+geom_bar(position = "stack" )+
                # facet_wrap(~get(groupby),ncol = 1,strip.position="top")+
                scale_fill_manual(values=ditto[14:c(lev+14)])+ theme_classic()+ theme_bw() +
                scale_y_continuous(expand = c(0,0))+ theme(axis.title = element_blank(),legend.title = element_blank())

                ggsave(file.path(output_dir,paste0('Grouped_',groupby,'_Column_barplot.png')), p,bg="white",height=4)
                ggsave(file.path(output_dir,paste0('Grouped_',groupby,'_Column_barplot.pdf')), p,bg="white",height=4)

#### ridgeplot

library(ggridges)
library(ggplot2)
library(monocle)

cds<-readRDS("CD8_T_cells/Pseuduotime/pseudotime_results.rds")
cds <- orderCells(cds, root_state = 2)
output_dir="CD8_T_cells/Pseuduotime/"
p <- ggplot(aes(x=Pseudotime, y = 0 ,fill=group), data = pData(cds)) +
                    geom_density_ridges(scale = 3) +
                    ylab("group") + facet_wrap(~group,ncol = 1) +
                    xlab("Pseudotime") + #xlim(0,100) +
                    theme_bw() +  scale_fill_manual(values = c("#7fc97f","#beaed4","#492d73","#cd4275")) +
                    theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 1, colour = "black"))
ggsave(file.path(output_dir, paste0("ridges_group_psecol.pdf")), plot = p, width = 15, height = 6)
ggsave(file.path(output_dir, paste0("ridges_group_psecol.png")), plot = p, width = 15, height = 6)

