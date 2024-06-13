casematch = function (search, match=rownames(seurat_ob),possibleOrigins=c(9606,10090),possibleTargets)
{
    search.match = sapply(X = search, FUN = function(s) {
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
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))

opt=list()
opt$RDS="../../rds/SNc_Neurons_clustering_res1.2.rds"
opt$output="../output/dotplot"
opt$extraGene="files/genelist2.xls"

seurat_ob = readRDS(opt$RDS)
DefaultAssay(seurat_ob) = "RNA"
metadata = seurat_ob@meta.data
seurat_ob = SetIdent( seurat_ob, value = "clusters")
seurat_ob@meta.data$clusters = factor(seurat_ob@meta.data$clusters, levels = sort( unique(as.numeric(as.vector( seurat_ob@meta.data$clusters)))))
groupby = "clusters"
splitby = NULL
root_dir = opt$output
dir.create(root_dir, recursive = T)
topn_markers = data.frame()
extra_gene = read.table(opt$extraGene, sep="\t", header = T, check.names=F)
if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
match = casematch(search = as.vector(formated_extra_gene$GENE),match = rownames(seurat_ob))
filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
if(length(filtered_gene)!=0){
    filtered_gene = as.data.frame(filtered_gene)
    colnames(filtered_gene) = "Gene"
    write.table(filtered_gene,file.path(root_dir,"genes_not_matched.xls"),quote = F,row.names=F)
}
formated_extra_gene = formated_extra_gene %>% dplyr::filter(GENE %in% names(match)) %>% dplyr::mutate(gene = match,folder_suffix = cluster) %>% dplyr::select(cluster, gene,folder_suffix)
topn_markers = rbind(topn_markers, formated_extra_gene)
pointsize = 0

clusterx = "extra"
topn = topn_markers %>% filter( folder_suffix == clusterx) %>% select(cluster,gene,folder_suffix)
topn_markers2vis = unique(as.vector(topn$gene))
path4vis = file.path(root_dir,paste0("markers_vis4",clusterx,collapse = ""))
if ( file.exists( path4vis ) ){
    output_dir = path4vis
}else{
    output_dir = path4vis
    dir.create(output_dir, recursive = T)
}
seurat_ob = SetIdent( seurat_ob, value = groupby )
# dot_plot = DotPlot(object = seurat_ob, features = topn_markers2vis ) + RotatedAxis() +sc
#             dot_plot + guides(color = guide_colorbar(order = 1, title = "Average Expression"))
# ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.pdf", collapse = "_")),width=(0.3*length(topn_markers2vis)+2.5+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10))
            



# 用dotplot函数获得绘图数据
a = DotPlot(object = seurat_ob, features = topn_markers2vis)
# > a$data
#                avg.exp     pct.exp features.plot id avg.exp.scaled
# Th        2.437918e-02  1.51821862            Th  2   -0.280170621
# Slc18a2   4.012648e-01 20.10796221       Slc18a2  2   -0.399079153
# Drd2      1.561609e-01  8.87314440          Drd2  2   -0.104784604
# Slc6a3    2.611161e-01 13.36032389        Slc6a3  2   -0.226140274
plot_data = a$data

# 基因分三组  DA_neuron Glut_neuron GABA_neuron -> red, blue, green
DA_neuron=c("Th","Slc18a2","Drd2","Slc6a3","Gch1","Aldh1a1","Sox6")
Glut_neuron=c("Slc17a6","Grm2","Pnoc","Oprm1")
GABA_neuron=c("Slc32a1","Gad1","Htr2c","Calb1","Dcx","Foxg1","Sox2","Gad2","Tph1","Sst","Pvalb","Vip","Cck","Calb2")
totalgene = c(DA_neuron,Glut_neuron,GABA_neuron)

install.packages("/public/scRNA_works/works/yangfan/packages/ggnewscale_0.4.8.tar.gz",
 lib = "/public/scRNA_works/works/yangfan/packages",
 repos = NULL,
 type="source")
library(ggnewscale, lib.loc = "/public/scRNA_works/works/yangfan/packages")
library(reshape2)


dot_plot = ggplot() + 
labs(color = "DA_neuron avg_exp") +
geom_point(data = plot_data[plot_data$features.plot %in% DA_neuron,], 
    aes(x = features.plot, y = id, size = pct.exp, color = pct.exp)) +
scale_color_gradient(low = "lightgrey", high = "cornflowerblue") +

new_scale("color") + 
labs(color = "Glut_neuron avg_exp") +
geom_point(data = plot_data[plot_data$features.plot %in% Glut_neuron,],
    aes(x = features.plot, y = id, size = pct.exp, color=pct.exp)) +
scale_color_gradient(low = "lightgrey", high = "pink") + 

new_scale("color") + 
labs(color = "GABA_neuron avg_exp") +
geom_point(data = plot_data[plot_data$features.plot %in% GABA_neuron,],
    aes(x = features.plot, y = id, size = pct.exp, color=pct.exp)) +
scale_color_gradient(low = "lightgrey", high = "darkolivegreen4") + 

scale_x_discrete(limits=totalgene) +
labs(x = "Features", y = "clusters") +
theme_classic() +
guides(size = guide_legend(title = "Percent Expressed")) +
theme(axis.text.x = element_text(angle = 45, vjust = 0.6))

ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.pdf", collapse = "_")),width=(0.3*length(topn_markers2vis)+6+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10), plot=dot_plot)
ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.png", collapse = "_")),width=(0.3*length(topn_markers2vis)+6+max(nchar(names(table(seurat_ob@meta.data[,groupby]))))/10), plot=dot_plot)
