casematch = function (search, match=rownames(seurat_ob),possibleOrigins=c(9606,10090),possibleTargets)
{
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
############
library(ggplot2)
library(Seurat)
library(OESingleCell)
library("RColorBrewer")

lev=c("ISC,TA_1,TA_2,RP,EP_1,EP_2,EP_3,EC_1,EC_2,EC_3,EC_4,EC_5,EC_6,EC_7,EC_8,EC_9,EC_10,Goblet cell,Paneth cell,EEC_1,EEC_2,TC")                                                                                                                                 
level = unlist(strsplit(lev,","))  
data_ob@meta.data$sub_celltype=factor(data_ob@meta.data$sub_celltype,levels=level)

SeuratDisk::UpdateH5Seurat(file = "/public/scRNA_works/works/guokaiqi/project/scRNA/DZOE2023040276/houxu-20230717/epi/seurat.h5seurat", object =  data_ob, verbose = FALSE )

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))
# 不同基因还需要不同细胞类型画
type = list(rp_ep = c("RP","EP_1","EP_2","EP_3"),
            isc_ta = c("ISC","TA_1","TA_2"),
            ec1_4 = c("EC_1","EC_2","EC_3","EC_4"),
            ec5_10 = c("EC_5","EC_6","EC_7","EC_8","EC_9","EC_10"))
groupby="sub_celltype2"

for (i in c("isc_ta","rp_ep","ec1_4","ec5_10")){
output_dir=i
df = OESingleCell::colData(data_ob)
desired_cells= subset(df, df$sub_celltype %in% type[[i]])
data_sub = data_ob[, rownames(desired_cells)]

extra_gene = read.delim(paste0(i,".txt"),sep="\t")
if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene"){colnames(extra_gene)[1]="extra"}
formated_extra_gene = as.data.frame(tidyr::gather(extra_gene,key = "cluster",value = "GENE"))
match = casematch(search = as.vector(formated_extra_gene$GENE),match = rownames(data_ob))
filtered_gene = formated_extra_gene$GENE[!formated_extra_gene$GENE %in% names(match )& formated_extra_gene$GENE != ""]
formated_extra_gene = formated_extra_gene %>% dplyr::filter(GENE %in% names(match)) %>% dplyr::mutate(gene = match,folder_suffix = cluster) %>% dplyr::select(cluster, gene,folder_suffix)
#topn_markers =  formated_extra_gene
#formated_extra_gene = unique(as.vector(formated_extra_gene$gene))
# 提取 GO ID（假设 GO ID 总为末尾 7 位数字）
formated_extra_gene$folder_suffix <- gsub(".*\\.go\\.(\\d{7})\\.", "GO.\\1", formated_extra_gene$folder_suffix)

mybool <- !duplicated(formated_extra_gene$gene)
# 提取不重复的数据子集
unique_formated_extra_gene <- formated_extra_gene[mybool, ]
# 根据不同term转换为list
mylist <- split(unique_formated_extra_gene[,c(2,3)], unique_formated_extra_gene$folder_suffix)
topn_markers <- lapply(mylist, function(subset) { 
    a = unique(as.vector(subset$gene)) 
    if (length(unique(as.vector(subset$gene))) > 100){ a = NULL} # 某一term基因数大于100个就不画了
    return(a)
    })

topn_markers2vis = unique(unlist(topn_markers))

data_sub = SetIdent( data_sub, value = groupby )
# DotPlot函数输入基因列表就可以自动进行分面
dot_plot = DotPlot(object = data_sub, features = topn_markers ) + RotatedAxis() + sc 
dot_plot + guides(color = guide_colorbar(order = 1, title = "Average Expression"))
ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.pdf", collapse = "_")),width=(0.3*length(topn_markers2vis)+2+max(nchar(names(table(data_sub@meta.data[,groupby]))))/10),limitsize = FALSE)
#ggsave(file.path(output_dir,paste0("top", "marker_gene_dotplot.png", collapse = "_")),dpi = 1000 ,limitsize = F ,width=(0.3*length(topn_markers2vis)+2.5+max(nchar(names(table(data_sub@meta.data[,groupby]))))/10))

}


