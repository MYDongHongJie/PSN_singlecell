# module purge 
# module load OESingleCell/3.0.d  

########################### 获取 cellchat 的配受体数据 ##########################
suppressPackageStartupMessages({
    # library("CellChat")
    library("Seurat")
    library("pheatmap")
    library("dplyr")
    library("ComplexHeatmap")
    library(tidyverse)

})

# sr_celldb <- readRDS("/public/scRNA/works/liuhongyan/Project/scRNA/DR2021/HT2020-14895_HT2020-16850_HT2021-10354-human/re12S/Further_analysis_2022-04-14/1.Cellchat/cellchat/cellchat_list.rds")

st_celldb <- readRDS("/public/scRNA/works/liuhongyan/Project/scST/SA/HT2020-17989-h/Further_analysis_2022-04-14/cellphoneDB/cellphonedb_results.rds")

output_dir="heatmap"
if ( ! file.exists( output_dir ) ) {
    dir.create(output_dir, recursive=T )
}

data = st_celldb$ligrec
data$receptor_cell <- as.numeric(data$receptor_cell)
data$ligand_cell <- as.numeric(data$ligand_cell)

## 2. 结合单细胞的cellchat结果，将EC和Myeloid的通讯与空转cluster3内部通讯做overlap分析，做如下展示
### 筛选显著的配受体对
net = data %>% dplyr::filter(pval < 0.05) %>% dplyr::filter( receptor_cell=="3") %>% dplyr::filter(ligand_cell=="3") 
write.table(net,file.path(output_dir, "st_cluster3_and_cluster3_significant_LR.xls"),quote=F,sep="\t",row.names=F)

## 单细胞 cellphone 分析结果 ## 
sc_NEGATIVE <- readRDS("/public/scRNA/works/liuhongyan/Project/scST/SA/HT2020-17989-h/Further_analysis_2022-04-14/singlecell_cellphoneDB/cellphoneDB_NEGATIVE/cellphonedb_results.rds")
data = sc_NEGATIVE$ligrec
data <- data %>% dplyr::filter(pval < 0.05) %>% dplyr::filter( (receptor_cell %in% c("Endothelial cells","Myeloid cells") ) & (ligand_cell  %in% c("Endothelial cells","Myeloid cells") )  ) 
write.table(data,file.path(output_dir, "sc_NEGATIVE_EC_and_Mye_significant_LR.xls"),quote=F,sep="\t",row.names=F)

sc<- data %>%select("ligand","receptor")%>%mutate(sc=1)%>%distinct() 
st<- net %>%select("ligand","receptor")%>%mutate(st=1)%>%distinct() 

# > dim(sc)
# [1] 187   3
# > dim(st)
# [1] 203   3

ligand_gene <- c(sc$ligand,st$ligand) %>% unique 
receptor_gene <- c(sc$receptor,st$receptor) %>% unique 

ligand_receptor <- matrix(data=NA, nrow = length(ligand_gene), ncol = length(receptor_gene))%>% data.frame()
colnames(ligand_receptor) <- receptor_gene
rownames(ligand_receptor) <-ligand_gene
ligand_receptor_matrix<- ligand_receptor%>%
                        tibble::rownames_to_column("ligand")%>%
                        tidyr::pivot_longer(!ligand, names_to = "receptor", values_to = "status")%>% 
                        left_join(sc,by= c("ligand","receptor"))%>%
                        left_join(st,by= c("ligand","receptor"))%>%
                        mutate(status=case_when(sc==1 & st==1 ~ "sig(both)",
                                                is.na(sc) & st==1 ~ "sig(st)",
                                                sc==1 & is.na(st)  ~ "sig(sc)",
                                                TRUE ~ " "
                                               ))%>%select(!c("sc","st"))


head(ligand_receptor_matrix)

ligand_receptor_matrix<- ligand_receptor_matrix %>%pivot_wider(names_from = receptor, values_from = status)
head(ligand_receptor_matrix)

write.table(as.data.frame(ligand_receptor_matrix) ,file.path(output_dir, "ligand_receptor_status_st_and_sc_in_NEGATIVE.xls"), sep="\t", quote=F, row.names=F )

# > ligand_receptor_matrix %>% dim()
# [1] 120 156

# options(repr.plot.width=20, repr.plot.height=12)
plot<-ComplexHeatmap::Heatmap(ligand_receptor_matrix%>%column_to_rownames("ligand")%>%as.matrix,
                        #name='pvalue',
                        column_title ="L-R pairs st and sc in NEGATIVE",
                        column_title_gp = grid::gpar(fontsize = 25, fontface = "bold"),
                        row_names_side = "left",
                        cluster_rows = FALSE,
                        column_names_side = "top",
                        cluster_columns=FALSE,
                        col=c(`sig(both)`='#C92B86',`sig(st)`="lightblue",`sig(sc)`="lightgreen"," "="white"),
                        rect_gp = grid::gpar(col = "grey", lwd = 2),
                        na_col = "white",
                        #left_annotation = h_lig_fc,
                        #bottom_annotation= h_recptor_exp_ord_new ,
                        heatmap_legend_param = list(title=c('','')),
                        layer_fun = function(j, i, x, y, width, height, fill) {
                         v = ligand_receptor_matrix%>%column_to_rownames("ligand")%>%as.matrix 
                         l = c(v == "sig(both)"  )
                         #if(v=="sig(both)") {grid::grid.points(x[l], y[l], pch =16, size = grid::unit(2, "mm"))} 
                         grid::grid.points(x[l], y[l], pch =16, size = grid::unit(4, "mm"))
                        # if(v=="sig(st)") {grid::grid.points(x[l], y[l], pch =12, size = grid::unit(4, "mm"))} 
                        }  
                        )

OESingleCell::save_ggplots(filename=file.path(output_dir,"heatmap_NEGATIVE"),
                           plot = ggplotify::as.ggplot(plot),
                           width = 40 ,
                           height =30  )


#############################################  Positive ##############################################
sc_POSITIVE <- readRDS("/public/scRNA/works/liuhongyan/Project/scST/SA/HT2020-17989-h/Further_analysis_2022-04-14/singlecell_cellphoneDB/cellphoneDB_POSITIVE/cellphonedb_results.rds")
data = sc_POSITIVE$ligrec
data <- data %>% dplyr::filter(pval < 0.05) %>% dplyr::filter( (receptor_cell %in% c("Endothelial cells","Myeloid cells") ) & (ligand_cell  %in% c("Endothelial cells","Myeloid cells") )  ) 
write.table(data,file.path(output_dir, "sc_POSITIVE_EC_and_Mye_significant_LR.xls"),quote=F,sep="\t",row.names=F)

sc<- data %>%select("ligand","receptor")%>%mutate(sc=1)%>%distinct() 
st<- net %>%select("ligand","receptor")%>%mutate(st=1)%>%distinct() 

ligand_gene <- c(sc$ligand,st$ligand) %>% unique 
receptor_gene <- c(sc$receptor,st$receptor) %>% unique 

ligand_receptor <- matrix(data=NA, nrow = length(ligand_gene), ncol = length(receptor_gene))%>% data.frame()
colnames(ligand_receptor) <- receptor_gene
rownames(ligand_receptor) <-ligand_gene
ligand_receptor_matrix<- ligand_receptor%>%
                        tibble::rownames_to_column("ligand")%>%
                        tidyr::pivot_longer(!ligand, names_to = "receptor", values_to = "status")%>% 
                        left_join(sc,by= c("ligand","receptor"))%>%
                        left_join(st,by= c("ligand","receptor"))%>%
                        mutate(status=case_when(sc==1 & st==1 ~ "sig(both)",
                                                is.na(sc) & st==1 ~ "sig(st)",
                                                sc==1 & is.na(st)  ~ "sig(sc)",
                                                TRUE ~ " "
                                               ))%>%select(!c("sc","st"))
head(ligand_receptor_matrix)

ligand_receptor_matrix<- ligand_receptor_matrix %>%pivot_wider(names_from = receptor, values_from = status)
head(ligand_receptor_matrix)
write.table(as.data.frame(ligand_receptor_matrix) ,file.path(output_dir, "ligand_receptor_status_st_and_sc_in_POSITIVE.xls"), sep="\t", quote=F, row.names=F )

# options(repr.plot.width=20, repr.plot.height=12)
plot<-ComplexHeatmap::Heatmap(ligand_receptor_matrix%>%column_to_rownames("ligand")%>%as.matrix,
                        #name='pvalue',
                        column_title ="L-R pairs st and sc in POSITIVE",
                        column_title_gp = grid::gpar(fontsize = 25, fontface = "bold"),
                        row_names_side = "left",
                        cluster_rows = FALSE,
                        column_names_side = "top",
                        cluster_columns=FALSE,
                        col=c(`sig(both)`='#C92B86',`sig(st)`="lightblue",`sig(sc)`="lightgreen"," "="white"),
                        rect_gp = grid::gpar(col = "grey", lwd = 2),
                        na_col = "white",
                        #left_annotation = h_lig_fc,
                        #bottom_annotation= h_recptor_exp_ord_new ,
                        heatmap_legend_param = list(title=c('','')),
                        layer_fun = function(j, i, x, y, width, height, fill) {
                         v = ligand_receptor_matrix%>%column_to_rownames("ligand")%>%as.matrix 
                         l = c(v == "sig(both)"  )
                         #if(v=="sig(both)") {grid::grid.points(x[l], y[l], pch =16, size = grid::unit(2, "mm"))} 
                         grid::grid.points(x[l], y[l], pch =16, size = grid::unit(4, "mm"))
                        # if(v=="sig(st)") {grid::grid.points(x[l], y[l], pch =12, size = grid::unit(4, "mm"))} 
                        }  
                        )

# dim(ligand_receptor_matrix)
# [1] 122 152

OESingleCell::save_ggplots(filename=file.path(output_dir,"heatmap_POSITIVE"),
                           plot = ggplotify::as.ggplot(plot),
                           width = 40 ,
                           height =30  )

