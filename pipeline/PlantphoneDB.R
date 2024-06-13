## usage :
## Rscript PlantphoneDB.R -i input_file -o output_file -c new_celltype -s ath -v cellphonedb -p 0.1 -m Product

# 载入PlantPhoneDB
# load functions from PlantPhoneDB R package. https://plantphonedb.readthedocs.io/en/latest/08-Software.html
RNAAnnotateCelltype <- function(RNA, genes, signatures = "human.immune.CIBERSORT", min.score = 0, orig.ident = NULL, outdir = "."){
    if(is.null(orig.ident)){
        if(class(signatures) == "character"){
            signatures_name = signatures
            data(list = signatures)
            signatures = get(signatures)
        } else {
            signatures_name = ""
        }
        celltypes <- as.character(unique(signatures[,1]))
        signature_list <- sapply(1:length(celltypes),function(x){
                          return(toupper(as.character(signatures[which(signatures[,1]==celltypes[x]),2])))})
        names(signature_list) <- celltypes
            
        cluster_celltype_score = sapply(as.integer(unique(RNA@meta.data$seurat_clusters))-1, function(x){
            idx = genes$cluster==x
            avglogFC = genes$avg_logFC[idx]
            names(avglogFC) = toupper(genes$gene[idx])
            score_cluster = sapply(signature_list, function(y){
                score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
                return(score)
            })
        })
        colnames(cluster_celltype_score) = as.character(as.integer(unique(RNA@meta.data$seurat_clusters))-1)
        
        cellscore_max = apply(cluster_celltype_score, 2, max, na.rm = TRUE)
        # aa=cluster_celltype_score %>% drop_na()
        A<-na.omit(cluster_celltype_score)
        cellscore_max_celltype = apply(A, 2, function(x){
            if (max(x) < min.score){
               return("Others")
            }else{
                return(rownames(cluster_celltype_score)[which.max(x)])
            }
        })
        
        RNA@meta.data$assign.ident = as.character(as.integer(RNA@meta.data$seurat_clusters)-1)
        current.cluster.ids = as.character(as.integer(unique(RNA@meta.data$seurat_clusters))-1)
        new.cluster.ids = cellscore_max_celltype
        
        RNA@meta.data$assign.score = cellscore_max[RNA@meta.data$assign.ident]
        RNA@meta.data$assign.ident = plyr::mapvalues(x = RNA@meta.data$assign.ident,
                                                           from = current.cluster.ids, to = new.cluster.ids)
        if(signatures_name == "human.immune.CIBERSORT") {
            current.cluster.ids = sort(celltypes)
            new.cluster.ids = c("DC", "Mast", "CD4Tconv", "NK", "pDC", "CD8T", 
                                "CD8Tex", "Endothelial", "Eosinophils", "Fibroblasts", "Mono/Macro", "Mono/Macro", 
                                "Mono/Macro", "B", "Mono/Macro", "Myofibroblasts", "B", "CD4Tconv", 
                                "Neutrophils", "Plasma", "DC", "Mast", "CD4Tconv", "NK", 
                                "pDC", "CD4Tconv", "TMKI67", "Treg")
            RNA@meta.data$assign.ident = plyr::mapvalues(x = RNA@meta.data$assign.ident,
                                                         from = current.cluster.ids, to = new.cluster.ids)
        }
        p = DimPlot(object = RNA, label = TRUE, pt.size = 0.2, group.by = "assign.ident", label.size = 3, repel = T)
        ggsave(file.path(outdir, paste0(RNA@project.name, "_annotated.png")), p, width=6, height=4)}
    else{
        RNA$assign.ident <- orig.ident
        p = DimPlot(object = RNA, label = TRUE, pt.size = 0.2, group.by = "assign.ident", label.size = 3, repel = T)
        ggsave(file.path(outdir, paste0(RNA@project.name, "_original.png")), p,  width=5.5, height=4)}
    
    return(RNA)    
}

LRTargetFisher <- function(Targets,GeneSets,TotalGene){
  a <- intersect(Targets,GeneSets)
  b <- setdiff(GeneSets,a)
  c <- setdiff(Targets,a)
  abc <- union(a,b)
  abc <- union(abc,c)
  d <- setdiff(TotalGene,abc)
  tbl <- matrix(c(length(a), length(b), length(c), length(d)), nrow = 2)
  fisher <- fisher.test(tbl, alternative = "two.sided", simulate.p.value=TRUE)
  pvalue <- fisher$p.value
  overlap <- paste(a,collapse = ", ")
  return(list(pvalue=pvalue,overlap=overlap))
}

LRscore <- function(expr, LRdb, cluster = NULL, min.pct = 0.1, method='LRscore',iterations=100, seed=123, ...){
    expr <- as.data.frame(expr)
    expr <- expr[rowSums(expr)>0,]
    u <- sum(expr)/(nrow(expr)*ncol(expr))
    expr_flt <- function(expr, LR_gene, ident='NK'){
        tmp <- expr[rownames(expr) %in% LR_gene, new_cluster==ident]
        pct <- apply(tmp, 1, function(x){
            sum(x>0)
        })
        pct <- pct/ncol(tmp)
        tmp <- tmp[pct>min.pct, ]
        tmp <- apply(tmp,1,mean)
        return(tmp)
    }
    score <- function(l,r,u=NULL,method='LRscore'){
        if(method=='LRscore'){
            if(!is.null(u)){
                s <- (l*r)^(1/2)/((l*r)^(1/2)+u)
            }else{
                print('u is NULL')
            }
        }else if(method=='WeightProduct'){
            s <- l*r
            z_scores <- (s-mean(s))/sd(s)
            ## max min normalization
            s <- (z_scores - min(z_scores))/(max(z_scores)-min(z_scores))
        }else if(method=='Average'){
            s <- (l+r)/2
        }else if(method=='Product'){
            s <- (l*r)
        }
      return(as.numeric(s))
    }

    LR_result <- NULL
    set.seed(seed)
    for(Ligands_cell in unique(cluster)){
        for(Receptors_cell in unique(cluster)){
            new_cluster <- cluster
            tmp.L <- expr_flt(expr, LRdb$Ligands, Ligands_cell)
            tmp.R <- expr_flt(expr, LRdb$Receptors, Receptors_cell)
            tmp.LR <- LRdb[LRdb$Receptors %in% names(tmp.R) & LRdb$Ligands %in% names(tmp.L),]
			if(nrow(tmp.LR)<=0){
				next
			}
            tmp.LR$Ligands_cell <- Ligands_cell
            tmp.LR$Receptors_cell <- Receptors_cell
            tmp.LR$Ligands_expr <- tmp.L[tmp.LR$Ligands]
            tmp.LR$Receptors_expr <- tmp.R[tmp.LR$Receptors]
            tmp.LR$Score <- score(l=tmp.L[tmp.LR$Ligands], r=tmp.R[tmp.LR$Receptors],u=u, method)
            Score <- matrix(0,nrow(tmp.LR),iterations)
            if(method %in% c('Product','Average')){
                for(i in 1:iterations){
                    cols <- sample(1:ncol(expr),ncol(expr))
                    new_cluster <- cluster[cols]
                    sample.L <- expr_flt(expr, LRdb$Ligands, Ligands_cell)
                    sample.R <- expr_flt(expr, LRdb$Receptors, Receptors_cell)
                    L <- sample.L[tmp.LR$Ligands]
                    L[is.na(L)] <- 0
                    R <- sample.R[tmp.LR$Receptors]
                    R[is.na(R)] <- 0
                    SS <- score(l=L, r=R,u=u, method)
                    Score[,i] <- SS
                }
                p.value <- lapply(1:nrow(tmp.LR),function(x){
                    random <- as.numeric(Score[x,])
                    ttest <- t.test(random,mu=tmp.LR$Score[x],alternative = 'less')
                    return(ttest$p.value)
                })
                tmp.LR$Pvalue <- unlist(p.value)
            }
            LR_result <- rbind(LR_result,tmp.LR)
        }
    }
    rownames(LR_result) <- 1:nrow(LR_result)
    LR_result$Type <- ifelse(LR_result$Ligands_cell==LR_result$Receptors_cell,'Autocrine','Paracrine')
    LR_result$LR_pair <- paste0(LR_result$Ligands,'->',LR_result$Receptors)
    LR_result$Cell_pair <- paste0(LR_result$Ligands_cell,'->',LR_result$Receptors_cell)
    return(LR_result)
}

LR_pathway <- function(lr, objs, CellA, CellB, neighbor, geneSet, cor_value=0.013, ...){
    if(CellA!=CellB){
        degs <- FindMarkers(objs, ident.1 = CellA, ident.2 = CellB, only.pos = TRUE, verbose = FALSE)
    }else{
        degs <- FindMarkers(objs, ident.1 = CellA, only.pos = TRUE, verbose = FALSE)
    }
    degs <- degs %>%
        filter(p_val_adj<0.05)

    expr <- objs@assays$RNA@data
    expr <- as.data.frame(expr)
    expr <- expr[rowSums(expr)>0,]
    cluster <- Idents(objs)
    re <- NULL
    for(i in 1:nrow(lr)){
        needGene <- c(lr$Ligands[i],lr$Receptors[i],rownames(degs))
        if(CellA==CellB){
            expr_flt <- expr[rownames(expr) %in% needGene, cluster %in% CellA]
        }else{
            expr_flt <- expr[rownames(expr) %in% needGene, cluster %in% c(CellA,CellB)]
        }
        pair <- t(data.frame(apply(expr_flt[c(lr$Ligands[i],lr$Receptors[i]),],2,mean)))
        rownames(pair) <- lr$LR_pair[i]
        expr_flt <- rbind(expr_flt,pair)
        expr_flt <- expr_flt[!rownames(expr_flt) %in% c(lr$Ligands[i],lr$Receptors[i]),]
        corr <- cor(t(expr_flt),method="spearman")
        corGene <- corr[,lr$LR_pair[i]]
        corGene <- corGene[corGene>cor_value]
        miGene <- rownames(expr_flt)[rownames(expr_flt) %in% names(corGene)]
        ## calculate mutual information
        dat <- discretize(t(expr_flt[miGene,]))
        mi <- mutinformation(dat,method= "mm")
        mi <- aracne.m(mi, 0.15)
        for(j in 1:ncol(mi)){
            mi[j,j] <- 0
        }
        g <- graph.adjacency(mi,mode="directed",weighted=T)
        #rank <- page.rank(g)$vector
        #rankGene <- miGene[order(rank,decreasing = TRUE)]
        #degs_flt <- rankGene[1:(length(rankGene)*0.5)]
		selegoV <- ego(g, order=neighbor, nodes=lr$LR_pair[i],mode = "all", mindist = 0)
		selegoG <- induced_subgraph(g,unlist(selegoV))
		degs_flt <- get.vertex.attribute(selegoG)$name
        for(k in unique(geneSet$`Gene set name`)){
            sets <- subset(geneSet,`Gene set name`==k)
            a <- intersect(degs_flt,sets$Gene)
            if(length(a)>0){
                fisher <- LRTargetFisher(Targets=degs_flt,GeneSets=sets$Gene,TotalGene=rownames(expr))
				GeneRatio <- paste(length(a),length(degs_flt),sep="/")
				BgRatio <- paste(length(sets$Gene),length(rownames(expr)),sep="/")
                tmp <- data.frame(LR_pair=lr$LR_pair[i],Cell_pair=lr$Cell_pair[i],Pathway=unique(sets$`Gene set name`),GeneRatio,BgRatio,
                                  Pvalue=fisher$pvalue,OverlapGene=fisher$overlap)
            }else{
                tmp <- NULL
            }
            re <- rbind(re,tmp)
        }
    }
    return(re)

}

CCI_circle <- function(interaction_count, mycolor){
    # parameters
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
    par(mar = rep(0, 4))
	arr.col <- interaction_count %>%
        select(Ligands_cell,Receptors_cell) %>%
        unique() %>%
        mutate(color='black') %>%
        as.data.frame()
    # color palette
    #mycolor <- viridis(num_cluster, alpha = 1, begin = 0, end = 1, option = "D")
    #mycolor <- mycolor[sample(1:num_cluster)]
    chordDiagram(
      x = interaction_count,
      grid.col = mycolor,
      transparency = 0.25,
      directional = 1,
      direction.type = c("arrows", "diffHeight"),
      diffHeight  = -0.04,
      annotationTrack = "grid",
      annotationTrackHeight = c(0.05, 0.1),
      #link.arr.type = "big.arrow",
        symmetric = TRUE,
      link.sort = TRUE,
      link.arr.col = arr.col, link.arr.length = 0.3,
      link.largest.ontop = TRUE)

    # Add text and axis
    circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        # Add names to the sector.
        circos.text(
          x = mean(xlim),
          y = 3.2,
          labels = sector.index,
          facing = "bending",
          cex = 1.5
          )
      }
    )
}

CCI_network <- function(interaction_count,mycolor, vertex.label.cex=1.5, edge.label.cex=1, title="",edgeLabel=TRUE){
    color <- data.frame(Ligands_cell=unique(interaction_count$Ligands_cell),color=mycolor[unique(interaction_count$Ligands_cell)])
    interaction_count <- interaction_count %>%
        inner_join(color)
    net <- graph_from_data_frame(interaction_count)
    karate_groups <- cluster_optimal(net)
    coords <- layout_in_circle(net, order= order(membership(karate_groups)))
    E(net)$width  <- E(net)$Number/10
    V(net)$color <- mycolor[get.vertex.attribute(net)$name]
    #E(net)$color <- mycolor
    if(edgeLabel){
        E(net)$label <- E(net)$Number
    }
    pic <- plot(net, edge.arrow.size=1,
     edge.curved=0.2,
     vertex.label.color="black",
     layout = coords,
	 edge.label.cex= edge.label.cex,
     vertex.label.cex=vertex.label.cex,main=title)
}

heatmap_count <- function(interaction_count,log10=FALSE,decimal=1,text_size=15,number_size=5,title=NULL){
    if(log10){
        interaction_count$Number <- log10(interaction_count$Number+1)
    }
    pic <- ggplot(interaction_count,aes(Ligands_cell,Receptors_cell, fill=Number)) +
      geom_tile(colour='white') +
      scale_fill_viridis_c(position="right")+
      theme(axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank())+
      theme_bw()+
      theme(axis.title=element_text(size=text_size),
            axis.text=element_text(size=text_size,color='black'),
            axis.text.x=element_text(size=text_size,angle=60,hjust=1),
            legend.text=element_text(size=text_size),
            plot.title = element_text(size = text_size,hjust=0.5))+
        coord_equal()+
        ggtitle(title)

    ### annotate number
    for(i in 1:nrow(interaction_count)){
        if(interaction_count$Number[i]<median(interaction_count$Number)){
            color <- 'white'
        }else{
            color <- 'black'
        }
    pic <- pic+
    annotate("text", x = interaction_count$Ligands_cell[i], y = interaction_count$Receptors_cell[i],
             label = round(interaction_count$Number[i],decimal),color=color,size=number_size)
    }
    return(pic)
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
# 载入需要的R包
# Load the required packages
.libPaths(c("/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/R/library","/home/chenhaoruo/R/x86_64-conda-linux-gnu-library/4.0"))
suppressPackageStartupMessages(library(discretization))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(lsa))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parmigene))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(infotheo))
suppressPackageStartupMessages(library(igraph))
# suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library("optparse"))

##装不上
# install.packages("muxViz")
# library(muxViz)
# install.packages("chorddiag")
# library(chorddiag)



#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
        help = "[REQUIRED]The filtered exprssion matrix in prefered format."),
    make_option( c("--celltype", "-c" ), type = "character",
        help = "[REQUIRED]The selected column to be compared."),
    make_option( c("--method", "-m" ), type = "character",default = "Product" ,
        help = "The selected column to be compared. Choose from 'Product','WeightProduct','lrscore','Average'."),
    make_option( c("--minpct", "-p" ), type = "numeric",default = 0.1 ,
        help = "Set the min pct to calculate the score."),
    make_option( c("--informat", "-f" ), type = "character", default = "h5seurat",
        help = "The indication of type of input expression matrix, the possible type can be:
    #                         seurat: the seurat object from the clustering results."),
    make_option( c("--spieces", "-s"), type = "character",
        help = "[REQUIRED]the spieces abstraction, the current options can be 'ath' for Arabidopsis thaliana and 'rice' for Oryza sativa."),
    make_option( c("--vis", "-v"), type = "character",default = "cellphonedb" ,
        help = "What kind of visualization you want to get, the current options can be 'cellphonedb' for cellphonedb style and 'origin' for original plantphonedb style."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="the output directory.", metavar="outputdir")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# objs <- readRDS("/public/scRNA/works/luyao/project/HT2020-20853-Arabidopsis_thaliana/mapping_raw_cluster/sub_C1-2_mapping_raw_cluster.rds")
if(is.null(opt$input)){
    stop("Please provide the input file.")
}else{
    if(opt$informat == "rds"){
        objs <- readRDS(opt$input)
    }else if(opt$informat == "h5seurat"){
        objs = OESingleCell::ReadX(input = opt$input, informat = 'h5seurat',verbose = F)
    }else{
        stop("Please provide the correct input format.")
    }
}

if(is.null(opt$output)){
    output_dir="./"
}else{
    if ( file.exists( opt$output ) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir, recursive = T)
    }
}



# 配体-受体互作数据集
# A dataset with ligand-receptor pairs from https://jasonxu.shinyapps.io/PlantPhoneDB/
if (opt$spieces == "ath"){
    a=load("/public/scRNA_works/works/chenhaoruo/test/plantphone20220824/LR_pair_ath_update.RDa")
    LR_pair=test
    LR_pair <- LR_pair %>%
        filter(source!="orthologs") %>%
        select(Ligands, Receptors) %>%
        unique()
}else if (opt$spieces == "rice"){
    a=load("/public/scRNA_works/works/chenhaoruo/test/plantphone20220824/LR_pair_osa_update.RDa")
    LR_pair=LR_pair_osa
    LR_pair <- LR_pair %>%
        select(Ligands, Receptors) %>%
        unique() 
}else{
    stop("The species is not supported, please check the input.")
}

if(is.null(opt$celltype)){
    stop("Please provide the celltype.")
}else{
    groupby=opt$celltype
}

Idents(objs)=objs@meta.data[,groupby]
# 计算配体-受体互作的分数
# Score ligand-receptor interactions
Score <- LRscore(objs@assays$RNA@data, LRdb=LR_pair, cluster = Idents(objs), min.pct = opt$minpct,iterations=100, method=opt$method)

write.table(Score,sep="\t",quote=F,row.names=F,file.path(output_dir,"PlantphoneDB_score.xls"))

if(opt$vis=="cellphonedb"){
    # 生成cellphonedb格式的数据
    # Generate data in cellphonedb format
    Score_cellphone = Score %>%
        rename("ligand" = "Ligands", 
            "receptor" = "Receptors",
            "ligand_cell" = "Ligands_cell",
            "receptor_cell" = "Receptors_cell",
            "ligand_expr" = "Ligands_expr",
            "receptor_expr" = "Receptors_expr",
            "expr" = "Score",
            "pval" = "Pvalue") %>%
        select(receptor,ligand,receptor_cell,ligand_cell,expr,pval,receptor_expr,ligand_expr)
    write.table(Score_cellphone,sep="\t",quote=F,row.names=F,file.path(output_dir,"PlantphoneDB_score.xls"))

    cellphone_vis_R = "/public/scRNA_works/works/chenhaoruo/test/plantphone20221124/visualize_cellcomm.R"
    path4vis = file.path(output_dir,"PlantphoneDB_vis")
    system(glue::glue("module purge && module load OESingleCell/2.0.0 &&
    Rscript {cellphone_vis_R}   -i {output_dir}/PlantphoneDB_score.xls -f matrix -d  network,circos,dotplot,chorddiagram,heatmap,barplot -o {path4vis}  -n 5 && module purge && module load OESingleCell/3.0.d" )) 
    system(glue::glue("rm -f {path4vis}/*.docx"))
    system(glue::glue("cp /public/scRNA_works/works/chenhaoruo/test/plantphone20221124/test_cellphone/PlantphoneDB_vis/PlantphoneDB细胞通讯分析说明.docx {path4vis}/"))

}else if(opt$vis=="origin"){
    interaction_count <- Score %>%
        group_by(Ligands_cell,Receptors_cell) %>%
        summarise(Number=n(),.groups = 'drop')
    mycolor = CustomCol2(1:length(unique(objs@meta.data[,groupby])))
    names(mycolor) = unique(objs@meta.data[,groupby])
    # interaction_count %>%
    #     mutate(Type=ifelse(Ligands_cell==Receptors_cell,"Autocrine","Paracrine")) %>%
    #     group_by(Type) %>%
    #     summarise(Number=sum(Number))
    Autocrine <- interaction_count[interaction_count$Ligands_cell==interaction_count$Receptors_cell,]
    Paracrine <- interaction_count[interaction_count$Ligands_cell!=interaction_count$Receptors_cell,]
    setwd(output_dir)
    pdf('no_self_CCI_circle.pdf',width=8,height=5)
    CCI_circle(Paracrine, mycolor)
    dev.off()
    pdf('Self_CCI_circle.pdf',width=8,height=5)
    CCI_circle(Autocrine, mycolor)
    dev.off()
    pdf('CCI_circle.pdf')
    CCI_circle(interaction_count, mycolor)
    dev.off()
    # Heatmap shows number of ligand-receptor interactions between pairwise cell type.
    pdf('number_of_ligand-receptor_interactions_heatmap.pdf',width=8,height=5)
    heatmap_count(interaction_count,text_size=15,number_size=5,decimal=4)
    dev.off()

# Top 10 ligand-receptor pairs with P value < 0.05 show different regulatory pattern. Columns are scaled by max ligand-receptor expression.
    Top10 <- Score %>%
        arrange(desc(Score)) %>%
        select(LR_pair) %>%
        unique() %>%
        head(10) %>%
        inner_join(Score) %>%
        select(LR_pair,Cell_pair,Score) %>%
        spread(.,Cell_pair,Score) %>%
        replace(is.na(.), 0)
    rownames(Top10) <- Top10$LR_pair
    Top10 <- Top10[,-1]
    Top10 <- t(Top10)
    Top10 <- apply(Top10,2,function(x){x/max(x)})
    options(repr.plot.width=5, repr.plot.height=10)
    Top10_pic <- pheatmap(Top10, scale="none",angle_col=45,fontsize_row=4,
            cluster_rows = T,cluster_cols = F,show_colnames=T)
    ggsave(Top10_pic,file="Top10_ligand-receptor_heatmap.pdf")
    print("Convert pdf to png...")
    system("for i in  `ls *.pdf`;do /data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -verbose -density 500 -trim  $i  -quality 100  -flatten  ${i/.pdf/.png}  ;done")
}else{
    stop("The visualization type is not supported.")
}




# geneSet <- fread('./PlantPhoneDB/plantGSAD/Ara_ALL.KEGG.txt')
# geneSet$Gene <- toupper(geneSet$Gene)

# CellA <- "Endodermis"
# CellB <- "Endosperm"
# lr <- subset(Heat_sig,Ligands_cell==CellA & Receptors_cell==CellB )



# # Construction of intracellular signaling pathway
# pathway_result2 <- LR_pathway(lr, objs_heat, CellA, CellB, neighbor=2, geneSet)

# pathway_result2$FDR <- p.adjust(pathway_result2$Pvalue,method='BH')

# pathway_result2 <- pathway_result2[order(pathway_result2$Pvalue),]


# dim(pathway_result2)



