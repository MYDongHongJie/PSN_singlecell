#plot fig
library(dplyr)

#该函数构建一个链接并返回
plotKEGGURL <- function(pathwayid,red_list,green_list){
  str1<-paste0('https://www.kegg.jp/kegg-bin/show_pathway?map=map',pathwayid,'&multi_query=')
  str2<-paste0(red_list,'+%23FFFFFF,red',collapse = '%0d%0a')
  str3<-paste0(green_list,'+%23FFFFFF,green',collapse = '%0d%0a')
  return(paste0(str1,str2,'%0d%0a',str3))
}

#该函数读取物种的Ko数据库，返回Ko号
getKObyentrezid <- function(entrezid,species){
  raw_db = read.delim(file.path("/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell_test/KeggDB",paste0(species,"_keggDB.txt")),header = F)
	db = raw_db$V1	
  names(db) = raw_db$V2
	return(gsub('^ko:','',db[paste0(species,':',entrezid)]))
}




perform_KEGG_enrichment <- function(species.org, diff.cluster_genes, updiff_genes, downdiff_genes, enrichment_res, species) {
  
  # 检查输入是否为空
  if (is.null(species.org)) {
    stop("species.org cannot be NULL")
  }
  
  # 提取所有基因的 ENTREZ ID
  all_entrez <- mapIds(species.org, diff.cluster_genes, 'ENTREZID', 'SYMBOL')
  all_entrez <- all_entrez[!is.na(all_entrez)]
  
  # 提前处理上调和下调基因的 ENTREZ ID
  up_entrez <- mapIds(species.org, updiff_genes, 'ENTREZID', 'SYMBOL')
	if(!is.null(downdiff_genes)){
		down_entrez <- mapIds(species.org, downdiff_genes, 'ENTREZID', 'SYMBOL')
	}else {
		 down_entrez=NULL
	}
  
  
  # 检查 KEGG 富集文件是否存在
  kegg_file_path <- file.path(enrichment_res, "KEGG_enrichment.xls")
  if (!file.exists(kegg_file_path)) {
    print(paste("File not found:", kegg_file_path))
		return( NULL)
  }else {
		 # 读取 KEGG 富集数据
		kegg_data <- read.table(kegg_file_path, header = TRUE, sep = "\t")
		
		# 处理每一行 KEGG 数据
		result_data <- list()
		for (num in 1:nrow(kegg_data)) {
			
			# 获取 pathway ID 并提取基因信息
			pathway_id <- gsub("\\D", "", kegg_data[num, 1])
			gene_symbols <- unlist(strsplit(kegg_data[num, 7], '/'))
			
			# 分别获取上调和下调基因
			up_gene <- intersect(gene_symbols, updiff_genes)
			down_gene <- intersect(gene_symbols, downdiff_genes)
			
			# 使用提前映射好的 ENTREZ ID 进行 KEGG pathway 处理
			up_id <- getKObyentrezid(up_entrez[up_gene], species)
			down_id <- getKObyentrezid(down_entrez[down_gene], species)
			
			# 存储结果
			result_data[[num]] <- list(pathway_id = pathway_id, up_id = up_id, down_id = down_id)
		}
		
		# 检查结果长度是否匹配
		if (length(result_data) != nrow(kegg_data)) {
			stop("KEGG table and result length do not match")
		}
		
		# 生成 KEGG 链接并更新表格
		for (i in 1:length(result_data)) {
			link <- plotKEGGURL(result_data[[i]]$pathway_id, result_data[[i]]$up_id, result_data[[i]]$down_id)
			kegg_data[i, "link"] <- link
		}
		write.table(kegg_data, kegg_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
		# 返回处理过的 kegg_data
		return(kegg_data)
	}
  
  
}
























groupDiffAuto <- function(seurat_obj,file_out,annocol,species,avg_log2FC,topn){
    seurat_diff_cluster_dir=paste(file_out,"Diff_Group",sep = "/")
    #find different gene between sample for each cluster test
    Idents(seurat_obj)<- seurat_obj@meta.data[,annocol]
    if(!file.exists(seurat_diff_cluster_dir)){
        dir.create(seurat_diff_cluster_dir)
    }

    seurat_obj$cluster_group <- paste(seurat_obj$group,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
		seurat_obj$cluster_sample <- paste(seurat_obj$sample,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
    seurat_obj$celltype <- Idents(seurat_obj)
    Idents(seurat_obj) <- "cluster_sample"
    sample_cluster_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
		Idents(seurat_obj) <- "cluster_group"
    cluster_group_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
    write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
		write.table(cluster_group_avg,paste(seurat_diff_cluster_dir,"group_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

    cluster_group=as.vector(unique(Idents(seurat_obj)))

    sample_group <- unique(seurat_obj$group)
		
    if(length(unique(sample_group))>1){
				
        sample_group_length=length(unique(sample_group))-1
        sample_group_length_1=length(unique(sample_group))
        for (idx in 1:sample_group_length) {
            idx_index=idx+1
            for(idx1 in idx_index:sample_group_length_1){
                if(idx != idx1){ 
                    compare=paste(unique(sample_group)[idx],unique(sample_group)[idx1],sep="_vs_")
                    compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
                    if(!file.exists(compare_dir)){
                        dir.create(compare_dir)
                    }
										all_diff = data.frame()
                    for(sub_cluster in unique(seurat_obj$celltype)){
                        cp1=paste(unique(sample_group)[idx],paste("cluster",sub_cluster,sep=""),sep="_")
                        cp2=paste(unique(sample_group)[idx1],paste("cluster",sub_cluster,sep=""),sep="_")
                        if(!cp1 %in% cluster_group || ! cp2 %in% cluster_group){
                            next;
                        }
                        ncell1=CellsByIdentities(seurat_obj,cp1)
                        ncell2=CellsByIdentities(seurat_obj,cp2)
                        if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
                            next;
                        }
          
                        cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
                       if(!file.exists(cluster_compare_dir)){
                           dir.create(cluster_compare_dir)
                       }
											 allcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/all",sep="/")
                       upcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/up",sep="/")
                       downcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/down",sep="/")
                       if(!file.exists(upcluster_compare_dir_enrich)){dir.create(upcluster_compare_dir_enrich,recursive =TRUE)}
                       if(!file.exists(downcluster_compare_dir_enrich)){dir.create(downcluster_compare_dir_enrich,recursive =TRUE)}
											 if(!file.exists(allcluster_compare_dir_enrich)){dir.create(allcluster_compare_dir_enrich,recursive =TRUE)}
                       diff.cluster=FindMarkers(seurat_obj, ident.1 = cp1, ident.2 = cp2, verbose = FALSE,min.pct = 0.10, logfc.threshold = avg_log2FC)
                       ###增加基因描述
											 if (!is.null(species.org)) {
											 diff.cluster$description= unname(mapIds(x = species.org,
                              keys = rownames(diff.cluster),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first"))}


                       write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
                      
											 diff.cluster$cluster = sub_cluster
											 diff.cluster <- diff.cluster %>% rownames_to_column(var = "gene")
	                     all_diff = rbind(all_diff,diff.cluster)
                       updiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC>0.25)
                       downdiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC< -0.25)
                       updiff_genes=updiff.cluster$gene
                       downdiff_genes=downdiff.cluster$gene
                       diff.cluster_genes=diff.cluster$gene
											 try(enrichment(species=species,outDir=upcluster_compare_dir_enrich,geneList=updiff_genes))
                       try(enrichment(species=species,outDir=downcluster_compare_dir_enrich,geneList=downdiff_genes))
											 try(enrichment(species=species,outDir=allcluster_compare_dir_enrich,geneList=diff.cluster_genes))
											 for (enrichment_res in c(upcluster_compare_dir_enrich, downcluster_compare_dir_enrich, allcluster_compare_dir_enrich)) {
											 		kegg_data = perform_KEGG_enrichment(species.org, diff.cluster_genes, updiff_genes, downdiff_genes, enrichment_res, species)
											 		if(is.null(kegg_data)){
														print(paste0(enrichment_res,"kegg file is empty"))
													}
											 }
											 
                   }
									 all_diff_table = file.path(compare_dir,"All_diff_gene.xls")
									 write.table(all_diff,all_diff_table,sep='\t',row.names = F,quote=F)
      	   }
        }

      }
			  
    }
}

groupDiffSpeci <- function(seurat_obj,file_out,annocol,cmpfile,species,avg_log2FC,topn){
    cmpdf <- read.table(cmpfile,header=T,sep="\t",stringsAsFactors=FALSE, colClasses = c("character"))
    #head(cmpdf)
    seurat_diff_cluster_dir=paste(file_out,"Diff_Group",sep = "/")
    #find different gene between sample for each cluster test
    Idents(seurat_obj)<- seurat_obj@meta.data[,annocol]
    if(!file.exists(seurat_diff_cluster_dir)){
        dir.create(seurat_diff_cluster_dir)
    }

    seurat_obj$cluster_sample <- paste(seurat_obj$sample,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
		seurat_obj$cluster_group <- paste(seurat_obj$group,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
    seurat_obj$celltype <- Idents(seurat_obj)
    Idents(seurat_obj) <- "cluster_sample"
    sample_cluster_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
		Idents(seurat_obj) <- "cluster_group"
    group_cluster_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
    write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
		write.table(group_cluster_avg,paste(seurat_diff_cluster_dir,"group_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    cluster_group=as.vector(unique(Idents(seurat_obj)))

    for(x in  c(1:nrow(cmpdf))){
        compare=paste(cmpdf[x,1],cmpdf[x,2],sep="_vs_")
        compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
        if(!file.exists(compare_dir)){
            dir.create(compare_dir)
        }
				all_diff = data.frame()
        for(sub_cluster in unique(seurat_obj$celltype)){
            print(sub_cluster)
            cp1=paste(cmpdf[x,1],paste("cluster",sub_cluster,sep=""),sep="_")
            cp2=paste(cmpdf[x,2],paste("cluster",sub_cluster,sep=""),sep="_")
            if(!cp1 %in% cluster_group || ! cp2 %in% cluster_group){
                next;
            }
            ncell1=CellsByIdentities(seurat_obj,cp1)
            ncell2=CellsByIdentities(seurat_obj,cp2)
            if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
                next;
            }

            cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
            if(!file.exists(cluster_compare_dir)){
                dir.create(cluster_compare_dir)
            }
            cluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment",sep="/")
						allcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/all",sep="/")
            upcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/up",sep="/")
            downcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/down",sep="/")
            if(!file.exists(upcluster_compare_dir_enrich)){dir.create(upcluster_compare_dir_enrich,recursive =TRUE)}
            if(!file.exists(downcluster_compare_dir_enrich)){dir.create(downcluster_compare_dir_enrich,recursive =TRUE)}
						if(!file.exists(allcluster_compare_dir_enrich)){dir.create(allcluster_compare_dir_enrich,recursive =TRUE)}
            diff.cluster=FindMarkers(seurat_obj, ident.1 = cp1, ident.2 = cp2, verbose = FALSE,min.pct = 0.10,logfc.threshold=avg_log2FC)
						###增加基因描述
						if (!is.null(species.org)) {
						diff.cluster$description= unname(mapIds(x = species.org,
									keys = rownames(diff.cluster),
									column = "GENENAME",
									keytype = "SYMBOL",
									multiVals = "first"))}
            write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            # sub_compare_cluster1=paste(cmpdf[x,1],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
            # sub_compare_cluster2=paste(cmpdf[x,2],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
            # diff.cluster$gene=row.names(diff.cluster)
						diff.cluster$cluster = sub_cluster
						diff.cluster <- diff.cluster %>% rownames_to_column(var = "gene")
	          all_diff = rbind(all_diff,diff.cluster)
            #diff.cluster=arrange(diff.cluster, desc(avg_log2FC))
            #top3_diff_gene=diff.cluster$gene[1:3]
            updiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC>0.25)
						downdiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC< -0.25)
            # write.table(updiff.cluster,paste(upcluster_compare_dir_enrich,"diff_gene.enrichment.up.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            # head(updiff.cluster)
            # downdiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC< -0.25)
            # write.table(downdiff.cluster,paste(downcluster_compare_dir_enrich,"diff_gene.enrichment.down.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            # head(downdiff.cluster)
						#筛选对应的上下调和所有基因
            updiff_genes=updiff.cluster$gene
            downdiff_genes=downdiff.cluster$gene
						diff.cluster_genes=diff.cluster$gene
						##富集
            try(enrichment(species=species,outDir=upcluster_compare_dir_enrich,geneList=updiff_genes))
            try(enrichment(species=species,outDir=downcluster_compare_dir_enrich,geneList=downdiff_genes))
						try(enrichment(species=species,outDir=allcluster_compare_dir_enrich,geneList=diff.cluster_genes))
						for (enrichment_res in c(upcluster_compare_dir_enrich, downcluster_compare_dir_enrich, allcluster_compare_dir_enrich)) {
							kegg_data = perform_KEGG_enrichment(species.org, diff.cluster_genes, updiff_genes, downdiff_genes, enrichment_res, species)
							if(is.null(kegg_data)){
								print(paste0(enrichment_res,"kegg file is empty"))
							}
						}
        }
				all_diff_table = file.path(compare_dir,"All_diff_gene.xls")
				write.table(all_diff,all_diff_table,sep='\t',row.names = F,quote=F)
				#system(glue::glue("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/PlotKEGGNet.R -m  {all_diff_table} -o {seurat_diff_cluster_dir} -s {species} -t {compare} -n {topn}" ))
    }

}

groupDiffSpeciAll <- function(seurat_obj,file_out,annocol,cmpfile,species,avg_log2FC,topn){
    cmpdf <- read.table(cmpfile,header=T,sep="\t",stringsAsFactors=FALSE, colClasses = c("character"))
 
    seurat_diff_cluster_dir=paste(file_out,"Diff_Group",sep = "/")
    #find different gene between sample for each cluster test
    Idents(seurat_obj)<- seurat_obj@meta.data[,annocol]
    if(!file.exists(seurat_diff_cluster_dir)){
        dir.create(seurat_diff_cluster_dir)
    }

    seurat_obj$cluster_sample <- paste(seurat_obj$sample,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
		seurat_obj$cluster_group <- paste(seurat_obj$group,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
    seurat_obj$celltype <- Idents(seurat_obj)
    Idents(seurat_obj) <- "cluster_sample"
    sample_cluster_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
		Idents(seurat_obj) <- "cluster_group"
		cluster_group_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
    write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    write.table(cluster_group_avg,paste(seurat_diff_cluster_dir,"group_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    cluster_group=as.vector(unique(Idents(seurat_obj)))

    for(x in  c(1:nrow(cmpdf))){
        compare=paste(cmpdf[x,1],cmpdf[x,2],sep="_vs_")
        compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
        if(!file.exists(compare_dir)){
            dir.create(compare_dir)
        }
        for(sub_cluster in unique(seurat_obj$celltype)){
            print(sub_cluster)
            cp1=paste(cmpdf[x,1],paste("cluster",sub_cluster,sep=""),sep="_")
            cp2=paste(cmpdf[x,2],paste("cluster",sub_cluster,sep=""),sep="_")
            if(!cp1 %in% cluster_group || ! cp2 %in% cluster_group){
                next;
            }
            ncell1=CellsByIdentities(seurat_obj,cp1)
            ncell2=CellsByIdentities(seurat_obj,cp2)
            if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
                next;
            }

            cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
            if(!file.exists(cluster_compare_dir)){
                dir.create(cluster_compare_dir)
             }
            allcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/all",sep="/")
						upcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/up",sep="/")
						downcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/down",sep="/")
						if(!file.exists(upcluster_compare_dir_enrich)){dir.create(upcluster_compare_dir_enrich,recursive =TRUE)}
						if(!file.exists(downcluster_compare_dir_enrich)){dir.create(downcluster_compare_dir_enrich,recursive =TRUE)}
						if(!file.exists(allcluster_compare_dir_enrich)){dir.create(allcluster_compare_dir_enrich,recursive =TRUE)}
						diff.cluster=FindMarkers(seurat_obj, ident.1 = cp1, ident.2 = cp2, verbose = FALSE,min.pct = 0.10, logfc.threshold = avg_log2FC)
						###增加基因描述
						if (!is.null(species.org)) {
						diff.cluster$description= unname(mapIds(x = species.org,
									keys = rownames(diff.cluster),
									column = "GENENAME",
									keytype = "SYMBOL",
									multiVals = "first"))}
						write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
						diff.cluster$cluster = sub_cluster
						diff.cluster <- diff.cluster %>% rownames_to_column(var = "gene")
						all_diff = rbind(all_diff,diff.cluster)
						updiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC>0.25)
						downdiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC< -0.25)
						updiff_genes=updiff.cluster$gene
						downdiff_genes=downdiff.cluster$gene
						diff.cluster_genes=diff.cluster$gene
						try(enrichment(species=species,outDir=upcluster_compare_dir_enrich,geneList=updiff_genes))
						try(enrichment(species=species,outDir=downcluster_compare_dir_enrich,geneList=downdiff_genes))
						try(enrichment(species=species,outDir=allcluster_compare_dir_enrich,geneList=diff.cluster_genes))
						for (enrichment_res in c(upcluster_compare_dir_enrich, downcluster_compare_dir_enrich, allcluster_compare_dir_enrich)) {
							kegg_data = perform_KEGG_enrichment(species.org, diff.cluster_genes, updiff_genes, downdiff_genes, enrichment_res, species)
							if(is.null(kegg_data)){
								print(paste0(enrichment_res,"kegg file is empty"))
							}
						}						
        }
				
				all_diff_table = file.path(compare_dir,"All_diff_gene.xls")
				write.table(all_diff,all_diff_table,sep='\t',row.names = F,quote=F)
				#system(glue::glue("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/PlotKEGGNet.R -m  {all_diff_table} -o {seurat_diff_cluster_dir} -s {species} -t {compare} -n {topn}" ))
    }
}


groupDiffSpeci_Auto <- function(seurat_obj,file_out,annocol,cmpfile,species,avg_log2FC,topn){
    
    #head(cmpdf)
    seurat_diff_cluster_dir=paste(file_out,"Diff_Group",sep = "/")
    #find different gene between sample for each cluster test
    Idents(seurat_obj)<- seurat_obj@meta.data[,annocol]
    if(!file.exists(seurat_diff_cluster_dir)){
        dir.create(seurat_diff_cluster_dir)
    }

    seurat_obj$cluster_sample <- paste(seurat_obj$sample,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
		seurat_obj$cluster_group <- paste(seurat_obj$group,paste("cluster",Idents(seurat_obj),sep="") , sep = "_")
    seurat_obj$celltype <- Idents(seurat_obj)
    Idents(seurat_obj) <- "cluster_sample"
    sample_cluster_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
		Idents(seurat_obj) <- "cluster_group"
    cluster_group_avg=AverageExpression(seurat_obj,"RNA",slot = "data")$RNA
    write.table(sample_cluster_avg,paste(seurat_diff_cluster_dir,"sample_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
		write.table(cluster_group_avg,paste(seurat_diff_cluster_dir,"group_avgExpression_cluster.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    cluster_group=as.vector(unique(Idents(seurat_obj)))

    for(x in  cmpfile){
				com1 = unlist(strsplit(x,'/'))[1]
				com1_str = unlist(strsplit(com1,"\\+"))
				com2 = unlist(strsplit(x,'/'))[2]
				com2_str = unlist(strsplit(com2,"\\+"))
        compare=paste(com1,com2,sep="_vs_")
        compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
        if(!file.exists(compare_dir)){
            dir.create(compare_dir)
        }
				all_diff = data.frame()
        for(sub_cluster in unique(seurat_obj$celltype)){
            print(sub_cluster)
            cp1=paste(com1_str,paste("cluster",sub_cluster,sep=""),sep="_")
            cp2=paste(com2_str,paste("cluster",sub_cluster,sep=""),sep="_")
            if(!any(cp1 %in% cluster_group) || !any(cp2 %in% cluster_group)){
                next;								
            }
            ncell1=unlist(CellsByIdentities(seurat_obj,cp1))
            ncell2=unlist(CellsByIdentities(seurat_obj,cp2))
            if (length(ncell1)<3 || length(ncell2)<3){
                next;
            }

            cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
            if(!file.exists(cluster_compare_dir)){
                dir.create(cluster_compare_dir)
            }
            cluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment",sep="/")
						allcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/all",sep="/")
            upcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/up",sep="/")
            downcluster_compare_dir_enrich=paste(cluster_compare_dir,"enrichment/down",sep="/")
            if(!file.exists(upcluster_compare_dir_enrich)){dir.create(upcluster_compare_dir_enrich,recursive =TRUE)}
            if(!file.exists(downcluster_compare_dir_enrich)){dir.create(downcluster_compare_dir_enrich,recursive =TRUE)}
						if(!file.exists(allcluster_compare_dir_enrich)){dir.create(allcluster_compare_dir_enrich,recursive =TRUE)}
            diff.cluster=FindMarkers(seurat_obj, ident.1 = cp1, ident.2 = cp2, verbose = FALSE,min.pct = 0.10,logfc.threshold=avg_log2FC)
					  if (!is.null(species.org)) {
						diff.cluster$description= unname(mapIds(x = species.org,
									keys = rownames(diff.cluster),
									column = "GENENAME",
									keytype = "SYMBOL",
									multiVals = "first"))}
            write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            # sub_compare_cluster1=paste(cmpdf[x,1],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
            # sub_compare_cluster2=paste(cmpdf[x,2],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
            # diff.cluster$gene=row.names(diff.cluster)
						diff.cluster$cluster = sub_cluster
						diff.cluster <- diff.cluster %>% rownames_to_column(var = "gene")
	          all_diff = rbind(all_diff,diff.cluster)
            #diff.cluster=arrange(diff.cluster, desc(avg_log2FC))
            #top3_diff_gene=diff.cluster$gene[1:3]
            updiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC>0.25)
						downdiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC< -0.25)
            # write.table(updiff.cluster,paste(upcluster_compare_dir_enrich,"diff_gene.enrichment.up.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            # head(updiff.cluster)
            # downdiff.cluster<-subset(diff.cluster,p_val<0.05 & avg_log2FC< -0.25)
            # write.table(downdiff.cluster,paste(downcluster_compare_dir_enrich,"diff_gene.enrichment.down.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
            # head(downdiff.cluster)
						#筛选对应的上下调和所有基因
            updiff_genes=updiff.cluster$gene
            downdiff_genes=downdiff.cluster$gene
						diff.cluster_genes=diff.cluster$gene
						##富集
            try(enrichment(species=species,outDir=upcluster_compare_dir_enrich,geneList=updiff_genes))
            try(enrichment(species=species,outDir=downcluster_compare_dir_enrich,geneList=downdiff_genes))
						try(enrichment(species=species,outDir=allcluster_compare_dir_enrich,geneList=diff.cluster_genes))
						for (enrichment_res in c(upcluster_compare_dir_enrich, downcluster_compare_dir_enrich, allcluster_compare_dir_enrich)) {
							kegg_data = perform_KEGG_enrichment(species.org, diff.cluster_genes, updiff_genes, downdiff_genes, enrichment_res, species)
							if(is.null(kegg_data)){
								print(paste0(enrichment_res,"kegg file is empty"))
							}
						}						
        }
				print("检查报错")
				all_diff_table = file.path(compare_dir,"All_diff_gene.xls")
				write.table(all_diff,all_diff_table,sep='\t',row.names = F,quote=F)
				#system(glue::glue("Rscript /PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/subcluster/PlotKEGGNet.R -m  {all_diff_table} -o {seurat_diff_cluster_dir} -s {species} -t {compare} -n {topn}" ))
    }

}