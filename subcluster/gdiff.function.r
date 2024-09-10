#plot fig
library(dplyr)
groupDiffAuto <- function(seurat_obj,file_out,annocol,species,avg_log2FC){
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
                       write.table(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
                       #sub_compare_cluster1=paste(unique(sample_group)[idx],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
                       #sub_compare_cluster2=paste(unique(sample_group)[idx1],paste("cluster",unique(seurat_obj$celltype),sep=""),sep="_")
                       #diff.cluster$gene=row.names(diff.cluster)
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
                   }
									 write.table(all_diff,file.path(compare_dir,"All_diff_gene.xls"),sep='\t',row.names = F,quote=F)
      	   }
        }

      }
			  
    }

}

groupDiffSpeci <- function(seurat_obj,file_out,annocol,cmpfile,species,avg_log2FC){
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
        }
				write.table(all_diff,file.path(compare_dir,"All_diff_gene.xls"),sep='\t',row.names = F,quote=F)
    }

}

groupDiffSpeciAll <- function(seurat_obj,file_out,annocol,cmpfile,species,avg_log2FC){
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
        }
				write.table(all_diff,file.path(compare_dir,"All_diff_gene.xls"),sep='\t',row.names = F,quote=F)
    }
}


groupDiffSpeci_Auto <- function(seurat_obj,file_out,annocol,cmpfile,species,avg_log2FC){
    
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
        }
				write.table(all_diff,file.path(compare_dir,"All_diff_gene.xls"),sep='\t',row.names = F,quote=F)
    }

}