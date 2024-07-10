DoubletFinder <- function(sample_list,ifnb.list,seurat_qc_dir){
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
												x <- Seurat::NormalizeData(x)
        							  x <- Seurat::FindVariableFeatures(x)
          							x <- Seurat::ScaleData(x)
                        x <- RunPCA(x)
                        x <- RunUMAP(x, dims = 1:10)
                        x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
                        x <- FindClusters(x)
                        ##doublet finder need test
                        annotations <- x@meta.data$seurat_clusters
                        sweep.res.list <- paramSweep_v3(x, PCs = 1:10,sct=T)
                        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
                        bcmvn <- find.pK(sweep.stats)
                        mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
                        ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
                        homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
                        double_rate <- 0.0005272+0.000007589*length(x$sample)
                        nExp_poi <- round(double_rate*length(x$sample))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
                        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
                        x <- my_doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE,sct = T)
                        x <- my_doubletFinder_v3_1(x, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = "DF.classifications_raw", sct = T)
                        x@meta.data[,"Double_status"] <- x@meta.data$DF.classifications_raw
                        sample_name=unique(x$sample) #get sample name
                        out_gene_umi=as.matrix(x@meta.data)
                        row_name=c("Cell",rownames(out_gene_umi))
                        out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
                        out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
                        write.table(out_gene_umi,file=paste0(seurat_qc_dir,"/",sample_name,".cells_doublet_info.tsv"),sep="\t")
                        p1<-DimPlot(x, group.by="seurat_clusters")
                        p2<-DimPlot(x, group.by="Double_status",order=c("Singlet","Doublet"),cols =c("gray","red"))
                        p1+p2
                        ggsave(paste0(seurat_qc_dir,"/",sample_name,".cells_qc_filter.pdf"),width =16,height = 7)
                        x<-subset(x,Double_status=="Singlet")
                        x<-x
                        }
                 )
    return(ifnb.list)
}

CellCycle <- function(immune.combined,species,seurat_qc_dir){
    flag <- FALSE
    if(file.exists(paste0("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/cyclegene/",species,".cycle.gene.txt"))){
        cycledf <- read.table(paste0("/PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/public/cyclegene/",species,".cycle.gene.txt"),header=T,sep="\t")
        immune.combined <- CellCycleScoring(immune.combined, s.features = cycledf$s.genes, g2m.features = cycledf$g2m.genes, set.ident = FALSE)
        #DimPlot(immune.combined,reduction="umap",group.by="Phase",cols=colors)
        immune.combined@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()
				preprocess = file.path(seurat_qc_dir,'1.preprocess')
				if(!file.exists(preprocess)){
        dir.create(preprocess)
    		} 
        ggsave(paste(preprocess,"cells_cell_cycle.pdf",sep="/"),width = 8,height = 7)
        ggsave(paste(preprocess,"cells_cell_cycle.png",sep="/"),width = 8,height = 7)
        p<-dittoBarPlot(immune.combined,"Phase", group.by = "sample",color.panel=colors)
        ggsave(p,filename=paste(preprocess,"cell_cycle_per_sample.pdf",sep="/"))
        ggsave(p,filename=paste(preprocess,"cell_cycle_per_sample.png",sep="/"))
				
				
        immune.combined$CC.Difference <- immune.combined$S.Score - immune.combined$G2M.Score
        #immune.combined <- ScaleData(immune.combined, vars.to.regress = "CC.Difference", features = rownames(immune.combined))
        flag <- TRUE
    }
    if(flag){	
						print('去除周期影响')
            immune.combined$CC.Difference <- immune.combined$S.Score - immune.combined$G2M.Score
            immune.combined <- ScaleData(immune.combined, vars.to.regress = "CC.Difference",features = VariableFeatures(immune.combined))#features = rownames(immune.combined),)
    }else{
				print('不去除周期影响')
        immune.combined <- ScaleData(immune.combined,features = VariableFeatures(immune.combined))
    }
    return(immune.combined)
}

BatchEffect <- function(method,species,ifnb.list,sample_list,seurat_qc_dir){
    
    print(sample_list)
    if(method=="harmony"){
        if(length(ifnb.list) > 1){
            immune.combined <- merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_list) %>% JoinLayers() %>% NormalizeData()
        }else{immune.combined <- ifnb.list[[1]]%>%NormalizeData()}
        #VariableFeatures(immune.combined) <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
				
				immune.combined=Seurat::FindVariableFeatures(immune.combined, loess.span = 0.3,
                              clip.max = "auto", mean.function = "FastExpMean",
                              dispersion.function = "FastLogVMR", num.bin = 20,
                              nfeature = 3000, binning.method = "equal_width" )

        immune.combined <- CellCycle(immune.combined,species,seurat_qc_dir)
        
        immune.combined <- immune.combined %>% RunPCA(verbose = FALSE) %>% RunHarmony( group.by.vars = "sample")
        immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:20)
        immune.combined <- RunTSNE(immune.combined, reduction = "harmony", dims = 1:20)
        immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:20) %>% FindClusters()
    }
    if(method=="CCA"){
        ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
        features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
        ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
        immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",anchor.features = features)
        immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
        DefaultAssay(immune.combined) <- "integrated"
        immune.combined <- CellCycle(immune.combined,species,seurat_qc_dir)
        immune.combined <- RunPCA(immune.combined,  features = VariableFeatures(object = immune.combined) ,verbose = FALSE)
        immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
        immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
        immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
        immune.combined <- FindClusters(immune.combined)
    }
    if(method=="xbatch"){
       if(length(ifnb.list) > 1){
            immune.combined <- merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_list) %>% JoinLayers() %>% NormalizeData()
        }else{immune.combined <- ifnb.list[[1]]%>%NormalizeData()}     
        immune.combined=Seurat::FindVariableFeatures(immune.combined, loess.span = 0.3,
                              clip.max = "auto", mean.function = "FastExpMean",
                              dispersion.function = "FastLogVMR", num.bin = 20,
                              nfeature = 3000, binning.method = "equal_width" )
        immune.combined <- CellCycle(immune.combined,species,seurat_qc_dir)
        immune.combined <- immune.combined %>% RunPCA(verbose = FALSE) %>% RunUMAP(reduction = "pca", dims = 1:20)
        immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
        immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20) %>% FindClusters()    
    }
    return(immune.combined)
}
