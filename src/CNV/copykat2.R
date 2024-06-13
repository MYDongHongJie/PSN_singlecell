#' copycat main_func.
#'
#' @param rawmat raw data matrix; genes in rows; cell names in columns.
#' @param  id.type gene id type: Symbol or Ensemble.
#' @param  cell.line if the data are from pure cell line,put "yes"; if cellline data are a mixture of tumor and normal cells, still put "no".
#' @param LOW.DR minimal population fractoins of genes for smoothing.
#' @param UP.DR minimal population fractoins of genes for segmentation.
#' @param win.size minimal window sizes for segmentation.
#' @param norm.cell.names a vector of normal cell names.
#' @param KS.cut segmentation parameters, input 0 to 1; larger looser criteria.
#' @param sam.name sample name.
#' @param n.cores number of cores for parallel computing.
#' @param ngene.chr minimal number of genes per chromosome for cell filtering.
#' @param distance  distance methods include euclidean, and correlation coverted distance include pearson and spearman.
#' @return 1) aneuploid/diploid prediction results; 2) CNA results in 220KB windows; 3) heatmap; 4) hclustering object.
#'
#' @examples
#' test.ck <- copykat(rawmat=rawdata, sam.name="test", n.cores=10)
#'
#' test.pred <- test.ck$prediction
#' @export


copykat <- function(rawmat=rawdata, id.type="S", cell.line="no", ngene.chr=5,LOW.DR=0.05, UP.DR=0.1, win.size=25, norm.cell.names="", KS.cut=0.1, sam.name="", distance="euclidean", n.cores=1){
  start_time <- Sys.time()
  set.seed(1)
  sample.name <- paste(sam.name,"_copykat_", sep="")

  print("running copykat v1.0.3")
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))

  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))

  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
  }
  ##
  der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)

  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
  }

  WNS1 <- "data quality is ok"
  if(nrow(rawmat) < 7000){
    WNS1 <- "low data quality"
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }

  print("step 2: annotations gene coordinates ...")
  anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  anno.mat <- anno.mat[order(anno.mat$abspos, decreasing = FALSE),]

# print(paste(nrow(anno.mat)," genes annotated", sep=""))

  ### module 3 removing genes that are involved in cell cycling
  HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
  toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
  if(length(toRev)>0){
    anno.mat <- anno.mat[-toRev, ]
  }

#  print(paste(nrow(anno.mat)," genes after rm cell cycle genes", sep=""))
  ### secondary filtering
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<23|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }

  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")

  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }

 # print(paste("filtered out ", length(ToRemov2), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat) <-  colnames(rawmat3)

  #print(paste("A total of ", ncol(norm.mat), " cells, ", nrow(norm.mat), " genes after preprocessing", sep=""))

  ##smooth data
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }

  test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)

  colnames(norm.mat.smooth) <- colnames(norm.mat)

  print("step 4: measuring baselines ...")
  if (cell.line=="yes"){
  	print("running pure cell line mode")
  	    relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
		norm.mat.relat <- relt$expr.relat
		CL <- relt$cl
        WNS <- "run with cell line mode"
    	preN <- NULL

  } else if(length(norm.cell.names)>1){

    #print(paste(length(norm.cell.names), "normal cells provided", sep=""))
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", sep=""))

    if (NNN==0) stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")

    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")

    d <- parallelDist::parDist(t(norm.mat.smooth),threads =n.cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells

    km <- 6
    fit <- hclust(d, method="ward.D2")
    CL <- cutree(fit, km)

    while(!all(table(CL)>5)){
      km <- km -1
      CL <- cutree(fit, k=km)
      if(km==2){
        break
      }
    }

    WNS <- "run with known normal"
    preN <- norm.cell.names
     ##relative expression using pred.normal cells
  	norm.mat.relat <- norm.mat.smooth-basel

  }else {
      basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
      basel <- basa$basel
      WNS <- basa$WNS
      preN <- basa$preN
      CL <- basa$cl
      if (WNS =="unclassified.prediction"){
        Tc <- colnames(rawmat)[which(as.numeric(apply(rawmat[which(rownames(rawmat) %in% c("PTPRC", "LYZ", "PECAM")),],2, mean)) >1)]; length(Tc)
        preN <- intersect(Tc, colnames(norm.mat.smooth))

        if(length(preN)> 5){
          print("start manual mode")
          WNS <- paste("copykat failed in locating normal cells; manual adjust performed with ", length(preN), " immune cells", sep="")
          print(WNS)
          basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% preN)], 1,mean)

            }else{
                    basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
                    basel <-basa$basel
                    WNS <- basa$WNS
                    preN <- basa$preN

            }
      }
  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.smooth-basel

  }

  ###use a smaller set of genes to perform segmentation
  DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)
  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.relat[which(DR2>=UP.DR),]

  ###filter cells
  anno.mat2 <- anno.mat[which(DR2>=UP.DR), ]

  ToRemov3 <- NULL
  for(i in 8:ncol(anno.mat2)){
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    } else if(length(rle(cell[,1])$length)<23|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i<- i+1
  }

  if(length(ToRemov3)==ncol(norm.mat.relat)) stop ("all cells are filtered")

  if(length(ToRemov3)>0){
    norm.mat.relat <-norm.mat.relat[, -which(colnames(norm.mat.relat) %in% ToRemov3)]
   # print(paste("filtered out ", length(ToRemov3), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  }

  #print(paste("final segmentation: ", nrow(norm.mat.relat), " genes; ", ncol(norm.mat.relat), " cells", sep=""))

  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]

  print("step 5: segmentation...")
  results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*0.5*KS.cut, n.cores=n.cores)
  }

  if(length(results$breaks)<25) stop ("too few segments; try to decrease KS.cut; or improve data")

  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)

  print("step 6: convert to genomic bins...") ###need multi-core
  Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat=RNA.copycat, n.cores = n.cores)

  uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])

  print("step 7: adjust baseline ...")

if(cell.line=="yes"){

  mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
  write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.xls", sep=""), sep="\t", row.names = FALSE, quote = F)

  if(distance=="euclidean"){
    hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
  }else {
    hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
  }


  saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

#plot heatmap
 print("step 8: ploting heatmap ...")
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)


  if (ncol(mat.adj)< 3000){
    h <- 10
  } else {
    h <- 15
  }

  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
library(parallelDist)

  # if(distance=="euclidean"){
  # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
   # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
            # ColSideColors=chr1,Colv=NA, Rowv=TRUE,
            # notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            # keysize=1, density.info="none", trace="none",
            # cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            # symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  # dev.off()
  # } else {
    # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
    # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                 # ColSideColors=chr1,Colv=NA, Rowv=TRUE,
              # notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              # keysize=1, density.info="none", trace="none",
              # cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              # symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

    # dev.off()
  # }
  end_time<- Sys.time()
  print(end_time -start_time)

  reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
  names(reslts) <- c("CNAmat","hclustering")
  return(reslts)
} else {
  ################removed baseline adjustment
  if(distance=="euclidean"){
  hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
  }else {
  hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
  }
  hc.umap <- cutree(hcc,2)
  names(hc.umap) <- colnames(results.com)

  cl.ID <- NULL
  for(i in 1:max(hc.umap)){
    cli <- names(hc.umap)[which(hc.umap==i)]
    pid <- length(intersect(cli, preN))/length(cli)
    cl.ID <- c(cl.ID, pid)
    i<- i+1
  }

  com.pred <- names(hc.umap)
  com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
  com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "nondiploid"
  names(com.pred) <- names(hc.umap)

  ################removed baseline adjustment
  results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
  results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
  results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)

  cf.h <- apply(results.com.rat.norm, 1, sd)
  base <- apply(results.com.rat.norm, 1, mean)

  adjN <- function(j){
    a <- results.com.rat[, j]
    a[abs(a-base) <= 0.25*cf.h] <- mean(a)
    a
  }


  mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
  adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
  colnames(adj.results) <- colnames(results.com.rat)

  rang <- 0.5*(max(adj.results)-min(adj.results))
  mat.adj <- adj.results/rang

  print("step 8: final prediction ...")

  if(distance=="euclidean"){
    hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
  }else {
    hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
  }

  hc.umap <- cutree(hcc,2)
  names(hc.umap) <- colnames(results.com)

  saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))

  cl.ID <- NULL
  for(i in 1:max(hc.umap)){
    cli <- names(hc.umap)[which(hc.umap==i)]
    pid <- length(intersect(cli, preN))/length(cli)
    cl.ID <- c(cl.ID, pid)
    i<- i+1
  }

  com.preN <- names(hc.umap)
  com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
  com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
  names(com.preN) <- names(hc.umap)

  if(WNS=="unclassified.prediction"){
    com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
    com.preN[which(com.preN == "nondiploid")] <- "c2:aneuploid:low.conf"
  }

  print("step 9: saving results...")
  res <- cbind(names(com.preN), com.preN)
  colnames(res) <- c("cell.names", "copykat.pred")

  write.table(res, paste(sample.name, "prediction.xls",sep=""), sep="\t", row.names = FALSE, quote = FALSE)

  ####save copycat CNA
  write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.xls", sep=""), sep="\t", row.names = FALSE, quote = F)


  ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
  # print("step 10: ploting heatmap ...")
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)

  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]

  cells <- rbind(compreN_pred,compreN_pred)

  if (ncol(mat.adj)< 3000){
    h <- 10
  } else {
    h <- 15
  }

  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

  # if(distance=="euclidean"){
  # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
   # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
            # ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            # notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            # keysize=1, density.info="none", trace="none",
            # cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            # symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

  # legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
  # dev.off()
  # } else {
    # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
    # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
                 # ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
              # notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
              # keysize=1, density.info="none", trace="none",
              # cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
              # symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))

    # legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)

    # dev.off()
  # }
  end_time<- Sys.time()
  print(end_time -start_time)

  reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
  names(reslts) <- c("prediction", "CNAmat","hclustering")
  return(reslts)
}

}


#' Creates a list of unique color values used for plotting
#'
#' @return A named vector of unique hexedecimal color values, either generated from a preselected
#'         vector of 20 unique colors, or from a sequence of colors in hsv colorspace.
#'
#' @param seurat.obj A singular preprocessed Seurat object
#' @param gradient Setting to TRUE will use a sequence of hsv colors instead of 20 unique colors,
#'                 useful for comparisons of more than 20 cell types.
#' @param value The Seurat metadata slot to generate colors for. Defaults to "celltype".
#'
#' @import SingleCellExperiment
#' @import Seurat
#'
#' @seealso \code{\link{as.SingleCellExperimentList}}
#' @seealso \code{\link{ExtractGenes}}
#' @seealso \code{\link{DecoderVariance}}
#' @seealso \code{\link{MeanDecoderVariance}}
#' @seealso \code{\link{GetCharMetadata}}
#'
#' @examples
#' DimPlot(object = seurat.obj,
#'         reduction = "tsne",
#'         cols = SelectColors(seurat.obj),
#'         group.by = "celltype",
#'         label = TRUE,
#'         repel = TRUE)
#'
#' @export
SelectColors <- function(
    object = NULL,
    palette = "blindless",
    value = "celltype",
    n = NULL
){
    if ( !is.null(object) ){
        if ( class(object) == "data.frame" ){
            colid <- ifelse( is.null(value), colnames(object)[1], value )
            if (is.factor(object[[colid]])) {
                names= levels(object[[colid]])
            } else {
                names <- unique(object[[colid]])
            }
        }
        if ( is.factor(object) ){
            names <- levels(object)
        }
        if ( is.vector(object) ){
            names <- unique(object)
        }
        n = length(names)
    }else if ( !is.null(n) ) {
        names = NULL
    }

    colors2pick = switch(palette,
    ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    ditto = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#666666", "#AD7700", "#1C91D4",
    "#007756", "#D5C711", "#005685", "#A04700", "#B14380",
    "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71",
    "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C"),
    CustomCol2 = c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
    "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
    "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
    "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17"),
    seurat = hcl( h = seq(15, 375, length = n+1), l = 65, c = 100),
    col50 =  c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
    "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
    "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
    "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
    "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
    "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
    "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
    "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
    "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
    "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"),
    paired = brewer.pal(n = n, 'Paired'),
    colx22 = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
    '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
    '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
    '#000075', '#808080', '#4f34ff', '#f340F0'),
    jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
    "#FF7F00", "red", "#7F0000" ),
    tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
    "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
    "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
    "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
    tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
    "#CDCC5D", "#6DCCDA"),
    colorblind10 = c("#006BA4", "#FF800E", "#ABABAB", "#595959",
    "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
    "#FFBC79", "#CFCFCF"),
    trafficlight = c("#B10318", "#DBA13A", "#309343", "#D82526",
    "#FFC156", "#69B764", "#F26C64", "#FFDD71",
    "#9FCD99"),
    purplegray12 = c("#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
    "#5F5A41", "#B4B19B", "#995688", "#D898BA",
    "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"),
    bluered12 = c("#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
    "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
    "#BD0A36", "#F4737A"),
    greenorange12 = c("#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
    "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
    "#39737C", "#86B4A9", "#82853B", "#CCC94D"),
    cyclic = c("#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
    "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
    "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
    "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"),
    CustomCol = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
    blindless = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
    "#BF5B17", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4",
    "#B3DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#BC80BD",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
    "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
    "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8",
    "#984EA3",  "#FFFF33", "#A65628", "#F781BF", "#999999","#FFED6F",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9","#666666"),
    )

    if ( is.null(n) ){
        colors_use <- colors2pick
    }else{
        colors_use <- colors2pick[1:n]
    }
    if ( !is.null(names) ){ names(colors_use) <- names }
    return(colors_use)
}


#=================================================================================
# library packages
#=================================================================================
suppressPackageStartupMessages( library("MCMCpack",lib.loc='/data/software/conda_envs/scrna_envs/copykat/lib/R/library') )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("Seurat") )
suppressPackageStartupMessages( library("ComplexHeatmap") )
suppressPackageStartupMessages( library("RColorBrewer") )
suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("tibble") )
suppressPackageStartupMessages( library("stringr") )
suppressPackageStartupMessages( library("circlize") )
suppressPackageStartupMessages( library("ggplot2") )
suppressPackageStartupMessages( library("copykat",lib.loc='/data/software/conda_envs/scrna_envs/copykat/lib/R/library') )
suppressPackageStartupMessages( library("dlm",lib.loc='/data/software/conda_envs/scrna_envs/copykat/lib/R/library') )
suppressPackageStartupMessages( library("matrixStats") )
#suppressPackageStartupMessages( library("MCMCpack",lib.loc='/data/software/conda_envs/scrna_envs/copykat/lib/R/library') )
#=================================================================================
# command line parameters setting
#=================================================================================
option_list = list(
    make_option( c("--RDS", "-i"), type = "character",
                 help = "[Required] The seurat object saved as R object in RDS format."),
    make_option( c("--celltype", "-l"), type = "character",
                 help = "[Required] The cell type annotation column name to use in seurat metadata."),
    make_option( c("--groupby", "-g"), type = "character",
                 help = "[Required] The groupping of cells in the metadata of seurat object to visualize."),
    make_option( c("--colormapping", "-m"), type = "character",
                 help = "[OPTIONAL] The color mapping for groupping column of cells set by the parameters '--groupby'.
                         The exmaple format is variable1:colorschema1,variable2:colorschema2.
                         The supported color schemas can be: blindless, col50, ditto, paired, CustomCol2."),
    make_option( c("--malignant", "-t"), type = "character",
                 help = "[OPTIONAL] The cell type name to be assumed as cancer cells"),
    make_option( c("--normalcells", "-r"), type = "character",
                 help = "[OPTIONAL] A comma seperated list containing the classifications of normal cells."),
    make_option( c("--assay", "-e"), type = "character", default = "RNA",
                 help = "[OPTIONAL] The array result to use in case of multimodal analysis."),
    make_option( c("--reduct"), type = "character", default = "tsne",
                 help = "[OPTIONAL] The reduction used in the DimPlot."),
    make_option( c("--pointsize"), type = "double", default = "0.8",
                 help = "[OPTIONAL] The pointsize used in the DimPlot."),
    make_option( c("--sample_ratio", "-s"), type = "double", default = 0.6,
                 help = "[OPTIONAL] The ratio of random subsample for each group. Only normal cells will be subseted.
                          The number of all cells must be less than 60,000 cells, or it will out of memory in step 4: measuring baselines..."),
    make_option( c("--outdir","-o"),type="character", default = "./",
                 help="[OPTIONAL] The output directory of results."),
    make_option( c("--ident2use", "-q" ), type = "character", default = NULL,
                 help = "[OPTIONAL] The column name in cell metadata used as identity of each cell combined with which_cell."),
    make_option( c("--which_cells", "-u" ), type = "character", default = NULL,
                 help = "[OPTIONAL] The subset of cluster ids used for analysis."),
    make_option( c("--ncores", "-j" ), type="integer", default = 10,
                 help="[OPTIONAL] the number of CPUs used to improve the performace.")
    );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$RDS) ){
    stop("the seurat object is NOT AVAILABLE!")
}else{
    seurat_ob = readRDS(opt$RDS)
    if ( seurat_ob@version < 3){
        seurat_ob = UpdateSeuratObject(seurat_ob) #make sure the seurat object match with the latest seurat package
    }
    # change the default assay for reduction if necessary
    if ( !is.null( opt$assay) ){
        DefaultAssay(seurat_ob) = opt$assay
    }else{
        DefaultAssay(seurat_ob) = "RNA"
    }
}

# get the subset of cells used for visualization if necessay
if ( !is.null(opt$which_cells)){
    cluster_list = unlist(strsplit( opt$which_cells,",",perl = T))
    seurat_ob = SubsetData(seurat_ob,  subset.name = opt$ident2use, accept.value = cluster_list)
}

# output directory setting
if ( is.null(opt$outdir) ){
    print("NO output directory specified, the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$outdir) ){
        output_dir = opt$outdir
    }else{
        output_dir = opt$outdir
        dir.create(output_dir, recursive = T)
    }
}
output_dir = normalizePath(output_dir)

if ( !is.null( opt$groupby) ){
    cell.annos = unlist( strsplit( opt$groupby, ",", perl =T))
    groupby = unlist( strsplit( opt$groupby, ",", perl =T))[1]
}else{
    stop("please provide the groupping of cells in the metadata of seurat object.")
}

# set the color schema of cell annotation bar
if ( !is.null( opt$colormapping) ){
    group_colors = list()
    for ( x in  unlist(strsplit(opt$colormapping, ",", perl =T))){
        m = unlist( strsplit( x, ":", perl =T))
        group_colors[[m[1]]] = m[-1]
    }
}else{
    group_colors= 1:(length(cell.annos)+1)
    names(group_colors) = c(cell.annos, "copykat.pred")
    group_colors = as.list(group_colors)
}

if ( !is.null( opt$normalcells ) ){
    normalcell = unlist(strsplit(opt$normalcells, ",", perl = T))
    normal.cells = rownames(seurat_ob@meta.data[which(seurat_ob@meta.data[, opt$celltype] %in% normalcell),] )
}else if ( !is.null( opt$malignant ) ){
    malignant = unlist(strsplit(opt$malignant, ",", perl = T))
    malignant.cells = rownames(seurat_ob@meta.data[which(seurat_ob@meta.data[, opt$celltype] %in% malignant),] )
    normal.cells = setdiff( Cells(seurat_ob), malignant.cells)
}else{
    print("NO normal cells or malignant cells are specified!")
    print("CopyKAT will predict normal cells automatically!")
    normal.cells = ""
}

## Needs to be less than 60,000 cells, or it will out of memory in step 4: measuring baselines...
if ( !is.null(opt$sample_ratio) ){
  if( normal.cells[1] != "" ){
    sampled_cellmeta = seurat_ob@meta.data[normal.cells,] %>%
                        rownames_to_column() %>%group_by( .dots= opt$celltype ) %>%
                        sample_frac( size = opt$sample_ratio, replace = F) %>% column_to_rownames()
    subset_cells = c(rownames(sampled_cellmeta), setdiff(Cells(seurat_ob), normal.cells))
    seurat_ob = SubsetData(seurat_ob, cells = subset_cells)
  }else{
    sampled_cellmeta = seurat_ob@meta.data %>% rownames_to_column() %>%
                        group_by( .dots= opt$celltype ) %>%
                        sample_frac( size = opt$sample_ratio, replace = F) %>% column_to_rownames()
    seurat_ob = SubsetData(seurat_ob, cells = rownames(sampled_cellmeta))
  }
}

##
setwd(output_dir)
copykat.res <- copykat(rawmat=GetAssayData(seurat_ob, slot="counts"), id.type="S", cell.line="no", ngene.chr=5,
                       win.size=25, KS.cut=0.1, sam.name="seurat", distance="euclidean",
                       norm.cell.names=normal.cells, n.cores=opt$ncores)

CNA_mat = data.frame(copykat.res$CNAmat)
cell.pred = data.frame(copykat.res$prediction)
cell.pred = cell.pred[,2, drop=F]
cell.pred = cell.pred %>% arrange(copykat.pred)

## TO DO (load existing results)
# CNA_mat= as.matrix(vroom::vroom("seurat_copykat_CNA_results.xls"))
# cell.pred = read.delim("seurat_copykat_prediction.xls",row.names=1)
# cell.pred <- cell.pred %>% tibble::rownames_to_column(var="barcodes") %>% arrange(copykat.pred) %>%  tibble::column_to_rownames("barcodes")

# cellanno = FetchData(seurat_ob, vars = opt$celltype ) %>% tibble::rownames_to_column(var = "cellbarcode")

meta.data = seurat_ob@meta.data[rownames(cell.pred),]
cell.anno = cbind(meta.data, cell.pred)
cell.anno = cell.anno[,c(cell.annos,"copykat.pred")]
cell.anno_data = tibble::rownames_to_column(as.data.frame(cell.anno),"barcode")
write.table(cell.anno_data, file.path(output_dir,"copykat_prediction.xls"), quote = F, row.names = F, sep = "\t")

gene.anno = as.data.frame(CNA_mat[,1, drop=F])
plot_data = t(CNA_mat[,4:ncol(CNA_mat)])

## TO DO (define subpopulations of aneuploid tumor cells)
# tumor.cells <- rownames(cell.pred[which(cell.pred$copykat.pred=="aneuploid"),])
# tumor.mat <- CNA_mat[, which(colnames(CNA_mat) %in% tumor.cells)]
# hcc <- hclust(parallelDist::parDist(t(tumor.mat), threads =10, method = "euclidean"), method = "ward.D2")
# hc.umap <- cutree(hcc,5)   ###

color_schema_row = list()
color_schema_col = list()

for ( x in colnames(cell.anno) ){
    nlevels <- length(unique(cell.anno[,x]))
    color_schema_row[[x]] <- SelectColors(cell.anno, value = x, n = nlevels, palette = group_colors[[x]])
}

color_schema_col[["chrom"]] <- SelectColors(object=gene.anno, value="chrom", n=length(unique(gene.anno[,"chrom"])), palette="ditto")

cellAnnotation = ComplexHeatmap::HeatmapAnnotation(
                                            df = cell.anno,
                                            col = color_schema_row,
                                            # annotation_legend_param = annotation_legend_params,
                                            which = "row",
                                            show_annotation_name = T)

geneAnnotation = ComplexHeatmap::HeatmapAnnotation(
                                            df = gene.anno, name = "chr",
                                            col = color_schema_col,
                                            # annotation_legend_param = annotation_legend_params,
                                            show_annotation_name = T)

col.split.by = "chrom"
coldata_split = CNA_mat[,col.split.by, drop =F]
plot_data = plot_data[rownames(cell.anno),]

pdf( file.path(output_dir,"copykat_heatmap.pdf"), width=10)
ComplexHeatmap::Heatmap(t(scale(t(plot_data))),name="CNV level",
        cluster_rows = FALSE, cluster_columns = FALSE,
        # clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
        col = colorRamp2(c(-2, 0, 2), c("#406AA8", "white", "#D91216")),
        # colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
                        # colorRamp2(c(-2, 0, 2), c("#053061", "white", "#67001F"))
        row_order = rownames(cell.anno),
        # row_split = rowdata_split,
        row_gap = unit(0, "mm"),
        # column_order = colnames(plot_data),
        column_split = coldata_split,
        column_gap = unit(0, "mm"),
        column_title = unique(CNA_mat[,"chrom"]),
        column_title_gp = gpar(fontsize = 7),
        column_title_side="bottom",
        border = T,
        # top_annotation = geneAnnotation ,
        left_annotation = cellAnnotation,
        use_raster=TRUE,
        raster_quality=4,
        show_heatmap_legend = TRUE,
        show_row_names =F,
        show_column_names = F)
dev.off()

## Dimplot
sub_ob = SubsetData(seurat_ob , cells=rownames(cell.anno))
sub_ob@meta.data$copykat.pred = cell.anno[rownames(sub_ob@meta.data),"copykat.pred"]
sub_ob = SetIdent(sub_ob, value="copykat.pred")
ggtsne = DimPlot(object = sub_ob, reduction = opt$reduct, pt.size = opt$pointsize ) + scale_colour_manual(values = c("#7fc97f","#beaed4"))
ggsave(file.path(output_dir,paste0("copykat_prediction_",opt$reduct,"_plot.png")), dpi=1000)
ggsave(file.path(output_dir,paste0("copykat_prediction_",opt$reduct,"_plot.pdf")) )

if(!file.exists(file.path(output_dir, "CopyKAT分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/CopyKAT分析说明.docx",
  file.path(output_dir, "CopyKAT分析说明.docx"))
}

