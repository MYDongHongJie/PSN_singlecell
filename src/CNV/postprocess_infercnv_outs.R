library( phylogram )
library( dendextend ) 
#library( gplots )
library( data.table )
#library( infercnv )
library( DECIPHER )

## standard output from InferCNV

infercnv.observations <- "infercnv.observations_dendrogram.txt"

# get and visualize the dendrogram
dend <- as.hclust(ReadDendrogram(file = infercnv.observations,))
plot(dend,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

length(dend$labels)

# cut the dendrogram 
NCLONES <- 4

tp <- color_branches(dend, k = NCLONES,groupLabels = T)
plot(tp,leaflab = 'none',horiz = F) # set horiz = TRUE to plot dendrogram as in infercnv heatmap

g <- cutree(dend,k = NCLONES)
#barcodes <- unlist(lapply(names(g),clean_barcode))

groups <- data.frame(BARCODE=names(g),#BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)


## custom analysis from InferCNV final matrix

# Get the the tumor cell matrix data values
input <- 'infercnv.observations.txt'
finalmat <- read.table(file = input,header = T,row.names = 1)
dim(finalmat)

genes <- read.delim('gene_ordering_file.txt',as.is = T,stringsAsFactors = F,header = F)
colnames(genes) <- c("gene_id","chr","start","end")

keep_genes <- which(rownames(finalmat) %in% unique(genes$gene_id[which(genes$chr %in% c('chr5','chr6','chr10','chr18'))]))

m <- t(as.matrix(finalmat[keep_genes,]))

# which(apply(m, 1, var)==0)
# m <- m[ , apply(m, 1, var) != 0]
# pca <- prcomp(t(m),scale. = T)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 10)

heatmap.2(m,
          col = my_palette,
          Rowv = TRUE,
          Colv = FALSE,
          trace="none",
          scale="none",
          dendrogram = "row",
          labCol=FALSE,labRow = FALSE)

hc.rows<- hclust(dist(m))

plot(hc.rows,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

# cut the dendrogram 
NCLONES <- 5

tp <- color_branches(hc.rows, k = NCLONES,groupLabels = T)
plot(tp,leaflab = 'none',horiz = F) # set horiz = TRUE to plot dendrogram as in infercnv heatmap

g <- cutree(hc.rows,k = NCLONES)
barcodes <- unlist(lapply(names(g),clean_barcode))

groups <- data.frame(BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)

#save(groups,file = 'TS6ZX9_Barcodes_Clones_custom_clustering.RData',compress = TRUE)
