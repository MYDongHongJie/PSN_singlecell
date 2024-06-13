#!/usr/bin/env Rscript

#' make a dotplot of the cellphonedb results using the pvalues.txt and means.txt file
#'
#' @param pval a data.frame of the pvalues.txt or just the file path
#' @param mean a data.frame of the means.txt or just the file path
#' @param features the desired interacting gene pairs in the above to plot,all will be used if NULL
#' @param cells the desired interacting cell pairs to plot, all will be used if NULL
LRDotplot <- function(
  data,
  is_onlySig = T,
  topn = 5,
  xangle = 90,
  xsize = 10,
  palette = c("black", "blue", "yellow", "red")
){
  if( class(data) == "data.frame" ){
      all_data = data[,c("receptor","ligand","receptor_cell","ligand_cell","pval","expr","receptor_expr","ligand_expr")]
  }else{
      stop("NO cell communication matrix is procided!")
  }

  all_data = all_data %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")

  if ( is_onlySig ){
    desired_pair = ""
    for (i in unique(all_data$clusters)) {
      desired_pair2 <- all_data %>% dplyr::filter( pval < opt$pvalue & receptor_expr != 0 & ligand_expr != 0 ) %>% arrange( desc(expr)) %>% filter( clusters == i ) %>% pull( pair )
      desired_pair <- c(desired_pair,unique(desired_pair2)[1:topn])
    }
    filter_data = all_data %>% dplyr::filter( pair %in% unique(desired_pair))
  }else{
    filter_data = all_data
  }

  filter_data$pval[filter_data$pval==0] = 0.0009
  while ( length(filter_data$expr[filter_data$expr==0]) >0 ) {
          filter_data$expr[filter_data$expr==0][1] = (2^filter_data[filter_data$expr==0,][1,]$receptor_expr+2^filter_data[filter_data$expr==0,][1,]$ligand_expr)/2 
  }
  filter_data$mean = as.numeric(log2(filter_data$expr))

  plot.data = filter_data[,c("pair","clusters","pval","mean")]
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(palette, alpha=TRUE)(n=399)
  dotplot <- ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=xsize, colour = "black"),
          axis.text.x = element_text(size = xsize, angle = xangle, hjust = 1),
          axis.text.y = element_text(size = xsize, colour = "black"),
          axis.title=element_blank(),
          text = element_text(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  out=list(dotplot,plot.data)
  return(out)
}


SelectColors <- function(
  palette = "blindless",
  n = NULL
){
  colors2pick = switch(palette,
         seurat = hcl( h = seq(15, 375, length = n+1), l = 65, c = 100),
        ##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
         blindless = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
                     "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                     "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
                     "#B2DF8A", "g#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                     "#CAB2D6", "g#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
                     "#CCEBC5", "g#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                     "#F2F2F2", "g#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
                     "#FFF2AE", "g#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A",
                     "#984EA3", "g#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                     "#66C2A5", "g#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                     "#E5C494", "g#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                     "#80B1D3", "g#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                     "#CCEBC5", "#FFED6F"),
         CustomCol2 = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
                      "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
                      "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
                      "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
                      "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
                      "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
                      "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
                      "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
                      "#3fdacb","#bf5b17"),
         paried = brewer.pal(n = n, 'Paired'),
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
                    "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB")
  )

  colors_use <- colors2pick[1:n]
  return(colors_use)
}

#' Network Viewing of cell-cell communication
#'
#' This function loads the significant interactions as a dataframe, and colors
#' represent different types of cells as a structure. The width of edges represent
#' the strength of the communication. Labels on the edges show exactly how many
#' interactions exist between two types of cells.
#'
#' @references Csardi G, Nepusz T: The igraph software package for complex network
#' research, InterJournal, Complex Systems 1695. 2006.
#' http://igraph.org
#' @param data A dataframe containing ligand-receptor pairs and corresponding
#' cell types used to do the plotting
#' @param col Colors used to represent different cell types
#' @param label Whether or not shows the label of edges (number of connections
#' between different cell types)
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”,
#' “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see
#' vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the
#' plot.igraph command. “none” does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout, e.g.layout=in_circle(),layout.circle
#' specification function.
#' @param vertex.size The size of vertex
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.max.width The maximum edge width
#' @param edge.arrow.size The arrow size
#' @import network
#' @import igraph
#' @return A network graph of the significant interactions
#' @export
LRNetwork<-function(
  data,
  col,
  label=TRUE,
  edge.curved=0.2,
  shape='circle',
  layout=nicely(),
  vertex.size=20,
  margin=0.4,
  vertex.label.cex=1.5,
  vertex.label.color='black',
  arrow.width=1.5,
  edge.label.color='black',
  edge.label.cex=1,
  edge.max.width=10,
  edge.arrow.size=1
){
  set.seed(1234)
  net<-data %>% group_by(ligand_cell, receptor_cell) %>% dplyr::summarize(n=n())
  net<-as.data.frame(net,stringsAsFactors=FALSE)
  g<-graph.data.frame(net,directed=TRUE)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  loop.angle<-ifelse(coords_scale[V(g),1]>0,-
                       atan(coords_scale[V(g),2]/coords_scale[V(g),1]),
                     pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  vertex.size <- data %>% group_by(ligand_cell) %>% dplyr::summarize(n=n())
  vertex.size <- vertex.size$n*0.07
  V(g)$size<-vertex.size
  V(g)$color<-col[V(g)]
  V(g)$label.color<-vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  V(g)$frame.color <- NA
  if(label){
    E(g)$label<-E(g)$n
  }
  if(max(E(g)$n)==min(E(g)$n)){
    E(g)$width<-2
  }else{
    E(g)$width<-1+edge.max.width/(max(E(g)$n)-min(E(g)$n))*(E(g)$n-min(E(g)$n))
  }
  E(g)$arrow.width<-arrow.width
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
#  E(g)$color<-V(g)$color[edge.start[,1]]
  E(g)$color<-"#BDBDBD"
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=layout.circle,margin=margin,edge.arrow.size=edge.arrow.size)
  return(g)
}


#' Plotting ligand-receptor pairs
#'
#' This function loads the significant interactions as a dataframe.
#' A circle plot will be generated using package circlize. The colors range
#' from green, black to red, which represents the expression level of the
#' ligand and receptor is from low to high. Users can select the colors
#' represent the cell type by their own or chose randomly by default.
#' width of arrow head represents the expression level/log fold change
#' of the receptor. Different line type of the link stands for whether
#' the interactions between the ligand and receptor are strong or week.
#' And different color of the arrow stands for the different types of
#' the interantions.
#'
#' @param data A datafram contains significant ligand-receptor pairs and related
#' information such as expression level and cell type.
#' @param screenvar The variable is that filtering data based on.
#' @param topn The topn numbers of the data to retain after filtering.
#' @param pval The p-value threshold value of the interaction of the receptor-ligand pairs.
#' @param gap.degree Gap between two neighbour sectors. It can be a single value
#'  or a vector.
#' @param color_expr The colors display of the gene expression.
#' @param expr_range The threshold range display of the gene expression.
#' @param cell_col Colors used to represent types of cells. If set NULL, it will
#' be generated randomly.
#' @param labels.cex The font size of the cell track labels.
#' @param arr.length The length of the single arrow head which is put in the
#' center of the belt.
#' @param link.lty Line type of the single line link which is put in the center
#' of the belt.
#' @param link.col Color of the single line link which is put in the center of
#' the belt.
#' @param link.lwd Line width of the single line link which is put in the center
#' of the belt.
#' @param link.arr.width Size of the single arrow head link which is put in the
#' center of the belt.
#' @param link.arr.type Type of the arrows, pass to Arrowhead. Default value is
#' triangle. There is an additional option big.arrow.
#' @param facing Facing of text.
#' @param track.height_1 Height of the gene and ligang or receptor notation track.
#' @param track.height_2 Height of the cell notation track.
#' @param text.vjust Adjustment on ’vertical’ (radical) direction. Besides to set it
#' as numeric values, the value can also be a string contain absoute unit, e.g.
#' "2.1mm", "-1 inche", but only "mm", "cm", "inches"/"inche" are allowed.
#' @import circlize
#' @import dplyr
#' @import ComplexHeatmap Legend
#' @import RColorBrewer brewer.pal
#' @improt OESingleCell CustomCol
#' @import gridBase
#' @example
#' pdf("cell_comm_circos.pdf", width = 20, height = 10)
#' LRCircos(cellb$ligrec,
#'          gap.degree = 0.05,
#'          cell_col = colx,topn = 5,
#'          labels.cex = 0.2,
#'          link.lwd=1,
#'          arr.length = 0.1)
#' dev.off()
#' @return A figure of the significant interactions
#' @export
LRCircos <- function(
  data,
  screenvar="expr",
  topn=5,
  pval=0.05,
  gap.degree=3,
  color_expr=c("green","yellow","red"),
#  expr_range=c(-15,-6,3),
  cell_col=NULL,
  link.lty=NULL,
  labels.cex=0.58,
  arr.length=0.2,
  link.col=NULL,
  link.lwd=2,
  link.arr.width=NULL,
  link.arr.type=NULL,
  facing='clockwise',
  track.height_1=uh(8,'mm'),
  track.height_2=uh(12,'mm'),
  text.vjust=0.5,
  ...
  ){
  if(class(data) == "data.frame"){
    data = data
  }else if(file.exists(data)){
    data <- read.csv(data, header = T, sep = '\t')
  }else{
    stop("No cell communication result is provided!")
  }

  if(is.null(screenvar)){
    screenvar = "expr"
  }

  if(is.null(topn)){
    topn = 5
  }

  data <- data %>% group_by(receptor_cell,ligand_cell) %>% filter( pval < opt$pvalue ) %>% top_n(topn, !!ensym(screenvar))
  data <- as.data.frame(data)
  data <- data %>% dplyr::mutate(lr = 'ligand', lr2 = 'receptor')
  ldf <- data %>% dplyr::select(ligand_cell, lr, ligand, ligand_expr,
                receptor_cell, lr2, receptor, receptor_expr, everything()) %>%
          dplyr::rename("cell" = "ligand_cell", "gene" = "ligand", "expr1" = "ligand_expr",
                "cell2" = "receptor_cell", "gene2" = "receptor", "expr2" = "receptor_expr")
  rdf <- data %>% dplyr::select(receptor_cell, lr2, receptor, receptor_expr,
                ligand_cell, lr, ligand, ligand_expr, everything()) %>%
          dplyr::rename("cell" = "receptor_cell", "lr2" = "lr", "gene" = "receptor", "expr1" = "receptor_expr",
                "cell2" = "ligand_cell", "lr" = "lr2", "gene2" = "ligand", "expr2" = "ligand_expr")
  df <- rbind(ldf,rdf)
  df <- df[order(df$cell),]
  df <- df %>% dplyr::mutate(gene_id = NA, gene_id2 = NA)
  df$gene_id <- paste("geneid", 1:nrow(df))
  d1 <- data.frame(index = paste(df$cell,df$lr,df$gene,df$cell2,df$lr2,df$gene2,sep = '_'), gene_id = df$gene_id)
  d2 <- data.frame(index = paste(df$cell2,df$lr2,df$gene2,df$cell,df$lr,df$gene,sep = '_'), gene_id = df$gene_id2)
  for(i in 1:nrow(d1)){
    index = which(d2$index %in% d1[i,1])
    for(ind in index){
      d2[ind,2] = as.character(d1[i,2])
    }
  }
  df$gene_id2 <- d2$gene_id
  df <- df %>% dplyr::rename("ligand_cell" = "cell", "ligand" = "gene", "ligand_expr" = "expr1",
               "receptor_cell" = "cell2", "receptor" = "gene2", "receptor_expr" = "expr2") %>%
        select(ligand_cell, lr, ligand, gene_id, ligand_expr, receptor_cell, lr2,
               receptor, gene_id2, receptor_expr, everything())
  lrid <- df[which(df$lr == 'ligand'),]
  
  max_range = max(ceiling(max(lrid$ligand_expr)),ceiling(max(lrid$receptor_expr)))
  min_range = min(floor(min(lrid$ligand_expr)), floor(min(lrid$receptor_expr)))

  expr_range = c(min_range, (max_range+min_range)/2, max_range)

  if(is.null(color_expr)){
    stop("No color set of gene expression in the plot is provide!")
  }else{
    color_expr = color_expr
  }

  if(is.null(pval)){
    pval = 0.05
  }else{
    pval = pval
  }
  if(is.null(link.lty)){
    link.lty = structure(ifelse(lrid$pval < pval, 'solid', 'dashed'),
                             names = paste(lrid$ligand_cell,lrid$receptor))
  }
  if(is.null(lrid$comm_type)){
    lrid <-lrid %>% dplyr::mutate(link_col = 'black')
  }else{
    if(is.null(link.col)){
      comm_col <- structure(brewer.pal(length(unique(lrid$comm_type)),'Paired'),
                            names = as.character(unique(lrid$comm_type)))
      lrid <-lrid %>% dplyr::mutate(link_col = as.character(comm_col[as.character(lrid$comm_type)]))
    }else{
      lrid$link_col = link.col
    }
  }
  if(is.null(cell_col)){
    cell_col<-structure(c("#7FC97F","#BEAED4","#FDC086","#FBB4AE","#1B9E77","#FFFF99","#386CB0","#F0027F",
                          "#666666","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666",
                          "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
                          "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#B3CDE3","#CCEBC5","#DECBE4","#BF5B17",
                          "#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC","#CBD5E8",
                          "#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC","#E41A1C","#377EB8","#4DAF4A",
                          "#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                          "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3",
                          "#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
                          "#CCEBC5", "#FFED6F"),names=as.character(unique(df$receptor_cell)))
#    cell_col <- structure(CustomCol(1:length(unique(df$receptor_cell))),
#                          names = as.character(unique(df$receptor_cell)))
  }
  if(is.null(link.lwd)){
    link.lwd = 2
  }else{
    link.lwd = link.lwd
  }
  if(is.null(labels.cex)){
    labels.cex = 0.58
  }else{
    labels.cex = labels.cex
  }
  if(is.null(arr.length)){
    arr.length = 0.2
  }else{
    arr.length = arr.length
  }

  circle_size = unit(1, "snpc")
  fa = df$gene_id
  fa = factor(fa,levels = fa)
  circos.par(canvas.xlim = c(-1.1,1.1), canvas.ylim = c(-1.1,1.1),
             cell.padding = c(0, 0, 0, 0), gap.degree = gap.degree)
  circos.initialize(factors = fa, xlim = c(0,1))


  col_fun <- colorRamp2(expr_range, color_expr)
  circos.trackPlotRegion(
    ylim = c(0,1), track.height = track.height_1, bg.border = 'black',
    bg.col = col_fun(df$ligand_expr),
    panel.fun = function(x, y){
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  for (i in 1:nrow(df)) {
    circos.axis(sector.index = df[i,4],direction = "outside",
                labels = df[i,3], labels.facing = "clockwise", labels.cex = labels.cex,
                col = 'black', labels.away.percentage = 0.1,
                minor.ticks = 0, major.at = seq(0.5, length(df$ligand)/2))
  }


  circos.trackPlotRegion(
    ylim = c(0,1), track.height = track.height_2, bg.border = NA,
    panel.fun = function(x, y){
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  for(i in unique(df$ligand_cell)){
    num1 = min(which(df$ligand_cell == i))
    num2 = max(which(df$ligand_cell == i))
    highlight.sector(as.character(df$gene_id[num1:num2]),
                     track.index = 2, text = i,
                     niceFacing = T, font = 2,cex=1,
                     col = cell_col[which(unique(df$ligand_cell) == i)])
  }


  circos.trackPlotRegion(
    ylim = c(0,1),track.height = track.height_1, bg.border = NA,
    panel.fun = function(x, y){
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  for(i in unique(df$ligand_cell)){
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "ligand"))) > 0 ){
      num1 = min(intersect(which(df$ligand_cell == i),which(df$lr == "ligand")))
    }else{
      num1 = NULL
    }
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "ligand"))) > 0 ){
      num2 = max(intersect(which(df$ligand_cell == i),which(df$lr == "ligand")))
    }else{
      num2 = NULL
    }
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "receptor"))) > 0 ){
      num3 = min(intersect(which(df$ligand_cell == i),which(df$lr == "receptor")))
    }else{
      num3 = NULL
    }
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "receptor"))) > 0 ){
      num4 = max(intersect(which(df$ligand_cell == i),which(df$lr == "receptor")))
    }else{
      num4 = NULL
    }
    if ( !(is.null(num1) | is.null(is.integer(num2))) ) {
      highlight.sector(as.character(df$gene_id[num1:num2]),
                       track.index = 3, text = 'L', cex = 0.85,
                       text.vjust = text.vjust,
                       text.col = 'white', niceFacing = T,
                       col = cell_col[which(unique(df$ligand_cell) == i)])
    }
    if ( !(is.null(num3) | is.null(is.integer(num4))) ) {
    highlight.sector(as.character(df$gene_id[num3:num4]),
                     track.index = 3, text = 'R',
                     text.vjust = text.vjust,
                     text.col = 'white', niceFacing = T, cex = 0.85,
                     col = cell_col[which(unique(df$ligand_cell) == i)])
    }
  }


  for(i in 1:nrow(lrid)){
    circos.link(sector.index1 = lrid[i,4], point1 = 0.5,
                sector.index2 = lrid[i,9], point2 = 0.5,
                directional = 1, h = 0.85, lwd = link.lwd,
                col = lrid[i,'link_col'],
                lty = ifelse(length(link.lty) == 1, link.lty, link.lty[i]),
                arr.length = 0.2 , arr.col = col)
  }


  lgd_expr <- Legend(title = 'gene expression', at = expr_range,
                     col_fun = col_fun, title_position = 'topleft')
#  lgd_pval = Legend(title = "ligand_receptor pval", at = 1:2,
#                    labels = c(paste('lr_pval < ', pval), paste('lr_pval >= ', pval)),
#                    type = "lines", legend_gp = gpar(lty = 1:2))
  lgd_pval = Legend(title = "ligand_receptor pval", at = 1:2,
                    labels = c(paste('lr_pval < ', pval)),
                    type = "lines", legend_gp = gpar(lty = 1:1))
  lgd_cell = Legend(title = "cell type", at = 1:length(unique(df$ligand_cell)),
                    labels = as.character(unique(df$ligand_cell)),
                    legend_gp = gpar(fill = cell_col[1:length(unique(df$ligand_cell))]),
                    title_position = "topleft")
  if(is.null(lrid$comm_type)){
    lgd_list_vertical = packLegend(lgd_expr, lgd_cell, lgd_pval)
  }else{
    lgd_comm <- Legend(title = "communication type", at = 1:length(unique(lrid$comm_type)),
                       labels = as.character(unique(lrid$comm_type)), type = 'lines',
                       legend_gp = gpar( col = comm_col[1:length(unique(lrid$comm_type))]),
                       title_position = "topleft")
    lgd_list_vertical = packLegend(lgd_expr, lgd_cell, lgd_pval, lgd_comm)
  }
  pushViewport(viewport(x = 0.75, y = 0.5))
  draw(lgd_list_vertical, x = circle_size, just = "left")
  upViewport()

  circos.clear()
}



#' make a chorddiagram of the cellphonedb results 
#'
#' @param data a data.frame of the cellphonedb results
#' @param grid_orbit_h setting the height of grid tracks
#' @param label_h settiing the height of the label
#' @param diffHeight chord height
LRChorddiagram <- function(
  data,
  grid_orbit_h=0.02,
  label_h=0.04,
  diffHeight=2,
  output_dir="./",
  ...
){
  ## significant summary   
  data <- data %>% group_by(ligand_cell,receptor_cell) %>%
      filter( pval < opt$pvalue ) %>% 
      dplyr::summarize(n=n()) %>% 
      dplyr::rename(significant_pairs =  n )
#  write.table(data,file.path(output_dir, "cell_comm_chorddiagram_summary.xls"),quote=F,sep="\t",row.names=F)

  ## col set
  celltypes = sort( unique(c(data$receptor_cell, data$ligand_cell)) )
  colx = SelectColors(palette = opt$colorschema, length(celltypes))
  names(colx) <- celltypes  

  ## plot 
  chordDiagram(data, 
              grid.col=colx,
              directional = 1,
              diffHeight = -uh(diffHeight, "mm"), 
              direction.type = c("diffHeight", "arrows"),
              link.arr.type = "big.arrow",
              annotationTrack = c("name", "grid"), 
              annotationTrackHeight = c(label_h, grid_orbit_h) 
  )
}

#' Heatmap Viewing of cell-cell communication
#'
#' This function loads the significant interactions as a dataframe, and colors
#' represent the number of relations between different types of cells.
#'
#' @param data A dataframe containing ligand-receptor pairs and corresponding
#' cell types used to do the plotting
#' 
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @return A heatmap graph of the significant interactions
#' @export
LRHeatmap <- function(
  data,
  ...
){
  ## significant summary
  sum_data = data %>% filter(pval < 0.05) %>% group_by(ligand_cell, receptor_cell) %>%
        dplyr::summarize(n=n()) %>% tidyr::spread(receptor_cell,n)
#  write.table(sum_data, file.path(output_dir, "cell_comm_heatmap_summary.xls"), quote=F,sep="\t", row.names=F)

  plot_data = sum_data %>% tibble::column_to_rownames(var="ligand_cell")

  ## plot
  Heatmap(as.matrix(plot_data),
          col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(50),
          rect_gp = gpar(col = "white"),
          row_title = "Cell types for ligands",
          row_title_side = c("left"),
          row_title_gp = gpar(fontsize = 14),
          column_title = "Cell types for receptors",
          column_title_side = c("top"),
          column_title_gp = gpar(fontsize = 14),
          cluster_rows = F,
          cluster_columns = F,
          row_labels = rownames(plot_data),
          row_names_side = c("right"),
          show_row_names = TRUE,
          column_labels = colnames(plot_data),
          column_names_side = c("bottom"),
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 12),
          column_names_gp = gpar(fontsize = 12),
          name = "Number",
          row_names_max_width = unit(8, "cm"),
          column_names_max_height = unit(8, "cm"))
#  return(p)
}

#' Stacked histogram Viewing of cell-cell communication
#'
#' This function loads the significant interactions as a dataframe, and colors
#' represent the number of relations between different types of cells.
#'
#' @param data A dataframe containing ligand-receptor pairs and corresponding
#' cell types used to do the plotting
#' 
#' @return A histogram graph of the significant interactions
#' @export

LRBarplot <- function(
  data,
  bar_width = 0.6,
  ...
) {

  ## significant summary
  plot_data = data %>% filter(pval < 0.05) %>% group_by(ligand_cell, receptor_cell) %>%
        dplyr::summarize(n=n()) %>% group_by(ligand_cell) %>% arrange(ligand_cell) %>% mutate(cumsum = sum(n))
  celltype = unique(plot_data[order(as.matrix(plot_data)[,4],decreasing=T),]$ligand_cell)

  ## sort the ligand_cell
  plot_data$ligand_cell = factor(plot_data$ligand_cell,levels=celltype)

  ## plot
  p = ggplot(plot_data,
              aes(fill = receptor_cell, y = n, x = ligand_cell)) +
      geom_bar(position="stack", stat="identity", width = bar_width) +
      labs(x="Cell types of ligands", y = "Number of interactions")+
#      scale_fill_discrete(name="Cell types of receptors") +
      scale_fill_manual("Cell types of receptors", values = colx) +
      theme(panel.background = element_rect(fill = "transparent", colour = NA), 
            panel.grid = element_blank(), 
            axis.title = element_text(size = 15, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(size = 12, colour = "black"),
            axis.line = element_line(size=0.5, colour = "black"),
            plot.margin = unit(c(3,3,3,3), "cm"))   
  return(p)
}


suppressPackageStartupMessages( library("dplyr") )
suppressPackageStartupMessages( library("optparse") )
suppressPackageStartupMessages( library("ggplot2") )
suppressPackageStartupMessages( library("tidyr") )
#=command line parameters setting=============================
option_list = list(
    make_option( c("--input", "-i" ), type = "character",
        help="[REQUIRED]The cellphonedb results in RDS."),
    make_option( c("--informat", "-f" ), type = "character",default = "rds",
        help="[OPTIONAL]the format of the input expression matrix format. The
              option can be rds(seurat object), matrix"),
    make_option( c("--plot", "-d"), type="character", default="dotplot",
        help="[REQUIRED]the comma seperated plot type list to visualize the cell communications, 
              currently dotplot, network, chorddiagram, heatmap, barplot and circos are supported."),
    make_option( c("--pvalue", "-p"), type="double", default=0.05,
        help="[OPTIONAL]the P value threshold for significant interacting pairs."),
    make_option( c("--topn", "-n"), type="integer", default= 5,
        help="[OPTIONAL]the top N interacting gene pairs for each cell pairs. This is setted for dotplot and circos plot."),
    make_option( c("--topby", "-t"), type="character", default="expr",
        help="[OPTIONAL]the variable used to order the interacting pairs for circos plot."),
    make_option( c("--onlySigs", "-v"), type="logical", default="true",
        help="[OPTIONAL]only keep the significant gene pairs in all interacting cell pairs for dot plot."),
    make_option( c("--colorschema", "-c"), type="character", default="CustomCol2",
        help="[OPTIONAL]the supported color schema currently:seurat,blindless,paried,colx22,jet,tableau20,tableau10medium,
                                                   colorblind10,cyclic."),
    make_option( c("--output","-o"),type="character", default = "./",
        help="[OPTIONAL]the output directory of Clustering results." )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#=================================================================================
#parse the command line parameters
#=================================================================================
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    output_dir = getwd()
}else{
    if ( file.exists(opt$output) ){
        output_dir = opt$output
    }else{
        output_dir = opt$output
        dir.create(output_dir)
    }
}

if ( tolower(opt$onlySigs) == "true" ){
    is_onlySig = TRUE
}else{
    is_onlySig = FALSE
}

output_dir = normalizePath(output_dir)

#=================================================================================
# the main workflow
#=================================================================================
if ( is.null(opt$informat) || opt$informat == "rds"){
    cellp = readRDS(opt$input)
    data = cellp$ligrec
}else if ( opt$informat == "matrix"){
    data = read.table( opt$input, header = T, sep = "\t", stringsAsFactors = F )
    if(is.numeric(sort(unique(data$ligand_cell))) && is.numeric(sort(unique(data$receptor_cell)))){
        data$ligand_cell = factor(data$ligand_cell,levels = sort(unique(data$ligand_cell)))
        data$receptor_cell = factor(data$receptor_cell,levels = sort(unique(data$receptor_cell)))
	  }
}

# summary
net = data %>% filter(pval < opt$pvalue) %>% group_by(ligand_cell, receptor_cell) %>%
    dplyr::summarize(n=n()) %>% rename(significant_pairs =  n )
write.table(net,file.path(output_dir, "cell_comm_summary.xls"),quote=F,sep="\t",row.names=F)

for ( plotx in unlist( strsplit(opt$plot, ",", perl = T) ) ){

  if ( plotx == "dotplot" ){

    out = LRDotplot(data = data, is_onlySig = is_onlySig, topn = opt$topn, xsize = 15, xangle = 45)
    ggdot <- out[[1]]
    plot_data <- out[[2]]
    png_width = length(unique(plot_data$clusters))
    png_height = length(unique(plot_data$pair))
    ggsave(file.path(output_dir, "cell_comm_dotplot.pdf"), limitsize = F, 
           plot = ggdot, width = png_width*0.4+6, height = png_height*0.3)
    # ggsave(file.path(output_dir, "cell_comm_dotplot.png"), limitsize = F, 
    #       plot = ggdot, width = png_width*0.4+6, height = png_height*0.3)
    system(paste0("convert ",file.path(output_dir, "cell_comm_dotplot.pdf "),file.path(output_dir, "cell_comm_dotplot.png")))
  }
  
  if ( plotx == "network" ){
    suppressPackageStartupMessages( library("network") )
    suppressPackageStartupMessages( library("igraph") )
    celltypes = unique(c(data$receptor_cell, data$ligand_cell))
    colx = SelectColors(palette = opt$colorschema, length(celltypes))
    pdf(file.path(output_dir, "cell_comm_network.pdf"))
    LRNetwork( data %>% filter(pval < opt$pvalue),
              col = colx, edge.label.cex = 0.3, 
              edge.max.width = 5, arrow.width = 0.5, 
              vertex.label.cex = 0.5, vertex.size = 10)
    dev.off()
    png(file.path(output_dir, "cell_comm_network.png"), width = 7, height = 7, res = 500, units = "in" )
    LRNetwork( data %>% filter(pval < opt$pvalue),
              col = colx, edge.label.cex = 0.3, 
              edge.max.width = 5, arrow.width = 0.5, 
              vertex.label.cex = 0.5, vertex.size = 10)
    dev.off()
  }
  
  if ( plotx == "circos" ){
    suppressPackageStartupMessages( library("circlize") )
    suppressPackageStartupMessages( library("ComplexHeatmap") )
    celltypes = unique(c(data$receptor_cell, data$ligand_cell))
    colx = SelectColors(palette = opt$colorschema, length(celltypes))
    pdf( file.path(output_dir, "cell_comm_circos_plot.pdf"), width = 20, height = 10)
    LRCircos( data, gap.degree = 0.05, cell_col = colx, screenvar = opt$topby,
             topn = opt$topn, labels.cex = 0.35,link.lwd=1, arr.length = 0.1)
    dev.off()
    png( file.path(output_dir, "cell_comm_circos_plot.png"), width = 20, height = 10, res = 500, units = "in" )
    LRCircos( data, gap.degree = 0.05, cell_col = colx, screenvar = opt$topby,
             topn = opt$topn, labels.cex = 0.35,link.lwd=1, arr.length = 0.1)
    dev.off()
  }

  if ( plotx == "chorddiagram" ){
    suppressPackageStartupMessages( library("circlize") ) 
    pdf( file.path(output_dir, "cell_comm_chorddiagram_plot.pdf") )
    LRChorddiagram(data, grid_orbit_h=0.02,label_h=0.04,diffHeight=4,output_dir=output_dir)
    dev.off()
    png( file.path(output_dir, "cell_comm_chorddiagram_plot.png"), width = 7, height = 7, res = 96, units = "in" )
    LRChorddiagram(data, grid_orbit_h=0.02,label_h=0.04,diffHeight=4,output_dir=output_dir)
    dev.off()
    circos.clear()
  }

  if ( plotx == "heatmap" ){
    suppressPackageStartupMessages( library("ComplexHeatmap") )
    suppressPackageStartupMessages( library("RColorBrewer") )
    pdf( file.path(output_dir, "cell_comm_heatmap_plot.pdf") )
    print(LRHeatmap(data))
    dev.off()
    png( file.path(output_dir, "cell_comm_heatmap_plot.png"), width = 7, height = 7, res = 96, units = "in" )
    print(LRHeatmap(data))
    dev.off()
  }

  if( plotx == "barplot"){
    celltypes = unique(c(data$receptor_cell, data$ligand_cell))
    colx = SelectColors(palette = opt$colorschema, length(celltypes))
    names(colx) = celltypes
    out = LRBarplot(data = data, bar_width = 0.6)
    ggsave(file.path(output_dir, "cell_comm_histogram_plot.pdf"), plot = out, width = ifelse( length(celltypes)<7,7,length(celltypes)))
    ggsave(file.path(output_dir, "cell_comm_histogram_plot.png"), plot = out, width = ifelse( length(celltypes)<7,7,length(celltypes)), dpi=1000)
  }

}
if(!file.exists(file.path(output_dir, "cellphoneDB细胞通讯分析说明.docx"))){
  file.copy("/public/scRNA_works/Documents/cellphoneDB细胞通讯分析说明.docx",
  file.path(output_dir, "cellphoneDB细胞通讯分析说明.docx"))
}
