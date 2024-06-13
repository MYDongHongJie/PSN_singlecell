library(psych)
library(qgraph)
library(igraph)
library(purrr)
netf<- "significant_pairs_network.txt"
mynet <- read.delim("cell_comm_summary.xls", check.names = FALSE)
head(mynet)
net<- graph_from_data_frame(mynet)

# 绘图代码如plot(x,y,...)
dev.off()
plot(net)
         CustomCol2 = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
                      "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
                      "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
                      "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
                      "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
                      "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
                      "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
                      "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
                      "#3fdacb","#bf5b17")

karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局

E(net)$width  <- E(net)$significant_pairs/10  # 边点权重（粗细）
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=CustomCol2 ,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
net2 <- net  # 复制一份备用

for (i in 1: length(unique(mynet$receptor_cell)) ){
  E(net)[map(unique(mynet$receptor_cell),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$receptor_cell)[i],x))
  })%>% unlist()]$color <- CustomCol2[i]
}  # 这波操作谁有更好的解决方案？ 

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=CustomCol2,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
plot(net, edge.arrow.size=.1, 
     edge.curved=0.2, # 只是调了这个参数
     vertex.color=CustomCol2,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
##
dev.off()
length(unique(mynet$receptor_cell)) # 查看需要绘制多少张图，以方便布局
par(mfrow=c(3,4), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$receptor_cell)) ){
  net1<-net2
  E(net1)[map(unique(mynet$receptor_cell),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$receptor_cell)[i],x))
  })%>% unlist()]$color <- CustomCol2[i]

  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       vertex.color=CustomCol2,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1) 

}
dev.off()
length(unique(mynet$receptor_cell))
par(mfrow=c(3,4), mar=c(.4,.4,.4,.4))

for (i in 1: length(unique(mynet$receptor_cell)) ){
  net1<-net2

  E(net1)$significant_pairs <- ""
  E(net1)[map(unique(mynet$receptor_cell),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$receptor_cell)[i],x))
  })%>% unlist()]$significant_pairs  <- E(net2)[map(unique(mynet$receptor_cell),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$receptor_cell)[i],x))
  })%>% unlist()]$significant_pairs  # 故技重施

  E(net1)[map(unique(mynet$receptor_cell),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$receptor_cell)[i],x))
  })%>% unlist()]$color <- CustomCol2[i]

  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       edge.label = E(net1)$significant_pairs, # 绘制边的权重
       vertex.color=CustomCol2,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1.4
  ) 

}
