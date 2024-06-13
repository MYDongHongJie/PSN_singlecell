library(rayshader)  # 在本地电脑上安装此包，linux终端可能出不来图
library(ggplot2)
library(Seurat)
data_ob=readRDS('sampleid.rds')
p=DimPlot(data_ob,group.by="clusters",reduction="SCT_umap")
p$data$clusters = factor(p$data$clusters,levels=c(1:length(unique(p$data$clusters))))
gg_tsne = p + stat_density_2d(geom="raster",contour = F,data=p$data,aes(x=UMAP_1,y=UMAP_2,fill=..density..),alpha=1)+ scale_fill_viridis_c()+theme_test() +theme(legend.position="none")

plot_gg(gg_tsne,multicore = T)

render_camera(zoom=0.6, phi=50,theta = 60,fov=80)  # zoom 图形放大缩小倍数；phi方位角；theta旋转角度。该步骤类似于调整好角度用照相机拍快照。
render_snapshot()  #保存或输出当前视图




