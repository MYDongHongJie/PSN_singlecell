library(Seurat)
library(ggplot2)
library(ggsci)
library(grid)
library(ggrepel)
library(tidyverse)
library(SeuratData)
library(tidydr)
library(patchwork)
library(pbmc3k.SeuratData)

pc12 <- Embeddings(object = seurat_ob,reduction = 'tsne') %>% data.frame()
lower <- floor(min(min(pc12$tSNE_1),min(pc12$tSNE_2))) - 2
linelen <- abs(0.3*lower) + lower
mid <- abs(0.3*lower)/2 + lower
label <- data.frame(lab = c('t-SNE2','t-SNE1'),angle = c(90,0),x = c(lower - 3,mid),y = c(mid,lower - 2.5))
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),group = c(1,1,2,2),label = rep(c('t-SNE2','t-SNE1'),each = 2))

p = DimPlot(seurat_ob,reduction = 'tsne',label = T) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +
  scale_colour_lancet() +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')