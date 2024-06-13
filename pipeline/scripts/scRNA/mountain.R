library(ggridges)
library(tidyverse)
library(gridExtra)
library(viridis)
library(ggplot2)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("OESingleCell"))
suppressPackageStartupMessages(library("future.apply"))
suppressPackageStartupMessages(library("magrittr"))

data_ob = OESingleCell::ReadX(input = 'new_celltype.h5seurat', informat = 'h5seurat', verbose = F)

a=data_ob@assays$RNA@data["G0s2",]
a=as.data.frame(a)
data_ob@meta.data$G0s2=""
for (i in unique(a$a) ){
    subset_ob= subset(a, a==i)
    data_ob@meta.data[as.character(rownames(subset_ob)) ,"G0s2"] = i
}
plotdata=data_ob@meta.data %>%
                               dplyr::rename( "Barcode" = "rawbc") %>%
                               dplyr::select( Barcode, new_celltype, G0s2)

plotdata2=plotdata
plotdata2=subset(plotdata2,G0s2!=0 )
plotdata2=rbind(plotdata2,c(1,"B cells",0))
plotdata2=rbind(plotdata2,c(1,"LSEC",0))
plotdata2$G0s2=as.numeric(plotdata2$G0s2)


gg = ggplot(plotdata2, aes(x=G0s2, y=new_celltype, fill=new_celltype))+
    geom_density_ridges()+
    scale_x_continuous(expand = c(0.01, 0))+
    scale_y_discrete(expand = c(0.01,0))+
    # scale_fill_viridis(name="Temp. [F]", option = "C")+
    # labs(title="Temperature in Lincoln NE",subtitle="Mean temperature (Fahrenheit) by month for 2016\nData:Orogin CSV from the Weather Underground ")+
    theme_ridges(font_size = 13, grid = FALSE)+
    theme(axis.title.y = element_blank())
ggsave(gg,file="mountain2.pdf")
ggsave(gg,file="mountain2.png",dpi=1000)