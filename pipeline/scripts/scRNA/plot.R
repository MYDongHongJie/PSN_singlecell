library(ggplot2)
library(hrbrthemes)
library(ggsci)

p = ggplot(data=data, aes(x=proportion, color=sampleid)) +
  geom_density(adjust=1.5, alpha=1,lwd=1.5,linetype = 1) +
  theme_ipsum()+
  guides(color = guide_legend(title = "sampleid"))+scale_colour_manual(values =c("#CC3366","#FF0099","#FF33CC","#FF66CC","#3333CC","#0000FF","#3366FF","#0099FF","#6699FF","#99CCFF",
    "#FFFF00","#FFFF66","#FFFF99","#FFFFCC"))+theme(panel.background = element_rect(fill ='white'))

g = ggplot(data=data, aes(x=proportion,fill=sampleid,color=sampleid)) +
  geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1)+scale_fill_manual(values=c("#CC3366","#FF0099","#FF33CC","#FF66CC","#3333CC","#0000FF","#3366FF","#0099FF","#6699FF","#99CCFF",
    "#FFFF00","#FFFF66","#FFFF99","#FFFFCC"))+
  theme_ipsum()+scale_colour_manual(values =c("#CC3366","#FF0099","#FF33CC","#FF66CC","#3333CC","#0000FF","#3366FF","#0099FF","#6699FF","#99CCFF",
    "#FFFF00","#FFFF66","#FFFF99","#FFFFCC"))+theme(panel.background = element_rect(fill ='white'))