install.packages("openxlsx")
library(openxlsx)
library(ggplot2)
library(dplyr)
setwd("D:\\Desktop\\temp\\20220329")
data<- read.xlsx("term.xlsx", sheet = 1)
data2<-read.xlsx("term.xlsx", sheet = 2)
data3=data %>% select(Term,term,Interested_group,pval,t) %>% filter(pval<0.05) %>% mutate(color="neg")
data3[which(data3$Interested_group=="Myofibroblast"),"Interested_group"]="PLN+myoFB"
data3[which(data3$Interested_group=="MyoCAF"),"Interested_group"]="PI15+myoFB"
data3[which(data3$Interested_group=="CAF_1"),"Interested_group"]="MMP1+FB"
data3[which(data3$Interested_group=="Fibroblast"),"Interested_group"]="PTGS2+FB"
data3[which(data3$Interested_group=="CAF_2"),"Interested_group"]="COMP+FB"
data3[which(data3$t>0),"color"]="pos"
myPalette <- colorRampPalette( RColorBrewer::brewer.pal(9,'GnBu'))

data3$Interested_group=factor(data3$Interested_group,levels = c("PI15+myoFB","PLN+myoFB","MMP1+FB","PTGS2+FB","COMP+FB"))
data3$term=factor(data3$term,levels = 
rev(c("Lipoprotein localization"," MSC proliferation",                             
"Toll signaling pathway","Cellular response to bacterial lipoprotein",     
"Pos-reg of ¦Ã¦ÄT cell activation","Antigen processing and presentation via MHCIb",    
"Epithelial structure maintenance","Epithelial cell-cell adhesion",                                       
"Pos-reg of epi-cell proliferation in wound healing","Pos-reg of epidermis development",                 
"Regulation of keratinocyte migration","Pos-reg of keratinocyte proliferation",            
"Positive regulation of EMT","Response to TGF¦Â",                               
"Extracellular matrix disassembly","Collagen metabolic process" ) ))

ggplot(data3,aes(term,Interested_group))+
  geom_point(aes(size=abs(t),color=factor(color)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=10),
        axis.title.y=element_text(size=10))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(title = "Tvalue"),
         color=guide_legend(title = "Tvalue"))+
  #scale_colour_gradient(low="red",high="blue")+
  #scale_colour_gradientn(colours = myPalette(256) )+
  #scale_fill_manual(values = c("hight"="#0000FF","low"="#FF0000"))+
  scale_color_manual(values =c("pos"='red',"neg"='black'))
#guides(fill="none")
#guides(color=guide_legend())

ggsave("X_Go_dotplot.pdf",width = 5,height=4.5)
ggsave("X_Go_dotplot.png",dpi = 1000 ,limitsize = F,width = 5,height=4.5)

data3$term=factor(data3$term,levels = 
                    (c("Lipoprotein localization"," MSC proliferation",                             
                          "Toll signaling pathway","Cellular response to bacterial lipoprotein",     
                          "Pos-reg of ¦Ã¦ÄT cell activation","Antigen processing and presentation via MHCIb",    
                          "Epithelial structure maintenance","Epithelial cell-cell adhesion",                                       
                          "Pos-reg of epi-cell proliferation in wound healing","Pos-reg of epidermis development",                 
                          "Regulation of keratinocyte migration","Pos-reg of keratinocyte proliferation",            
                          "Positive regulation of EMT","Response to TGF¦Â",                               
                          "Extracellular matrix disassembly","Collagen metabolic process" ) ))
ggplot(data3,aes(Interested_group,term))+
  geom_point(aes(size=abs(t),color=factor(color)))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=10),
        axis.title.y=element_text(size=10))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(title = "Tvalue"),
         color=guide_legend(title = "Tvalue"))+
  #scale_colour_gradient(low="red",high="blue")+
  #scale_colour_gradientn(colours = myPalette(256) )+
  #scale_fill_manual(values = c("hight"="#0000FF","low"="#FF0000"))+
  scale_color_manual(values =c("pos"='red',"neg"='black'))

ggsave("Y_Go_dotplot.pdf",width = 5,height=4.5)
ggsave("Y_Go_dotplot.png",dpi = 1000 ,limitsize = F,width = 5,height=4.5)
