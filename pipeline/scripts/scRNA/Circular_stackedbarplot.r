# module purge 
# module load OESingleCell/2.0.0 

library(tidyverse)

inputFile="3-20.txt"             
data=read.table(inputFile,sep="\t",header=T,check.names=F)
data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
head(data)
empty_bar <- 2
nObsType <- nlevels(as.factor(data$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
head(data)
                                
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))

number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.2) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)),id=ceiling(mean(c(start, end))))

base_data_angle <- left_join(base_data,label_data,by='id') %>% subset(,c('group','start','end','title','id','angle'))
base_data_angle$angle <- base_data_angle$angle + 90


# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
p <- ggplot(data) +      
  
  # Add the stacked bar
# geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
#scale_fill_brewer(palette = "Paired") +
geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=1) +
scale_fill_brewer(palette = "RdBu") +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 8, xend = start, yend = 8), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 16, xend = start, yend = 16), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 24, xend = start, yend = 24), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 32, xend = start, yend = 32), colour = "black", alpha=1, size=1 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 8, 16, 24, 32), label = c("0", "8", "16", "24", "32") , color="black", size=2 , angle=0,fontface="bold", hjust=1) +
  
 ylim(-40,max(label_data$tot, na.rm=T)) +
 theme_minimal() +
 theme(
    #legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
#  geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", family="Times",fontface="bold", alpha=0.6, size=2, angle= label_data$angle, inherit.aes = FALSE ) +
  
  geom_text(data=label_data, aes(x=id, y=tot+5, label=individual, hjust=hjust), size=2, angle= label_data$angle, inherit.aes = FALSE ) +

 # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -8, label=group), hjust=0.5,vjust=0.5, angle=base_data_angle$angle,colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE) 



ggsave("circos.pdf",p, width=13, height=11)
ggsave("circos.png",p, dpi=1000, width=13, height=11)
