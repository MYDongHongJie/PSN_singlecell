df <- mtcars

## 2.相关系数计算
res <- cor(df)

## 3.安装并加载corrplot包
# install.packages("corrplot")
library(corrplot)

#简单相关性系数可视化
corrplot(res)  

# # 不同method的相关性系数图表
# corrplot(res,method="pie",tl.col="black",tl.srt=45, title = "method=pie的饼图", cex.main = 1, mar = c(2,2,3,2))  #饼图
# corrplot(res,method="ellipse",shade.col=NA,tl.col="black",tl.srt=45,  title = "method=ellipse的饼图", cex.main = 1, mar = c(2,2,3,2)) #椭圆
# corrplot(res, method="number",shade.col=NA,tl.col="black",tl.srt=45,  title = "method=number的饼图", cex.main = 1, mar = c(2,2,3,2))#数字




#相关系数可视化混合矩阵图
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45, 
         title = "type = upper的数字+饼图", mar = c(2,2,3,2))  #上三角


corrplot.mixed(res, title = "图形和数值混合矩阵", mar = c(2,2,3,2)) #图形和数值混合矩阵

corrplot.mixed(res, lower.col = "black", number.cex = .7, 
               title = "文字看不清，可以设置文字为黑色lower.col和大小number.cex", mar = c(2,2,3,2))


corrplot(res, order = "hclust", addrect = 2, 
         title = "按hclust聚类方式排序", mar = c(2,2,3,2))  
#按hclust聚类方式排序，addrect是添加分组矩形，可自定义分组类
#类似于平时热图的kmean分组方式。