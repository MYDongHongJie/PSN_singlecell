library(ggplot2)
library(tidyverse)
library(ggrepel)

df <- read.delim('data.xls')
# 添加显著性标签：
df$label <- ifelse(df$p_val_adj < 0.05, "adjust P-val<0.05", "adjust P-val>=0.05")

# 获取每个cluster中表达差异最显著的10个基因；
top10sig1 <- filter(df, cluster == "1") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))
top10sig2 <- filter(df, cluster == "2") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))
top10sig3 <- filter(df, cluster == "3") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))
top10sig4 <- filter(df, cluster == "4") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))
top10sig5 <- filter(df, cluster == "5") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))
top10sig6 <- filter(df, cluster == "6") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))
top10sig7 <- filter(df, cluster == "7") %>%
    distinct(geneID, .keep_all = T) %>%
    top_n(10, abs(log2FC))


# 将提取所有cluster的Top10基因表格合并：
top10sig <- rbind(top10sig1, top10sig2, top10sig3, top10sig4, top10sig5, top10sig6, top10sig7)

# 新增一列，将Top10的差异基因标记为2，其他的标记为1；
df$size <- case_when(
    !(df$geneID %in% top10sig$geneID) ~ 1,
    df$geneID %in% top10sig$geneID ~ 2
)

# 提取非Top10的基因表格；
dt <- filter(df, size == 1)
# 绘制每个Cluster Top10以外基因的散点火山图：
p <- ggplot() +
    geom_jitter(
        data = dt,
        aes(x = cluster, y = log2FC, color = label),
        size = 0.85,
        width = 0.4
    )
# 根据图p中log2FC区间确定背景柱长度：
dfbar <- data.frame(
    x = c(1, 2, 3, 4, 5, 6, 7),   ###
    y = c(2.8, 2.2, 3.6, 2, 2, 1.9, 2.3)
)
dfbar1 <- data.frame(
    x = c(1, 2, 3, 4, 5, 6, 7),
    y = c(-3.5, -3.2, -3.2, -2.2, -3.9, -2.5, -2.6)
)
# # 绘制背景柱：
# p1 <- ggplot() +
#     geom_col(
#         data = dfbar,
#         mapping = aes(x = x, y = y),
#         fill = "#dcdcdc", alpha = 0.6
#     ) +
#     geom_col(
#         data = dfbar1,
#         mapping = aes(x = x, y = y),
#         fill = "#dcdcdc", alpha = 0.6
#     )

# 把散点火山图叠加到背景柱上：
p2 <- ggplot() +
    geom_col(
        data = dfbar,
        mapping = aes(x = x, y = y),
        fill = "#dcdcdc", alpha = 0.6
    ) +
    geom_col(
        data = dfbar1,
        mapping = aes(x = x, y = y),
        fill = "#dcdcdc", alpha = 0.6
    ) +
    geom_jitter(
        data = dt,
        aes(x = cluster, y = log2FC, color = label),
        size = 0.85,
         position=position_jitter(seed=1)
    ) +
    geom_jitter(
        data = top10sig,
        aes(x = cluster, y = log2FC, color = label),
        size = 1,
         position=position_jitter(seed=1)
    )
# 添加X轴的cluster色块标签：
dfcol <- data.frame(
    x = c(1:7),
    y = 0,
    label = c(1:7)
)
mycol <- c("#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F")
p3 <- p2 + geom_tile(
    data = dfcol,
    aes(x = x, y = y),
    height = 0.6,
    color = "black",
    fill = mycol,
    alpha = 0.6,
    show.legend = F
)


# 给每个Cluster差异表达前Top10基因加上标签：
p4 <- p3 +
    geom_text_repel(
        data = top10sig,
        aes(x = cluster, y = log2FC, label = geneID),
        force = 1.2, position=position_jitter(seed=1) ,
  max.overlaps=50,
        arrow = arrow(
            length = unit(0.008, "npc"),
            type = "open", ends = "last"
        )
    )
# 散点颜色调整：
p5 <- p4 +
    scale_color_manual(
        name = NULL,
        values = c("red", "black")
    )
# 修改X/Y轴标题和添加cluster数字：
p6 <- p5 +
    labs(x = "Cluster", y = "average logFC") +
    geom_text(
        data = dfcol,
        aes(x = x, y = y, label = label),
        size = 6,
        color = "white"
    )
# 自定义主题美化：
p7 <- p6 +
    theme_minimal() +
    theme(
        axis.title = element_text(
            size = 13,
            color = "black",
            face = "bold"
        ),
        axis.line.y = element_line(
            color = "black",
            size = 1.2
        ),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1, 0),
        legend.text = element_text(size = 15)
    )
ggsave('plot.pdf')