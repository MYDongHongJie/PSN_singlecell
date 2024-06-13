## »æÖÆÄâÊ±ĞòÉ½ÂÍÍ¼
library(ggridges)
Library(ggplot2)
Library(monocle)

new_cds <- buildBranchCellDataSet(gbm_cds, branch_point = branchpoint, progenitor_method = "duplicate")
            cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
            cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
            branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))
            new_cds <- buildBranchCellDataSet(gbm_cds, branch_point = branchpoint, branch_labels = branch_labels, progenitor_method = "duplicate")

p <- ggplot(aes(x=Pseudotime, y = B_new_celltype , fill = B_new_celltype), data = pData(new_cds)) +
                    geom_density_ridges() +
                    ylab("B_new_celltype") + facet_wrap(~Branch) +
                    xlab("Pseudotime") + xlim(100, 0) +
                    theme_bw() +  scale_fill_manual(values = CustomCol2(1:50)) +
                    theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(size = 1, colour = "black"))
                ggsave(file.path(output_dir, paste0("ridges_B_new_celltype_branchtime2.pdf")), plot = p, width = 15, height = 3)
                ggsave(file.path(output_dir, paste0("ridges_B_new_celltype_branchtime2.png")), plot = p, width = 15, height = 3)