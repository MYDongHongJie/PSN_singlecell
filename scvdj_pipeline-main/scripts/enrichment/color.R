#!/usr/bin/env Rscript

fix_color <- function(names) {
    k <- c(
        "Metabolism", "Human Diseases",
        "Organismal Systems",
        "Environmental Information Processing",
        "Genetic Information Processing",
        "Cellular Processes",
        "Brite Hierarchies",
        "Not Included in Pathway or Brite"
    )

    g <- c(
        "biological_process",
        "molecular_function",
        "cellular_component"
    )

    base_color <- c(
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
        "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
        "#CCEBC5", "#FFED6F", "#71A99F", "#CCCC8F", "#9895AE",
        "#C9665B", "#668EA9", "#CA904E", "#8FB254", "#CAA4B7"
    )

    k_color <-  base_color[seq_len(length(k))]
    g_color <-  base_color[seq_len(length(g))]
    fixed_color <- c(k_color, g_color)
    names(fixed_color) <- c(k, g)

    return(fixed_color[names])
}
