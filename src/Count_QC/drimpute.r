
## Impute dropouts with the DrImpute package
drimpute_dropouts <- function(count) {
    suppressPackageStartupMessages(library(DrImpute))
    ####################################################################################################################
    ###Define ks to use ，聚类分组，默认设置为ks=10:15
    ####################################################################################################################
    n <- ncol(count) / 2
    ncl <- 6
    ks <- unique(round(2 : n))
    if (length(ks) <= ncl) {
        ks <- ks
    } else if (length(ks) <= (3 * ncl)) {
        ks <- ks[ks <= (4 * ncl)]
        ks <- ks[((length(ks) - ncl) / 2 + 1) : (length(ks) - (n - ncl) / 2)]
    } else {
        ks <- ks[ks <= (4 * ncl)]
        ks <- ks[seq(1, length(ks), 2)]
        ks <- ks[((length(ks) - ncl) / 2 + 1) : (length(ks) - (length(ks) - ncl) / 2)]
    }
    print(paste0("KS：",ks)
    count_orig <- as.matrix(count)
    ####################################################################################################################
    ## Normalize counts                                                                                             ####
    ####################################################################################################################
    count <- as.matrix(count)
    sf <- colSums(count)
    sf[sf == 0] <- 1
    sf <- sf / exp(mean(log(sf)))
    lognormcount <- log(sweep(count, 2, sf, "/") + 1)
    lognormcount_imp <- DrImpute(lognormcount, ks = ks, dists = c("spearman", "pearson"), method = "mean", cls = NULL,
    seed = 123, zerop = 0)
    ####################################################################################################################
    ## Go back to original count scale                                                                               ###
    ####################################################################################################################
    normcount_imp <- exp(lognormcount_imp) - 1
    count_imp <- sweep(normcount_imp, 2, sf, "*")
    colnames(count_imp) <- colnames(count_orig)
}
