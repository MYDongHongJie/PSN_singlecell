#!/usr/bin/env Rscript
# by The Coder, 20160726
# Packing script(add optparse), 20180111
suppressMessages(library(clusterProfiler))
suppressMessages(library(docopt))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessag es(library(stringr))
suppressMessages(library(grid))
suppressMessages(library(oebio))
suppressMessages(library(Cairo))

"
Usage: enrich
    enrich ora -i <genelist> -g <gmt> 
"