

# 重新自己安装，所有工具和包都安装最新版，能不用conda就不用，决不能用mamba

mp 
source /public/dev_scRNA/software/miniconda3/etc/profile.d/conda.sh

conda create --prefix /public/dev_scRNA/software/miniconda3/envs/OESingleCell_alternativesplicing
conda activate /public/dev_scRNA/software/miniconda3/envs/OESingleCell_alternativesplicing
conda install -c conda-forge r-base
conda install -c conda-forge cmake
conda install -c conda-forge make
conda install -c bioconda bioconductor-rhdf5filters


# 需要用BiocManager安装的工具
install.packages("BiocManager")
# 由于目前最新的R版本为4.3及以上，因此在使用BiocManager安装包时需要使用3.18版本
install.packages(version = '3.18')  # 是否更新包？否
BiocManager::install("AnnotationDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("GenomicRanges")
install.packages("ggpubr")
install.packages("gtools")
install.packages("hash")
# install.packages("devtools")  # 安装不成功，改用conda 
# devtools::install_github("wenweixiong/MARVEL")
install.packages("MARVEL")
install.packages("optparse")

BiocManager::install("rhdf5filters") # 安装不成功，改用conda  
BiocManager::install("BSgenome")   
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("phastCons100way.UCSC.hg38")

install.packages("scCustomize")
install.packages("factoextra")
install.packages("kSamples")
install.packages("pheatmap")
install.packages("textclean")
install.packages("twosamples")

install.packages("hdf5r")  # 安装不成功，改用conda  
remotes::install_github("mojaveazure/seurat-disk")

BiocManager::install("wiggleplotr")