#!/usr/bin/env python
# coding: utf-8
import click
import numpy as np
import os
import pandas as pd
import scanpy as sc
import stlearn as st
import sys

sys.path.append(
    '/public/dev_scRNA/xfyang/spatial_transcriptome/software_testing/stLearn/stLearn-master/stlearn/image_preprocessing/')
sys.path.append('/public/dev_scRNA/xfyang/spatial_transcriptome/software_testing/stLearn/stLearn-master/stlearn/')
from oebio.utils.log import getLogger
from model_zoo import encode, Model
from typing import Optional, Union
from anndata import AnnData
from pathlib import Path
from _compat import Literal
from PIL import Image
from pathlib import Path
from tqdm import tqdm
from sklearn.decomposition import PCA

logger = getLogger('oesc.stlearn_SMEclust')


# from stlearn import *
# from stlearn.image_preprocessing import *
###deal adata
def pre_process(anndata, mincell):
    ##过滤基因，至少1个细胞中表达
    st.pp.filter_genes(anndata, min_cells=1)
    ##标准化
    st.pp.normalize_total(anndata)
    ##归一化
    st.pp.log1p(anndata)
    return anndata


_CNN_BASE = Literal["resnet50", "vgg16", "inception_v3", "xception"]


def extract_feature_update(
    anndata: AnnData,
    cnn_base: _CNN_BASE = "resnet50",
    n_components: int = 50,
    verbose: bool = False,
    copy: bool = False,
    seeds: int = 1,
) -> Optional[AnnData]:
    """    Extract latent morphological features from H&E images using pre-trained
    convolutional neural network base
    Parameters
    ----------
    anndata
        Annotated data matrix.
    cnn_base
        Established convolutional neural network bases
        choose one from ['resnet50', 'vgg16', 'inception_v3', 'xception']
    n_components
        Number of principal components to compute for latent morphological features
    verbose
        Verbose output
    copy
        Return a copy instead of writing to anndata.
    seeds
        Fix random state
    Returns
    -------
    Depending on `copy`, returns or updates `anndata` with the following fields.
    **X_morphology** : `anndata.obsm` field
        Dimension reduced latent morphological features.
    """
    feature_df = pd.DataFrame()
    model = Model_update(cnn_base)
    if "tile_path" not in anndata.obs:
        raise ValueError("Please run the function stlearn.pp.tiling")
    with tqdm(
        total=len(anndata),
        desc="Extract feature",
        bar_format="{l_bar}{bar} [ time left: {remaining} ]",
    ) as pbar:
        for spot, tile_path in anndata.obs["tile_path"].items():
            tile = Image.open(tile_path)
            tile = np.asarray(tile, dtype="int32")
            tile = tile.astype(np.float32)
            tile = np.stack([tile])
            if verbose:
                print("extract feature for spot: {}".format(str(spot)))
            features = encode(tile, model)
            feature_df[spot] = features
            pbar.update(1)

    anndata.obsm["X_tile_feature"] = feature_df.transpose().to_numpy()

    pca = PCA(n_components=n_components, random_state=seeds)
    pca.fit(feature_df.transpose().to_numpy())
    anndata.obsm["X_morphology"] = pca.transform(feature_df.transpose().to_numpy())
    print("The morphology feature is added to anndata.obsm['X_morphology']!")
    return anndata if copy else None


###函数修改，用于支持本地预训练模型导入
script_dir = os.path.split(os.path.realpath(__file__))[0]
resne50_path = f"/public/dev_scRNA/xfyang/devlopment_project/python_script/cnn_model/resnet50_weights_tf_dim_ordering_tf_kernels_notop.h5"
vgg16_path = f"/public/dev_scRNA/xfyang/devlopment_project/python_script/cnn_model/vgg16_weights_tf_dim_ordering_tf_kernels_notop.h5"
inception_v3 = f"/public/dev_scRNA/xfyang/devlopment_project/python_script/cnn_model/inception_v3_weights_tf_dim_ordering_tf_kernels_notop.h5"
xception = f"/public/dev_scRNA/xfyang/devlopment_project/python_script/cnn_model/xception_weights_tf_dim_ordering_tf_kernels_notop.h5"


def encode(tiles, model):
    features = model.predict(tiles)
    features = features.ravel()
    return features


class Model_update:
    __name__ = "CNN base model"

    def __init__(self, base, batch_size=1):
        from tensorflow.keras import backend as K

        self.base = base
        self.model, self.preprocess = self.load_model()
        self.batch_size = batch_size
        self.data_format = K.image_data_format()

    def load_model(self):
        if self.base == "resnet50":
            from tensorflow.keras.applications.resnet50 import (
                ResNet50,
                preprocess_input,
            )

            cnn_base_model = ResNet50(
                include_top=False, weights=resne50_path, pooling="avg"
            )
        elif self.base == "vgg16":
            from tensorflow.keras.applications.vgg16 import VGG16, preprocess_input

            cnn_base_model = VGG16(include_top=False, weights=vgg16_path, pooling="avg")
        elif self.base == "inception_v3":
            from tensorflow.keras.applications.inception_v3 import (
                InceptionV3,
                preprocess_input,
            )

            cnn_base_model = InceptionV3(
                include_top=False, weights=inception_v3, pooling="avg"
            )
        elif self.base == "xception":
            from tensorflow.keras.applications.xception import (
                Xception,
                preprocess_input,
            )

            cnn_base_model = Xception(
                include_top=False, weights=xception, pooling="avg"
            )
        else:
            raise ValueError("{} is not a valid model".format(self.base))
        return cnn_base_model, preprocess_input

    def predict(self, x):
        from tensorflow.keras import backend as K

        if self.data_format == "channels_first":
            x = x.transpose(0, 3, 1, 2)
        x = self.preprocess(x.astype(K.floatx()))
        return self.model.predict(x, batch_size=self.batch_size)


##输入TILE_PATH
def tiling_and_extract(anndata, TILE_PATH, cnn_base):  ##cnn_base:["resnet50", "vgg16", "inception_v3", "xception"]
    st.pp.tiling(anndata, TILE_PATH)
    extract_feature_update(anndata, cnn_base)
    return anndata


# apply stSME to normalise log transformed data
def SME_cluster_louvain(anndata, n_clusters, n_neighbors, resolution):
    ##1.pca
    st.em.run_pca(anndata, n_comps=50)
    ##2.SME_normalize
    st.spatial.SME.SME_normalize(anndata, use_data="raw")
    anndata.X = anndata.obsm['raw_SME_normalized']
    st.pp.scale(anndata)
    ##3.pca ,kmeans,louvain
    st.em.run_pca(anndata, n_comps=50)
    st.tl.clustering.kmeans(anndata,
                            n_clusters=n_clusters,
                            use_data="X_pca",
                            key_added="X_pca_kmeans")
    st.pp.neighbors(anndata,
                    n_neighbors=n_neighbors,
                    use_rep='X_pca',
                    random_state=0)
    st.tl.clustering.louvain(anndata,
                             random_state=0,
                             resolution=resolution)
    # 修改clusters
    # anndata.obs["X_pca_kmeans"] = anndata.obs["X_pca_kmeans"].astype(int) + 1
    # anndata.obs["X_pca_kmeans"] = anndata.obs["X_pca_kmeans"].astype("str").astype("category")
    # orderKmeans = [str(x) for x in range(1, len(anndata.obs.X_pca_kmeans.unique()) + 1)]
    # anndata.obs.X_pca_kmeans.cat.reorder_categories(orderKmeans, inplace=True)
    # anndata.obs["louvain"] = anndata.obs["louvain"].astype(int) + 1
    # anndata.obs["louvain"] = anndata.obs["louvain"].astype("str").astype("category")
    # orderLouvain = [str(x) for x in range(1, len(anndata.obs.louvain.unique()) + 1)]
    # anndata.obs.louvain.cat.reorder_categories(orderLouvain, inplace=True)
    return anndata


###st.spatial.morphology.adjust与st.spatial.SME.SME_normalize(data_SME, use_data="raw")二选一
####spatial.morphology.adjust需要已有pca信息用于后续adjust
def adjust_and_louvain(anndata, n_neighbors, resolution):
    st.spatial.morphology.adjust(anndata,
                                 use_data="X_pca",
                                 radius=50,
                                 method="mean")  ####使用了其他标准化方法，因此使用st.spatial.morphology.adjust调整PCA
    st.pp.neighbors(anndata,
                    n_neighbors=n_neighbors,
                    use_rep='X_pca_morphology',
                    random_state=0)
    st.tl.clustering.louvain(anndata,
                             random_state=0,
                             resolution=resolution)
    return anndata


####h5ad subset
def h5ad_subset(anndata, colname, sample):
    anndata_sub = anndata[anndata.obs[colname] == sample]
    anndata_sub.uns['spatial'] = {sample: anndata_sub.uns['spatial'][sample]}
    return anndata_sub


def h5ad_adjust_cluster(anndata, TILE_PATH, OUT_PATH, cnn_base, samplename, n_neighbors, resolution):
    TILE_PATH_samplename = Path(f"{TILE_PATH}/{samplename}")
    TILE_PATH_samplename.mkdir(parents=True, exist_ok=True)
    tiling_and_extract(anndata, TILE_PATH_samplename, cnn_base)
    adjust_and_louvain(anndata, n_neighbors=n_neighbors, resolution=resolution)
    anndata.write(f'{OUT_PATH}/{samplename}_SME.h5ad')
    st.pl.cluster_plot(anndata, use_label='louvain',
                       show_cluster_labels=True,
                       crop=False,
                       fname=f'{OUT_PATH}/{samplename}_louvain_resolution_{resolution}.pdf', dpi=600, size=10)


def file_SME_cluster(anndata, TILE_PATH, OUT_PATH, cnn_base, samplename, n_clusters, n_neighbors, resolution):
    pre_process(anndata, 1)
    tiling_and_extract(anndata, TILE_PATH, cnn_base="resnet50")
    SME_cluster_louvain(anndata, n_clusters, n_neighbors, resolution)
    anndata.write(f'{OUT_PATH}/adata_SME.h5ad')
    st.pl.cluster_plot(anndata, use_label='X_pca_kmeans', show_cluster_labels=True, crop=False,
                       fname=f'{OUT_PATH}/adata_kmeans_{n_clusters}.pdf', dpi=600, size=10)
    st.pl.cluster_plot(anndata, use_label='louvain', show_cluster_labels=True, crop=False,
                       fname=f'{OUT_PATH}/adata_louvain_resolution_{resolution}.pdf', dpi=600, size=10)


@click.command()
@click.option('--analy_path', default="./", help="the file containing .h5ad file or a file containing .h5")
@click.option('--analy_file', help="h5ad file or a file containing h5 file and the spatial file")
@click.option('--n_clusters', default=15,
              help="the number of clusters to form as well as the number of centroids to generate,for kmeans")
@click.option('--n_neighbors', default=25,
              help=" The size of local neighborhood (in terms of number of neighboring data points) used for manifold "
                   "approximation. Larger values result in more global views of the manifold, while smaller values result"
                   "in more local data being preserved. In general values should be in the range 2 to 100")
@click.option('--resolution', default=1.0,
              help="resolution for louvain,higher resolution means finding more and smaller clusters")
@click.option('--out_path', help="output pathlib")
@click.option('--barcodes', default=None, help="to subset the adata")
@click.option('--size', default=16, help="point size of spots")
@click.option('--crop', default=True, help="whether to crop the image (CytAssist project should be True)")
@click.option('--cluster_labels', default=False, help="whether to show cluster labels")
@click.option('--image_alpha', default=0.5, help="the color transparency of image")
#####analy_path=存放filtered_feature_bc_matrix.h5的文件夹的上一级文件夹或.h5ad上一级, 如"/public/dev_scRNA/yfang/Stlearn_test/"
####analy_file=存放filtered_feature_bc_matrix.h5与spatial的文件夹或.h5ad, 如"Mouse_Brain_SpS1_normalized.h5ad" 或者/public/dev_scRNA/yfang/Stlearn_test/10x_visium/filter_h5_data/
#####1) 多个合并的样本必须有sample(样本标签)or library_id（样本标签）
#####2) h5ad文件必须要有.uns['spatial']信息
def st_SMEclust(analy_path, analy_file, n_clusters, n_neighbors, resolution, out_path, barcodes, size, cluster_labels,
                image_alpha, crop):
    # step1: get Base PATH
    BASE_PATH = Path(f"{analy_path}/{analy_file}")
    # step2: spot tile is the intermediate result of image pre-processing
    TILE_PATH = Path(f"{out_path}/tmp/tiles/")
    TILE_PATH.mkdir(parents=True, exist_ok=True)
    # output path
    OUT_PATH = Path(f"{out_path}/1_SMEnormalize/")
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    ###输入格式为h5ad时
    if "h5ad" in analy_file:
        print("read from h5ad")
        adata = sc.read_h5ad(BASE_PATH)
        if barcodes != None:
            bar = pd.read_csv(barcodes, sep='\t')
            intbarcodes = adata.obs_names & bar.iloc[:, 0]
            adata = adata[intbarcodes].copy()
        if "imagecol" in adata.obs.columns:
            sample_id = list(adata.uns['spatial'].keys())
            if "library_id" in adata.obs.columns:
                for i in sample_id:
                    print(f"from imagecol data wiht library_id subset {i}")
                    adata_copy = adata.copy()
                    adata_sub = h5ad_subset(anndata=adata_copy, colname='library_id', sample=i)
                    print(f"till {i}")
                    h5ad_adjust_cluster(anndata=adata_sub,
                                        TILE_PATH=TILE_PATH,
                                        OUT_PATH=OUT_PATH,
                                        cnn_base="resnet50",
                                        n_neighbors=n_neighbors,
                                        resolution=resolution,
                                        samplename=i)
                print("data with imagecol ")
            elif "sample" in adata.obs.columns:
                for i in sample_id:
                    print(f"from imagecol data wiht sample subset {i}")
                    adata_copy = adata.copy()
                    adata_sub = h5ad_subset(anndata=adata_copy, colname='sample', sample=i)
                    print(f"till {i}")
                    h5ad_adjust_cluster(anndata=adata_sub,
                                        TILE_PATH=TILE_PATH,
                                        OUT_PATH=OUT_PATH,
                                        cnn_base="resnet50",
                                        n_neighbors=n_neighbors,
                                        resolution=resolution,
                                        samplename=i)
                print("data with imagecol ")
            else:
                print(f"from imagecol data SME_cluster")
                tiling_and_extract(anndata=adata,
                                   TILE_PATH=TILE_PATH,
                                   cnn_base="resnet50")
                adjust_and_louvain(anndata=adata, n_neighbors=n_neighbors, resolution=resolution)
                adata.write(f'{OUT_PATH}/adata_SME.h5ad')
                st.pl.cluster_plot(adata,
                                   use_label='louvain',
                                   show_cluster_labels=cluster_labels,
                                   crop=crop,
                                   fname=f'{OUT_PATH}/adata_SME_louvain_resolution_{resolution}.pdf',
                                   dpi=600,
                                   size=size,
                                   image_alpha=image_alpha,
                                   figsize=(10, 10))
                print("data with imagecol")
            # return adata
        else:
            sample_id = list(adata.uns['spatial'].keys())
            if "library_id" in adata.obs.columns:
                for i in sample_id:
                    print(f"from noimagecol data with library_id subset {i}")
                    adata.uns['spatial'][i]["use_quality"] = "hires"
                    scale = adata.uns["spatial"][i]["scalefactors"]["tissue_hires_scalef"]
                    image_coor = adata.obsm["spatial"][adata.obs['library_id'] == i,] * scale
                    adata.obs.loc[adata.obs['library_id'] == i, "imagecol"] = adata.obsm["spatial"][
                        adata.obs['library_id'] == i, 0]
                    adata.obs.loc[adata.obs['library_id'] == i, "imagerow"] = adata.obsm["spatial"][
                        adata.obs['library_id'] == i, 1]
                    adata_copy = adata.copy()
                    adata_sub = h5ad_subset(anndata=adata_copy, colname='library_id', sample=i)
                    h5ad_adjust_cluster(anndata=adata_sub,
                                        TILE_PATH=TILE_PATH,
                                        OUT_PATH=OUT_PATH,
                                        cnn_base="resnet50",
                                        n_neighbors=n_neighbors,
                                        resolution=resolution,
                                        samplename=i)
                    print(f'till {i}')
                print("from adata.obs['library_id'] and adata.obsm['spatial'] add 'imagecol' and 'imagerow'")
            elif "sample" in adata.obs.columns:
                for i in sample_id:
                    print(f"from noimagecol data with sample subset {i}")
                    adata.uns['spatial'][i]["use_quality"] = "hires"
                    scale = adata.uns["spatial"][i]["scalefactors"]["tissue_hires_scalef"]
                    image_coor = adata.obsm["spatial"][adata.obs['library_id'] == i,] * scale
                    adata.obs.loc[adata.obs['sample'] == i, "imagecol"] = adata.obsm["spatial"][
                        adata.obs['sample'] == i, 0]
                    adata.obs.loc[adata.obs['sample'] == i, "imagerow"] = adata.obsm["spatial"][
                        adata.obs['sample'] == i, 1]
                    adata_copy = adata.copy()
                    adata_sub = h5ad_subset(anndata=adata_copy, colname='sample', sample=i)
                    h5ad_adjust_cluster(anndata=adata_sub,
                                        TILE_PATH=TILE_PATH,
                                        OUT_PATH=OUT_PATH,
                                        cnn_base="resnet50",
                                        n_neighbors=n_neighbors,
                                        resolution=resolution,
                                        samplename=i)
                    print(f"till {i}")
                print("from adata.obs['sample'] and adata.obsm['spatial'] add 'imagecol' and 'imagerow'")
            else:
                print(f"from noimagecol data ")
                adata.uns['spatial'][sample_id[0]]["use_quality"] = "hires"
                scale = adata.uns["spatial"][list(adata.uns['spatial'].keys())[0]]["scalefactors"][
                    "tissue_hires_scalef"]
                image_coor = adata.obsm["spatial"] * scale
                adata.obs["imagecol"] = image_coor[:, 0]
                adata.obs["imagerow"] = image_coor[:, 1]
                tiling_and_extract(anndata=adata, TILE_PATH=TILE_PATH, cnn_base="resnet50")
                adjust_and_louvain(anndata=adata, n_neighbors=n_neighbors, resolution=resolution)
                adata.write(f'{OUT_PATH}/adata_SME.h5ad')
                st.pl.cluster_plot(adata,
                                   use_label='louvain',
                                   show_cluster_labels=cluster_labels,
                                   crop=crop,
                                   fname=f'{OUT_PATH}/adata_SME_louvain.pdf',
                                   dpi=600,
                                   size=size,
                                   image_alpha=image_alpha,
                                   figsize=(10, 10))
                print(f"till {analy_file}")
            print("from adata.obsm['spatial'] add 'imagecol' and 'imagerow'")

    else:
        print("read from h5")
        ##1.读入10x cellranger数据
        adata = st.Read10X(BASE_PATH)
        if barcodes != None:
            bar = pd.read_csv(barcodes, sep='\t')
            intbarcodes = adata.obs_names & bar.iloc[:, 0]
            adata = adata[intbarcodes].copy()
        ##2.数据预处理
        pre_process(adata, 1)
        ##3.图像训练
        tiling_and_extract(anndata=adata, TILE_PATH=TILE_PATH, cnn_base="resnet50")
        ##4.降维聚类
        SME_cluster_louvain(anndata=adata, n_clusters=n_clusters, n_neighbors=n_neighbors, resolution=resolution)
        ##5.可视化及结果保存
        st.pl.cluster_plot(adata,
                           use_label='X_pca_kmeans',
                           show_cluster_labels=cluster_labels,
                           crop=crop,
                           fname=f'{OUT_PATH}/adata_SME_kmeans_{n_clusters}.pdf',
                           dpi=600,
                           size=size,
                           image_alpha=image_alpha,
                           figsize=(10, 10))
        st.pl.cluster_plot(adata,
                           use_label='louvain',
                           show_cluster_labels=cluster_labels,
                           crop=crop,
                           fname=f'{OUT_PATH}/adata_SME_louvain_resolution_{resolution}.pdf',
                           dpi=600,
                           size=size,
                           image_alpha=image_alpha,
                           figsize=(10, 10))
        adata.write(f'{OUT_PATH}/adata_SME.h5ad')
    return adata


# tiling_and_extract(adata,TILE_PATH,cnn_base="resnet50") #cnn_base=['resnet50', 'vgg16', 'inception_v3', 'xception']
if __name__ == '__main__':
    st_SMEclust()
