#!/usr/bin/env python
# coding: utf-8
import os, click
import stlearn as st
import scanpy as sc
import numpy as np
import pandas as pd
import warnings
import networkx as nx
import matplotlib.pyplot as plt
from oebio.utils.log import getLogger

from PIL import Image
from pathlib import Path
from scipy.stats import spearmanr
from fa2.fa2util import Node

logger = getLogger('oesc.stlearn_psd')


def _read_graph_local(anndata, graph_type):
    if graph_type == "PTS_graph":
        graph = nx.from_scipy_sparse_array(anndata.uns[graph_type]["graph"], create_using=nx.DiGraph)
    else:
        graph = nx.from_scipy_sparse_array(anndata.uns[graph_type]["graph"])
    node_dict = anndata.uns[graph_type]["node_dict"]
    node_dict = {int(k): int(v) for k, v in node_dict.items()}
    relabel_graph = nx.relabel_nodes(graph, node_dict)
    return relabel_graph


@click.command()
@click.option('--file', help="the absolute path of h5ad file")
@click.option('--use_label', default="louvain", help="the absolute path of h5ad file")
@click.option('--iroot', type=click.INT, default=0, help="the start of paga")
@click.option('--output', default="./", help="h5ad file or a file containing h5 file and the spatial file")
@click.option('--globlist', type=click.STRING, help="cluster list for glob trajectory")
@click.option('--barcodefile', type=click.STRING, default=None, help="to subset barcodes")
@click.option('--genelist', type=click.STRING, default=None, help="to plot your gene list in clade_X")
def st_trajectory(file, use_label, iroot, globlist, output, barcodefile, genelist):
    if os.path.exists(output):
        output = Path(output)
    else:
        os.makedirs(output)
        output = Path(output)
    if genelist is not None:
        import sys
        sys.path.append('/public/dev_scRNA/yfang/project/')
        import tra_gene_plot as tra
        geneplot = pd.read_csv(genelist)
    logger.info("step1.输出文件夹创建")
    sc.settings.verbosity = 3
    sc.settings.figdir = f"{output}/1.no_spatial_pseudotime_plot"
    logger.info("step2.读入h5ad文件")
    anndata = sc.read_h5ad(file)
    if barcodefile is not None:
        barcodes = pd.read_csv(barcodefile, sep="\t")
        intbarcodes = anndata.obs_names & barcodes.iloc[:, 0]
        anndata = anndata[intbarcodes]
    logger.info(f"step3.PAGA及Diffusion map:设置PAGA的iroot为cluster{iroot}: 第1个细胞")
    anndata.uns['iroot'] = np.flatnonzero(anndata.obs[use_label] == str(iroot))[1]
    st.spatial.trajectory.pseudotime(anndata,
                                     eps=50,
                                     use_rep="X_pca",
                                     use_label=use_label,
                                     use_sme=True)
    sc.pl.paga(anndata,
               color=use_label,
               save=f"_{use_label}.pdf",
               title="PAGA plot")
    sc.tl.draw_graph(anndata, init_pos="paga")
    logger.info("Gene expression (reduced dimension) plot:")
    sc.pl.draw_graph(anndata,
                     color=use_label,
                     title="Reduced dimension plot",
                     size=8,
                     legend_loc="on data",
                     save=f"_{use_label}.pdf")

    logger.info("Diffusion pseudotime plot:")
    sc.pl.draw_graph(anndata,
                     color="dpt_pseudotime",
                     title="Diffusion pseudotime plot",
                     size=8,
                     save=f"_diffusion_pseudotime_plot.pdf")

    print("Trajectory pseudotime")
    st.pl.trajectory.pseudotime_plot(anndata,
                                     use_label=use_label,
                                     pseudotime_key="dpt_pseudotime",
                                     show_graph=True,
                                     # show_trajectory=False,
                                     node_alpha=1,
                                     edge_alpha=0.4,
                                     node_size=6,
                                     cropped=False,
                                     margin=200,
                                     name="diffusion_pseudotime_plot.pdf",
                                     dpi=600,
                                     output=f"{output}/1.no_spatial_pseudotime_plot")
    ##4.局部轨迹分析：针对各个cluster, 部分cluster可能不适用，如果报错将自动跳过
    logger.info("step4: excute trajectory local: 针对各个cluster, 部分cluster可能不适用，如果报错将自动跳过")
    list_clusters = list(range(1, len(anndata.obs[use_label].unique())))
    if not os.path.exists(f"{output}/2.spatial_trajectory_plot/trajectory_for_each_cluster"):
        os.makedirs(f"{output}/2.spatial_trajectory_plot/trajectory_for_each_cluster")
    for cluster in list_clusters:
        try:
            st.spatial.trajectory.local_level(anndata,
                                              use_label=use_label,
                                              cluster=cluster)
            st.pl.subcluster_plot(anndata,
                                  use_label=use_label,
                                  cluster=cluster,
                                  crop=False,
                                  dpi=600,
                                  size=8,
                                  fname=f"{output}/2.spatial_trajectory_plot/trajectory_for_each_cluster"
                                        f"/cluster_{cluster}_{use_label}_subcluster.pdf")
            st.spatial.trajectory.pseudotimespace_local(anndata,
                                                        use_label=use_label,
                                                        cluster=cluster)
            st.pl.trajectory.local_plot(anndata,
                                        use_cluster=cluster,
                                        branch_alpha=1,
                                        dpi=600,
                                        spot_size=1,
                                        reverse=False,
                                        output=f"{output}/2.spatial_trajectory_plot/trajectory_for_each_cluster/",
                                        name=f"cluster_{cluster}_trajectory.local_plot.pdf",
                                        show_plot=True)
        except Exception as e:
            pass
        continue

    ##5.全局轨迹分析：选定clusters,至少2个cluster
    if globlist:
        logger.info("step5: excute trajectory global：针对选定clusters进行全局轨迹分析,至少包含2个cluster")
        logger.info(f"step5.1.全局轨迹分析:{globlist}")
        if not os.path.exists(
            f"{output}/2.spatial_trajectory_plot/for_clusters_{globlist.replace(',', '_')}"):
            os.makedirs(f"{output}/2.spatial_trajectory_plot/for_clusters_{globlist.replace(',', '_')}")
        list_clusters = globlist.split(",")
        st.spatial.trajectory.pseudotimespace_global(anndata,
                                                     use_label=use_label,
                                                     list_clusters=list_clusters)
        st.pl.cluster_plot(anndata,
                           use_label=use_label,
                           show_trajectories=True,
                           list_clusters=list_clusters,
                           show_subcluster=True,
                           crop=False,
                           dpi=600,
                           size=10,
                           fname=f"{output}/2.spatial_trajectory_plot"
                                 f"/for_clusters_{globlist.replace(',', '_')}/{use_label}_trajectories.pdf")
        st.pl.trajectory.tree_plot(anndata,
                                   output=f"{output}/2.spatial_trajectory_plot"
                                          f"/for_clusters_{globlist.replace(',', '_')}",
                                   name=f"{globlist.replace(',', '_')}_tree_plot.pdf")

        logger.info("step5.2.全局轨迹分析:不同轨迹的渐变基因分析")
        tmp = _read_graph_local(anndata=anndata, graph_type="PTS_graph")
        G = tmp.copy()
        remove = [edge for edge in G.edges if 9999 in edge]
        G.remove_edges_from(remove)
        G.remove_node(9999)
        print(int(list_clusters[0]))
        print(type(int(list_clusters[0])))
        for i in range(len(anndata.uns['split_node'][int(list_clusters[0])])):
            node = int(anndata.uns['split_node'][int(list_clusters[0])][i])
            nodes = []
            for edge in G.edges([node]):
                nodes.append(edge[0])
                nodes.append(edge[1])
            nodes = list(set(nodes))
            if len(nodes) > 0:
                st.spatial.trajectory.detect_transition_markers_clades(anndata,
                                                                       clade=node,
                                                                       use_raw_count=False,
                                                                       cutoff_spearman=0.3)
                st.pl.trajectory.transition_markers_plot(anndata,
                                                         top_genes=30,
                                                         trajectory=f"clade_{node}",
                                                         dpi=600,
                                                         output=f"{output}/2.spatial_trajectory_plot"
                                                                f"/for_clusters_{globlist.replace(',', '_')}",
                                                         name=f"subcluster_{node}.pdf")
                anndata.uns[f"clade_{node}"].to_csv(
                    f"{output}/2.spatial_trajectory_plot"
                    f"/for_clusters_{globlist.replace(',', '_')}/subcluster_{node}_gene.csv",
                    sep="\t", index=False)
                if genelist is not None:
                    st.pl.trajectory.transition_markers_plot(adata=anndata,
                                                             defgenelist=geneplot,
                                                             trajectory=f"clade_{node}",
                                                             dpi=600,
                                                             output=f"{output}/2.spatial_trajectory_plot"
                                                                    f"/for_clusters_{globlist.replace(',', '_')}",
                                                             name=f"clade_{node}_withdefgene.pdf")
    # anndata.write(f"{output}/anndata_sj.h5ad") #h5py版本过高导致，由于其它包需要该版本无法降级，暂时无法解决，因此不再保存
    return anndata


if __name__ == '__main__':
    st_trajectory()
