from sklearn.preprocessing import LabelEncoder as le
from stlearn.plotting import palettes_st
from pathlib import Path
from natsort import natsorted
from oebio.utils.log import getLogger
import argparse
import os
import pandas as pd
import stlearn as st
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import os
import stlearn.tools.microenv.cci.r_helpers as rhs
import sys


import io
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt
from anndata import AnnData
from scanpy.plotting import palettes
from stlearn.plotting import palettes_st
from typing import Optional, Union, Mapping  # Special
from typing import Sequence, Iterable  # ABCs
from typing import Tuple  # Classes
from enum import Enum
from matplotlib import rcParams, ticker, gridspec, axes
from matplotlib.axes import Axes
from abc import ABC

sys.path.append(
    '/public/dev_scRNA/software/conda_envs/stlearn/lib/python3.8/site-packages/stlearn-0.4.11-py3.8.egg/stlearn/tools/microenv/cci/')
from go import run_GO

logger = getLogger('oesc.stlearn_cci')

def rest_index(x):
    x = x[-thelen:]
    return x


def run_lr_go_updata(
    adata,
    r_path,
    n_top=100,
    bg_genes=None,
    min_sig_spots=1,
    species="human",
    p_cutoff=0.01,
    q_cutoff=0.5,
    onts="BP",
    verbose=True
):
    """Runs a basic GO analysis on the genes in the top ranked LR pairs.
        Only supported for human and mouse species.
    Parameters
    ----------
    adata: AnnData
        Must have had st.tl.cci_rank.run() called prior.
    r_path: str
        Path to R, must have clusterProfiler, org.Mm.eg.db, and org.Hs.eg.db
        installed.
    bg_genes: np.array
        Genes to be used as the background. If None, defaults to all genes in
        lr database: 'connectomeDB2020_put'.
    n_top: int
        The top number of LR pairs to use.
    min_sig_spots: int
        Minimum no. of significant spots pairs must have to be considered.
    species: str
        Species to perform the GO testing for.
    p_cutoff: float
        P-value & P-adj cutoff below which results will be returned.
    q_cutoff: float
        Q-value cutoff below which results will be returned.
    onts: str
        As per clusterProfiler; One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
    Returns
    -------
    adata: AnnData
        Relevant information stored in adata.uns['lr_go']
    """
    #### Making sure inputted correct species ####
    all_species = ["human", "mouse"]
    if species not in all_species:
        raise Exception(f"Got {species} for species, must be one of " f"{all_species}")

    #### Getting the genes from the top LR pairs ####
    if "lr_summary" not in adata.uns:
        raise Exception("Need to run st.tl.cci.run first.")
    lrs = adata.uns["lr_summary"].index.values.astype(str)
    n_sig = adata.uns["lr_summary"].loc[:, "n_spots_sig"].values.astype(int)
    top_lrs = lrs[n_sig > min_sig_spots][0:n_top]
    top_genes = np.unique([lr.split("_") for lr in top_lrs])

    ## Determining the background genes if not inputted ##
    if type(bg_genes) == type(None):
        all_lrs = st.tl.cci.load_lrs("connectomeDB2020_put", species=species)
        bg_genes = np.unique([lr_.split("_") for lr_ in all_lrs])

    #### Running the GO analysis ####
    go_results = run_GO(
        top_genes,
        bg_genes,
        species,
        r_path,
        p_cutoff=p_cutoff,
        q_cutoff=q_cutoff,
        onts=onts,
    )
    adata.uns["lr_go"] = go_results
    if verbose:
        print("GO results saved to adata.uns['lr_go']")


def get_cmap(cmap):
    """Checks inputted cmap string."""
    if cmap == "vega_10_scanpy":
        cmap = palettes.vega_10_scanpy
    elif cmap == "vega_20_scanpy":
        cmap = palettes.vega_20_scanpy
    elif cmap == "default_102":
        cmap = palettes.default_102
    elif cmap == "default_28":
        cmap = palettes.default_28
    elif cmap == "jana_40":
        cmap = palettes_st.jana_40
    elif cmap == "default":
        cmap = palettes_st.default
    elif type(cmap) == str:  # If refers to matplotlib cmap
        cmap_n = plt.get_cmap(cmap).N
        return plt.get_cmap(cmap), cmap_n
    elif type(cmap) == matplotlib.colors.LinearSegmentedColormap:  # already cmap
        cmap_n = cmap.N
        return cmap, cmap_n

    cmap_n = len(cmap)
    cmaps = matplotlib.colors.LinearSegmentedColormap.from_list("", cmap)

    cmap_ = plt.cm.get_cmap(cmaps)

    return cmap_, cmap_n


def check_cmap2(cmap):
    """Initialize cmap"""
    scanpy_cmap = ["vega_10_scanpy", "vega_20_scanpy", "default_102", "default_28"]
    stlearn_cmap = ["jana_40", "default"]
    cmap_available = plt.colormaps() + scanpy_cmap + stlearn_cmap
    error_msg = (
        "cmap must be a matplotlib.colors.LinearSegmentedColormap OR"
        "one of these: " + str(cmap_available)
    )
    if type(cmap) == str:
        assert cmap in cmap_available, error_msg
    elif type(cmap) != matplotlib.colors.LinearSegmentedColormap:
        raise Exception(error_msg)

    return cmap


def get_colors2(adata, obs_key, cmap="default", label_set=None):
    """Retrieves colors if present in adata.uns, if not present then will set
    them as per scanpy & return in order requested.
    """
    # Checking if colors are already set #
    col_key = f"{obs_key}_colors"
    if col_key in adata.uns:
        labels_ordered = adata.obs[obs_key].cat.categories
        colors_ordered = adata.uns[col_key]
    else:  # Colors not already present
        check_cmap2(cmap)
        cmap, cmap_n = get_cmap(cmap)

        if not hasattr(adata.obs[obs_key], "cat"):  # Ensure categorical
            adata.obs[obs_key] = adata.obs[obs_key].astype("category")
        labels_ordered = adata.obs[obs_key].cat.categories
        colors_ordered = [
            matplotlib.colors.rgb2hex(cmap(i / (len(labels_ordered) - 1)))
            for i in range(len(labels_ordered))
        ]
        adata.uns[col_key] = colors_ordered

    # Returning the colors of the desired labels in indicated order #
    if type(label_set) != type(None):
        colors_ordered = [
            colors_ordered[np.where(labels_ordered == label)[0][0]]
            for label in label_set
        ]
    return colors_ordered


def cci_check2(
    adata,
    use_label,
    figsize=(16, 10),
    cell_label_size=20,
    axis_text_size=18,
    tick_size=14,
    show=True,
):
    """Checks relationship between no. of significant CCI-LR interactions and cell type frequency.

    Parameters
    ----------
    adata: AnnData
        Data on which st.tl.cci.run & st.tl.cci.run_cci has been performed.
    use_label: str
        The cell type label information used when running st.tl.cci.run_cci
    figsize: tuple
        Size of outputted figure.
    cell_label_size: int
        Size of the cell labels put on top of the bar chart.
    axis_text_size: int
        Size of the axis text.
    tick_size: int
        Size of the ticks displayed at bottom of chart.
    show: bool
        Whether to show the plot or not; if false returns figure & axes.
    Returns
    -------
    Figure, Ax1, Ax2
        The figure, axes for the barchart, and twin axes for the lineplot.
    """
    labels = adata.obs[use_label].values.astype(str)
    label_set = np.array(list(adata.obs[use_label].cat.categories))
    colors = get_colors2(adata, use_label)
    xs = np.array(list(range(len(label_set))))
    int_dfs = adata.uns[f"per_lr_cci_{use_label}"]

    # Counting!!! #
    cell_counts = []  # Cell type frequencies
    cell_sigs = []  # Cell type significant interactions
    for j, label in enumerate(label_set):
        counts = sum(labels == label)
        cell_counts.append(counts)

        int_count = 0
        for lr in int_dfs:
            int_df = int_dfs[lr]
            label_index = np.where(int_df.index.values == label)[0][0]
            int_bool = int_df.values > 0
            int_count += sum(int_bool[label_index, :])
            int_count += sum(int_bool[:, label_index])
            # prevent double counts
            int_count -= int_bool[label_index, label_index]

        cell_sigs.append(int_count)

    cell_counts = np.array(cell_counts)
    cell_sigs = np.array(cell_sigs)
    order = np.argsort(cell_counts)
    cell_counts = cell_counts[order]
    cell_sigs = cell_sigs[order]
    colors = np.array(colors)[order]
    label_set = label_set[order]

    # Plotting bar plot #
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(xs, cell_counts, color=colors)
    text_dist = max(cell_counts) * 0.015
    fontdict = {"fontweight": "bold", "fontsize": cell_label_size}
    for j in range(len(xs)):
        ax.text(
            xs[j],
            cell_counts[j] + text_dist,
            label_set[j],
            rotation=90,
            fontdict=fontdict,
        )
    axis_text_fp = {"fontweight": "bold", "fontsize": axis_text_size}
    ax.set_ylabel("Cell type frequency", color="black", **axis_text_fp)
    ax.spines["top"].set_visible(False)
    ax.tick_params(labelsize=tick_size)
    ax.set_xlabel("Cell type rank", **axis_text_fp)

    # Line-plot of the interaction counts #
    ax2 = ax.twinx()
    ax2.set_ylabel("CCI-LR interactions", color="blue", **axis_text_fp)
    ax2.plot(xs, cell_sigs, color="blue", linewidth=2)
    ax2.tick_params(axis="y", labelcolor="blue", labelsize=tick_size)
    ax2.spines["top"].set_visible(False)
    ax2.tick_params(labelsize=tick_size)
    fig.tight_layout()

    if show:
        plt.show()
    else:
        return fig, ax, ax2


def savefig2(filepath, bbox_inches, dpi):
    plt.savefig(f'{filepath}.pdf', bbox_inches=bbox_inches, dpi=dpi)
    plt.savefig(f'{filepath}.png', bbox_inches=bbox_inches, dpi=dpi)


stcci = argparse.ArgumentParser(description='stlearn cci')
stcci.add_argument('--anndatafile', type=str, default=None,
                   help="the output of sapacerange")
stcci.add_argument('--n_pairs', type=int, default=10000,
                   help="Number of random pairs to generate; low as example, recommend to 10000")
stcci.add_argument('--labpath', type=str, default=None,
                   help="the result of spotlight or RCTD")
stcci.add_argument('--l_r_database', type=str, default='connectomeDB2020_lit',
                   help="input the database that can be  built by customer(default: connectomeDB2020_lit)")
stcci.add_argument('--top_n', type=int, default='10',
                   help="choose some top LRs to visualization")
stcci.add_argument('--specie', type=str, default='human')
stcci.add_argument('--outpath', type=str, default=None)

args = stcci.parse_args()

# 数据预处理##############################################################################################################
logger.info("step1.数据预处理")
# 数据集读取和处理
outpath = Path(args.outpath)
if not os.path.exists(outpath):
    os.makedirs(outpath)

anndata = st.Read10X(args.anndatafile)
# st.pp.filter_genes(anndata, min_cells=3)
st.pp.normalize_total(anndata)
top_n = args.top_n

# 细胞类型表格读取，支持两种分隔符 "\t" ","
lab = pd.read_csv(args.labpath, index_col=0, sep="\t")
thelen = len(anndata.obs_names[1])
if lab.shape[1] > 0:
    print("sep=\\\t")
else:
    del lab
    lab = pd.read_csv(args.labpath, index_col=0, sep=",")

# 使数据集和细胞类型表格中barcodes保持一致
if lab.index[0][-2:] == "-1":
    res_ind = list(map(rest_index, lab.index))
    lab.index = res_ind
else:
    lab.index = lab.index + "-1"
    res_ind = list(map(rest_index, lab.index))
    lab.index = res_ind
intbarcode = list(set(anndata.obs_names).intersection(set(lab.index)))
anndata = anndata[intbarcode, :]
lab = lab.loc[intbarcode, :]

logger.info("step2.导入细胞类型表格信息")
# 导入细胞类型表格信息
if lab.shape[1] == 1:  # 只有细胞类型注释信息没有占比信息的情况
    anndata.obs['cell_type'] = lab.iloc[:, 0].astype('category')
else:
    if "predicted.id" in lab.columns:  # 如果有多余信息则去除之后再导入
        anndata.obsm['cell_type'] = lab.drop(["predicted.id", "prediction.score.max"], axis=1)
        anndata.obsm["deconvolution"] = lab.drop(["predicted.id", "prediction.score.max"], axis=1)
        if "prediction.score.max" in lab.columns:
            anndata.uns['cell_type'] = lab.drop(["predicted.id", "prediction.score.max"], axis=1)
        else:
            anndata.uns['cell_type'] = lab.drop(["predicted.id"], axis=1)
    else:
        anndata.obsm['cell_type'] = lab
        anndata.obsm["deconvolution"] = lab
        anndata.uns['cell_type'] = lab
        anndata.uns['cell_type'] = anndata.uns['cell_type'].apply(pd.to_numeric, axis=1)

    anndata.obs['cell_type'] = lab.idxmax(axis=1).astype('category')

# 获取显著LR对############################################################################################################
logger.info("step3.获取显著LR对")
LR_outpath = Path(f'{outpath}/1.LR_analysis_results')
if not os.path.exists(LR_outpath):
    os.makedirs(LR_outpath)

# 选取database中的文件名称
lrs = st.tl.cci.load_lrs([args.l_r_database], species=args.specie)

# 根据一个中心点，扫描出它所有的近邻点。
# 对于这个近邻点，它的近邻点中表达高于阈值的R的点占比，与近邻点中L表达高于阈值的点的占比，相加为count，构建矩阵——行是L-R类型，列是spot的坐标，值为L-R对应于该中心点的count.
st.tl.cci.run(anndata, lrs,
              min_spots=20,  # Filter out any LR pairs with no scores for less than min_spots
              distance=None,  # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
              n_pairs=args.n_pairs,  # Number of random pairs to generate; low as example, recommend ~10,000
              n_cpus=4,  # Number of CPUs for parallel. If None, detects & use all available.
              pval_adj_cutoff=0.05
              )

lr_info = anndata.uns['lr_summary']
lr_info.to_csv(f'{LR_outpath}/LRs_summary.csv', sep="\t")

# LR统计学指标空间映射图
st.tl.cci.adj_pvals(anndata, correct_axis='spot',
                    pval_adj_cutoff=0.05, adj_method='fdr_bh')

# LR显著性排序图
st.pl.lr_summary(anndata, n_top=50, figsize=(10, 3))
savefig2(f'{LR_outpath}/LRs_rank', bbox_inches="tight", dpi=300)

# LR分析诊断图，如果相关性较强需要调高npairs参数
st.pl.lr_diagnostics(anndata, figsize=(10, 2.5))
savefig2(f'{LR_outpath}/LRs_diagnostics', bbox_inches="tight", dpi=300)

logger.info("step4.显著LR对GO作图")
# GO作图
r_path = '/public/dev_scRNA/software/conda_envs/stlearn/lib/R/'
run_lr_go_updata(anndata, r_path, n_top=100, min_sig_spots=1, species=args.specie, onts='ALL')

# GO Rank 根据p.adjust排序结果展示
st.pl.lr_go(anndata, lr_text_fp={'weight': 'bold', 'size': 10}, rot=15,
            figsize=(12, 9), n_top=15, show=False)
savefig2(f'{LR_outpath}/LRs_GO_top100', bbox_inches="tight", dpi=300)
LR_GO = anndata.uns['lr_go']
LR_GO.to_csv(f'{LR_outpath}/LRs_GO_top100.csv', sep="\t")

# 选取显著性前top_n的LR对进一步可视化===========================================================================================
logger.info(f'step5.选取显著性前{top_n}的LR对进一步可视化')
# 所有spots可视化
output_topL_R = Path(f'{LR_outpath}/top{top_n}_LRs_visualization_all')
if not os.path.exists(output_topL_R):
    os.makedirs(output_topL_R)

best_lr = anndata.uns['lr_summary'].iloc[0:top_n, :].index

stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
fig, axes = plt.subplots(ncols=len(stats), figsize=(16, 6))

for j in range(top_n):
    for i, stat in enumerate(stats):
        st.pl.lr_result_plot(anndata, use_result=stat, use_lr=best_lr[j], show_color_bar=False, ax=axes[i])
        axes[i].set_title(f'{best_lr[j]} {stat}')
        savefig2(f'{output_topL_R}/{best_lr[j]}', bbox_inches="tight", dpi=300)

# 显著spots可视化
output_topL_R_sig = Path(f'{LR_outpath}/top{top_n}_LRs_visualization_sig')
if not os.path.exists(output_topL_R_sig):
    os.makedirs(output_topL_R_sig)

for j in range(top_n):
    fig, axes = plt.subplots(ncols=2, figsize=(8, 6))
    st.pl.lr_result_plot(anndata, use_result='-log10(p_adjs)', use_lr=best_lr[j], show_color_bar=False, ax=axes[0])
    st.pl.lr_result_plot(anndata, use_result='lr_sig_scores', use_lr=best_lr[j], show_color_bar=False, ax=axes[1])
    axes[0].set_title(f'{best_lr[j]} -log10(p_adjs)')
    axes[1].set_title(f'{best_lr[j]} lr_sig_scores')
    savefig2(f'{output_topL_R_sig}/{best_lr[j]}', bbox_inches="tight", dpi=300)

# 加入细胞类型信息，获取细胞通讯分析结果#######################################################################################
logger.info("step6.获取细胞通讯分析结果")
CCI_outpath = Path(f'{outpath}/2.CCI_analysis_results')
if not os.path.exists(CCI_outpath):
    os.makedirs(CCI_outpath)

lrs = anndata.uns['lr_summary'].index.values[0:top_n]
st.tl.cci.run_cci(anndata, 'cell_type',  # Spot cell information either in data.obs or data.uns
                  min_spots=3,  # Minimum number of spots for LR to be tested.
                  spot_mixtures=True,  # If True will use the label transfer scores,
                  # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.2,  # Spot considered to have cell type if score>0.2
                  sig_spots=True,  # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=args.n_pairs  # Permutations of cell information to get background, recommend ~1000
                  )

logger.info("step7.细胞通讯分析结果可视化")
# CCI分析诊断图，如果相关性较强需要调高n_perms参数
if 'cell_type' in anndata.uns:
    cci_check2(anndata, 'cell_type')
    savefig2(f'{CCI_outpath}/CCI_diagnostics', bbox_inches="tight", dpi=300)

# 网状图=================================================================================================================

# 所有LR对
pos_1 = st.pl.ccinet_plot(anndata, 'cell_type', return_pos=True)
savefig2(f'{CCI_outpath}/ALL_LRs_ccinet_plot', bbox_inches="tight", dpi=300)

# topLR对
output_ccinet = Path(Path(f'{CCI_outpath}/top{top_n}_LRs_ccinet_plot'))
if not os.path.exists(output_ccinet):
    os.makedirs(output_ccinet)

# fig, axes = plt.subplots(ncols=4,nrows=5, figsize=(16,20))
for lr in lrs:
    st.pl.ccinet_plot(anndata, 'cell_type', lr, min_counts=2, figsize=(10, 7.5), pos=pos_1)
    savefig2(f'{output_ccinet}/CCI_{lr}_ccinet_plot', bbox_inches="tight", dpi=300)

# 弦状图=================================================================================================================

# 所有LR对
st.pl.lr_chord_plot(anndata, 'cell_type', label_size=8)
savefig2(f'{CCI_outpath}/ALL_LRs_chord_plot', bbox_inches="tight", dpi=300)

# topLR对
output_chord = Path(Path(f'{CCI_outpath}/top{top_n}_LRs_chord_plot'))
if not os.path.exists(output_chord):
    os.makedirs(output_chord)

for lr in lrs:
    st.pl.lr_chord_plot(anndata, 'cell_type', lr)
    savefig2(f'{output_chord}/CCI_{lr}_chord_plot', bbox_inches="tight", dpi=300)

# 热图===================================================================================================================

# topLR对的细胞类型互作热图
st.pl.lr_cci_map(anndata, 'cell_type', lrs=lrs, min_total=100, figsize=(20, 4))
savefig2(f'{CCI_outpath}/top{top_n}_LRs_heatmap', bbox_inches="tight", dpi=300)

# 考虑所有LR对的细胞类型热图
st.pl.cci_map(anndata, 'cell_type')
savefig2(f'{CCI_outpath}/ALL_LRs_heatmap', bbox_inches="tight", dpi=300)

# topLR对
output_heatmap = Path(f'{CCI_outpath}/top{top_n}_LRs_heatmap')
if not os.path.exists(output_heatmap):
    os.makedirs(output_heatmap)

for lr in lrs:
    st.pl.cci_map(anndata, 'cell_type', lr)
    savefig2(f'{output_heatmap}/CCI_{lr}_heatmap', bbox_inches="tight", dpi=300)

logger.info("step8.保存数据集")
# 保存数据集
lib = list(anndata.uns['spatial'].keys())[0]

anndata.write(f'{outpath}/{lib}_anndata_with_cci.h5ad')
