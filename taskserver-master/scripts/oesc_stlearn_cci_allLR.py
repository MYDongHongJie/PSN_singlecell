#!/usr/bin/env python
# coding: utf-8




import click


import os
import sys
import stlearn as st
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
from typing import Optional, Union
from anndata import AnnData
from natsort import natsorted

sys.path.append('/data/software/conda_envs/scanpy/lib/python3.8/site-packages/stlearn/tools/microenv/cci/')
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.spatial as spatial
from anndata import AnnData
from het import create_grids


def lr_test(
    adata: AnnData,
    use_lr: str = "cci_lr",
    distance: float = None,
    verbose: bool = True,
    L_R_pirs:str = "lr"
) -> AnnData:

    """Calculate the proportion of known ligand-receptor co-expression among the neighbouring spots or within spots
    Parameters
    ----------
    adata: AnnData          The data object to scan
    use_lr: str             object to keep the result (default: adata.uns['cci_lr'])
    distance: float         Distance to determine the neighbours (default: closest), distance=0 means within spot

    Returns
    -------
    adata: AnnData          The data object including the results
    """

    # automatically calculate distance if not given, won't overwrite distance=0 which is within-spot
    if not distance and distance != 0:
        # for arranged-spots
        scalefactors = next(iter(adata.uns["spatial"].values()))["scalefactors"]
        library_id = list(adata.uns["spatial"].keys())[0]
        distance = (
            scalefactors["spot_diameter_fullres"]
            * scalefactors[
                "tissue_" + adata.uns["spatial"][library_id]["use_quality"] + "_scalef"
            ]
            * 2
        )

    df = adata.to_df()

    # expand the LR pairs list by swapping ligand-receptor positions
    lr_pairs = adata.uns[L_R_pirs].copy()
    lr_pairs += [item.split("_")[1] + "_" + item.split("_")[0] for item in lr_pairs]
    print(lr_pairs)

    # get neighbour spots for each spot according to the specified distance
    coor = adata.obs[["imagerow", "imagecol"]]
    point_tree = spatial.cKDTree(coor)
    neighbours = []
    for spot in adata.obs_names:
        if distance == 0:
            neighbours.append([spot])
        else:
            n_index = point_tree.query_ball_point(
                np.array(
                    [adata.obs["imagerow"].loc[spot], adata.obs["imagecol"].loc[spot]]
                ),
                distance,
            )
            neighbours.append(
                [item for item in df.index[n_index] if not (item == spot)]
            )

    # filter out those LR pairs that do not exist in the dataset
    lr1 = [item.split("_")[0] for item in lr_pairs]
    lr2 = [item.split("_")[1] for item in lr_pairs]
    avail = [
        i for i, x in enumerate(lr1) if lr1[i] in df.columns and lr2[i] in df.columns
    ]
    spot_lr1 = df[[lr1[i] for i in avail]]
    spot_lr2 = df[[lr2[i] for i in avail]]
    if verbose:
        print("Altogether " + str(len(avail)) + " valid L-R pairs")

    # function to calculate mean of lr2 expression between neighbours or within spot (distance==0) for each spot
    def mean_lr2(x):
        # get lr2 expressions from the neighbour(s)
        nbs = spot_lr2.loc[neighbours[df.index.tolist().index(x.name)], :]
        if nbs.shape[0] > 0:  # if neighbour exists
            return (nbs > 0).sum() / nbs.shape[0]
        else:
            return 0

    # mean of lr2 expressions from neighbours of each spot
    nb_lr2 = spot_lr2.apply(mean_lr2, axis=1)

    # check whether neighbours exist
    try:
        nb_lr2.shape[1]
    except:
        raise ValueError("No neighbours found within given distance.")

    # keep value of nb_lr2 only when lr1 is also expressed on the spots
    spot_lr = pd.DataFrame(
        spot_lr1.values * (nb_lr2.values > 0) + (spot_lr1.values > 0) * nb_lr2.values,
        index=df.index,
        columns=[lr_pairs[i] for i in avail],
    )
    num= spot_lr.shape[1]/2
    spot_lr_sum=spot_lr.iloc[:,0:int(num)].values+spot_lr.iloc[:,int(num):(spot_lr.shape[1]+1)].values
    spot_lr2=pd.DataFrame(spot_lr_sum/2,
                          index=df.index,
                          columns=spot_lr.columns[0:int(num)],)

    
    return spot_lr2



# In[555]:


import sys, os, random, scipy
sys.path.append('/data/software/conda_envs/scanpy/lib/python3.8/site-packages/stlearn/tools/microenv/cci/')
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from anndata import AnnData
import numpy as np
import pandas as pd
from anndata import AnnData





import scanpy as sc
import stlearn as st
def read_and_addlabs(anndatafile,labpath):
    
    
    """
        anndatafile:输入的anndata文件绝对路径,或者保存.h5和spatial文件夹的文件夹
        labpath对应的细胞标签。example:表头是barecode,行名为细胞类型。值是没给barcode中细胞比例
        
    """
    if "h5ad" in anndatafile:
        anndata=sc.read_h5ad(anndatafile)
        outlist=anndatafile.split("/")
        outpath="/".join(outlist[0:len(outlist)-1])
        outpath=Path(f"{outpath}/cci")
    else:
        anndata = st.Read10X(anndatafile)
        st.pp.filter_genes(anndata,min_cells=3)
        st.pp.normalize_total(anndata)
        st.pp.scale(anndata)
        outpath=Path(f"{Path(anndatafile)}/cci")
    outpath.mkdir(parents=True, exist_ok=True)
    lab=pd.read_csv(labpath,index_col=0,sep="\t")
    
    if lab.shape[1]>0:
        print("sep=\\\t")
    else:
        del lab
        lab=pd.read_csv(labpath,index_col=0,sep=",")
    def rest_index(x):
        x=x[-18:]
        return x

    if lab.index[0][-2:]=="-1":
        res_ind=list(map(rest_index,lab.index))
        lab.index=res_ind        
    else:
        lab.index=lab.index+"-1"
        res_ind=list(map(rest_index,lab.index))
        lab.index=res_ind
    intbarcode=list(set(anndata.obs_names).intersection(set(lab.index))) 
    anndata=anndata[intbarcode,:]
    lab=lab.loc[intbarcode,:]
    if "predicted.id" in lab.columns :
        anndata.obsm["deconvolution"]=lab.drop(["predicted.id", "prediction.score.max"], axis=1)
    else:
        anndata.obsm["deconvolution"] = lab
    st.pl.deconvolution_plot(anndata,threshold=0.1,output=str(outpath),name='deconvolution_plot.pdf',dpi=300)
    
    if "predicted.id" in lab.columns:
        if "prediction.score.max" in lab.columns:
            anndata.uns["label_transfer"] = lab.drop(["predicted.id", "prediction.score.max"], axis=1)
        else:
            anndata.uns["label_transfer"] = lab.drop(["predicted.id"], axis=1)

        key_add = "predictions"
        key_source = "predicted.id"
        anndata.obs[key_add] = pd.Categorical(
            values=np.array(lab[key_source]).astype("U"),
            categories=natsorted(lab[key_source].unique().astype("U")),
        )
    else:
        anndata.uns["label_transfer"] = lab
        anndata.uns["label_transfer"]=anndata.uns["label_transfer"].apply(pd.to_numeric,axis=1)

        lab["predicted.id"]=lab.idxmax(axis=1)
        key_add = "predictions"
        key_source = "predicted.id"
        anndata.obs[key_add] = pd.Categorical(values=np.array(lab[key_source]).astype("U"),categories=natsorted(lab[key_source].unique().astype("U")),)
        #st.pl.cluster_plot(adata,use_label="predictions",fname=f'{outpath}/celltype_cluster.pdf',dpi=300)
    return anndata
            
    


# In[ ]:


@click.command()
@click.option('--anndatafile',default="data_SME.h5ad",type=click.STRING,help="输入的anndata文件绝对路径,或者保存.h5和spatial文件夹的文件夹")
@click.option('--labpath',default="label_transfer_bc.csv",type=click.STRING,help="labpath对应的细胞标签。example:表头是barecode,行名为细胞类型。值是没给barcode中细胞比例")
@click.option('--sig_mean',default='./cellphonedb_results/out/significant_means.txt',type=click.STRING,help="the output file of cellphonedb called significant_means.txt")
@click.option('--save',default='FALSE',type=click.STRING,help="是否保存h5ad文件")


# In[597]:


def spatial_cci(anndatafile,labpath,sig_mean,save):
    anndata=read_and_addlabs(anndatafile=anndatafile,labpath=labpath)
    if "h5ad" in anndatafile:
        outlist=anndatafile.split("/")
        outpath="/".join(outlist[0:len(outlist)-1])
        outpath=Path(f"{outpath}/cci")
    else:
        outpath=Path(f"{Path(anndatafile)}/cci")


    st.add.lr(anndata,db_filepath=sig_mean,source="cellphonedb", copy=False)
    anndata_df=lr_test(anndata)


    st.pl.cluster_plot(anndata,use_label="predictions",fname=f'{outpath}/celltype_cluster.pdf',dpi=300)
    print("点间细胞类型热图")
    st.tl.cci.het.count(anndata, use_label='label_transfer',use_het="cci_het_betweenspots")
    st.pl.het_plot(anndata, use_het='cci_het_betweenspots',fname=f'{outpath}/cci_het_betweenspots.pdf',dpi=300)

    def mut(x):
        return np.multiply(anndata.obsm['cci_het_betweenspots'],x)

    merge_d=anndata_df.apply(mut,axis=0)

    merge_d.to_csv(f'{outpath}/adata_cci_result.csv',sep=",")

    sig_mean_filter=anndata.uns['cpdb'].iloc[list(anndata.uns['cpdb']['interacting_pair'].isin(merge_d.columns)),:]
    inter_cell=list(sig_mean_filter.columns[12:sig_mean_filter.shape[1]])
    data_select=pd.DataFrame()
    for i in inter_cell:
        paires=sig_mean_filter.sort_values(by=['rank',i],ascending=[True,False])['interacting_pair'][0:5]
        paires_d=merge_d[list(paires)]
        paires_d.columns=i+"_"+(merge_d[list(paires)].columns)
        data_select=pd.concat([data_select, paires_d],axis=1)

    anndata.obs=pd.concat([anndata.obs, data_select],axis=1)
    os.chdir(outpath)
    for x,i in enumerate(data_select.columns):
        sc.pl.spatial(anndata, img_key="lowres", color=[i],alpha_img =0.5,size=1.5,save=f"{x}")
    
    anndata.uns['significant_means_merged']=merge_d
    if save == "TRUE" :
        anndata.write(f'{outpath}/adata_cci_result.h5ad')
        return anndata


    anndata.uns['significant_means_merged']=merge_d
    if save == "TRUE" :
        anndata.write(f'{outpath}/adata_cci_result.h5ad')
        return anndata

    else:
        return anndata
    


if __name__ =="__main__":
    spatial_cci()
    

