#!/usr/bin/env python
# coding: utf-8

#加载
from pathlib import Path
import argparse
import re
import scvelo as scv
import pandas as pd
import scanpy as sc
import numpy as np
import anndata2ri
import os
import difflib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import loompy

scv.set_figure_params()
anndata2ri.activate()
# get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
sc.settings.verbosity = 3




#输入
scvelo = argparse.ArgumentParser(description='RNAvelo')
scvelo.add_argument('--input', type=str, default = None)
scvelo.add_argument('--loom_dir', type=str, default = None)
scvelo.add_argument('--metadata', type=str, default = None)
scvelo.add_argument('--output', type=str, default = None)
scvelo.add_argument('--groupby', type=str, default ='clusters')
scvelo.add_argument('--base', type=str, default ='umap')
scvelo.add_argument('--onlygeneplot', type=str, default ='FALSE')
scvelo.add_argument('--genelist', type=str,default='NULL')
scvelo.add_argument('--order', type=str, default = None)

args = scvelo.parse_args()
# print(args.gpus, type(args.gpus))
# print(args.batch_size, type(args.batch_size))



output=Path(args.output)
if not os.path.exists(output):
    os.makedirs(output)

#颜色
Colorss=["#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b","#666666","#1b9e77","#d95f02","#7570b3",
        "#d01b2a","#43acde","#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff","#ecb9e5","#813139",
        "#743fd2","#434b7e","#e6908e","#214a00","#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
        "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7","#006874","#d2ad2c","#b5d7a5","#9e8442",
        "#4e1737","#e482a7","#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99","#3fdacb","#bf5b17"]

#颜色order
def get_colors(object, groupby, order) :
    if groupby+'_col' in object.obs.columns:
        groupby_col = groupby+'_col'
        temp_df = object.obs[[groupby,groupby_col]].drop_duplicates() 
        groupby_pal = pd.Series(list(temp_df[groupby_col]), index=temp_df[groupby])
        user_color_pal = list(groupby_pal[order[groupby]])
    elif groupby == 'clusters' :
        my_Colorss=pd.Series(Colorss)
        user_color_pal = list(my_Colorss.iloc[[int(i)-1 for i in order[groupby]]])
    else :
        user_color_pal = Colorss
    return user_color_pal


#如果基因不存在于高可变基因中，该基因无法画图
if args.onlygeneplot == 'TRUE':
    adata = scv.read(args.input)
    gene = pd.read_csv(args.genelist, sep="\t", na_values='NA')
    order={}
    order["clusters"]=[int(x) for x in adata.obs["clusters"].sort_values().unique() ]

    if gene.shape[1]<=1:
        gene2 = gene.loc[list(gene['gene'].isin(adata.var_names)), ['gene']] #或者adata.var_names
        for i in range(0, len(gene2['gene']), 10):
            genes = gene2['gene'][i:i + 10]
            scv.pl.scatter(adata,
                           ncols=2,
                           nrows=5,
                           var_names=genes,
                           use_raw=False,
                           ylabel=1,
                           color='clusters',
                           frameon=False,
                           legend_loc='right', palette=get_colors(adata,"clusters",order),
                           save=f'{output}/range{i}to{i + 10}genes.pdf')

            scv.pl.scatter(adata,
                           x='latent_time',
                           y=list(genes),
                           ncols=5,
                           use_raw=False,
                           frameon=False, palette=get_colors(adata,"clusters",order), color='clusters', 
                           legend_loc='right margin',
                           save=f'{output}/range{i}to{i + 10}genes_alonetime_groupby_clusters.pdf')


    else:
        gene2 = gene.loc[list(gene['gene'].isin(adata.var_names)), ['cluster', 'gene']]
        cluster = set(gene2['cluster'].unique())
        gene_dic = {}
        for i in cluster:
            genes = gene2.loc[gene2['cluster'] == i, 'gene']
            gene_dic[i] = genes

        for i in cluster:
            scv.pl.scatter(adata,
                           ncols=2,
                           nrows=5,
                           var_names=gene_dic[i],
                           use_raw=False,
                           ylabel=1,
                           color='clusters',
                           frameon=False,
                           legend_loc='right',palette=get_colors(adata,"clusters",order),
                           save=f'{output}/cluster{i}_top10genes.pdf')
            scv.pl.scatter(adata,
                           x='latent_time',
                           y=list(gene_dic[i]),
                           ncols=5,
                           use_raw=False,
                           frameon=False, palette=get_colors(adata,"clusters",order), color='clusters', 
                           legend_loc='right margin',
                           save=f'{output}/cluster{i}_top10genes_alonetime_groupby_clusters.pdf')

    genelist = set(gene2['gene']) & set(adata.var_names)
    high = len(genelist) / 5
    scv.pl.heatmap(adata, yticklabels=True, sort=True, colorbar=True, show=True, layer="MS",
                   var_names=genelist,
                   sortby='latent_time',
                   col_color='clusters',n_convolve=300,figsize=(8,high),
                   save=f'{output}/inputgene_latent_time_heatmap.pdf')



else:
    #读入loom,做为data
    loom_dir=Path(args.loom_dir)
    loom_file=os.listdir(loom_dir)
    files=[]
    loom_data={}
    index_sample=[]
    for file in loom_file:
        if '.loom' in file :
            print(str(loom_dir)+'/'+file)
            data=scv.read(str(loom_dir)+'/'+file, cache=True)
            data.var_names_make_unique()
            data.obs_names = data.obs_names.str.replace(":","_")
            data.obs_names = data.obs_names.str.replace("x","")
            index_sample=index_sample+list(data.obs_names)
            loom_data[file]=data

    #读取metadata以及后h5d,做为adata
    metadata=pd.read_table(args.metadata,sep='\t')
    order_data=pd.read_table(args.order,sep='\t')
    adata = scv.read(args.input)
    if ('X_pca' in adata.obsm.keys()) == False and 'X_mnn' in adata.obsm.keys():
        adata.obsm['X_pca'] = adata.obsm['X_mnn']
    adata.obs=metadata
    adata.obs_names=adata.obs_names.str.replace("-", "_")
    # adata.X=adata.raw.X #测试
    ldata=loom_data[list(loom_data.keys())[0]]
    for i in range(1,len(loom_data)):
        ldata=ldata.concatenate(loom_data[list(loom_data.keys())[i]])
    ldata.obs_names=index_sample
    
    #merge数据，后续正式分析使用adata2
    adata2 = scv.utils.merge(adata, ldata)

    print('step 1')
    #程序会自动判断这个X是否看起来像标准化后的样子，如果是，则自动选择进行非log
    # Preprocessing requisites consist of gene selection by detection( with a minimum number of counts)
    # and high variability (dispersion), normalizing every cell by its total size and logarithmizing X.
    # Filtering and normalization is applied in the same vein to spliced / unspliced counts and X.
    # Logarithmizing is only applied to X.If X is already preprocessed from former analysis, it will not be touched
    scv.pp.filter_and_normalize(adata2)
    scv.pp.moments(adata2, n_pcs=30, n_neighbors=30)

    print('step 2')
    #proportions画图
    groupby=args.groupby.split(',')
    order={}
    for i in groupby:
        if i == 'clusters':
            adata2.obs[i]=adata2.obs[i].astype('str')
            adata2.obs[i]=adata2.obs[i].astype('category')
            scv.pl.proportions(adata2, groupby=i,save=f'{output}/1.proportions_spliced_unspliced_counts_groupby_{i}.pdf')
            order[i]=[x for x in order_data[i][0].split(",")]
            adata2.obs[i].cat.reorder_categories(order[i],ordered=True,inplace=True)
        else:
            adata2.obs[i]=adata2.obs[i].astype('category')
            order[i]=order_data[i][0].split(",")
            adata2.obs[i].cat.reorder_categories(order[i],ordered=True,inplace=True)
            scv.pl.proportions(adata2, groupby=i,save=f'{output}/1.proportions_spliced_unspliced_counts_groupby_{i}.pdf')
    
    print('step 3')
    ##使用动态模型，比稳态与随机稳态捕捉更多信息与细胞类型（参考nature文章）
    scv.tl.recover_dynamics(adata2)
    scv.tl.velocity(adata2,vkey='dynamical_velocity', mode='dynamical')
    scv.tl.velocity_graph(adata2,vkey='dynamical_velocity')
    
    #dynamical_velocity画图
    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.pl.velocity_embedding_stream(adata2, basis=args.base,
                                         vkey='dynamical_velocity',
                                         title='dynamical_velocity',
                                         color=i,
                                         legend_fontsize=14, legend_loc='right margin',palette=my_palette,
                                         dpi=120,save=f'{output}/2.1.dynamical_velocity_groupby_{i}.svg')

    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.pl.velocity_embedding(adata2, basis=args.base,
                                  vkey='dynamical_velocity',
                                  arrow_length=4,color=i,
                                  arrow_size=2, legend_loc='right margin',palette=my_palette,
                                  dpi=200,save=f'{output}/2.3.single-cell_level_dynamical_velocity_groupby_{i}.svg')

    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.pl.velocity_embedding_grid(adata2, basis=args.base,
                                  vkey='dynamical_velocity',
                                  color=i,
                                  dpi=200,
                                  scale=0.25, legend_loc='right margin',palette=my_palette,
                                  save=f'{output}/2.2.scvelo_embedding_grid_groupby_{i}.pdf')




    #动态模型：基因的似然用于确定反应rates的参数（基因的transcriptional state 和 cell-internal latent time）和潜在细胞特异性变量
    #The dynamical model recovers the latent time of the underlying cellular processes.
    #This latent time represents the cell’s internal clock and approximates the real time experienced by cells as they differentiate,
    #based only on its transcriptional dynamics.

    #潜伏时间(latent_time)仅基于其转录动力学并代表细胞的内部时钟。它比基于相似性的扩散伪时间更好地捕捉实际时间的各个方面
    scv.tl.recover_latent_time(adata2, vkey='dynamical_velocity')
    scv.pl.scatter(adata2,basis=args.base,
                   color='latent_time',
                   fontsize=24,
                   size=100,
                   color_map='gnuplot',
                   perc=[2, 98],
                   colorbar=True,
                   rescale_color=[0,1],palette=Colorss,
                   save=f'{output}/3.1.latent_time_by_dynamical.pdf')


    #Gene expression dynamics resolved along latent time shows a clear cascade of transcription in the top 300 likelihood-ranked genes.
    ##高似然基于，更具latent time排序，看转录状态
    top_genes = adata2.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    df = pd.DataFrame(top_genes)
    df.to_csv(f'{output}/4.3.top_300_likelihood_ranked_genes.csv')
    #heatmap提取的df，列为基因，行为barcodes。画heatmap时，会做转置，行为基因。做scale参数standard_scale 默认为0（行）。如果基因表达值都为0，则整行为空
    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.pl.heatmap(adata2,
                    var_names=top_genes,
                    sortby='latent_time',
                    col_color=i,
                    n_convolve=300,palette=my_palette,
                    save=f'{output}/3.2.latent_time_heatmap.pdf')


    #动态模型中的高似然基因可以认为是驱动基因。这些基因展现出明显的动态特征，并且可以明显的地被检测到
    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.pl.scatter(adata2,
                       basis=top_genes[:15],
                       ncols=5,
                       use_raw=False,
                       frameon=False,
                       color=i,palette=my_palette,legend_loc='right margin',
                       save=f'{output}/4.1.Top-likelihood_genes_groupby_{i}.pdf')

    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.pl.scatter(adata2,
                       x='latent_time',
                       y=top_genes[:15],
                       ncols=5,
                       use_raw=False,
                       frameon=False,palette=my_palette,color=i,legend_loc='right margin',
                       save=f'{output}/4.2.Top-likelihood_genes_alonetime_groupby_{i}.pdf')

    ###根据动态模型去识别簇特异性的驱动基因
    def pl_spatial_dyn_gene(anndata,group):
        scv.tl.rank_dynamical_genes(anndata,groupby=group)
        df=scv.get_df(anndata, 'rank_dynamical_genes/names')
        return df


    for i in groupby:
        my_palette=get_colors(adata2,i,order) 
        output_topgenet=Path(f'{output}/{i}_specific_top_gene')
        if not os.path.exists(output_topgenet):
            os.makedirs(output_topgenet)
        df=pl_spatial_dyn_gene(adata2,group=i)
        for j in df.columns:
            scv.pl.scatter(adata2,
                           df[j][:5],
                           use_raw=False,
                           ylabel= j,
                           color=i,
                           frameon=False,palette=my_palette,
                           save=f'{output_topgenet}/{i}{j}_specific_top-likelihood_genes.pdf')


    ##系统性地识别基因，这些基因可以解释结果的vector field 和谱系分化：t检验簇特异性的差异速度基因
    #The module scv.tl.rank_velocity_genes runs a differential velocity t-test and outpus a gene ranking for each cluster

    for i in groupby:
        my_palette=get_colors(adata2,i,order)
        scv.tl.rank_velocity_genes(adata2, vkey='dynamical_velocity',groupby=i, min_corr=.3)
        df_velocity_genes = scv.DataFrame(adata2.uns['rank_velocity_genes']['names'])
        output_topgenet2=Path(f'{output}/{i}_specific_velocity_genes')
        if not os.path.exists(output_topgenet2):
            os.makedirs(output_topgenet2)
        for j in df_velocity_genes:
            scv.pl.scatter(adata2,
                           df_velocity_genes[j][:5],
                           ylabel=j,
                           color=i,
                           frameon=False,palette=my_palette,
                          save=f'{output_topgenet2}/{i}{j}specific_velocity_genes.pdf')



    #The cell cycle detected by RNA velocity, is biologically affirmed by cell cycle scores (standardized scores of mean expression levels of phase marker genes).
    # scv.tl.score_genes_cell_cycle(adata2)
    # scv.pl.scatter(adata2, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95],save=f'{output}/cellcycle.pdf')

    scv.tl.velocity_confidence(adata2,vkey='dynamical_velocity')

    ##Two more useful stats: - The speed or rate of differentiation is given by the length of the velocity vector.
    ##每一个细胞都有一个长度为g(g个基因)的速度向量。以此向量的长度来替代该细胞分化的速率（速度）。但这是无方向的
    ##The coherence of the vector field (i.e., how a velocity vector correlates with its neighboring velocities) provides a measure of confidence.
    ##细胞速度与邻近细胞的速度是否共表达？（余弦距离）以此来评估细胞与细胞间的相干性。此想干性为置信度（理论上，近邻的细胞同簇的细胞是有相同的干性的）
    keys = 'dynamical_velocity_length', 'dynamical_velocity_confidence'
    scv.pl.scatter(adata2,
                   c=keys,
                   cmap='coolwarm',
                   basis=args.base,
                   perc=[5, 95],palette=Colorss,
                   save=f'{output}/5.Speed_and_coherence.pdf')

    for i in groupby:
        df_v_c = adata2.obs.groupby(i)[keys].mean().T
        df_v_c.style.background_gradient(cmap='coolwarm', axis=1)


    ##paga 画图
    adata2.uns['neighbors']['distances'] = adata2.obsp['distances']
    adata2.uns['neighbors']['connectivities'] = adata2.obsp['connectivities']

    # for i in groupby:
        # my_palette=get_colors(adata2,i,order)
        # scv.tl.paga(adata2, groups=i,vkey='dynamical_velocity')
        # df_page = scv.get_df(adata2, 'paga/transitions_confidence', precision=2).T
        # df_page.style.background_gradient(cmap='Blues').format('{:.2g}')
        # scv.pl.paga(adata2, basis=args.base, size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,legend_loc='right',palette=my_palette,save=f'{output}/paga_velocity-graph.pdf')
        #保存数据
    adata2.__dict__['_raw'].__dict__['_var'] = adata2.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    adata2.write(f'{output}/adata_with_scvelo.h5ad')
    latent_time=adata2.obs[['sampleid']+groupby+['dynamical_velocity_pseudotime','latent_time']]
    latent_time.to_csv(f'{output}/3.3.latent_time.xls',sep="\t")
    


    #输入基因列表画图
    if args.genelist == 'None':
        print("no geneplot")
    else:
        print('plot gene list')
        gene = pd.read_csv(args.genelist, sep="\t", na_values='NA')
        geneplot_output=Path(str(output)+'/'+"geneplot")
        if not os.path.exists(geneplot_output):
            os.makedirs(geneplot_output)
        if gene.shape[1] <= 1:
            gene2 = gene.loc[list(gene['gene'].isin(adata2.var_names)), ['gene']]
            for i in range(0, len(gene2['gene']), 10):
                genes = gene2['gene'][i:i + 10]
                scv.pl.scatter(adata2,
                               ncols=2,
                               nrows=5,
                               var_names=genes,
                               use_raw=False,
                               ylabel=1,
                               color='clusters',
                               frameon=False,
                               legend_loc='right', palette=get_colors(adata2,"clusters",order),
                               save=f'{geneplot_output}/range{i}to{i + 10}genes.pdf')
        else:
            gene2 = gene.loc[list(gene['gene'].isin(adata2.var_names)), ['cluster', 'gene']]
            cluster = set(gene2['cluster'].unique())
            gene_dic = {}
            for i in cluster:
                genes = gene2.loc[gene2['cluster'] == i, 'gene']
                gene_dic[i] = genes
            geneplot_output = Path(str(output) + '/' + "geneplot")
            for i in cluster:
                scv.pl.scatter(adata2,
                               ncols=2,
                               nrows=5,
                               var_names=gene_dic[i],
                               use_raw=False,
                               ylabel=1,
                               color='clusters',
                               frameon=False,
                               legend_loc='right', palette=get_colors(adata2,"clusters",order),
                               save=f'{geneplot_output}/cluster{i}_top10genes.pdf')

        # gene2 = gene.loc[list(gene['gene_name'].isin(adata2.var_names)), ['cluster', 'gene_name']]
        # cluster = set(gene2['cluster'].unique())
        # gene_dic = {}
        # for i in cluster:
        #     genes = gene2.loc[gene2['cluster'] == i, 'gene_name']
        #     gene_dic[i] = genes
        # geneplot_output=Path(str(output)+'/'+"geneplot")
        # if not os.path.exists(geneplot_output):
        #     os.makedirs(geneplot_output)
        # for i in cluster:
        #     scv.pl.scatter(adata2,
        #                    ncols=2,
        #                    nrows=5,
        #                    var_names=gene_dic[i],
        #                    use_raw=False,
        #                    ylabel=1,
        #                    color='clusters',
        #                    frameon=False,
        #                    legend_loc='right',palette=Colorss,
        #                    save=f'{geneplot_output}/cluster{i}_top10genes.pdf')

