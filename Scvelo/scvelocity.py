#!/usr/bin/env python
# coding: utf-8
# author: donghongjie
# time: 2024/04/07
#	email: 979853020@qq.com
import scvelo as scv
import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
import argparse
import re
#import anndata2ri
import os
import difflib
import matplotlib as matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import loompy
import shutil
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


#读取参数
scvelo = argparse.ArgumentParser(description='scvelocity')
scvelo.add_argument('--input', type=str, default = None)
scvelo.add_argument('--loom_dir', type=str, default = None)
scvelo.add_argument('--metadata', type=str, default = None)
scvelo.add_argument('--output', type=str, default = None)
scvelo.add_argument('--groupby', type=str, default ='clusters')
scvelo.add_argument('--basis', type=str, default ='umap')

args = scvelo.parse_args()

output=Path(args.output)
if not os.path.exists(output):
    os.makedirs(output)

Colorss=[
         "#E54B34","#4CBAD4","#009F86","#3B5387","#F29A7F","#8491B3","#91D1C1","#DC0000","#7E6047","#CCCCCC","#BC8B83","#33ADAD","#347988","#9F7685","#C1969A","#8BB0BB","#CE8662","#B04929","#A59487","#E3907E","#D46F5B","#41B4C1","#278C87","#726486","#DA988C","#88A0B7","#B9AC91","#C63517","#927A66","#DBAEA4","#97A4AB","#21A69A","#3A6688","#C98882","#A593A7","#8EC0BE","#D85935","#985738","#B9AFA9","#E67059","#E5BFB9","#B2CED4","#779F99","#747A87","#F2DCD5","#A7ABB3","#C1D1CD","#DCA5A5","#7E7770","#CCCCCC"]

loom_dir=Path(args.loom_dir)
loom_file=os.listdir(loom_dir)


# 初始化 loom 数据字典和样本索引列表
loom_data = {}
index_sample = []

# 遍历 loom 文件列表
for file in loom_file:
    # 检查文件是否为 loom 文件
    if '.loom' in file:
        # 打印文件路径
        print(f"{loom_dir}/{file}")
        
        # 读取 loom 文件数据
        data = sc.read(f"{loom_dir}/{file}", cache=True)
        
        # 确保基因名称唯一
        data.var_names_make_unique()
        
        # 替换观测名称中的特殊字符
        data.obs_names = data.obs_names.str.rstrip('x').str.replace(":","_")
        
        
        # 将样本索引添加到列表中
        index_sample += list(data.obs_names)
        
        # 将数据存储在 loom 数据字典中
        loom_data[file] = data


#metadata=pd.read_table(args.metadata,sep='\t')
adata = scv.read(args.input)

adata.obs_names = [re.sub(r'-\d+$', '', name) for name in adata.obs_names]
processed_names = [name.rsplit('_', 1)[-1] for name in adata.obs_names]
adata.obs_names = adata.obs['sample'] + '_' + processed_names


ldata=loom_data[list(loom_data.keys())[0]]
for i in range(1,len(loom_data)):
		ldata=ldata.concatenate(loom_data[list(loom_data.keys())[i]])

ldata.obs_names=index_sample

adata2 = scv.utils.merge(adata, ldata)



groupby=args.groupby.split(',')
for i in groupby:
		adata2.obs[i]=adata2.obs[i].astype('category')
#数据质量展示
spliced_unspliced_Statistics=Path(f'{output}/1.spliced_unspliced_Statistics')
if not os.path.exists(spliced_unspliced_Statistics):
				os.makedirs(spliced_unspliced_Statistics)

for i in groupby:
		scv.pl.proportions(adata2, groupby=i,save=f'{spliced_unspliced_Statistics}/proportions_spliced_unspliced_by_{i}.pdf')



print('step 1:进行过滤和标准化，如果该数据已经进行过标准化，则跳过！')

scv.pp.filter_and_normalize(adata2)
scv.pp.moments(adata2, n_pcs=30, n_neighbors=30)



print('step2:建立动态模型')
##使用动态模型，比稳态与随机稳态捕捉更多信息与细胞类型
scv.tl.recover_dynamics(adata2,n_jobs=10)
scv.tl.velocity(adata2,vkey='dynamical_velocity', mode='dynamical')
scv.tl.velocity_graph(adata2,vkey='dynamical_velocity',n_jobs=10)


print('step 3：动态模型绘图')
scvelo_stream=Path(f'{output}/2.scvelo_stream')
if not os.path.exists(scvelo_stream):
				os.makedirs(scvelo_stream)


for i in groupby:
		scv.pl.velocity_embedding_stream(adata2,
                                      #basis=args.basis,
																			vkey='dynamical_velocity',
																			title='dynamical_velocity',
																			color=i,
																			legend_fontsize=14, legend_loc='right margin',palette=Colorss,
																			dpi=1000,save=f'{scvelo_stream}/scvelo_stream_by{i}.pdf')


for i in groupby:
		scv.pl.velocity_embedding(adata2, basis=args.basis,
															vkey='dynamical_velocity',
															arrow_length=2,color=i,
															arrow_size=1, legend_loc='right margin',palette=Colorss,
															dpi=1000,save=f'{scvelo_stream}/scvelo_embedding_by{i}.pdf')


for i in groupby:
		scv.pl.velocity_embedding_grid(adata2, basis=args.basis,
															vkey='dynamical_velocity',
															color=i,
															dpi=1000,
															scale=0.25, legend_loc='right margin',palette=Colorss,
															save=f'{scvelo_stream}/scvelo_grid_by{i}.pdf')




##潜伏时间

#The dynamical model recovers the latent time of the underlying cellular processes.
#This latent time represents the cell’s internal clock and approximates the real time experienced by cells as they differentiate,
#based only on its transcriptional dynamics.

#潜伏时间(latent_time)仅基于其转录动力学并代表细胞的内部时钟。它比基于相似性的扩散伪时间更好地捕捉实际时间的各个方面
latent_time=Path(f'{output}/3.latent_time')
if not os.path.exists(latent_time):
				os.makedirs(latent_time)
scv.tl.recover_latent_time(adata2, vkey='dynamical_velocity')
scv.pl.scatter(adata2,
                basis=args.basis,
								color='latent_time',
								fontsize=12,
								size=100,
								color_map='gnuplot',
								perc=[2, 98],
								colorbar=True,
								rescale_color=[0,1])
 #       ,save=f'{latent_time}/latent_time_by_dynamical.pdf')

legend = plt.gca().get_legend()
for ax in plt.gcf().get_axes():
    for item in ax.get_children():
        # 如果子对象是图例对象，则设置其边框属性
        if isinstance(item, matplotlib.legend.Legend):
            item.set_frame_on(False)
plt.tight_layout()            
legend = plt.gca().get_legend()

plt.subplots_adjust(right=0.9)
plt.gcf().set_size_inches(7,5)
plt.savefig(f'{latent_time}/latent_time_by_dynamical.pdf')

#新增一个小提琴图
for i in groupby:
	sc.pl.violin(adata2, keys='latent_time',groupby=i,palette=Colorss,
              rotation=90)
plt.savefig(f'{latent_time}/scVelo-violin-latent_time_by{i}.pdf', dpi=1000, bbox_inches='tight')

#Gene expression dynamics resolved along latent time shows a clear cascade of transcription in the top 300 likelihood-ranked genes.

top_genes = adata2.var['fit_likelihood'].sort_values(ascending=False).index[:300]
#heatmap提取的df，列为基因，行为barcodes。画heatmap时，会做转置，行为基因。做scale参数standard_scale 默认为0（行）。如果基因表达值都为0，则整行为空
for i in groupby:
	scv.pl.heatmap(adata2,
									var_names=top_genes,
									sortby='latent_time',
									col_color=groupby,
									n_convolve=300,
									colorbar=True,
									save=f'{latent_time}/latent_time_heatmap_by{i}.pdf')

#动态模型中的高似然基因可以认为是驱动基因。这些基因展现出明显的动态特征，并且可以明显的地被检测到
likelihood_gene=Path(f'{output}/4.likelihood_gene')
if not os.path.exists(likelihood_gene):
				os.makedirs(likelihood_gene)

for i in groupby:
		scv.pl.scatter(adata2,
										basis=top_genes[:15],
										ncols=5,
										use_raw=False,
										frameon=False,
										color=i,palette=Colorss,legend_loc='right margin',
										save=f'{likelihood_gene}/Top15-likelihood_genes_by_{i}.pdf')

for i in groupby:
		scv.pl.scatter(adata2,
										x='latent_time',
										y=top_genes[:15],
										ncols=5,
										use_raw=False,
										frameon=False,palette=Colorss,color=i,legend_loc='right margin',
										save=f'{likelihood_gene}/Top15-likelihood_genes_latent_time_by_{i}.pdf')
###根据动态模型去识别簇特异性的驱动基因
def Calculate_genes(anndata,group):
		scv.tl.rank_dynamical_genes(anndata,groupby=group, n_genes=10)
		df=scv.get_df(anndata, 'rank_dynamical_genes/names')
		return df

for i in groupby:
		output_specific_gene=Path(f'{likelihood_gene}/Top10_specific_likelihood_gene_by{i}')
		if not os.path.exists(output_specific_gene):
				os.makedirs(output_specific_gene)
		df=Calculate_genes(adata2,group=i)
		for j in df.columns:
				scv.pl.scatter(adata2,
												df[j],
												use_raw=False,
												color=i,
												ncols=5,
												frameon=False,palette=Colorss,
												save=f'{output_specific_gene}/{j}_Top10_specific_likelihood_gene.pdf')


#应用了差异表达检验（Welch t检验，高估方差是保守的），以在一个聚类中找到与所有其他聚类相比表现出不同转录调控动态的基因（例如，该聚类中的诱导和剩余种群中的稳态）,这些基因可以解释结果的vector field 和谱系分化：t检验簇特异性的差异速度基因
for i in groupby:
		scv.tl.rank_velocity_genes(adata2, vkey='dynamical_velocity',groupby=i, min_corr=.3,n_genes=10)
		df_velocity_genes = pd.DataFrame(adata2.uns['rank_velocity_genes']['names'])
		specific_velocity_genes=Path(f'{likelihood_gene}/Top10_specific_velocity_genes_by{i}')
		if not os.path.exists(specific_velocity_genes):
				os.makedirs(specific_velocity_genes)
		for j in df_velocity_genes:
				scv.pl.scatter(adata2,
												df_velocity_genes[j],
												ncols=5,
												color=i,
												frameon=False,palette=Colorss,
											save=f'{specific_velocity_genes}/{j}_Top10_specific_velocity_genes_by{j}.pdf')

scv.tl.velocity_confidence(adata2,vkey='dynamical_velocity')

##Two more useful stats: - The speed or rate of differentiation is given by the length of the velocity vector.
##每一个细胞都有一个长度为g(g个基因)的速度向量。以此向量的长度来替代该细胞分化的速率（速度）。但这是无方向的
##The coherence of the vector field (i.e., how a velocity vector correlates with its neighboring velocities) provides a measure of confidence.
##细胞速度与邻近细胞的速度是否共表达？（余弦距离）以此来评估细胞与细胞间的相干性。此想干性为置信度（理论上，近邻的细胞同簇的细胞是有相同的干性的）
#keys = 'dynamical_velocity_length', 'dynamical_velocity_confidence'
#TODO:改到这里啦
speed_and_coherence=Path(f'{output}/5.speed_and_coherence')
if not os.path.exists(speed_and_coherence):
		os.makedirs(speed_and_coherence)
scv.pl.scatter(adata2,
								color='dynamical_velocity_confidence',
								cmap='coolwarm',
								basis=args.basis,legend_align_text="x",
								perc=[5, 95],palette=Colorss
                )
legend = plt.gca().get_legend()
for ax in plt.gcf().get_axes():
    for item in ax.get_children():
        # 如果子对象是图例对象，则设置其边框属性
        if isinstance(item, matplotlib.legend.Legend):
            item.set_frame_on(False)
plt.tight_layout()            
legend = plt.gca().get_legend()

plt.subplots_adjust(right=0.9)
plt.gcf().set_size_inches(7,5)

plt.savefig(f'{speed_and_coherence}/dynamical_velocity_confidence.pdf')




scv.pl.scatter(adata2,
								color='dynamical_velocity_length',
								cmap='coolwarm',legend_align_text="x",
								basis=args.basis,
								perc=[5, 95],palette=Colorss
								)

legend = plt.gca().get_legend()
for ax in plt.gcf().get_axes():
    for item in ax.get_children():
        # 如果子对象是图例对象，则设置其边框属性
        if isinstance(item, matplotlib.legend.Legend):
            item.set_frame_on(False)
plt.tight_layout()            
legend = plt.gca().get_legend()

plt.subplots_adjust(right=0.9)
plt.gcf().set_size_inches(7,5)

plt.savefig(f'{speed_and_coherence}/dynamical_velocity_length.pdf')


adata2.obs.columns = adata2.obs.columns.str.replace('sample','sampleid')
adata2.write(f'{output}/adata_with_scvelo.h5ad')
latent_time=adata2.obs[['sampleid','group']+groupby+['dynamical_velocity_pseudotime','latent_time','dynamical_velocity_length','dynamical_velocity_confidence']]

latent_time.to_csv(f'{output}/velocity_data.xls',sep="\t")

shutil.copy('/PERSONALBIO/work/singlecell/s04/Test/donghongjie/PSN_singlecell/Scvelo/Scvelocity 结果说明.docx',f'{output}')

# ##PAGA 分析
# for i in groupby:
#   sc.tl.paga(adata2, groups=i)
#   sc.pl.paga(adata2, save='paga_plot.png')
  