import os
from oebio.report import Report
from oebio.app import *
from configparser import ConfigParser

#定义函数批量获取路径
def get_url(target,tag):
    urls = []
    for f in Path(target).glob(tag):
        text = f'[{f}]({f})\n\n'
        urls.append(text)
    return ''.join(urls)

##########################读取项目配置文件信息##########################
config = ConfigParser()
config.read('config.txt',encoding="utf-8")

# Define header info
report = Report('10X单细胞转录组测序分析结题报告', title = '10X单细胞转录组测序分析结题报告', header_info = dict(config['par']), oe_welcome='''
####  感谢您选择欧易生物！

**关于网页报告使用有以下几点提示：**

1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节，便于您快速查看相应内容。
2. **报告正文**中的**图片**均可以点击后进行**放大查看**，且图片点击后可以包含更多的细节信息，比如左上角会显示具体的**图片数量**，右上角有相关的**图片工具**，**键盘上的左
右方向键/鼠标滚轮**可以对图片进行**切换**。
3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张，更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**，当表格多于15行的时候，表格会嵌入到网页之中，可以利用滚动栏进行调整设置，同时会在表格上方增加搜索设置，便于快速
查询表格信息。
5. 在报告中浏览，需要返回顶部的时候，可以使用网页**右下角的白色箭头标识**快速返回顶部。
6. 本提示可以点击右上方的不再显示，则后续打开网页均不会显示本提示。
''')
report.add_yaml_config('/home/luyao/oebio_report/report.yaml')
os.environ['OEBIO'] = "/home/luyao/oebio_report/pic"


##############################项目概况##############################
project_module = report.add_section('项目概况')

##########################项目摘要##########################
summ = pd.read_csv('2.Count_QC/cell_statitics_before_after_QC.xls',index_col=0,sep='\t')
num = len(summ)   #样本数量
min_beforeQC_cell = min(summ.Total_cells_beforeQC)
max_beforeQC_cell = max(summ.Total_cells_beforeQC)
min_afterQC_cell = min(summ.Total_cells_afterQC)
max_afterQC_cell = max(summ.Total_cells_afterQC)
min_afterQC_umi = format(min(summ.Mean_nUMI_afterQC.astype(float)),'.2f')
max_afterQC_umi = format(max(summ.Mean_nUMI_afterQC.astype(float)),'.2f')
min_afterQC_gene = format(min(summ.Mean_nGene_afterQC.astype(float)),'.2f')
max_afterQC_gene = format(max(summ.Mean_nGene_afterQC.astype(float)),'.2f')
min_afterQC_mito = format(min(summ['Mean_mito.percent_afterQC'].astype(float)),'.4f')
max_afterQC_mito = format(max(summ['Mean_mito.percent_afterQC'].astype(float)),'.4f')
celltype_file =  list(glob('5.Celltyping/*_simplified_celltyping_results.xls'))
cell_results = pd.read_csv(celltype_file[0], sep='\t')    ####
cluster_num = max(cell_results['clusters'])         ####
celltype_tmp = [ str(i) for i in set(cell_results['celltype']) ]      ####
celltype = ', '.join(celltype_tmp)

# diff = pd.read_csv('6.Diffexp/diffexp_results_stat.xls',sep='\t')
# num_diff = len(diff)
# diff_num_tmp = [ str(i) for i in diff.iloc[:,4]]  ### diff_num_tmp = [ str(i) for i in diff.iloc[:,4].to_list()]
# if len(diff_num_tmp)>0:
	# diff_num =diff_num_tmp[0]
# else:
	# diff_num = ','.join(diff_num_tmp)

abstract_module = project_module.add_section('项目摘要',num=num,min_beforeQC_cell=min_beforeQC_cell,max_beforeQC_cell=max_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,max_afterQC_cell=max_afterQC_cell,min_afterQC_umi=min_afterQC_umi,max_afterQC_umi=max_afterQC_umi,min_afterQC_gene=min_afterQC_gene,max_afterQC_gene=max_afterQC_gene,min_afterQC_mito=min_afterQC_mito,max_afterQC_mito=max_afterQC_mito,cluster_num=cluster_num,celltype=celltype)

# abstract_module = project_module.add_section('项目摘要',num=num,min_beforeQC_cell=min_beforeQC_cell,max_beforeQC_cell=max_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,max_afterQC_cell=max_afterQC_cell,min_afterQC_umi=min_afterQC_umi,max_afterQC_umi=max_afterQC_umi,min_afterQC_gene=min_afterQC_gene,max_afterQC_gene=max_afterQC_gene,min_afterQC_mito=min_afterQC_mito,max_afterQC_mito=max_afterQC_mito,cluster_num=cluster_num,celltype=celltype,num_diff=num_diff,diff_num=diff_num)

##########################项目快捷链接##########################
link_refdiff = 'supplemental_material/单细胞转录组常见问题(FAQ).pdf'
linkage = project_module.add_section('项目快捷链接',link_refdiff=link_refdiff)


##########################建库测序流程##########################
library_module = report.add_section('实验测序流程')
library_module_c = library_module.add_comment('实验建库测序流程图如下：')
library_module_plot = library_module.add_plot('Experiment.png',caption ='建库测序流程图',content='建库测序流程图')

##########################生物信息分析流程##########################
analysis_module = report.add_section('生物信息分析流程')
analysis_module_c = analysis_module.add_comment('生物信息分析流程图如下：')
analysis_module_plot = analysis_module.add_plot('Pipeline.png', caption = '生物信息分析流程',content='生物信息分析流程')

##########################项目分析结果##########################
result_module = report.add_section('项目分析结果')

# ##########################测序数据质量评估##########################
# fastqc_module = result_module.add_section('测序数据质量评估')

# fastqc1_module = fastqc_module.add_section('原始测序数据')
# fastqc1_module_c1 = fastqc1_module.add_comment('高通量测序得到的原始图像数据文件经碱基识别（Base Calling）分析转化为原始测序序列，我们称之为 Raw Data 或 Raw Reads，结果以 FASTQ（简称为fq）文件格式存储，其中包含测序序列（reads）的序列信息以及其对应的测序质量信息。FASTQ 格式文件中每个 read 由四行描述，如下图所示：')

# fastqc1_plot1 = fastqc1_module.add_plot('FASTQ.png', caption = '原始数据示例',content='原始数据 FASTQ 格式示例')
# fastqc1_module_c2 = fastqc1_module.add_comment("每个序列共有 4 行，第 1 行和第 3 行是序列名称（有的 fq 文件为了节省存储空间会省略第三行“＋”后面的序列名称），由测序仪产生；第 2 行是序列；第 4 行是序列的测序质量，每个字符对应第 2 行每个碱基。第 4 行中每个字符对应的 ASCII 值减去 33，即为对应第二行碱基的测序质量值。如果测序错误率用 E 表示，Illumina 的碱基质量值用 Qphred 表示，对应关系如下所示 ：")
# fastqc1_plot2 = fastqc1_module.add_plot('Qphred.png', caption = '碱基质量值换算公式', content='碱基质量值换算公式')
# quality_table = fastqc1_module.add_table('Illumina_quality.txt', caption = 'Illumina 测序错误率与测序质量值简明对应关系表')

# ##########################碱基质量分析##########################
# fastqc2_module = fastqc_module.add_section('测序碱基质量分析')

# fastqc2_plots = sorted(list(glob('1.Fastqc/*/Images/per_base_quality.png')),reverse=False)
# fastqc2_contents = ['测序碱基质量分布图。图片说明：横坐标为碱基在 reads 中的位置，纵坐标为碱基质量值，红线为质量值中位数，蓝线为质量值平均值，箱线的上下线为质量值的极值。由于 R1 端由 Barcode+UMI（共28pbs）和转录本3’端 polyT 组成，故在28bps之后碱基质量下降明显，一般只使用前 28bps 进行下游分析。']*len(fastqc2_plots)
# fastqc2_plots = list(zip(fastqc2_plots, fastqc2_contents))
# fastqc2_module_plot = fastqc2_module.add_plot(fastqc2_plots, caption = '测序碱基质量分布图',description = '图片说明：左右两张图片分别为reads的R1和R2碱基质量分布图，横坐标为碱基在 reads 中的位置，纵坐标为碱基质量值，红线为质量值中位数，蓝线为质量值平均值，箱线的上下线为质量值的极值。由于 R1 端由 Barcode+UMI（共28pbs）和转录本3’端 polyT 组成，故在28bps之后碱基质量下降明显，一般只使用前 28bps 进行下游分析。' )

# ##########################碱基含量分布检查##########################
# fastqc3_module = fastqc_module.add_section('碱基含量分布检查')
# fastqc3_plots = sorted(list(glob('1.Fastqc/*/Images/per_base_sequence_content.png')),reverse=False)
# fastqc3_contents = ['测序碱基含量分布图。横坐标是 reads 碱基坐标，纵坐标是所有reads在该测序位置 A、C、G、T 碱基分别占的百分比，不同碱基用不同颜色表示。序列的起始位置与测序的引物接头相连，因此 A、C、G、T 在起始端会有所波动，后面会趋于稳定。R1 端由于在 UMI 序列后存在用于捕获转录本的 polyT，因此在28bp后 T 碱基的比例会相应升高。']*len(fastqc3_plots)
# fastqc3_plots = list(zip(fastqc3_plots, fastqc3_contents))
# fastqc3_module_plot= fastqc3_module.add_plot(fastqc3_plots, caption = '测序碱基含量分布图', description = '图片说明：横坐标是 reads 碱基坐标，纵坐标是所有reads在该测序位置 A、C、G、T 碱基分别占的百分比，不同碱基用不同颜色表示。序列的起始位置与测序的引物接头相连，因此 A、C、G、T 在起始端会有所波动，后面会趋于稳定。R1 端由于在 UMI 序列后存在用于捕获转录本的 polyT，因此在 28bp 后 T 碱基的比例会相应升高。')

##########################基因定量质控##########################
QC_1 = result_module.add_section('基因定量质控')
QC_1_plots = sorted(list(glob('1.Cellranger/*summary.png')),reverse=False)
QC_1_contents = ['样本细胞质量统计结果']*len(QC_1_plots)
QC_1_plots = list(zip(QC_1_plots, QC_1_contents))
QC1_plot = QC_1.add_plot(QC_1_plots, caption = '样本细胞质量统计结果')
QC1_table = QC_1.add_table('scRNA_qc_stat.txt', caption = '样本细胞质量统计结果说明')

# ## 测序饱和度分析
# QC_2 = QC_module.add_section('测序饱和度分析')
# QC2_plot = QC_2.add_plot('1.Cellranger/*rarefraction_curve.png', caption = '样本测序饱和度分析结果' , description = '图片说明：左图为测序饱和度分布图，对 reads 抽样，观察不同抽样条件下检测到的转录本数量占检测到所有转录本的比例，并绘制曲线。如果曲线末端区域接近平缓，则表明增大测序数据量与检测到更多的转录本之间的获益不明显，说明测序接近饱和。右图为每个细胞内的基因中位数分布，随着 reads 数目增加，检测到每个细胞内的基因中位值不发生显著增加，则测序深度已达饱和。' )


##########################基因定量质控##########################
Quantification_module = result_module.add_section('定量后质控')

Quantification1_module = Quantification_module.add_section('过滤低质量细胞')
Quantification1_plots = list(glob('2.Count_QC/outliers.png'))
Quantification1_contents = ['拟合广义线性模型过滤离域细胞图。图片说明：横轴为每个细胞内的 UMI 数，纵轴为每个细胞内的基因数，根据二者的线性关系拟合分布模型，着色的点表示离域细胞，在下游分析中会进行剔除。']
Quantification1_plots = list(zip(Quantification1_plots, Quantification1_contents))

Quantification1_plot = Quantification1_module.add_plot(Quantification1_plots, caption = '拟合广义线性模型过滤离域细胞图',description = '图片说明：横轴为每个细胞内的 UMI 数，纵轴为每个细胞内的基因数，根据二者的线性关系拟合分布模型，着色的点表示离域细胞，在下游分析中会剔除。')

mito_plots = sorted(list(glob('2.Count_QC/*QC_mitochondroin_transcript_ratio_in_each_cell_violin_plot.png')),reverse=True)
mito_contents = ['质控前每个细胞中线粒体基因比例的小提琴分布图。图片说明：纵轴表示线粒体基因占单个细胞所有基因的比例（%），图中每个点代表一个油包水微滴中细胞的线粒体基因比例，每个小提琴图反映对应样本中所有细胞的线粒体基因在细胞所有基因中所占的比例，一般要求大部分细胞的线粒体基因比例越低越好（特殊样本除外）。', '质控后每个细胞中线粒体基因比例的小提琴分布图。图片说明：纵轴表示线粒体基因占单个细胞所有基因的比例（%），图中每个点代表一个油包水微滴中细胞的线粒体基因比例，每个小提琴图反映对应样本中所有细胞的线粒体基因在细胞所有基因中所占的比例，一般要求大部分细胞的线粒体基因比例越低越好（特殊样本除外）。']
mito_plots = list(zip(mito_plots, mito_contents))
mito_plot = Quantification1_module.add_plot(mito_plots, caption = '质控前后每个细胞中线粒体基因比例的小提琴分布图' , description = '图片说明：纵轴表示线粒体基因占单个细胞所有基因的比例（%），图中每个点代表一个油包水微滴中细胞的线粒体基因比例，每个小提琴图反映对应样本中所有细胞的线粒体基因在细胞所有基因中所占的比例，一般要求大部分细胞的线粒体基因比例越低越好（特殊样本除外）。左图为质控前，右图为质控后。')

ngene_plots = sorted(list(glob('2.Count_QC/*QC_total_genes4each_cell_on_violin_plot.png')),reverse=True)
ngene_contents = ['质控前每个细胞中基因数目的分布图。图片说明：纵轴表示细胞中有表达的基因数目，图中每个点代表一个油包水微滴中细胞的基因数目。该图反映样本中的每一个细胞表达基因的数目，基因数目异常过多的点很有可能是由于对应的油包水微滴中包含多个细胞，需要根据实际需求设置合理的阈值过滤掉。左图为质控前，右图为质控后。', '质控后每个细胞中基因数目的分布图。图片说明：纵轴表示细胞中有表达的基因数目，图中每个点代表一个油包水微滴中细胞的基因数目。该图反映样本中的每一个细胞表达基因的数目，基因数目异常过多的点很有可能是由于对应的油包水微滴中包含多个细胞，需要根据实际需求设置合理的阈值过滤掉。左图为质控前，右图为质控后。']
ngene_plots = list(zip(ngene_plots, ngene_contents))
ngene_plot = Quantification1_module.add_plot(ngene_plots, caption = '质控前后每个细胞中基因表达数目的小提琴分布图', description = '图片说明：纵轴表示细胞中有表达的基因数目，图中每个点代表一个油包水微滴中细胞的基因数目。该图反映样本中的每一个细胞表达基因的数目，基因数目异常过多的点很有可能是由于对应的油包水微滴中包含多个细胞，需要根据实际需求设置合理的阈值过滤掉。左图为质控前，右图为质控后。')

nUMI_plots = sorted(list(glob('2.Count_QC/*QC_total_UMIs4each_cell_on_violin_plot.png')),reverse=True)
nUMI_contents = ['质控前每个细胞中 UMI 数目的分布图。图片说明：纵轴表示 UMI 数，图中的每个点代表一个油包水微滴中细胞的 UMI 数目，即转录本的数目，该图反映样本中每一个细胞的转录本数目，转录本数目异常过多的细胞需要通过设置合理的阈值将其过滤。', '质控后每个细胞中 UMI 数目的分布图。图片说明：纵轴表示 UMI 数，图中的每个点代表一个油包水微滴中细胞的 UMI 数目，即转录本的数目，该图反映样本中每一个细胞的转录本数目，转录本数目异常过多的细胞需要通过设置合理的阈值将其过滤。']
nUMI_plots = list(zip(nUMI_plots, nUMI_contents))
nUMI_plot = Quantification1_module.add_plot(nUMI_plots, caption = '质控前后每个细胞中 UMI 数目的小提琴分布图', description = '图片说明：纵轴表示 UMI 数，图中的每个点代表一个油包水微滴中细胞的 UMI 数目，即转录本的数目，该图反映样本中每一个细胞的转录本数目，转录本数目异常过多的细胞需要通过设置合理的阈值将其过滤。左图为质控前，右图为质控后。')

Quantification1_c = Quantification1_module.add_comment(" ")
Quantification1_c2 = Quantification1_module.add_comment("定量质控前后细胞统计情况如下表所示：")
Quantification1_table = Quantification1_module.add_table('2.Count_QC/cell_statitics_before_after_QC.xls', caption = '质控前后细胞数目统计表', headerdict={'sampleid':'样本名称','Mean_nUMI_beforeQC':'质控前细胞中的平均UMI数','Mean_nGene_beforeQC':'质控前细胞中的平均基因数','Mean_mito.percent_beforeQC':'质控前细胞中的平均线粒体比例','Total_cells_beforeQC':'质控前的总细胞数','Mean_nUMI_afterQC':'质控后细胞中的平均UMI数','Mean_nGene_afterQC':'质控后细胞中的平均基因数','Mean_mito.percent_afterQC':'质控后细胞中的平均线粒体比例','Total_cells_afterQC':'质控后的总细胞数'})


##########################降维与聚类分析##########################
Clustering_module = result_module.add_section('降维与聚类分析')
tsne_module = Clustering_module.add_section('t-SNE 降维聚类分析')
tsne_plots = list(glob('3.Clustering/tsne_Dimension_Reduction/tsne_groupby_cluster_resolution*_plot.png'))
tsne_contents = ['t-SNE 降维聚类结果图。图片说明：横纵坐标分别代表t-SNE降维第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。']
tsne_plots = list(zip(tsne_plots, tsne_contents))

tsne_plot = tsne_module.add_plot(tsne_plots, caption = 't-SNE 降维聚类结果图', description ='图片说明：横纵坐标分别代表t-SNE降维第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')

# Clustering_link = get_url('3.Clustering/tSNE_Dimension_Reduction/', 'tsne_resolution*_output.xls')   ####
# tsne_table = tsne_module.add_comment('''t-SNE 二维降维聚类坐标表格：

# {Clustering_link}''',Clustering_link=Clustering_link)


##########################Marker基因鉴定##########################
Marker_module = result_module.add_section('Marker基因鉴定')

Marker_table = Marker_module.add_table('4.Marker/top10_for_each_clsuters.xls', caption = '每个细胞群 Top10 Marker 基因列表', show_rows=10 , headerdict={'gene':'基因名','p_val':'p值','avg_logFC':'对数转化后的平均Foldchange值','pct.1':'表达Marker基因的细胞在当前群中的占比','pct.2':'表达Marker基因的细胞在其余群中的占比','p_val_adj':'校正后的p值','cluster':'细胞群编号','gene_diff':'pct.1与pct.2的比值。Top10 Marker以gene_diff列为筛选依据。'})

Marker_link1 = '4.Marker/all_DEGs_for_all_clsuters.xls'
Marker_link2 = '4.Marker/top10_for_each_clsuters.xls'
Marker_table1 = Marker_module.add_comment('''每个细胞群中所有 Marker 基因结果表格：

[{Marker_link1}]({Marker_link1})
	''',Marker_link1=Marker_link1)
Marker_table2 = Marker_module.add_comment('''每个细胞群中 Top10 Marker 基因结果表格：
	
[{Marker_link2}]({Marker_link2})''',Marker_link2=Marker_link2)

heatmap_plots = list(glob('4.Marker/topmarker_gene_heatmap.png'))
heatmap_contents = ['Top10 Marker 基因表达热图。图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。']
heatmap_plots = list(zip(heatmap_plots, heatmap_contents))


heatmap_plot = Marker_module.add_plot(heatmap_plots, caption = 'Top10 Marker 基因表达热图', description = '图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')

Marker_c = Marker_module.add_comment('每个细胞群中 Top10 Marker 基因的可视化展示如下：')
feature_plots = ["4.Marker/markers_vis4cluster"+str(i+1)+"/topmarker_gene_featureplot.png" for i in range(cluster_num) ]
#### feature_plots = sorted(list(glob('4.Marker/markers_vis4cluster*/*marker_gene_featureplot.png')))
feature_contents = ["第 "+str(i+1)+" 群细胞的 Top10 Marker 基因在 t-SNE 聚类结果中的可视化图。图片说明：红色越深表示该细胞中对应基因的表达量越高。" for i in range(len(feature_plots)) ]    ### str(i+1)
feature_plots = list(zip(feature_plots, feature_contents))
feature_plot = Marker_module.add_plot(feature_plots, caption = 'Top10 Marker 基因在 t-SNE 聚类结果中的可视化图', description = '图片说明：红色越深表示该细胞中对应基因的表达量越高。')
vln_plots = ["4.Marker/markers_vis4cluster"+str(i+1)+"/topmarker_gene_violin_plot.png" for i in range(cluster_num) ]
#### vln_plots = sorted(list(glob('4.Marker/markers_vis4cluster'+i+'/*marker_gene_violin_plot.png' for i in range(12)  )))
vln_contents = ["第 "+str(i+1)+" 群细胞的 Top10 Marker基因表达量小提琴图。图片说明：横坐标为细胞群的编号，纵坐标为标准化后的基因表达值，图中每个点表示每个细胞中 Marker 基因的表达值。" for i in range(len(vln_plots)) ]    ### str(i+1)
vln_plots = list(zip(vln_plots, vln_contents))
vln_plot = Marker_module.add_plot(vln_plots, caption = 'Top10 Marker 基因表达量小提琴图', description = '图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值，图中每个点表示每个细胞中 Marker 基因的表达值。')

##########################细胞类型鉴定##########################
Celltype_module = result_module.add_section('细胞类型鉴定')


Celltype_link = list(glob('5.Celltyping/*top.*celltyping_on_tsne*.pdf'))[0] 
Celltype_tsne = Celltype_module.add_comment('''细胞类型鉴定图结果见：

[{Celltype_link}]({Celltype_link})
	''',Celltype_link=Celltype_link)

Celltype_plots = list(glob('5.Celltyping/*top.*celltyping_on_tsne*.png'))
Celltype_contents = ['细胞类型鉴定图。图片说明：细胞类型注释结果在 t-SNE 图上的展示，每种细胞类型以不同颜色区分。']
Celltype_plots = list(zip(Celltype_plots, Celltype_contents))
Celltype_plot = Celltype_module.add_plot(Celltype_plots, caption = '细胞类型鉴定图',description = '图片说明：细胞类型注释结果在 t-SNE 图上的展示，每种细胞类型以不同颜色区分。')

Celltype_heatmap_link = list(glob('5.Celltyping/*celltyping_heatmap.pdf'))[0]
Celltype_heatmap = Celltype_module.add_comment('''细胞类型鉴定相关性热图结果见：

[{Celltype_heatmap_link}]({Celltype_heatmap_link})
	''',Celltype_heatmap_link=Celltype_heatmap_link)

Celltype_heatmap_plots = list(glob('5.Celltyping/*celltyping_heatmap.png'))
Celltype_heatmap_contents = ['细胞类型鉴定相关性热图。图片说明：纵轴表示每一个待鉴定的细胞，横轴表示参考数据集中的细胞类型注释名称。颜色越红代表相关性值越大，表明待鉴定的细胞类型与参考数据集中的该细胞类型注释最为相似。']
Celltype_heatmap_plots = list(zip(Celltype_heatmap_plots, Celltype_heatmap_contents))
Celltype_heatmap_plot = Celltype_module.add_plot(Celltype_heatmap_plots, caption = '细胞类型鉴定相关性热图', description = '图片说明：纵轴表示每一个待鉴定的细胞，横轴表示参考数据集中的细胞类型注释名称。颜色越红表示相关性值越大，表明待鉴定的细胞类型与参考数据集中的该种细胞类型最为相似。')

Celltype_table_link = list(glob('5.Celltyping/*_simplified_celltyping_results.xls'))[0]
Celltype_table = Celltype_module.add_comment('''细胞类型注释表格见：

[{Celltype_table_link}]({Celltype_table_link})
	''',Celltype_table_link=Celltype_table_link)
Celltype_table = Celltype_module.add_table('5.Celltyping/*_simplified_celltyping_results.xls', caption = '细胞类型注释结果表格', show_rows=10 , headerdict={'cell_barcode':'细胞Barcode','sampleid':'样本名称','celltype':'细胞类型注释','clusters':'细胞群'})

celltype_link1 = get_url('5.Celltyping/', '*_celltyping_statistics.xls')
celltype_c = Celltype_module.add_comment('''各细胞群中的细胞类型数目统计表格见：
	
{celltype_link1} ''',celltype_link1=celltype_link1)


# ##########################差异表达基因及富集##########################
# Diff_module = result_module.add_section('差异表达基因筛选')
# Diffexp_t = Diff_module.add_comment("各分组差异基因数目统计表如下：")
# Diffexp_table = Diff_module.add_table('6.Diffexp/diffexp_results_stat.xls', caption = '差异表达基因统计表',  headerdict={'Case':'实验组名','Control':'对照组名','Up_diff':'显著上调差异基因数量','Down_diff':'显著下调差异基因数量','Total_diff':'显著性差异基因总数量'})

# Diffexp_c = Diff_module.add_comment("差异显著基因表格结果示例：")
# difffiles = list(glob('6.Diffexp/*diff*pval*xls'))
# Diffexp_table = Diff_module.add_table(difffiles[0], caption = '差异表达基因结果文件',show_rows =10 , headerdict={'GeneID':'基因名','pvalue':'p值','pct.1':'表达该基因的细胞在比较组细胞群中的占比','pct.2':'表达该基因的细胞在对照组细胞群中的占比','padj':'校正后的p值','FoldChange':'差异倍数','log2FoldChange':'log2转化后的差异倍数','up_down':'基因上下调描述'})


# ##########################差异基因 GO 富集分析##########################

# enrich_module = result_module.add_section('差异基因富集分析')



# #定义函数获取比较组名称
# def get_diffgroup_name(target,tag):
    # diff_group = []
    # for f in Path(target).glob(tag):
        # t = str(f).split('/')[-1]
        # t = t.split('-diff')[0]
        # diff_group.append(t)
    # return diff_group

# link = '7.enrichment/GO_enrichment'
# diff_group = get_diffgroup_name('6.Diffexp', '*diff*FC*.xls')

# diff_gene_GO = enrich_module.add_section('差异基因 GO 富集分析')
# diff_gene_GO_plot2 = diff_gene_GO.add_plot('enrich.png',caption='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式',content='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式。其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。',description='其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。')
# diff_gene_GO_c = diff_gene_GO.add_comment('''输出文件结果目录：
# [{link}]({link})

# GO 富集性分析结果示例：''',link=link)
# diff_gene_GO_table1 = diff_gene_GO.add_table('7.enrichment/GO_enrichment/%s/enrichment-go-*-Total.xls' %(diff_group[0]),caption='富集分析结果表格',show_rows=10,headerdict={'id':'条目在Gene Ontology的登录号','Term':'该条目的描述','Gene':'属于该条目的差异gene'})
# diff_gene_GO_plot3_c = diff_gene_GO.add_comment('''GO 富集分析 top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log<sub>10</sub>pValue 由大到小排序的各 10 条）条形图展示如下：''')

# diff_gene_GO3_plots = list(glob('7.enrichment/GO_enrichment/*/GO.top.*.png'))
# diff_gene_GO3_contents = ['GO 富集分析结果展示。图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>pValue。']*len(diff_gene_GO3_plots)
# diff_gene_GO3_plots = list(zip(diff_gene_GO3_plots, diff_gene_GO3_contents))
# diff_gene_GO_plot3 = diff_gene_GO.add_plot(diff_gene_GO3_plots,caption='GO 富集分析结果展示',description='图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>pValue。')
# diff_gene_GO_fish = diff_gene_GO.add_comment('''使用 fisher 算法分别对样本间差异基因进行 CC，BP，MF 富集分析，并使用 topGO[^topGO] 对富集到的 Term 绘制有向无环图。topGO 有向无环图能直观展示差异表达基因富集的 GO 节点（Term）及其层级关系，是差异表达基因 GO 富集分析的结果图形化展示，分支代表的包含关系，从上至下所定义的功能描述范围越来越具体。
# [^topGO]:Alexa A, Rahnenfuhrer J. topGO: enrichment analysis for gene ontology. R package version 2.8,2010.''')
# diff_gene_GO_plot4 = diff_gene_GO.add_plot('topGO.png',caption='差异基因topGO有向无环示例图展示',content='差异基因topGO有向无环示例图展示。图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。',description='''图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。''')
# diff_gene_GO_plot5_c = diff_gene_GO.add_comment('''根据功能分级，一般将 GO 分为三个层级，level1 包含三个条目：biological process、cellular component和molecular function，level2 包含 biological adhesion、cell 和 binding 等 64 个条目，level3 即为常规富集使用的数万个条目。从 level1到 level3 功能更具体，反之，更概括。

# 差异基因和所有基因在 GO Level2 水平分布比较图如下：''')

# diff_gene_GO5_plots = list(glob('7.enrichment/GO_enrichment/*/ALL_vs_DEG.GO.level2.stat.png'))
# diff_gene_GO5_contents = ['差异表达基因及所有基因在 GO Level2 水平分布比较图。图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。']*len(diff_gene_GO5_plots)
# diff_gene_GO5_plots = list(zip(diff_gene_GO5_plots, diff_gene_GO5_contents))
# diff_gene_GO_plot5 = diff_gene_GO.add_plot(diff_gene_GO5_plots,caption='差异表达基因及所有基因在 GO Level2 水平分布比较图',description='图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')
# diff_gene_GO_plot6_c = diff_gene_GO.add_comment('''上调差异基因和下调差异基因在 GO Level2 水平分布比较图如下：''')
# diff_gene_GO6_plots = list(glob('7.enrichment/GO_enrichment/*/Up_vs_Down.GO.level2.stat.png'))
# diff_gene_GO6_contents = ['上调差异基因和下调差异基因在 GO Level2 水平分布比较图。图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比']*len(diff_gene_GO6_plots)
# diff_gene_GO6_plots = list(zip(diff_gene_GO6_plots, diff_gene_GO5_contents))
# diff_gene_GO_plot6 = diff_gene_GO.add_plot(diff_gene_GO6_plots,caption='上调差异基因和下调差异基因在 GO Level2 水平分布比较图',description='图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')



# ##########################差异基因 KEGG 富集分析##########################
# link1 = '7.enrichment/KEGG_enrichment'
# diff_gene_KEGG = enrich_module.add_section('差异基因 KEGG 富集分析')
# diff_gene_KEGG_c = diff_gene_KEGG.add_comment('''差异基因 KEGG 富集分析结果：[{link1}]({link1})

# KEGG 富集分析 top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序）气泡图如下：''',link1=link1)

# diff_gene_KEGG_top20_plots = list(glob('7.enrichment/KEGG_enrichment/*/KEGG.top.*.png'))
# diff_gene_KEGG_top20_contents = ['KEGG富集 top20 气泡图。图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。']*len(diff_gene_KEGG_top20_plots)
# diff_gene_KEGG_top20_plots = list(zip(diff_gene_KEGG_top20_plots, diff_gene_KEGG_top20_contents))
# diff_gene_KEGG_top20 = diff_gene_KEGG.add_plot(diff_gene_KEGG_top20_plots,caption='KEGG富集 top20 气泡图',description='图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。')
# diff_gene_KEGG_levle2_c = diff_gene_KEGG.add_comment('''根据功能分级，通常将 KEGG 分为三个层级，level1 包含六个分类：Metabolism、Genetic Information Processing、Environmental Information Processing、Cellular Processes、Organismal Systems 和 Human Diseases（具体物种注释可能有删减）。level2 包含 Cell growth and death、Transcription 和 Development 等 44 个分类（具体物种注释可能有删减），level3 即为常规富集使用的数百个 Pathway，从 level1 到 level3 功能更具体，反之，更概括。

# 差异表达基因及所有基因在 KEGG Level2 水平分布比较图如下：''')

# diff_gene_KEGG_level2_plots = list(glob('7.enrichment/KEGG_enrichment/*/ALL_vs_DEG.KEGG_Classification.png'))
# diff_gene_KEGG_level2_contents = ['差异表达基因及所有基因在 KEGG Level2 水平分布比较图。图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。']*len(diff_gene_KEGG_level2_plots)
# diff_gene_KEGG_level2_plots = list(zip(diff_gene_KEGG_level2_plots, diff_gene_KEGG_level2_contents))

# diff_gene_KEGG_level2 = diff_gene_KEGG.add_plot(diff_gene_KEGG_level2_plots,caption='差异表达基因及所有基因在 KEGG Level2 水平分布比较图',description='图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。')
# diff_gene_KEGG_updown_level2_c = diff_gene_KEGG.add_comment('上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图如下：')
# diff_gene_KEGG_updown_level2_plots = list(glob('7.enrichment/KEGG_enrichment/*/Up_vs_Down.KEGG_Classification.png'))
# diff_gene_KEGG_updown_level2_contents = ['上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图。图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。']*len(diff_gene_KEGG_updown_level2_plots)
# diff_gene_KEGG_updown_level2_plots = list(zip(diff_gene_KEGG_updown_level2_plots, diff_gene_KEGG_updown_level2_contents))
# diff_gene_KEGG_updown_level2 = diff_gene_KEGG.add_plot(diff_gene_KEGG_updown_level2_plots,caption='上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图',description='图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。')




##########################高级分析##########################
# further_module = result_module.add_section('高级分析')

# ##细胞亚群分析
# subtype_module = further_module.add_section('细胞亚群分析')

# ##亚群分析
# subcluster_module = subtype_module.add_section('亚群聚类')
# subcluster_plot = subcluster_module.add_plot('subtype/2D_tSNE_groupby_cluster_resolution*_plot.png', caption = '亚群聚类结果')
# submarker_module = subtype_module.add_section('亚群Marker基因鉴定')
# submarker_plot = submarker_module.add_plot('subtype/topmarker_gene_heatmap.png', caption = '每个cluster亚群的Marker基因热图')
# subcelltype_module = subtype_module.add_section('亚群细胞类型注释')
# subcelltype_plot = subcelltype_module.add_plot('subtype/*_top.subtype_celltyping_on_tsne_resolution*_plot.png', caption = '亚群细胞类型鉴定结果')


# ##拟时序分析
# pseudo_module = further_module.add_section('拟时序分析')
# pseudotime_plot = pseudo_module.add_plot('pseudo/cell_trajectory_color_by_Pseudotime.png', caption = '拟时间轨迹结果', description = '图片说明：颜色由深到浅，表示分化的时间由早及晚。')
# pseudocluster_plot = pseudo_module.add_plot('pseudo/cell_trajectory_color_by_clusters.png', caption = '细胞分化拟时间轨迹结果', description = '图片说明：对应上一张图的细胞群分布情况，可初步推断各细胞群之间的分化关系。')

# pseudogene_plot = pseudo_module.add_plot('pseudo/gene.png', caption = '拟时序基因表达量图', description = '图片说明：横轴由左到右表示时间从早到晚，纵轴为基因的表达量，不同颜色的点代表不同细胞群，可以看到随着时间变化在不同细胞群中基因表达趋势的变化。')

# pseudoheatmap_plot = pseudo_module.add_plot('pseudo/heatmap.png', caption = '拟时序基因表达谱热图', description = '图片说明：蓝色表示低表达，红色表示高表达，并且对时间发育轨迹中具有相似表达模式的基因进行聚类。')


# ##细胞周期分析
# cellcycle_module = further_module.add_section('细胞周期分析')
# cellcycle_plot = cellcycle_module.add_plot('Cell_cycle/tSNE_resolution1_visualize_CellCycle_plot.png', caption = '细胞周期鉴定tSNE展示图', description = '图片说明：5种颜色分别代表5个分裂时期，展现细胞分裂时期在tSNE聚类图中的分布情况。')
# pseudocluster_plot = cellcycle_module.add_plot('Cell_cycle/cell_cycle_annotation_heatmap.png', caption = '细胞周期鉴定热图', description = '图片说明：列为细胞，行为五个细胞分裂时期：G1.S，S，G2.M，M，M.G1。')





##########################附录##########################
appendix_module = report.add_section('附录')

##########################实验技术方法说明##########################
link1 = 'supplemental_material/欧易生物单细胞转录组实验技术方法说明_中文.pdf'
link2 = 'supplemental_material/欧易生物单细胞转录组实验技术方法说明_英文.pdf'
supply_section_bio = appendix_module.add_section('实验技术方法说明',link1=link1,link2=link2)


##########################生信分析方法说明##########################
link3 = 'supplemental_material/欧易生物单细胞转录组生信分析方法_中文.pdf'
link4 = 'supplemental_material/欧易生物单细胞转录组生信分析方法_英文.pdf'
supply_section_analysis = appendix_module.add_section('生信分析方法说明',link3=link3,link4=link4)

##########################FAQ##########################
faq = appendix_module.add_section('常见问题FAQ',link_refdiff=link_refdiff)

##########################软件及数据库信息##########################
supply_section_database = appendix_module.add_section('数据库信息')
database_table = appendix_module.add_table('scRNA_database.txt', caption = '数据库信息')
supply_section_software = appendix_module.add_section('数据分析软件')
software_table = appendix_module.add_table('scRNA_software.txt', caption = '使用软件及版本')

other_module = report.add_section('申明')

# Generate Report HTML
report.write_to('report.html')

