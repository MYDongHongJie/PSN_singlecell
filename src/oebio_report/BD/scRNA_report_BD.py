# encoding:utf-8

import shutil
import os
from oebio.report import oeweb_register,Report
from oebio.app import *
from configparser import ConfigParser
import click
import subprocess
import sys
import re

@click.command()
@click.option('-i','--input', prompt='the program path',
			 help='the directory of your program. ')
@click.option('-c','--configfile', prompt='config file', default='./config.txt',
			 help='the config file which contain your program information.')
@click.option('-r','--rds',
			 help='the rds file of the report.')


def scrna_report(input, configfile, rds):

	""" python3 scRNA_report_BD.py  -i ./ -c ./config.txt"""
	config = ConfigParser()
	config.read(configfile, encoding="utf-8")
	program_num = dict(config['par'])['项目编号']
	program_path = os.path.abspath(input)
	orig_path = os.getcwd()

	# Report
	os.chdir(program_path)
	log_file = open("%s/report_log.txt" % (program_path), "w")
	if os.path.exists("%s/%s_Report" % (program_path, program_num)):
		shutil.rmtree("%s/%s_Report" % (program_path, program_num), ignore_errors=True)
	else:
		os.makedirs("%s/%s_Report" % (program_path, program_num))
		print("Create Report...")
		log_file.write("Create Report..." + "\n")

	num = 0
	outdir = "%s/%s_Report" % (program_path, program_num)
	# QC
	if os.path.exists("%s/Analysis"  % (program_path)):
		num += 1
		name = [str(i).split('/')[-1] for i in glob("Analysis/*" )]
		for j in name:
			os.makedirs("%s/%s.QC/%s" % (outdir, num, j))
			subprocess.call('ln -s %s/Analysis/%s/* %s/%s.QC/%s ' % (program_path, j, outdir, num, j), shell=True)
			subprocess.call('rm %s/%s.QC/%s/*.BAM*' % (outdir, num, j), shell=True)
			subprocess.call('rm -r %s/%s.QC/%s/Logs' % (outdir, num, j), shell=True)
	else:
		print("Can not find the QC Analysis results!")
		log_file.write("Can not find the QC Analysis results!" + "\n")

	# Count_QC
	if os.path.exists("%s/Count_QC"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Count_QC" % (outdir, num))
		subprocess.call('cp  %s/Count_QC/* %s/%s.Count_QC/' % (program_path, outdir, num), shell=True) ########################
		
	else:
		print("Can not find the Count_QC results!!!!!!!!!!!!!!!!")
		log_file.write("Can not find the Count_QC results!!!!!!!!!!!!!!!!!" + "\n")
		
	# Clustering
	if os.path.exists("%s/tsne_Dimension_Reduction"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Clustering" % (outdir, num))
		subprocess.call('cp %s/tsne_Dimension_Reduction/* %s/%s.Clustering/' % (program_path, outdir, num), shell=True)
	elif os.path.exists("%s/umap_Dimension_Reduction"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Clustering" % (outdir, num))
		subprocess.call('cp %s/umap_Dimension_Reduction/* %s/%s.Clustering/' % (program_path, outdir, num), shell=True)
	else:
		print("Can not find the Clustering results!!!!!!!!!!!!!!!!")
		log_file.write("Can not find the Clustering results!!!!!!!!!!!!!!!!" + "\n")
	if os.path.exists("%s/visualize_cluster_by_clusters" %(program_path)):
		subprocess.call('cp -r  %s/visualize_cluster_by_* %s/%s.Clustering/' % (program_path, outdir, num), shell=True)
	else:
		print("There is no visualize_cluster_by_clusters result.")
		log_file.write("There is no visualize_cluster_by_clusters result." + "\n")
	# Clusters_correlation
	if os.path.exists("%s/Clusters_correlation"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Clusters_correlation" % (outdir, num))
		subprocess.call('cp %s/Clusters_correlation/* %s/%s.Clusters_correlation/' % (program_path, outdir, num), shell=True)
	else:
		print("Can not find the Clusters_correlation results!!!!!!!!!!!!!!!!")
		log_file.write("Can not find the Clusters_correlation results!!!!!!!!!!!!!!!!" + "\n")
		
	# Marker
	if os.path.exists("%s/Marker"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Marker/top10_markers_visualize_for_each_cluster" % (outdir, num))
		subprocess.call('cp -r  %s/Marker/markers_vis4cluster* %s/%s.Marker/top10_markers_visualize_for_each_cluster/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp -r  %s/Marker/topmarker_gene_heatmap.p* %s/%s.Marker/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp -r  %s/Marker/top10_markers_for_each_cluster_anno.xls %s/%s.Marker/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp -r  %s/Marker/all_markers_for_each_cluster_anno.xls %s/%s.Marker/' % (program_path, outdir, num), shell=True)
	
	else:
		print("Can not find the Marker results!!!!!!!!!!!!!!!!")
		log_file.write("Can not find the Marker results!!!!!!!!!!!!!!!!" + "\n")
		
	# Reference_celltype
	if os.path.exists("%s/Reference_celltype"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Reference_celltype" % (outdir, num))
		subprocess.call('cp %s/Reference_celltype/* %s/%s.Reference_celltype/' % (program_path, outdir, num), shell=True)
		subprocess.call('rm %s/%s.Reference_celltype/*.rds' % (outdir, num), shell=True)
	else:
		print("Can not find the Reference_celltype results!!!!!!!!!!!!!!!!")
		log_file.write("Can not find the Reference_celltype results!!!!!!!!!!!!!!!!" + "\n")
		
	# Diffexp
	if os.path.exists("%s/Diffexp"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.Diffexp" % (outdir, num))
		subprocess.call('cp %s/Diffexp/diffexp_results_stat.xls %s/%s.Diffexp/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp %s/Diffexp/*-vs-*-diff-*_anno.xls %s/%s.Diffexp/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp %s/Diffexp/*-vs-*-all_diffexp_genes_anno.xls %s/%s.Diffexp/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp %s/Diffexp/top*-vs-*_genes.xls %s/%s.Diffexp/' % (program_path, outdir, num), shell=True)
		subprocess.call('cp %s/Diffexp/*_heatmap.p* %s/%s.Diffexp/' % (program_path, outdir, num), shell=True)
	else:
		print("There is no Diffexp result.")
		log_file.write("There is no Diffexp result." + "\n")

	# enrichment
	if os.path.exists("%s/enrichment"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.enrichment" % (outdir, num))
		subprocess.call('cp -r %s/enrichment/* %s/%s.enrichment/' % (program_path, outdir, num), shell=True)
	else:
		print("There is no enrichment result.")
		log_file.write("There is no enrichment result." + "\n")
	
	# PPI
	if os.path.exists("%s/PPI"  % (program_path)):
		num += 1
		os.makedirs("%s/%s.PPI" % (outdir, num))
		subprocess.call('cp -r %s/PPI/*string_protein-protein-interaction.* %s/%s.PPI/' % (program_path, outdir, num), shell=True)
	else:
		print("There is no PPI result.")
		log_file.write("There is no PPI result." + "\n")
	
	# subprocess.call('cp -r /public/hstore1/works/scRNA_works/luyao/test/git/scRNA-seq/oebio_report/BD/supplemental_material %s/' % (outdir), shell=True)
	
	# rds transfer
	if not (rds is None):
		rds_file = os.path.join(orig_path,rds)
		if os.path.isfile(rds_file):
			os.system("ssh scrna@192.168.10.24 mkdir /oecloud/userfiles/%s" %(program_num))
			if ( os.system("scp %s scrna@192.168.10.24:/oecloud/userfiles/%s/" %(rds_file, program_num)) != 0 ):
				print("Fails to transfer rds file. Please check your network and transfer again.")
				log_file.write("Fails to transfer rds file. Please check your network and transfer again." + " \n")
			else:
				print("%s rds uploaded successfully." % rds )
				log_file.write("%s rds uploaded successfully."  % rds + "\n ")
		else:
			print("Couldn't find rds input file!!!!! Please check the file path.")
			log_file.write("Couldn't find rds input file!!! Please check the file path." + "\n")
	

	log_file.close()
	
	
	#######################################################Generate report############################################################
	os.chdir(outdir)
	#定义函数批量获取路径
	def get_url(target,tag):
		urls = []
		for f in Path(target).glob(tag):
			text = f'[{f}]({f})\n\n'
			urls.append(text)
		return ''.join(urls)
	
	# Define header info
	report = Report('BD单细胞转录组项目结题报告', title = 'BD单细胞转录组项目结题报告', header_info = dict(config['par']), oe_welcome='''
	
####  感谢您选择欧易生物！

**关于网页报告使用有以下几点提示:**

1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节，便于您快速查看相应内容。
2. **报告正文**中的**图片**均可以点击后进行**放大查看**，且图片点击后可以包含更多的细节信息，比如左上角会显示具体的**图片数量**，右上角有相关的**图片工具**，**键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张，更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**，当表格多于20行的时候，表格会嵌入到网页之中，可以利用滚动栏进行调整设置，同时会在表格上方增加搜索设置，便于快速查询表格信息。
5. 在报告中浏览，需要返回顶部的时候，可以使用网页**右下角的白色箭头标识**快速返回顶部。
6. 本提示可以点击右上方的不再显示，则后续打开网页均不会显示本提示。
	''')
	report.add_yaml_config('/home/luyao/10X_scRNAseq_v3/src/oebio_report/BD/report_BD.yaml')
	os.environ['OEBIO'] = "/home/luyao/10X_scRNAseq_v3/src/oebio_report/BD/pic"
	os.environ['oeweb_register_token'] = 'oearray'
	#######################################项目概况#######################################
	project_module = report.add_section('项目概况')
	register_info = oeweb_register(project_id=program_num, target_url='https://cloud.oebiotech.cn/task/category/scrna', note='本报告包含项目基本分析内容，如要快速做其它分析调整，如细胞表达可视化等')
	if register_info:
		project_module.add_comment(register_info)

	##########################项目摘要、快捷链接##########################
	QC_folder = list(glob(r"[0-9]*.Count_QC"))[0]
	summ = pd.read_csv('%s/cell_statitics_before_after_QC.xls' %(QC_folder),index_col=0,sep='\t')
	sample_num = len(summ)   #样本数量
	min_beforeQC_cell = min(summ.Total_cells_beforeQC)
	max_beforeQC_cell = max(summ.Total_cells_beforeQC)
	min_afterQC_cell = min(summ.Total_cells_afterQC)
	max_afterQC_cell = max(summ.Total_cells_afterQC)
	min_afterQC_umi = format(min(summ.Mean_nUMI_afterQC.astype(float)),'.0f')
	max_afterQC_umi = format(max(summ.Mean_nUMI_afterQC.astype(float)),'.0f')
	min_afterQC_gene = format(min(summ.Mean_nGene_afterQC.astype(float)),'.0f')
	max_afterQC_gene = format(max(summ.Mean_nGene_afterQC.astype(float)),'.0f')
	min_afterQC_mito = format(min(summ['Mean_mito.percent_afterQC'].astype(float)),'.4f')
	max_afterQC_mito = format(max(summ['Mean_mito.percent_afterQC'].astype(float)),'.4f')
	
	QC0_folder = list(glob(r"[0-9]*.QC"))[0]
	Clustering_folder = list(glob(r"[0-9]*.Clustering"))[0]
	Clusters_correlation_folder = list(glob(r"[0-9]*.Clusters_correlation"))[0]
	Marker_folder = list(glob(r"[0-9]*.Marker"))[0]
	Marker_file =  list(glob('%s/top10_markers_for_each_cluster_anno.xls' %(Marker_folder)))[0]
	Marker_results = pd.read_csv(Marker_file, sep='\t')
	cluster_num = max(Marker_results['cluster'])       #获取细胞群数
	
	# link_refdiff = 'supplemental_material/单细胞转录组常见问题(FAQ).pdf'
	# link0 = 'supplemental_material/Loupe-Cell-Browser使用教程.html'
	
	celltype_folder = list(glob(r"[0-9]*.Reference_celltype"))
	
	if sample_num == 1:        #1个样本
		if len(celltype_folder) == 0:        #特殊物种
			abstract_module = project_module.add_section('项目摘要_单样本_无细胞类型',sample_num=sample_num,min_beforeQC_cell=min_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,min_afterQC_umi=min_afterQC_umi,min_afterQC_gene=min_afterQC_gene,min_afterQC_mito=min_afterQC_mito,cluster_num=cluster_num)              #摘要
			linkage = project_module.add_section('项目快捷链接_无细胞类型_无差异')        #快捷链接
		else:
			celltype_file = list(glob('%s/*_simplified_celltyping_results.csv' %(celltype_folder[0])))
			cell_results = pd.read_csv(celltype_file[0], sep=',')
			celltype_tmp = [ str(i) for i in set(cell_results['celltype']) ]
			celltype = ', '.join(celltype_tmp)
			Diffexp_folder = list(glob(r"[0-9]*.Diffexp"))
			abstract_module = project_module.add_section('项目摘要_单样本',sample_num=sample_num,min_beforeQC_cell=min_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,min_afterQC_umi=min_afterQC_umi,min_afterQC_gene=min_afterQC_gene,min_afterQC_mito=min_afterQC_mito,cluster_num=cluster_num,celltype=celltype)
			inkage = project_module.add_section('项目快捷链接_无差异')        #快捷链接
	
	else:
		if len(celltype_folder) == 0:        #特殊物种
			Diffexp_folder = list(glob(r"[0-9]*.Diffexp"))
			if len(Diffexp_folder) == 0:
				abstract_module = project_module.add_section('项目摘要_无细胞类型_无差异',sample_num=sample_num,min_beforeQC_cell=min_beforeQC_cell,max_beforeQC_cell=max_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,max_afterQC_cell=max_afterQC_cell,min_afterQC_umi=min_afterQC_umi,max_afterQC_umi=max_afterQC_umi,min_afterQC_gene=min_afterQC_gene,max_afterQC_gene=max_afterQC_gene,min_afterQC_mito=min_afterQC_mito,max_afterQC_mito=max_afterQC_mito,cluster_num=cluster_num)              #摘要
				linkage = project_module.add_section('项目快捷链接_无细胞类型_无差异')        #快捷链接
			else:
				diff = pd.read_csv('%s/diffexp_results_stat.xls' %(Diffexp_folder[0]),sep='\t')
				num_diff = len(diff)
				diff_num_tmp = [ str(i) for i in diff.iloc[:,4]]  ### diff_num_tmp = [ str(i) for i in diff.iloc[:,4].to_list()]
				diff_num = ', '.join(diff_num_tmp)
				abstract_module = project_module.add_section('项目摘要_无细胞类型_有差异',sample_num=sample_num,min_beforeQC_cell=min_beforeQC_cell,max_beforeQC_cell=max_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,max_afterQC_cell=max_afterQC_cell,min_afterQC_umi=min_afterQC_umi,max_afterQC_umi=max_afterQC_umi,min_afterQC_gene=min_afterQC_gene,max_afterQC_gene=max_afterQC_gene,min_afterQC_mito=min_afterQC_mito,max_afterQC_mito=max_afterQC_mito,cluster_num=cluster_num,num_diff=num_diff,diff_num=diff_num)
				linkage = project_module.add_section('项目快捷链接_无细胞类型_有差异')        #快捷链接
			
		else:
			celltype_file = list(glob('%s/*_simplified_celltyping_results.csv' %(celltype_folder[0])))
			cell_results = pd.read_csv(celltype_file[0], sep=',')
			celltype_tmp = [ str(i) for i in set(cell_results['celltype']) ]
			celltype = ', '.join(celltype_tmp)
			
			Diffexp_folder = list(glob(r"[0-9]*.Diffexp"))
			if len(Diffexp_folder) == 0:
				abstract_module = project_module.add_section('项目摘要_有细胞类型_无差异',sample_num=sample_num,min_beforeQC_cell=min_beforeQC_cell,max_beforeQC_cell=max_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,max_afterQC_cell=max_afterQC_cell,min_afterQC_umi=min_afterQC_umi,max_afterQC_umi=max_afterQC_umi,min_afterQC_gene=min_afterQC_gene,max_afterQC_gene=max_afterQC_gene,min_afterQC_mito=min_afterQC_mito,max_afterQC_mito=max_afterQC_mito,cluster_num=cluster_num,celltype=celltype)
				linkage = project_module.add_section('项目快捷链接_无差异')        #快捷链接
				
			else:
				diff = pd.read_csv('%s/diffexp_results_stat.xls' %(Diffexp_folder[0]),sep='\t')
				num_diff = len(diff)
				diff_num_tmp = [ str(i) for i in diff.iloc[:,4]]  ### diff_num_tmp = [ str(i) for i in diff.iloc[:,4].to_list()]
				diff_num = ', '.join(diff_num_tmp)
				abstract_module = project_module.add_section('项目摘要_有细胞类型_有差异',sample_num=sample_num,min_beforeQC_cell=min_beforeQC_cell,max_beforeQC_cell=max_beforeQC_cell,min_afterQC_cell=min_afterQC_cell,max_afterQC_cell=max_afterQC_cell,min_afterQC_umi=min_afterQC_umi,max_afterQC_umi=max_afterQC_umi,min_afterQC_gene=min_afterQC_gene,max_afterQC_gene=max_afterQC_gene,min_afterQC_mito=min_afterQC_mito,max_afterQC_mito=max_afterQC_mito,cluster_num=cluster_num,celltype=celltype,num_diff=num_diff,diff_num=diff_num)
				linkage = project_module.add_section('项目快捷链接_差异')        #快捷链接
			
		
	##########################建库测序流程##########################
	library_module = report.add_section('实验技术流程')
	library_module_c = library_module.add_comment('实验建库测序流程图如下：')
	library_module_plot = library_module.add_plot('Experiment.png',caption ='建库测序流程图',content='建库测序流程图')

	##########################生信分析流程##########################
	analysis_module = report.add_section('生信分析流程')
	analysis_module_c = analysis_module.add_comment('生物信息分析流程图如下：')
	analysis_module_plot = analysis_module.add_plot('Pipeline.png', caption = '生物信息分析流程',content='生物信息分析流程图')

	##########################项目分析结果##########################
	result_module = report.add_section('项目分析结果')
	
	##########################基因定量质控##########################
	QC_1 = result_module.add_section('基因定量质控')
	QC_link = get_url(QC0_folder,'*/*Metrics_Summary.csv')
	QC_table0 = QC_1.add_comment('基因定量质控统计表格见：')
	QC_table = QC_1.add_comment("{QC_link}",QC_link = QC_link) 

	QC1_table = QC_1.add_table('scRNA_qc_stat.txt', caption = '基因定量质控统计结果说明')
	
	##########################基因定量质控##########################
	Quantification_module = result_module.add_section('定量后质控')
	nGene_nUMI_fold2sd = dict(config['parameter'])['ngene_numi_fold2sd']
	percent_mito = int(float(dict(config['parameter'])['percent_mito'])*100)
	Quantification1_module = Quantification_module.add_section('过滤低质量细胞', nGene_nUMI_fold2sd=nGene_nUMI_fold2sd, percent_mito=percent_mito)
	
	Quantification1_plots = list(glob('*.Count_QC/outliers.png'))
	Quantification1_contents = ['拟合广义线性模型过滤离域细胞图。图片说明：横轴为每个细胞内的 UMI 数，纵轴为每个细胞内的基因数，根据二者的线性关系拟合分布模型，着色的点表示离域细胞，在下游分析中会进行剔除。']
	Quantification1_plots = list(zip(Quantification1_plots, Quantification1_contents))

	Quantification1_plot = Quantification1_module.add_plot(Quantification1_plots, caption = '拟合广义线性模型过滤离域细胞图',description = '图片说明：横轴为每个细胞内的 UMI 数，纵轴为每个细胞内的基因数，根据二者的线性关系拟合分布模型，着色的点表示离域细胞，在下游分析中会剔除。')

	mito_plots = sorted(list(glob('*.Count_QC/*QC_mitochondroin_transcript_ratio_in_each_cell_violin_plot.png')),reverse=True)
	mito_contents = ['质控前每个细胞中线粒体基因比例的小提琴分布图。图片说明：纵轴表示线粒体基因占单个细胞所有基因的比例，图中每个点代表一个油包水微滴中细胞的线粒体基因比例，每个小提琴图反映对应样本中所有细胞的线粒体基因在细胞所有基因中所占的比例，一般要求大部分细胞的线粒体基因比例越低越好（特殊样本除外）。', '质控后每个细胞中线粒体基因比例的小提琴分布图。图片说明：纵轴表示线粒体基因占单个细胞所有基因的比例，图中每个点代表一个油包水微滴中细胞的线粒体基因比例，每个小提琴图反映对应样本中所有细胞的线粒体基因在细胞所有基因中所占的比例，一般要求大部分细胞的线粒体基因比例越低越好（特殊样本除外）。']
	mito_plots = list(zip(mito_plots, mito_contents))
	mito_plot = Quantification1_module.add_plot(mito_plots, caption = '质控前后每个细胞中线粒体基因比例的小提琴分布图' , description = '图片说明：纵轴表示线粒体基因占单个细胞所有基因的比例，图中每个点代表一个油包水微滴中细胞的线粒体基因比例，每个小提琴图反映对应样本中所有细胞的线粒体基因在细胞所有基因中所占的比例，一般要求大部分细胞的线粒体基因比例越低越好（特殊样本除外）。左图为质控前，右图为质控后。')

	ngene_plots = sorted(list(glob('*.Count_QC/*QC_total_genes4each_cell_on_violin_plot.png')),reverse=True)
	ngene_contents = ['质控前每个细胞中基因数目的分布图。图片说明：纵轴表示细胞中有表达的基因数目，图中每个点代表一个油包水微滴中细胞的基因数目。该图反映样本中的每一个细胞表达基因的数目，基因数目异常过多的点很可能是由于对应的油包水微滴中包含多个细胞，需要通过设置合理的阈值将其过滤。左图为质控前，右图为质控后。', '质控后每个细胞中基因数目的分布图。图片说明：纵轴表示细胞中有表达的基因数目，图中每个点代表一个油包水微滴中细胞的基因数目。该图反映样本中的每一个细胞表达基因的数目，基因数目异常过多的点很可能是由于对应的油包水微滴中包含多个细胞，需要通过设置合理的阈值将其过滤。左图为质控前，右图为质控后。']
	ngene_plots = list(zip(ngene_plots, ngene_contents))
	ngene_plot = Quantification1_module.add_plot(ngene_plots, caption = '质控前后每个细胞中基因表达数目的小提琴分布图', description = '图片说明：纵轴表示细胞中有表达的基因数目，图中每个点代表一个油包水微滴中细胞的基因数目。该图反映样本中的每一个细胞表达基因的数目，基因数目异常过多的点很可能是由于对应的油包水微滴中包含多个细胞，需要通过设置合理的阈值将其过滤。左图为质控前，右图为质控后。')

	nUMI_plots = sorted(list(glob('*.Count_QC/*QC_total_UMIs4each_cell_on_violin_plot.png')),reverse=True)
	nUMI_contents = ['质控前每个细胞中 UMI 数目的分布图。图片说明：纵轴表示 UMI 数，图中的每个点代表一个油包水微滴中细胞的 UMI 数目，即转录本的数目，该图反映样本中每一个细胞的转录本数目，转录本数目异常过多的点很可能是由于对应的油包水微滴中包含多个细胞，需要通过设置合理的阈值将其过滤。', '质控后每个细胞中 UMI 数目的分布图。图片说明：纵轴表示 UMI 数，图中的每个点代表一个油包水微滴中细胞的 UMI 数目，即转录本的数目，该图反映样本中每一个细胞的转录本数目，转录本数目异常过多的点很可能是由于对应的油包水微滴中包含多个细胞，需要通过设置合理的阈值将其过滤。']
	nUMI_plots = list(zip(nUMI_plots, nUMI_contents))
	nUMI_plot = Quantification1_module.add_plot(nUMI_plots, caption = '质控前后每个细胞中 UMI 数目的小提琴分布图', description = '图片说明：纵轴表示 UMI 数，图中的每个点代表一个油包水微滴中细胞的 UMI 数目，即转录本的数目，该图反映样本中每一个细胞的转录本数目，转录本数目异常过多的点很可能是由于对应的油包水微滴中包含多个细胞，需要通过设置合理的阈值将其过滤。左图为质控前，右图为质控后。')

	Quantification1_c2 = Quantification1_module.add_comment("定量质控前后的细胞数统计情况如下表所示：")
	Quantification1_table = Quantification1_module.add_table('*.Count_QC/cell_statitics_before_after_QC.xls', caption = '质控前后细胞数目统计表', headerdict={'sampleid':'样本名称','Mean_nUMI_beforeQC':'质控前细胞中的平均UMI数','Mean_nGene_beforeQC':'质控前细胞中的平均基因数','Mean_mito.percent_beforeQC':'质控前细胞中的平均线粒体比例','Total_cells_beforeQC':'质控前的总细胞数','Mean_nUMI_afterQC':'质控后细胞中的平均UMI数','Mean_nGene_afterQC':'质控后细胞中的平均基因数','Mean_mito.percent_afterQC':'质控后细胞中的平均线粒体比例','Total_cells_afterQC':'质控后的总细胞数'})
	
	##########################降维与聚类分析##########################
	if os.path.exists("%s/pca_Dimension_Reduction"  % (program_path)):
		reduct1 = "PCA(Principal Components Analysis, 主成分分析)"
	elif os.path.exists("%s/mnn_Dimension_Reduction"  % (program_path)):
		reduct1 = "MNN(mutual nearest neighbors, 互享最近邻)"
	if os.path.exists("%s/tsne_Dimension_Reduction"  % (program_path)):
		reduct2 = "t-SNE"
	elif os.path.exists("%s/umap_Dimension_Reduction"  % (program_path)):
		reduct2 = "UMAP"
	
	Clustering_module = result_module.add_section('降维与聚类分析', nav_emphasize = True, nav_title = "项目核心内容")
	tsne_module = Clustering_module.add_section('降维聚类结果')
	if reduct2 == "t-SNE":
		tsne_c1 = tsne_module.add_comment("本项目中采用的降维算法为 %s 和 t-SNE(t-distributed Stochastic Neighbor Embedding, t-分布邻域嵌入)算法。基于 %s 的降维结果通过 t-SNE 对单细胞群聚类进行可视化，聚类算法采用SNN，最终获得最优细胞分群。" %(reduct1, reduct1))
	elif reduct2 == "UMAP":
		tsne_c1 = tsne_module.add_comment("本项目中采用的降维算法为 %s 和 UMAP(统一流形逼近与投影)算法。UMAP 主要基于流形理论和拓扑算法的理论，对高维数据进行降维，从而能够保留更多数据的全局结构，并且具有优越的运行性能。基于 %s 降维结果通过 UMAP 对单细胞群聚类进行可视化，聚类算法采用SNN，最终获得最优细胞分群。" %(reduct1, reduct1))

	tsne_plots = list(glob('*.Clustering/*_groupby_cluster_resolution*_plot.png'))
	tsne_contents = ['降维聚类结果图。图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。']
	tsne_plots = list(zip(tsne_plots, tsne_contents))

	tsne_plot = tsne_module.add_plot(tsne_plots, caption = '降维聚类结果图', description ='图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')
	
	coordination_table = list(glob('*.Clustering/*_Dimension_Reduction_coordination.csv'))[0]
	Clustering_link = get_url(Clustering_folder, str(coordination_table).split('/')[-1])   ####
	tsne_table0 = tsne_module.add_comment('降维聚类二维坐标表格见：')
	tsne_table = tsne_module.add_comment("{Clustering_link}",Clustering_link=Clustering_link) 
	
	tsne_3dplot = list(glob('*.Clustering/*_resolution*_3D_plot.html'))[0]
	tsne_3dplot_link= get_url(Clustering_folder, str(tsne_3dplot).split('/')[-1])
	tsne_3dplot_c = tsne_module.add_comment('降维聚类 3D 展示图见：')
	tsne_3dplot_c2 = tsne_module.add_comment("{tsne_3dplot_link}",tsne_3dplot_link=tsne_3dplot_link) 
	tsne_3dplot_2 = tsne_module.add_plot('3D_plot.png', caption = '降维聚类 3D 示例图',content='降维聚类 3D 示例图', description ='图片说明：降维聚类 3D 展示图为交互式网页界面，可以使用鼠标自由拖动、缩放，查看三维状态下任一角度的降维聚类结果，点击右上角的相机按钮即可保存图片。')

	if os.path.exists("%s/visualize_cluster_by_clusters" %(Clustering_folder) ):
		groupby_module = Clustering_module.add_section('样本间降维聚类分组展示')
		groupby_plots = list(glob('*.Clustering/visualize_cluster_by_clusters/groupby-sampleid_resolution*_contrast_plot.png' ))
		groupby_contents = ['多样本降维聚类分组展示图。图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。']
		groupby_plots = list(zip(groupby_plots, groupby_contents))
		groupby_plot = groupby_module.add_plot(groupby_plots, caption = '样本间降维聚类分组展示图', description ='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。')
		
		groupby_c = groupby_module.add_comment('每个细胞群中样本占比的柱状统计图如下：')
		groupby_plots3 = list(glob('*.Clustering/visualize_cluster_by_clusters/groupby-clusters_resolution*_summary_plot.png'))
		groupby_contents = ['细胞群样本占比柱状图。图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。']
		groupby_plots3 = list(zip(groupby_plots3, groupby_contents))
		groupby_plot3 = groupby_module.add_plot(groupby_plots3, caption = '细胞群样本占比柱状图', description ='图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。')
		
		groupby_c = groupby_module.add_comment('多个样本的降维聚类分面展示图如下：')
		groupby_plots2 = list(glob('*.Clustering/visualize_cluster_by_clusters/splitby-sampleid_resolution*_split_plot.png' ))
		groupby_contents = ['多样本降维聚类分面展示图。图片说明：横纵坐标分别代表降维的第一和第二主成分，不同群的细胞以不同颜色区分。']
		groupby_plots2 = list(zip(groupby_plots2, groupby_contents))
		groupby_plot2 = groupby_module.add_plot(groupby_plots2, caption = '多样本降维聚类分面展示图', description ='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同群的细胞以不同颜色区分。')
		
		groupby_c = groupby_module.add_comment('每个样本中细胞群占比的柱状统计图如下：')
		groupby_plots4 = list(glob('*.Clustering/visualize_cluster_by_clusters/groupby-sampleid_resolution*_summary_plot.png' ))
		groupby_contents = ['样本间细胞群占比柱状图。图片说明：横坐标表示不同样本，纵坐标表示不同组别中细胞数目所占的百分比。']
		groupby_plots4 = list(zip(groupby_plots4, groupby_contents))
		groupby_plot4 = groupby_module.add_plot(groupby_plots4, caption = '样本间细胞群占比柱状图', description ='图片说明：横坐标表示不同样本，纵坐标表示不同组别中细胞数目所占的百分比。')

		Quantification1_c2 = groupby_module.add_comment("各样本在不同细胞群中的细胞数目统计表如下：")
		Quantification1_table = groupby_module.add_table('*.Clustering/visualize_cluster_by_clusters/clust_cond_freq_info.xls', caption = '样本间细胞群数目统计表', headerdict={'sampleid':'样本名称','clusters':'细胞群编号','cell_number':'细胞数目','freq':'该样本在当前细胞群的细胞数占该样本所有细胞数的百分比'})
		
		groupby_link = get_url(Clustering_folder,'visualize_cluster_by_clusters')
		groupby_c2 = groupby_module.add_comment("详细结果目录见 : ")
		groupby_link2 = groupby_module.add_comment("{groupby_link}",groupby_link=groupby_link) 
	
	##########################细胞群间相似性分析##########################
	Clusters_correlation_module = result_module.add_section('细胞群间相似性分析')
	Clusters_correlation_plots = list(glob('*.Clusters_correlation/coefficient_heatmap.png'))
	Clusters_correlatio_contents = ['细胞群间相关性热图。图片说明：横纵坐标为不同细胞群，热图中的数字为皮尔森相关性系数值。该值越大，热图颜色越红，表示细胞群之间的相关性程度越高。']
	Clusters_correlatio_plots2 = list(zip(Clusters_correlation_plots, Clusters_correlatio_contents))
	Clusters_correlation_plot = Clusters_correlation_module.add_plot(Clusters_correlatio_plots2, caption = '细胞群间相关性热图', description = '图片说明：横纵坐标为不同细胞群，图中的数字为皮尔森相关性系数。该值越大，热图颜色越红，表示细胞群之间的相关性程度越高。')
	
	Clusters_correlation_c = Clusters_correlation_module.add_comment("各细胞群中的基因平均表达量表格见：")
	Clusters_correlation_link = get_url(Clusters_correlation_folder,'normalized_data_groupby_clusters.xls')
	Clusters_correlation_link2 = Clusters_correlation_module.add_comment("{Clusters_correlation_link}",Clusters_correlation_link=Clusters_correlation_link) 
	
	##########################Marker基因鉴定##########################
	marker_test = dict(config['parameter'])['marker_test']
	Marker_module = result_module.add_section('Marker基因鉴定', marker_test=marker_test, nav_emphasize = True, nav_title = "项目核心内容")

	Marker_table = Marker_module.add_table('*.Marker/top10_markers_for_each_cluster_anno.xls', caption = '每个细胞群 Top10 Marker 基因列表', show_rows=10 , headerdict={'gene':'基因名','p_val':'p值','avg_logFC':'对数转化后的平均Foldchange值','pct.1':'表达Marker基因的细胞在当前群中的占比','pct.2':'表达Marker基因的细胞在其余群中的占比','p_val_adj':'校正后的p值','cluster':'细胞群编号','gene_diff':'pct.1与pct.2的比值。Top10 Marker以gene_diff列为筛选依据。','ensemble_id':'ensemble基因ID','Dbxref_GeneID':' NCBI基因ID索引号','gene_type':'基因的类型','gene_description':'基因描述','GO_id':'Gene Ontology登录号','GO_term':'GO条目描述','pathway':'KEGG通路号','pathway_description':'KEGG通路描述'})

	Marker_link1 = list(glob('%s/all_markers_for_each_cluster_anno.xls' %(Marker_folder)))[0]
	Marker_link2 = list(glob('%s/top10_markers_for_each_cluster_anno.xls' %(Marker_folder)))[0]
	Marker_table0 = Marker_module.add_comment("每个细胞群中所有 Marker 基因结果表格：")
	Marker_table1 = Marker_module.add_comment("[{Marker_link1}]({Marker_link1})",Marker_link1=Marker_link1)
	Marker_table2 = Marker_module.add_comment("每个细胞群中 Top10 Marker 基因结果表格：")
	Marker_table3 = Marker_module.add_comment("[{Marker_link2}]({Marker_link2})",Marker_link2=Marker_link2)

	heatmap_plots = list(glob('*.Marker/topmarker_gene_heatmap.png'))
	heatmap_contents = ['Top10 Marker 基因表达热图。图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。']
	heatmap_plots = list(zip(heatmap_plots, heatmap_contents))

	heatmap_plot = Marker_module.add_plot(heatmap_plots, caption = 'Top10 Marker 基因表达热图', description = '图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')

	Marker_c = Marker_module.add_comment('每个细胞群中 Top10 Marker 基因的可视化展示如下：')
	feature_plots = ["*.Marker/top10_markers_visualize_for_each_cluster/markers_vis4cluster"+str(i+1)+"/topmarker_gene_featureplot.png" for i in range(cluster_num) ]
	feature_contents = ["第 "+str(i+1)+" 群细胞的 Top10 Marker 基因在 t-SNE 聚类结果中的可视化图。图片说明：红色越深表示该细胞中对应基因的表达量越高。" for i in range(len(feature_plots)) ]    ### str(i+1)
	feature_plots = list(zip(feature_plots, feature_contents))
	feature_plot = Marker_module.add_plot(feature_plots, caption = 'Top10 Marker 基因在 t-SNE 聚类结果中的可视化图', description = '图片说明：红色越深表示该细胞中对应基因的表达量越高。')
	vln_plots = ["*.Marker/top10_markers_visualize_for_each_cluster/markers_vis4cluster"+str(i+1)+"/topmarker_gene_violin_plot.png" for i in range(cluster_num) ]
	vln_contents = ["第 "+str(i+1)+" 群细胞的 Top10 Marker基因表达量小提琴图。图片说明：横坐标为细胞群的编号，纵坐标为标准化后的基因表达值。" for i in range(len(vln_plots)) ]    ### str(i+1)
	vln_plots = list(zip(vln_plots, vln_contents))
	vln_plot = Marker_module.add_plot(vln_plots, caption = 'Top10 Marker 基因表达量小提琴图', description = '图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')
	
	
	##########################细胞类型鉴定##########################
	if len(celltype_folder) > 0:
		Celltype_module = result_module.add_section('细胞类型鉴定', nav_emphasize = True, nav_title = "项目核心内容")
		Celltype_link = list(glob('%s/*top.*celltyping_on_tsne*.pdf' %(celltype_folder[0])))[0]
		refcelltype = str(Celltype_link).split('/')[-1].split("_")[1]
		if refcelltype.lower() == "hpca":
			Celltype_c1 = Celltype_module.add_comment('''本项目采用 SingleR[^SingleR] 包，基于 HPCA[^hpca]参考数据集进行细胞类型注释。该方法通过计算单细胞参考表达谱数据集与待鉴定的细胞表达谱之间的相关性，将待鉴定细胞注释为与参考数据集中相关性最高的一种细胞类型。

报告中的数据集鉴定结果供参考，后续可根据文献已有相关基因对细胞群的特征加以描述和验证。
[^SingleR]: Aran D, Looney A P, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage[J]. Nature immunology, 2019, 20(2): 163-172.
[^hpca]: Mabbott N A, Baillie J K, Brown H, et al. An expression atlas of human primary cells: inference of gene function from coexpression networks[J]. BMC genomics, 2013, 14(1): 632.''')
		elif refcelltype.lower() == "blueprint_encode":
			Celltype_c1 = Celltype_module.add_comment('''本项目采用 SingleR[^SingleR] 包，基于 Blueprint+Encode[^Blueprint][^Encode]参考数据集进行细胞类型注释。该方法通过计算单细胞参考表达谱数据集与待鉴定的细胞表达谱之间的相关性，将待鉴定细胞注释为与参考数据集中相关性最高的一种细胞类型。

报告中的数据集鉴定结果供参考，后续可根据文献已有相关基因对细胞群的特征加以描述和验证。
[^SingleR]: Aran D, Looney A P, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage[J]. Nature immunology, 2019, 20(2): 163-172.
[^Blueprint]: Martens, J. H. A, Stunnenberg, H. G. BLUEPRINT: mapping human blood cell epigenomes[J]. Haematologica, 98(10):1487-1489.
[^Encode]: Bernstein B E, Birney E, Dunham I, et al. An integrated encyclopedia of DNA elements in the human genome[J]. Nature, 2012, 489(7414): 57-74.
			''')
		elif refcelltype.lower() == "schcl":
			Celltype_c1 = Celltype_module.add_comment('''本项目采用 SingleR[^SingleR] 包，基于 schcl[^schcl]考数据集进行细胞类型注释。该方法通过计算单细胞参考表达谱数据集与待鉴定的细胞表达谱之间的相关性，将待鉴定细胞注释为与参考数据集中相关性最高的一种细胞类型。

报告中的数据集鉴定结果供参考，后续可根据文献已有相关基因对细胞群的特征加以描述和验证。
[^SingleR]: Aran D, Looney A P, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage[J]. Nature immunology, 2019, 20(2): 163-172.
[^schcl]: http://bis.zju.edu.cn/HCL/index.html
			''')
		elif refcelltype.lower() == "immgen":
			Celltype_c1 = Celltype_module.add_comment('''本项目采用 SingleR[^SingleR] 包，基于 immgen[^immgen]参考数据集进行细胞类型注释。该方法通过计算单细胞参考表达谱数据集与待鉴定的细胞表达谱之间的相关性，将待鉴定细胞注释为与参考数据集中相关性最高的一种细胞类型。

报告中的数据集鉴定结果供参考，后续可根据文献已有相关基因对细胞群的特征加以描述和验证。
[^SingleR]: Aran D, Looney A P, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage[J]. Nature immunology, 2019, 20(2): 163-172.
[^immgen]: Heng T S P , Painter M W , Elpek K , et al. The Immunological Genome Project: networks of gene expression in immune cells[J]. Nature Immunology, 2008, 9(10):1091-1094.
			''')
		elif refcelltype.lower() == "mouse.rnaseq":
			Celltype_c1 = Celltype_module.add_comment('''本项目采用 SingleR[^SingleR] 包，基于 mouse.rnaseq[^mouse.rnaseq]参考数据集进行细胞类型注释。该方法通过计算单细胞参考表达谱数据集与待鉴定的细胞表达谱之间的相关性，将待鉴定细胞注释为与参考数据集中相关性最高的一种细胞类型。

报告中的数据集鉴定结果供参考，后续可根据文献已有相关基因对细胞群的特征加以描述和验证。
[^SingleR]: Aran D, Looney A P, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage[J]. Nature immunology, 2019, 20(2): 163-172.
[^mouse.rnaseq]: Benayoun B A, Pollina E A, Singh P P, et al. Remodeling of epigenome and transcriptome landscapes with aging in mice reveals widespread induction of inflammatory responses[J]. Genome research, 2019, 29(4): 697-709.
			''')
		elif refcelltype.lower() == "scmca":
			Celltype_c1 = Celltype_module.add_comment('''本项目采用 SingleR[^SingleR] 包，基于 scmca[^scmca]参考数据集进行细胞类型注释。该方法通过计算单细胞参考表达谱数据集与待鉴定的细胞表达谱之间的相关性，将待鉴定细胞注释为与参考数据集中相关性最高的一种细胞类型。

报告中的数据集鉴定结果供参考，后续可根据文献已有相关基因对细胞群的特征加以描述和验证。
[^SingleR]: Aran D, Looney A P, Liu L, et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage[J]. Nature immunology, 2019, 20(2): 163-172.
[^scmca]: Xiaoping Han, Renying Wang, Yincong Zhou,et al. Mapping the Mouse Cell Atlas by Microwell-Seq[J]. Cell, 2018, 172(5):1091-1107.e17.
			''')

		Celltype_tsne0 = Celltype_module.add_comment("数据集方法细胞类型鉴定参考结果见：")
		Celltype_tsne1 = Celltype_module.add_comment("[{Celltype_link}]({Celltype_link})",Celltype_link=Celltype_link)
		
		Celltype_plots = list(glob('%s/*top.*celltyping_on_tsne*.png' %(celltype_folder[0])))
		Celltype_contents = ['细胞类型参考结果图。图片说明：细胞类型注释结果在 t-SNE 图上的展示，每种细胞类型以不同颜色区分。']
		Celltype_plots = list(zip(Celltype_plots, Celltype_contents))
		Celltype_plot = Celltype_module.add_plot(Celltype_plots, caption = '细胞类型参考结果图',description = '图片说明：细胞类型注释结果在 t-SNE 图上的展示，每种细胞类型以不同颜色区分。')
		Celltype_heatmap_link = list(glob('%s/*celltyping_heatmap.pdf' %(celltype_folder[0])))[0]
		Celltype_heatmap0 = Celltype_module.add_comment("细胞类型鉴定相关性热图结果见：")
		Celltype_heatmap1 = Celltype_module.add_comment("[{Celltype_heatmap_link}]({Celltype_heatmap_link})",Celltype_heatmap_link=Celltype_heatmap_link)
		Celltype_heatmap_plots = list(glob('%s/*celltyping_heatmap.png' %(celltype_folder[0])))
		Celltype_heatmap_contents = ['细胞类型鉴定相关性热图。图片说明：纵轴表示每一个待鉴定的细胞，横轴表示参考数据集中的细胞类型注释名称。颜色越红代表相关性值越大，表明待鉴定的细胞类型与参考数据集中的该细胞类型注释最为相似。']
		Celltype_heatmap_plots = list(zip(Celltype_heatmap_plots, Celltype_heatmap_contents))
		Celltype_heatmap_plot = Celltype_module.add_plot(Celltype_heatmap_plots, caption = '细胞类型鉴定相关性热图', description = '图片说明：每一行表示参考数据集中的细胞类型注释名称，每一列表示待鉴定的细胞。颜色越红表示相关性值越大，表明待鉴定的细胞类型与参考数据集中的该种细胞类型最为相似。')

		Celltype_table_link = list(glob('%s/*_simplified_celltyping_results.csv' %(celltype_folder[0])))[0]
		Celltype_table0 = Celltype_module.add_comment("细胞类型注释表格见：")
		Celltype_table1 = Celltype_module.add_comment("[{Celltype_table_link}]({Celltype_table_link})",Celltype_table_link=Celltype_table_link)
		Celltype_table = Celltype_module.add_table('%s/*_simplified_celltyping_results.csv' %(celltype_folder[0]), caption = '细胞类型注释结果表格', show_rows=10 , headerdict={'Barcode':'细胞Barcode','sampleid':'样本名称','celltype':'细胞类型注释','clusters':'细胞群','group':'样本分组'})

		celltype_link1 = get_url('%s/' %(celltype_folder[0]), '*_celltyping_statistics.xls')
		celltype_c0 = Celltype_module.add_comment("各细胞群中的原始细胞类型数目统计表格见：")
		celltype_c0 = Celltype_module.add_comment("{celltype_link1}",celltype_link1=celltype_link1)
		
	##########################差异表达基因及富集##########################
	if len(Diffexp_folder) > 0:
		diffexp = dict(config['parameter'])['diffexp']
		Diff_module = result_module.add_section('差异表达基因筛选', diffexp=diffexp)
		link = '%s' %(Diffexp_folder[0])
		diff_gene_c = Diff_module.add_comment('''详细结果见目录：
		[{link}]({link})''',link=link)
		Diffexp_t = Diff_module.add_comment("各分组差异基因数目统计表如下：")
		Diffexp_table = Diff_module.add_table('%s/diffexp_results_stat.xls' %(Diffexp_folder[0]), caption = '差异表达基因统计表',  headerdict={'Case':'实验组名','Control':'对照组名','Up_diff':'显著上调差异基因数量','Down_diff':'显著下调差异基因数量','Total_diff(pvalue<0.05&FoldChange>1.5)':'显著性差异基因总数量'})
		Diffexp_c = Diff_module.add_comment("差异显著基因结果示例：")
		difffiles = list(glob('%s/*diff*pval*xls' %(Diffexp_folder[0])))
		Diffexp_table = Diff_module.add_table(difffiles[0], caption = '差异显著基因结果表格',show_rows =10 , headerdict={'GeneID':'基因名','pvalue':'p值','pct.1':'表达该基因的细胞在比较组细胞群中的占比','pct.2':'表达该基因的细胞在对照组细胞群中的占比','padj':'校正后的p值','FoldChange':'差异倍数','log2FoldChange':'log2转化后的差异倍数','up_down':'基因上下调描述','ensemble_id':'ensemble基因ID','Dbxref_GeneID':' NCBI基因ID索引号','gene_type':'基因的类型','gene_description':'基因描述','GO_id':'Gene Ontology登录号','GO_term':'GO条目描述','pathway':'KEGG通路号','pathway_description':'KEGG通路描述'})

		topdiff_c = Diff_module.add_comment("将差异基因的差异倍数（FoldChange）从大到小排列，上下调各选取 25 个基因绘制热图：")
		topdiff_heatmap_plots = list(glob('%s/*_heatmap.png' %(Diffexp_folder[0])))
		topdiff_heatmap_plots = [ str(i) for i in topdiff_heatmap_plots ]
		topdiff_heatmap_prefix = [re.findall("top25_(.*)_heatmap.png",i)[0]  for i in topdiff_heatmap_plots]
		topdiff_heatmap_contents = [j+' 上下调 Top25 差异基因热图。图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足25个，则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。' for j in topdiff_heatmap_prefix]
		topdiff_heatmap_plots = list(zip(topdiff_heatmap_plots, topdiff_heatmap_contents))
		topdiff_heatmap_plot = Diff_module.add_plot(topdiff_heatmap_plots, caption = '上下调 Top25 差异基因热图', description = '图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足25个，则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。')

		enrichment_folder = list(glob(r"[0-9]*.enrichment"))
		ppi_folder = list(glob(r"[0-9]*.PPI"))
		if len(enrichment_folder)>0:
			##########################差异基因 GO 富集分析##########################
			enrich_module = result_module.add_section('差异基因富集分析')
			#定义函数获取比较组名称
			def get_diffgroup_name(target,tag):
				diff_group = []
				for f in Path(target).glob(tag):
					t = str(f).split('/')[-1]
					t = t.split('-diff')[0]
					diff_group.append(t)
				return diff_group
			
			
			link = '%s/GO_enrichment' %(enrichment_folder[0])
			diff_group = get_diffgroup_name('%s' %(Diffexp_folder[0]), '*diff*FC*.xls')

			diff_gene_GO = enrich_module.add_section('差异基因 GO 富集分析')
			diff_gene_GO_plot2 = diff_gene_GO.add_plot('enrich.png',caption='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式',content='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式。其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。',description='其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。')
			diff_gene_GO_c = diff_gene_GO.add_comment('''输出文件结果目录：
			[{link}]({link})''',link=link)
			
			diff_gene_GO_c2 = diff_gene_GO.add_comment("GO 富集分析结果示例：")
			diff_gene_GO_table1 = diff_gene_GO.add_table('%s/GO_enrichment/%s/enrichment-go-*-Total.xls' %(enrichment_folder[0], diff_group[0]),caption='富集分析结果表格',show_rows=10,headerdict={'id':'条目在Gene Ontology的登录号','term':'该条目的描述','category':'GO分类','ListHits':'该GO条目中差异基因数','ListTotal':'注释到GO的总差异基因数','PopHits':'注释到该条目中的所有基因数目','PopTotal':'注释到GO的总基因数','pval':'富集显著性p值','padj':'校正后的p值','Enrichment_score':'富集打分','Gene':'属于该条目的差异gene'})
			
			diff_gene_GO3_plots = list(glob('%s/GO_enrichment/*/GO.top.*.png'  %(enrichment_folder[0])))
			if len(diff_gene_GO3_plots) >0:
				diff_gene_GO_plot3_c = diff_gene_GO.add_comment('''GO 富集分析 top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log<sub>10</sub>pValue 由大到小排序的各 10 条）条形图展示如下：''')
				diff_gene_GO3_contents = ['GO 富集分析结果展示。图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>pValue。']*len(diff_gene_GO3_plots)
				diff_gene_GO3_plots = list(zip(diff_gene_GO3_plots, diff_gene_GO3_contents))
				diff_gene_GO_plot3 = diff_gene_GO.add_plot(diff_gene_GO3_plots,caption='GO 富集分析结果展示',description='图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>pValue。')
			
			if len(list(glob('%s/GO_enrichment/*/topGO_*.png'  %(enrichment_folder[0]))))>0:
				diff_gene_GO_fish = diff_gene_GO.add_comment('''使用 fisher 算法分别对样本间差异基因进行 CC，BP，MF 富集分析，并使用 topGO[^topGO] 对富集到的 Term 绘制有向无环图。topGO 有向无环图能直观展示差异表达基因富集的 GO 节点（Term）及其层级关系，是差异表达基因 GO 富集分析的结果图形化展示，分支代表的包含关系，从上至下所定义的功能描述范围越来越具体。
[^topGO]:Alexa A, Rahnenfuhrer J. topGO: enrichment analysis for gene ontology. R package version 2.8,2010.''')
				diff_gene_GO_plot4 = diff_gene_GO.add_plot('topGO.png',caption='差异基因topGO有向无环示例图展示',content='差异基因topGO有向无环示例图展示。图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。',description='''图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。''')
				
			
			diff_gene_GO5_plots = list(glob('%s/GO_enrichment/*/ALL_vs_DEG.GO.level2.stat.png' %(enrichment_folder[0])))
			if len(diff_gene_GO5_plots)>0:

				diff_gene_GO_plot5_c0 = diff_gene_GO.add_comment("根据功能分级，一般将 GO 分为三个层级，level1 包含三个条目：biological process、cellular component和molecular function，level2 包含 biological adhesion、cell 和 binding 等 64 个条目，level3 即为常规富集使用的数万个条目。从 level1到 level3 功能更具体，反之，更概括。")

				diff_gene_GO_plot5_c1 = diff_gene_GO.add_comment("差异基因和所有基因在 GO Level2 水平分布比较图如下：")

				diff_gene_GO5_contents = ['差异表达基因及所有基因在 GO Level2 水平分布比较图。图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。']*len(diff_gene_GO5_plots)
				diff_gene_GO5_plots = list(zip(diff_gene_GO5_plots, diff_gene_GO5_contents))
				diff_gene_GO_plot5 = diff_gene_GO.add_plot(diff_gene_GO5_plots,caption='差异表达基因及所有基因在 GO Level2 水平分布比较图',description='图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')
				
			diff_gene_GO6_plots = list(glob('%s/GO_enrichment/*/Up_vs_Down.GO.level2.stat.png' %(enrichment_folder[0])))
			if len(diff_gene_GO6_plots)>0:
				diff_gene_GO_plot6_c = diff_gene_GO.add_comment('''上调差异基因和下调差异基因在 GO Level2 水平分布比较图如下：''')
				diff_gene_GO6_plots = list(glob('%s/GO_enrichment/*/Up_vs_Down.GO.level2.stat.png' %(enrichment_folder[0])))
				diff_gene_GO6_contents = ['上调差异基因和下调差异基因在 GO Level2 水平分布比较图。图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比']*len(diff_gene_GO6_plots)
				diff_gene_GO6_plots = list(zip(diff_gene_GO6_plots, diff_gene_GO5_contents))
				diff_gene_GO_plot6 = diff_gene_GO.add_plot(diff_gene_GO6_plots,caption='上调差异基因和下调差异基因在 GO Level2 水平分布比较图',description='图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')



			##########################差异基因 KEGG 富集分析##########################
			link1 = '%s/KEGG_enrichment' %(enrichment_folder[0])
			diff_gene_KEGG = enrich_module.add_section('差异基因 KEGG 富集分析')
			diff_gene_KEGG_c0 = diff_gene_KEGG.add_comment('''差异基因 KEGG 富集分析结果：[{link1}]({link1})''',link1=link1)

			
			diff_gene_KEGG_top20_plots = list(glob('%s/KEGG_enrichment/*/KEGG.top.*.png' %(enrichment_folder[0])))
			if len(diff_gene_KEGG_top20_plots)>0:
				diff_gene_KEGG_c0 = diff_gene_KEGG.add_comment("KEGG 富集分析 top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序）气泡图如下：")
				diff_gene_KEGG_top20_contents = ['KEGG富集 top20 气泡图。图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。']*len(diff_gene_KEGG_top20_plots)
				diff_gene_KEGG_top20_plots = list(zip(diff_gene_KEGG_top20_plots, diff_gene_KEGG_top20_contents))
				diff_gene_KEGG_top20 = diff_gene_KEGG.add_plot(diff_gene_KEGG_top20_plots,caption='KEGG富集 top20 气泡图',description='图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。')
				
			diff_gene_KEGG_level2_plots = list(glob('%s/KEGG_enrichment/*/ALL_vs_DEG.KEGG_Classification.png' %(enrichment_folder[0])))
			if len(diff_gene_KEGG_level2_plots)>0:
				diff_gene_KEGG_levle2_c0 = diff_gene_KEGG.add_comment('''根据功能分级，通常将 KEGG 分为三个层级，level1 包含六个分类：Metabolism、Genetic Information Processing、Environmental Information Processing、Cellular Processes、Organismal Systems 和 Human Diseases（具体物种注释可能有删减）。level2 包含 Cell growth and death、Transcription 和 Development 等 44 个分类（具体物种注释可能有删减），level3 即为常规富集使用的数百个 Pathway，从 level1 到 level3 功能更具体，反之，更概括。''')

				diff_gene_KEGG_levle2_c1 = diff_gene_KEGG.add_comment("差异表达基因及所有基因在 KEGG Level2 水平分布比较图如下：")

				diff_gene_KEGG_level2_contents = ['差异表达基因及所有基因在 KEGG Level2 水平分布比较图。图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。']*len(diff_gene_KEGG_level2_plots)
				diff_gene_KEGG_level2_plots = list(zip(diff_gene_KEGG_level2_plots, diff_gene_KEGG_level2_contents))

				diff_gene_KEGG_level2 = diff_gene_KEGG.add_plot(diff_gene_KEGG_level2_plots,caption='差异表达基因及所有基因在 KEGG Level2 水平分布比较图',description='图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。')
				
			diff_gene_KEGG_updown_level2_plots = list(glob('%s/KEGG_enrichment/*/Up_vs_Down.KEGG_Classification.png' %(enrichment_folder[0])))
			if len(diff_gene_KEGG_updown_level2_plots)>0:
				diff_gene_KEGG_updown_level2_c = diff_gene_KEGG.add_comment('上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图如下：')
				diff_gene_KEGG_updown_level2_contents = ['上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图。图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。']*len(diff_gene_KEGG_updown_level2_plots)
				diff_gene_KEGG_updown_level2_plots = list(zip(diff_gene_KEGG_updown_level2_plots, diff_gene_KEGG_updown_level2_contents))
				diff_gene_KEGG_updown_level2 = diff_gene_KEGG.add_plot(diff_gene_KEGG_updown_level2_plots,caption='上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图',description='图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。')
				
				
				
		if len(ppi_folder)>0:
			##########################差异基因 PPI 富集分析##########################
			ppi_module = result_module.add_section('差异基因蛋白网络互作分析')
			ppi_link1 = '%s' %(ppi_folder[0])
			ppi_c0 = ppi_module.add_comment('''结果目录：[{ppi_link1}]({ppi_link1})''',ppi_link1=ppi_link1)
			ppi_c1 = ppi_module.add_comment('''差异基因互作关系结果文件：''')
			ppi_table = list(glob('%s/*_protein-protein-interaction.tsv' %(ppi_folder[0])))[0]
			ppi_table1 = ppi_module.add_table(ppi_table, caption='差异基因互作关系表',show_rows=10, headerdict={'stringId_A':'基因A的stringId','stringId_B':'基因B的stringId','preferredName_A':'A基因名','preferredName_B':'B基因名','ncbiTaxonId':'NCBI的物种分类标识符','score':'结合不同渠道证据的概率矫正后的综合得分','nscore':'邻近得分（从基因间核苷酸计数计算）','fscore':'融合得分（来自其他物种的融合蛋白质）','pscore':'系统谱系的共存得分（源自相似的基因缺失/存在模式）','ascore':'共表达得分（源自由DNA阵列和类似技术测量相似的mRNA表达模式）','escore':'实验得分（来自实验数据，如亲和色谱）','dscore':'数据库得分（派生自各种数据库的精选数据）','tscore':'文本挖掘得分（源自摘要中同时存在的基因/蛋白质名称）'})
			
			ppi_c2 = ppi_module.add_comment("上下调 Top25 差异表达基因互作网络图展示如下")
			ppi_plots = list(glob('%s/*string_protein-protein-interaction.new_colors.png' %(ppi_folder[0])))
			ppi_plots = [ str(i) for i in ppi_plots ] 
			ppi_contents = ['差异基因互作网络图。图片说明：红色圆圈表示上调表达基因，绿色圆圈表示下调表达基因，节点之间的连线（或称为边）表示两蛋白之间具有相互作用，线的粗细表示相互作用关系的可靠性。']*len(ppi_plots)
			ppi_plots = list(zip(ppi_plots, ppi_contents))
			ppi_plot = ppi_module.add_plot(ppi_plots, caption = '差异基因互作网络图', description = '图片说明：红色圆圈表示上调表达基因，绿色圆圈表示下调表达基因，节点之间的连线（或称为边）表示两蛋白之间具有相互作用，线的粗细表示相互作用关系的可靠性。')

	# else:
		# ##########################差异表达基因及富集##########################
		# Diffexp_folder = list(glob(r"[0-9]*.Diffexp"))
		# if len(Diffexp_folder) > 0:
			# Diff_module = result_module.add_section('差异表达基因筛选')
			# link = '%s' %(Diffexp_folder[0])
			# diff_gene_c = Diff_module.add_comment('''详细结果见目录：
			# [{link}]({link})''',link=link)
			# Diffexp_t = Diff_module.add_comment("各分组差异基因数目统计表如下：")
			# Diffexp_table = Diff_module.add_table('%s/diffexp_results_stat.xls' %(Diffexp_folder[0]), caption = '差异表达基因统计表',  headerdict={'Case':'实验组名','Control':'对照组名','Up_diff':'显著上调差异基因数量','Down_diff':'显著下调差异基因数量','Total_diff(pvalue<0.05&FoldChange>1.5)':'显著性差异基因总数量'})
			# Diffexp_c = Diff_module.add_comment("差异显著基因结果示例：")
			# difffiles = list(glob('%s/*diff*pval*xls' %(Diffexp_folder[0])))
			# Diffexp_table = Diff_module.add_table(difffiles[0], caption = '差异显著基因结果表格',show_rows =10 , headerdict={'GeneID':'基因名','pvalue':'p值','pct.1':'表达该基因的细胞在比较组细胞群中的占比','pct.2':'表达该基因的细胞在对照组细胞群中的占比','padj':'校正后的p值','FoldChange':'差异倍数','log2FoldChange':'log2转化后的差异倍数','up_down':'基因上下调描述','ensemble_id':'ensemble基因ID','Dbxref_GeneID':' NCBI基因ID索引号','gene_type':'基因的类型','gene_description':'基因描述','GO_id':'Gene Ontology登录号','GO_term':'GO条目描述','pathway':'KEGG通路号','pathway_description':'KEGG通路描述'})

			# enrichment_folder = list(glob(r"[0-9]*.enrichment"))
			# if len(enrichment_folder)>0:
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
				
				
				# link = '%s/GO_enrichment' %(enrichment_folder[0])
				# diff_group = get_diffgroup_name('%s' %(Diffexp_folder[0]), '*diff*FC*.xls')

				# diff_gene_GO = enrich_module.add_section('差异基因 GO 富集分析')
				# diff_gene_GO_plot2 = diff_gene_GO.add_plot('enrich.png',caption='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式',content='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式。其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。',description='其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。')
				# diff_gene_GO_c = diff_gene_GO.add_comment('''输出文件结果目录：
				# [{link}]({link})''',link=link)
				
				# diff_gene_GO_c2 = diff_gene_GO.add_comment("GO 富集分析结果示例：")
				# diff_gene_GO_table1 = diff_gene_GO.add_table('%s/GO_enrichment/%s/enrichment-go-*-Total.xls' %(enrichment_folder[0], diff_group[0]),caption='富集分析结果表格',show_rows=10,headerdict={'id':'条目在Gene Ontology的登录号','term':'该条目的描述','category':'GO分类','ListHits':'该GO条目中差异基因数','ListTotal':'注释到GO的总差异基因数','PopHits':'注释到该条目中的所有基因数目','PopTotal':'注释到GO的总基因数','pval':'富集显著性p值','padj':'校正后的p值','Enrichment_score':'富集打分','Gene':'属于该条目的差异gene'})
				
				# diff_gene_GO3_plots = list(glob('%s/GO_enrichment/*/GO.top.*.png'  %(enrichment_folder[0])))
				# if len(diff_gene_GO3_plots) >0:
					# diff_gene_GO_plot3_c = diff_gene_GO.add_comment('''GO 富集分析 top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log<sub>10</sub>pValue 由大到小排序的各 10 条）条形图展示如下：''')

					# diff_gene_GO3_contents = ['GO 富集分析结果展示。图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>pValue。']*len(diff_gene_GO3_plots)
					# diff_gene_GO3_plots = list(zip(diff_gene_GO3_plots, diff_gene_GO3_contents))
					# diff_gene_GO_plot3 = diff_gene_GO.add_plot(diff_gene_GO3_plots,caption='GO 富集分析结果展示',description='图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>pValue。')
					
				# if len(list(glob('%s/GO_enrichment/*/topGO_*.png'  %(enrichment_folder[0]))))>0:
					# diff_gene_GO_fish = diff_gene_GO.add_comment('''使用 fisher 算法分别对样本间差异基因进行 CC，BP，MF 富集分析，并使用 topGO[^topGO] 对富集到的 Term 绘制有向无环图。topGO 有向无环图能直观展示差异表达基因富集的 GO 节点（Term）及其层级关系，是差异表达基因 GO 富集分析的结果图形化展示，分支代表的包含关系，从上至下所定义的功能描述范围越来越具体。
# [^topGO]:Alexa A, Rahnenfuhrer J. topGO: enrichment analysis for gene ontology. R package version 2.8,2010.''')
					# diff_gene_GO_plot4 = diff_gene_GO.add_plot('topGO.png',caption='差异基因topGO有向无环示例图展示',content='差异基因topGO有向无环示例图展示。图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。',description='''图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。''')
					
					
				
				# if len(diff_gene_GO_plot5_c1)>0:
					# diff_gene_GO_plot5_c0 = diff_gene_GO.add_comment("根据功能分级，一般将 GO 分为三个层级，level1 包含三个条目：biological process、cellular component和molecular function，level2 包含 biological adhesion、cell 和 binding 等 64 个条目，level3 即为常规富集使用的数万个条目。从 level1到 level3 功能更具体，反之，更概括。")

					# diff_gene_GO_plot5_c1 = diff_gene_GO.add_comment("差异基因和所有基因在 GO Level2 水平分布比较图如下：")
					# diff_gene_GO5_plots = list(glob('%s/GO_enrichment/*/ALL_vs_DEG.GO.level2.stat.png' %(enrichment_folder[0])))
					# diff_gene_GO5_contents = ['差异表达基因及所有基因在 GO Level2 水平分布比较图。图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。']*len(diff_gene_GO5_plots)
					# diff_gene_GO5_plots = list(zip(diff_gene_GO5_plots, diff_gene_GO5_contents))
					# diff_gene_GO_plot5 = diff_gene_GO.add_plot(diff_gene_GO5_plots,caption='差异表达基因及所有基因在 GO Level2 水平分布比较图',description='图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')
					# diff_gene_GO_plot6_c = diff_gene_GO.add_comment('''上调差异基因和下调差异基因在 GO Level2 水平分布比较图如下：''')
					# diff_gene_GO6_plots = list(glob('%s/GO_enrichment/*/Up_vs_Down.GO.level2.stat.png' %(enrichment_folder[0])))
					# diff_gene_GO6_contents = ['上调差异基因和下调差异基因在 GO Level2 水平分布比较图。图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比']*len(diff_gene_GO6_plots)
					# diff_gene_GO6_plots = list(zip(diff_gene_GO6_plots, diff_gene_GO5_contents))
					# diff_gene_GO_plot6 = diff_gene_GO.add_plot(diff_gene_GO6_plots,caption='上调差异基因和下调差异基因在 GO Level2 水平分布比较图',description='图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')



				# ##########################差异基因 KEGG 富集分析##########################
				# link1 = '%s/KEGG_enrichment' %(enrichment_folder[0])
				# diff_gene_KEGG = enrich_module.add_section('差异基因 KEGG 富集分析')
				# diff_gene_KEGG_c0 = diff_gene_KEGG.add_comment('''差异基因 KEGG 富集分析结果：[{link1}]({link1})''',link1=link1)
				
				# diff_gene_KEGG_top20_plots = list(glob('%s/KEGG_enrichment/*/KEGG.top.*.png' %(enrichment_folder[0])))
				# if len(diff_gene_KEGG_top20_plots)>0:
					# diff_gene_KEGG_c0 = diff_gene_KEGG.add_comment("KEGG 富集分析 top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序）气泡图如下：")

					# diff_gene_KEGG_top20_contents = ['KEGG富集 top20 气泡图。图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。']*len(diff_gene_KEGG_top20_plots)
					# diff_gene_KEGG_top20_plots = list(zip(diff_gene_KEGG_top20_plots, diff_gene_KEGG_top20_contents))
					# diff_gene_KEGG_top20 = diff_gene_KEGG.add_plot(diff_gene_KEGG_top20_plots,caption='KEGG富集 top20 气泡图',description='图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。')
					
				# diff_gene_KEGG_level2_plots = list(glob('%s/KEGG_enrichment/*/ALL_vs_DEG.KEGG_Classification.png' %(enrichment_folder[0])))
				# if len(diff_gene_KEGG_level2_plots)>0:
					# diff_gene_KEGG_levle2_c0 = diff_gene_KEGG.add_comment('''根据功能分级，通常将 KEGG 分为三个层级，level1 包含六个分类：Metabolism、Genetic Information Processing、Environmental Information Processing、Cellular Processes、Organismal Systems 和 Human Diseases（具体物种注释可能有删减）。level2 包含 Cell growth and death、Transcription 和 Development 等 44 个分类（具体物种注释可能有删减），level3 即为常规富集使用的数百个 Pathway，从 level1 到 level3 功能更具体，反之，更概括。''')

					# diff_gene_KEGG_levle2_c1 = diff_gene_KEGG.add_comment("差异表达基因及所有基因在 KEGG Level2 水平分布比较图如下：")

					# diff_gene_KEGG_level2_contents = ['差异表达基因及所有基因在 KEGG Level2 水平分布比较图。图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。']*len(diff_gene_KEGG_level2_plots)
					# diff_gene_KEGG_level2_plots = list(zip(diff_gene_KEGG_level2_plots, diff_gene_KEGG_level2_contents))

					# diff_gene_KEGG_level2 = diff_gene_KEGG.add_plot(diff_gene_KEGG_level2_plots,caption='差异表达基因及所有基因在 KEGG Level2 水平分布比较图',description='图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。')
					# diff_gene_KEGG_updown_level2_c = diff_gene_KEGG.add_comment('上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图如下：')
					# diff_gene_KEGG_updown_level2_plots = list(glob('%s/KEGG_enrichment/*/Up_vs_Down.KEGG_Classification.png' %(enrichment_folder[0])))
					# diff_gene_KEGG_updown_level2_contents = ['上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图。图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。']*len(diff_gene_KEGG_updown_level2_plots)
					# diff_gene_KEGG_updown_level2_plots = list(zip(diff_gene_KEGG_updown_level2_plots, diff_gene_KEGG_updown_level2_contents))
					# diff_gene_KEGG_updown_level2 = diff_gene_KEGG.add_plot(diff_gene_KEGG_updown_level2_plots,caption='上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图',description='图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。')


	##########################附录##########################
	appendix_module = report.add_section('附录')

	# ##########################loupe cell browser##########################
	# link0 = 'supplemental_material/Loupe-Cell-Browser使用教程.html'
	# supply_section_loupe = appendix_module.add_section('Loupe-Cell-Browser使用教程',link0=link0)
	
	
	# ##########################实验技术方法说明##########################
	# link1 = 'supplemental_material/欧易生物单细胞转录组实验技术方法说明_中文.pdf'
	# link2 = 'supplemental_material/欧易生物单细胞转录组实验技术方法说明_英文.pdf'
	# supply_section_bio = appendix_module.add_section('实验技术方法说明',link1=link1,link2=link2)

	# ##########################生信分析方法说明##########################
	# link3 = 'supplemental_material/欧易生物单细胞转录组生信分析方法_中文.pdf'
	# link4 = 'supplemental_material/欧易生物单细胞转录组生信分析方法_英文.pdf'
	# supply_section_analysis = appendix_module.add_section('生信分析方法说明',link3=link3,link4=link4)

	# ##########################FAQ##########################
	# faq = appendix_module.add_section('常见问题FAQ',link_refdiff=link_refdiff)

	##########################软件及数据库信息##########################
	database = pd.DataFrame([['Genome Database',dict(config['db'])['genome']],['GO Database','http://geneontology.org/'],['KEGG Database','http://www.genome.jp/kegg/']])
	database.columns = ['使用数据库', '网页链接']
	supply_section_database = appendix_module.add_section('数据库信息',description='')
	supply_section_database_table = supply_section_database.add_table(database,caption='数据库信息')
	supply_section_software = appendix_module.add_section('数据分析软件')
	software_table = appendix_module.add_table('scRNA_software.txt', caption = '使用软件及版本')

	other_module = report.add_section('申明')
	
	
	
	# Generate Report HTML
	report.write_to('report.html')



if __name__ == "__main__":
	scrna_report()
