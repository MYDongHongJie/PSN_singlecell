import shutil
import argparse
import os
from oebio.report import Report
from oebio.app import *
from configparser import ConfigParser
import click
import subprocess
import sys
import re,yaml
import time 
import glob
args  = argparse.ArgumentParser(description=' vdj report')
args.add_argument('-i','--input', type=str,default = "result/VDJ_aggr",
			help='the directory of your program. ')
args.add_argument('-c','--configfile', type=str, default='config/config.yaml',
			help='the config file which contain your program information.')
args.add_argument('-t','--vdjtype',type=str,
			help='Immune repertoire type: TCR or BCR')
args  = args.parse_args()

input = os.path.abspath(args.input)
configfile = os.path.abspath(args.configfile)
vdjtype = args.vdjtype


def vdj_report(input, configfile, vdjtype):

    """ python /public/scRNA_works/works/donghongjie/project/test/vdj/scripts/VDJ_scripts/VDJ_report.py -i result/VDJ_aggr -c config/config.yaml -t TCR"""
    cfg = open(configfile, 'r', encoding='utf-8').read()
    config = yaml.load(cfg,Loader=yaml.FullLoader)
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.abspath(input + "/" +  vdjtype)
    program_num = dict(config['report'])['Task_Num'].split("-")[0]
    orig_path = os.getcwd()
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    header_info = {key: config["report"][key] for key in config["report"] if
            key in ["Customer", "Species", "Executor", "Sample_Num", "Task_Num"]}
    header_info["项目编号"] = program_num
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["实验物种"] = header_info.pop("Species")
    header_info["执行编号"] = header_info.pop("Executor")
    header_info["样本数目"] = header_info.pop("Sample_Num")
    header_info["任务单号"] = header_info.pop("Task_Num")


	# Report
    os.chdir(program_path)
    log_file = open("%s/report_log.txt" % (program_path), "w")
    report_time = time.strftime("%Y_%m_%d")
    outdir = f"{program_num}_{vdjtype}_Report_{report_time}" 
    if os.path.exists( f'{outdir}' ):
        shutil.rmtree(f'{outdir}', ignore_errors=True)
    else:
        os.makedirs(f'{outdir}')
        print("Create Report...")
        log_file.write("Create Report..." + "\n")

    num = 0
	# Cellranger
    if os.path.exists("%s/../../cellranger"  % (program_path)):
        num += 1
        os.makedirs(f"{outdir}/{num}.CellRanger/")
        name = [str(i).split('/')[-2] for i in glob.glob("../../cellranger/*/outs")]
		
        # for j in name:
        #     os.makedirs("%s/%s.CellRanger/%s" % (outdir, num, j))
        #     subprocess.call('ln -s %s/Cellranger/%s/outs/* %s/%s.CellRanger/%s ' % (program_path, j, outdir, num, j), shell=True)
        #     if os.path.exists("%s/%s.CellRanger/%s/*.bam*" % (outdir, num, j)):
        #         subprocess.call('rm %s/%s.CellRanger/%s/*.bam*' % (outdir, num, j), shell=True)
        #if len(name) >1:
        #    if os.path.exists("%s/Cellranger/aggr"  % (program_path)):
        #        os.makedirs("%s/%s.CellRanger/aggr" % (outdir, num))
        #        subprocess.call('ln -s %s/Cellranger/aggr/*/outs/* %s/%s.CellRanger/aggr ' % (program_path, outdir, num), shell=True)
        #        subprocess.call('rm  %s/%s.CellRanger/aggr/aggregation.csv ' % ( outdir, num), shell=True)
        #    else:
        #        print("Can not find CellRanger Aggr results!!!!!!!!!!!!!!!!")
        #        log_file.write("Can not find CellRanger Aggr results!!!!!!!!!!!!!!!!" + "\n")
        #        sys.exit()
        subprocess.call('cp %s/../../cellranger/*%s.png %s/%s.CellRanger/ ' % (program_path,vdjtype, outdir, num), shell=True)
    else:
        print("Can not find CellRanger results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find CellRanger results!!!!!!!!!!!!!!!!" + "\n")
        name = [str(i).split('/')[-1].split('.')[-3] for i in glob.glob(os.path.join("Clonity", "*qstat.pdf"))]

	# clonotypes
    if os.path.exists("%s/Clonotypes"  % (program_path)):
        num += 1
        os.makedirs("%s/%s.Clonotypes" % (outdir, num))
        subprocess.call('cp  %s/Clonotypes/* %s/%s.Clonotypes/' % (program_path, outdir, num), shell=True) ########################
		
    else:
        print("Can not find clonotypes results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find clonotypes results!!!!!!!!!!!!!!!!!" + "\n")
		
	# diversity
    if len(name) >1:
        if os.path.exists("%s/Diversity"  % (program_path)):
            num += 1
            os.makedirs("%s/%s.Diversity" % (outdir, num))
            subprocess.call('cp %s/Diversity/* %s/%s.Diversity/' % (program_path, outdir, num), shell=True)
        else:
            print("Can not find diversity results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find diversity results!!!!!!!!!!!!!!!!" + "\n")
      # summary
        if os.path.exists("%s/summary"  % (program_path)):
            num += 1
            os.makedirs("%s/%s.summary" % (outdir, num))
            subprocess.call('cp %s/summary/* %s/%s.summary/' % (program_path, outdir, num), shell=True)
        else:
            print("Can not find summary results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find summary results!!!!!!!!!!!!!!!!" + "\n")
		
	# Gene_usage
    if os.path.exists("%s/Gene_usage"  % (program_path)):
        num += 1
        os.makedirs("%s/%s.Gene_usage" % (outdir, num))
        subprocess.call('cp %s/Gene_usage/* %s/%s.Gene_usage/' % (program_path, outdir, num), shell=True)
    else:
        print("Can not find Gene_usage results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find Gene_usage results!!!!!!!!!!!!!!!!" + "\n")
		
	# Clonity
    if os.path.exists("%s/Clonity"  % (program_path)):
        num += 1
        os.makedirs("%s/%s.Clonity" % (outdir, num))
        subprocess.call('cp %s/Clonity/* %s/%s.Clonity/' % (program_path, outdir, num), shell=True)
    else:
        print("Can not find Clonity results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find Clonity results!!!!!!!!!!!!!!!!" + "\n")

	# Spectratyping
    if os.path.exists("%s/Spectratyping"  % (program_path)):
        num += 1
        os.makedirs("%s/%s.Spectratyping" % (outdir, num))
        subprocess.call('cp %s/Spectratyping/* %s/%s.Spectratyping/' % (program_path, outdir, num), shell=True)
    else:
        print("Can not find Spectratyping results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find Spectratyping results!!!!!!!!!!!!!!!!" + "\n")
  
    #conjoint
    if os.path.exists("%s/vdj_scRNA_Conjoint_analysis"  % (program_path)):
        num += 1
        os.makedirs("%s/%s.vdj_scRNA_Conjoint_analysis" % (outdir, num))
        subprocess.call('cp  -r %s/vdj_scRNA_Conjoint_analysis/* %s/%s.vdj_scRNA_Conjoint_analysis/' % (program_path, outdir, num), shell=True)
        #subprocess.call('rm %s/*.vdj_scRNA_Conjoint_analysis/*h5seurat'  % (outdir), shell=True)
    else:
        print("Can not find vdj_scRNA_Conjoint_analysis results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find vdj_scRNA_Conjoint_analysis results!!!!!!!!!!!!!!!!" + "\n")

  ### supplemental_material
    num += 1
    os.makedirs("%s/%s.supplemental_material" % (outdir,num))
    subprocess.call('cp -r %s/supplemental_material/Loupe_V\(D\)J_Browser使用手册.doc %s/%s.supplemental_material' %(scriptdir,outdir, num), shell=True)

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
    report = Report(f"10X 单细胞 V(D)J 测序结题报告-{vdjtype}", title = f"10X 单细胞 V(D)J 测序结题报告-{vdjtype}", header_info = dict(header_info), oe_welcome='''
####  感谢您选择欧易生物！

**关于网页报告使用有以下几点提示:**

1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节，便于您快速查看相应内容。
2. **报告正文**中的**图片**均可以点击后进行**放大查看**，且图片点击后可以包含更多的细节信息，比如左上角会显示具体的**图片数量**，右上角有相关的**图片工具**，**键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张，更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**，当表格多于20行的时候，表格会嵌入到网页之中，可以利用滚动栏进行调整设置，同时会在表格上方增加搜索设置，便于快速查询表格信息。
5. 在报告中浏览，需要返回顶部的时候，可以使用网页**右下角的白色箭头标识**快速返回顶部。
6. 本提示可以点击右上方的不再显示，则后续打开网页均不会显示本提示。
    ''')
    report.add_yaml_config('scriptdir,report.yaml')
    report.add_yaml_config(os.path.join(scriptdir, 'report.yaml'))
    os.environ['OEBIO'] = os.path.join(scriptdir,"pic")

    #os.environ['oeweb_register_token'] = 'oearray'
	#######################################项目概况#######################################
    #统计信息
    summ = pd.read_csv(glob.glob("*Clonotypes/merged_vdj_contig_annotation.xls")[0], index_col=0, sep='\t')
    sample_num= len(summ["sampleid"].value_counts())
    max_barcodes=max(summ["sampleid"].value_counts())
    min_barcodes=min(summ["sampleid"].value_counts())
    project_module = report.add_section('项目概况')
    #register_info = oeweb_register(project_id=program_num, target_url='https://cloud.oebiotech.cn/task/category/scrna', note='本报告包含项目基本分析内容，如要快速做其它分析调整，如细胞表达可视化等')
    #if register_info:
    #    project_module.add_comment(register_info)

	##########################项目摘要、快捷链接##########################	
    #link_refdiff = f"{num}.supplemental_material/10X 单细胞 V(D)J 测序常见问题(FAQ).pdf"
    if sample_num == 1 :
        project_module.add_section('项目摘要_样本数目为1',sample_num=sample_num,min_barcodes=min_barcodes)
    else:
        project_module.add_section('项目摘要',sample_num=sample_num,max_barcodes=max_barcodes,min_barcodes=min_barcodes)

    link0 = f"{num}.supplemental_material/Loupe-Browser使用教程.html"

    if os.path.exists("%s/vdj_scRNA_Conjoint_analysis"  % (program_path)):
        project_module.add_section('项目快捷链接_有联合')
    else: 
        project_module.add_section('项目快捷链接_无联合')

	
	##########################建库测序流程##########################
    library_module = report.add_section('实验技术流程')
    library_module_c = library_module.add_comment('实验建库测序流程图如下：')
    library_module_plot = library_module.add_plot('vdj_Experiment.png',caption ="10x Genomics 单细胞 5' 端转录组及免疫组库测序实验技术流程",content="10x Genomics 单细胞 5' 端转录组及免疫组库测序实验技术流程")
    library_module_plot2 = library_module.add_plot('GEMs_cDNA.png',caption ="5' 端转录组文库构建",content="5' 端转录组文库构建")
    library_module_plot3 = library_module.add_plot('VDJ_cDNA.png',caption ='TCR/BCR 文库构建',content='TCR/BCR 文库构建')
 
	##########################生信分析流程##########################
    analysis_module = report.add_section('生信分析流程')
    analysis_module_c = analysis_module.add_comment('生物信息分析流程图如下：')
    analysis_module_plot = analysis_module.add_plot('vdj_Pipeline.png', caption = '生物信息分析流程',content='生物信息分析流程图')
    	
	##########################项目分析结果##########################
    result_module = report.add_section('项目分析结果')

	##########################V(D)J 序列质控##########################
    QC_1 = result_module.add_section('V(D)J 序列质控')
    QC_1_plots = sorted(list(glob.glob('1.CellRanger/*.png')),reverse=False)
    QC_1_contents = ['样本细胞质量统计结果']*len(QC_1_plots)
    QC_1_plots = list(zip(QC_1_plots, QC_1_contents))
    QC1_plot = QC_1.add_plot(QC_1_plots, caption = '样本细胞质量统计结果',content= '样本细胞质量统计结果')
    #QC1_table = QC_1.add_table('vdj_qc_stat.txt', caption = '样本细胞质量统计结果说明')
	
	##########################clonotypes 鉴定##########################
    clonotypes = result_module.add_section('clonotypes 鉴定')
    clonotypes.add_table("*Clonotypes/merged_vdj_contig_annotation.xls", 
           headerdict={'barcode':'细胞名称',
                          'sampleid':'样品名称',
                          'cdr3s_aa':'cdr3s 氨基酸序列',
                          'cdr3s_nt':'cdr3s 核苷酸序列',
                          'segments':'V(D)J 基因片段位置(基因:起始位置:终止位置)',
                          'vdj_read_count':'轻链和重链分别比对到 V(D)J 基因 reads 数目',
                          'vdj_umi_count':'轻链和重链分别比对到 V(D)J 基因 umi 数目',
                          'umis4cell':'每个细胞中所有 V(D)J 基因 umi 数',
                          'reads4cell':'每个细胞中所有 V(D)J 基因 reads 数',
                          'n':'同种克隆型细胞数目',
                          'clonotype_id':'克隆型 ID ,数值为按 n (同种克隆型细胞数目)从大到小排序的 rank 值',
                          'is_paired':'是否是配对链'
                        },
           show_rows=10 ,caption = "样本 clonotype 注释表")
    clonotypes.add_section("clonotypes 丰度分析")
    clonotypes.add_plot("clonotypes_abundances.png",caption = "clonotypes 丰度示例图",description = "图片说明：横坐标为按照细胞数从大到小排序的 clonotype 类型，纵坐标可以是每种 clonotype 所对应的细胞数。图中从左至右为依此为各的 clonotypes 丰度分布。该图片可以通过 Loupe V(D)J Browser 软件导出，Loupe V(D)J Browser 使用教程见附录。")
    
	
    ######################### 免疫组库多样性分析 ###########################
    diversity = result_module.add_section('免疫组库多样性分析')
    if len(name) >1:
        diversity_rare = diversity.add_section("稀释曲线")
        diversity_rare.add_plot("*Diversity/sample.rarefaction.strict.png",caption = "样本稀释曲线",description = "图片说明：图中每一条曲线代表一个样本，横坐标为随机抽样的深度（即抽样的序列数），纵坐标为指数数值。当曲线随着抽取序列数的增加而趋于平缓时，说明样本的测序数据量合理，更多的数据量只会产生少量新的 clonotypes，反之则表明继续测序还可能产生较多新的 clonotypes。实线和虚线分别标记稀疏曲线的插值和外推区域，点标记精确的样本大小和多样性。阴影区域标记 95％ 置信区间。")
        diversity_rank = diversity.add_section("Rank-abundance 分析")
        diversity_rank.add_plot("*summary/clonotypes_count_rank_aboundance_groupby_Sample.png",caption = "Rank-abundance 曲线",description = "图片说明： 各样本的 Rank-abundance 曲线，横坐标 clonotypes 按照其所包含的数目从多到少进行排序，例如“100”代表样本中丰度排名第 100 位的 clonotypes；纵坐标为该 clonotypes 的丰度。")
        if len(name) <= 5 :
            diversity_venn = diversity.add_section("免疫组库多样本 clonotypes 分析")
            diversity_venn.add_plot("*Diversity/sample.join.aa.venn.png",caption = "免疫组库多样本 clonotypes 的 Venn 图")
            
        diversity2 = diversity.add_section('多样性指数')
        diversity_d50 = diversity2.add_section("D50 指数")
        diversity_d50.add_plot("*Diversity/diversity.d50_groupby_Sample.png",caption = "D50 指数柱状图")
        diversity_chao1 = diversity2.add_section("Chao1 指数")
        diversity_chao1.add_plot("*Diversity/diversity.chao1_index_groupby_Sample.png",caption = "Chao1 指数柱状图")
        diversity_hill = diversity2.add_section("Hill 指数")
        diversity_hill.add_plot("*Diversity/diversity.hill_index_groupby_Sample.png",caption = "Hill指数曲线图")

    diversity3 = diversity.add_section('V/J 基因使用率丰度分析')
    diversity3.add_plot("vj_gene_usage_frequency.png",caption="V/J 基因使用率丰度示例图",description = "图片说明：横坐标可根据需求选择V、J、V/J paired，纵坐标可根据需求选择细胞支持丰度或该基因丰度占比进行展示。该图片可以通过 Loupe V(D)J Browser 软件导出，Loupe V(D)J Browser 使用教程见附录。")
    diversity3.add_comment("通过研究不同样本的 V-J 基因组合的丰度可以进一步关注免疫组库中特异基因对的变化，在不同样本免疫组库中 V-J 基因对可直接反应出 CDR3 或免疫组库的变化。通过对 V-J 基因组合的统计可进一步进行后续的免疫学研究，例如，通过查看得到某一免疫时期样本高丰度（高细胞支持）的 V-J 基因组合或比较两个不同免疫时期高丰度 V-J 基因组合对得到特异性表达的免疫基因。使用 Loupe V(D)J Browser 也对每个样本 Cell ranger 分析结果进一步进行 V-J 基因组合进行解析。解析结果可通过 Loupe V(D)J Browser 进行导出。")
    diversity3.add_plot("vj_gene_usage_frequency_heatmap.png",caption="V-J 基因组合使用频率示例热图",description = "图片说明:横坐标为 V 基因，纵坐标为 J 基因，横纵坐标对应的值表示 V-J 基因对对应的 clonotypes 数,颜色越深，代表该 V-J 基因组合在该样本中所占比例越高。该图片可以通过 Loupe V(D)J Browser 软件导出，Loupe V(D)J Browser 使用教程见附录。")
    diversity3.add_comment("V-J 基因组合频率 Circos 图：")
    diversity3.add_plot("*Gene_usage/*fancyvj.wt.png",caption=" V/J 基因组合 Circos 图",description = "图片说明：上图是样本的 V-J 基因对使用率 Circos 图展示结果。每个圈图中上部分为 V 基因类型，下部分为 J 基因类型。每个颜色块代表一种基因，颜色块弧度越宽，代表该基因对应的 clonotypes 数目越多。色块间的弧线代表一种 V-J 基因组合方式，弧线越宽代表该种 V-J 组合的 clonotypes 的数目越多。")

    if len(name) >1:
        diversity3.add_comment("V-J 基因组合使用频率的三维柱状图:")
        diversity3.add_plot("*Gene_usage/*_3D_barplot.png",caption = "V-J 基因组合三维柱状图",description = "图片说明：各样本 V-J 基因组合 clonotypes 的频率统计三维柱状图，X轴为 V 基因，Y 轴为 J 基因，Z 轴为对应的 clonotype 的数目。")

    diversity4 = diversity.add_section('单链 CDR3 多样性分析')
    diversity4.add_plot("*Clonity/*qstat.png",caption = "甜甜圈图", description = """图片说明：各样本的甜甜圈图展示。从内到外，第一层 ("set") 展示了根据 Clonotypes/*CR.xls 第一列 count 数目分别统计单链 CDR3 umi 数目 =1、=2、>=3 三个区间内所有 frequency 的和；第二层 ("quantile") 展示了umi 数目>= 3 的单链 CDR3 按 frequency 从大到小排列，按顺序分成 5 等份 (Q1-Q5) ，每个等份内所有 frequency 的和；第三层展示了根据 frequency 从大到小排取 top5 单链 CDR3 氨基酸序列。""")
    
    diversity5 = diversity.add_section('CDR3 特征分析')
    diversity_type = diversity5.add_section('spectratype 分析')
    diversity_type.add_plot('*Spectratyping/*fancyspectra.png',caption = "CDR3 氨基酸序列长度分布图",description= "图片说明：各样本的在不同 CDR3 氨基酸序列长度下 CDR3 氨基酸序列类型的分布直方图。")
    diversity_type.add_plot('*Spectratyping/*spectraV.wt.png',caption = "不同 CDR3 长度下 V 基因组成分布图",description= "图片说明：各样本的在不同 V 基因对应的 CDR3 氨基酸序列类型的分布直方图。")
    diversity_type.add_comment('CDR3 区的 indel 序列也可通过 Loupe V(D)J Browser 进行查看：')
    diversity_type.add_plot('cd3_indel.png',caption = "不同 CDR3 长度下 V 基因组成分布图",description= "图片说明：各样本的在不同 V 基因对应的 CDR3 氨基酸序列类型的分布直方图。该图片可以通过 Loupe V(D)J Browser 软件导出，Loupe V(D)J Browser 使用教程见附录。")

    #########################免疫组库和单细胞转录组联合分析 ###########################
    if glob.glob("*vdj_scRNA_Conjoint_analysis"):
        conjoin = result_module.add_section('免疫组库和单细胞转录组联合分析')
        clonotype = conjoin.add_section("clonotype 分布情况")
        clonotype.add_plot("*vdj_scRNA_Conjoint_analysis/clonotype.png",caption = "clonotype 映射结果图",description = "图片说明：横纵坐标分别代表降维的第一和第二主成分。")

        top10_clonotype = conjoin.add_section("top10 clonotype 分布情况")
        top10_clonotype.add_plot("*vdj_scRNA_Conjoint_analysis/topn_clonotype_id/top_10_clonotype_id.png",caption = "top10 clonotype 映射结果图",description = "图片说明：横纵坐标分别代表降维的第一和第二主成分，不同 clonotype 以不同颜色区分。")
        if list(glob.glob("*vdj_scRNA_Conjoint_analysis/topn_clonotype_id/top_10_splitby_sampleid_clonotype_id.png")):
            top10_clonotype.add_plot("*vdj_scRNA_Conjoint_analysis/topn_clonotype_id/top_10_splitby_sampleid_clonotype_id.png",caption = "多样本 clonotype 分面展示图",description = "图片说明：横纵坐标分别代表降维的第一和第二主成分，不同 clonotype 以不同颜色区分。")

        stat=conjoin.add_section("clonotype 柱状统计图") 
        stat.add_plot("*vdj_scRNA_Conjoint_analysis/stat/groupby_clusters_splitby_sampleid_summary_clonotype_id_freq_plot.png",caption = "clonotype 数目柱状图",description = "图片说明：每个样品中每个细胞群中去重后 clonotype 的数目。")

        def figure_type(n):
            if n <= 5 :
                return "韦恩图"
            elif 5 < n <= 15:
                return "UpSet图"
            else:
                return "花瓣图"
  
        if list(glob.glob("*vdj_scRNA_Conjoint_analysis/venn/sampleid/VennGraph.png")) :
            venn=conjoin.add_section("clonotype 共享情况")    
            venn.add_plot("*vdj_scRNA_Conjoint_analysis/venn/sampleid/VennGraph.png",caption = f"clonotype 的共享{figure_type(len(name))}",description = f"图片说明：{figure_type(len(name))}展示了每个样品中 clonotype 的共享情况。")


        clonotype_id_num=conjoin.add_section("clonotype 丰度映射情况")
        clonotype_id_num.add_plot("*vdj_scRNA_Conjoint_analysis/clonotype_id_num/clonotype_id_num.png",caption = "clonotype 丰度映射图",description = "图片说明：横纵坐标分别代表降维的第一和第二主成分，点的颜色代表不同 clonotype 丰度情况。>=3表示有不少于三个细胞为同一种 clonotype。")

	##########################附录##########################
    appendix_module = report.add_section('附录')

	##########################loupe browser##########################
    link0 = f"{num}.supplemental_material/Loupe_V(D)J_Browser使用手册.doc"
    supply_section_loupe = appendix_module.add_section('Loupe V(D)J Browser 使用教程',link0=link0)
	
	##########################实验技术方法说明##########################
	#link1 = 'supplemental_material/欧易生物单细胞转录组实验技术方法说明_中文.pdf'
	#link2 = 'supplemental_material/欧易生物单细胞转录组实验技术方法说明_英文.pdf'
	#supply_section_bio = appendix_module.add_section('实验技术方法说明',link1=link1,link2=link2)

	##########################生信分析方法说明##########################
	#link3 = 'supplemental_material/欧易生物单细胞转录组生信分析方法_中文.pdf'
	#link4 = 'supplemental_material/欧易生物单细胞转录组生信分析方法_英文.pdf'
	#supply_section_analysis = appendix_module.add_section('生信分析方法说明',link3=link3,link4=link4)

	##########################FAQ##########################
	#faq = appendix_module.add_section('常见问题FAQ',link_refdiff=link_refdiff)

	##########################AI修图使用说明##########################
	#ai_link = 'supplemental_material/AI使用说明.html'
	#ai = appendix_module.add_section('AI修图使用说明',ai_link=ai_link)

	##########################软件及数据库信息##########################
	#database = pd.DataFrame([['Genome Database',dict(config['db'])['genome']],['GO Database','http://geneontology.org/'],['KEGG Database','http://www.genome.jp/kegg/']])
	#database.columns = ['使用数据库', '网页链接']
	#supply_section_database = appendix_module.add_section('数据库信息',description='')
	#supply_section_database_table = supply_section_database.add_table(database,caption='数据库信息')
	#supply_section_software = appendix_module.add_section('数据分析软件')
	#software_table = appendix_module.add_table('scRNA_software.txt', caption = '使用软件及版本')

    other_module = report.add_section('申明')
	
	
	
	# Generate Report HTML
    report.write_to('分析说明.html',  zip_report_name=f'{outdir}.zip')
    log_file.close()

if __name__ == "__main__":
    vdj_report(input, configfile, vdjtype)
