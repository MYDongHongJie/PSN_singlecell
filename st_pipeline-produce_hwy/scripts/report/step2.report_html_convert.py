#  -*- coding: UTF-8 -*
"""
Project :convert yaml to html
Author  :Xiufeng.Yang
Contact : xiufeng.yang@oebiotech.com
File   : report_oebio.py
IDE    : PyCharm
Time   : 2020-04-20
Desc   : produce 10x visium html report by using oebio
"""

import click
from oebio.report import oeweb_register,Report
from oebio.app import *
from collections import OrderedDict
import sys,shutil,os,click,subprocess,sys,re,time,yaml,re


def ordered_yaml_load(stream, Loader=yaml.SafeLoader,
                      object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def _construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        _construct_mapping)
    return yaml.load(stream, OrderedLoader)

@click.command()
@click.option('-d', '--outdir', prompt='the report directory',
              help='the directory of your report files ')
@click.option('-c', '--configfile', prompt='config yaml file. [default:config/config.yaml]',
              default='config/config.yaml',
              help='the config file which contain your program information.')
def convert_html(outdir, configfile):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    # 项目摘要统计数据====================================================================================================
    with open(configfile) as f:
        config = ordered_yaml_load(f)
    # file = open(configfile, 'r', encoding='utf-8')
    os.chdir(outdir)
    # cfg = file.read()
    # config = yaml.load(cfg)
    QC_folder = list(glob(r"[0-9]*.Count_QC"))[0]
    summ = pd.read_csv('%s/statitics_for_QC.xls' % (QC_folder), index_col=0, sep='\t')
    sample_num = config['report']['Sample_Num']  # 样本数量
    slices_num = len(summ)  # 切片数量
    min_spot = min(summ.Total_Spots_QC)
    max_spot = max(summ.Total_Spots_QC)
    min_umi = format(min(summ.mean_nCount_Spatial_QC.astype(float)), '.0f')
    max_umi = format(max(summ.mean_nCount_Spatial_QC.astype(float)), '.0f')
    min_gene = format(min(summ.mean_nFeature_Spatial_QC.astype(float)), '.0f')
    max_gene = format(max(summ.mean_nFeature_Spatial_QC.astype(float)), '.0f')
    library_type = config['params']['library_type']
    MT_flie =  os.path.isfile(os.path.join(config['database']['reference'],'MT_genelist.gmt'))
    if (library_type == "fresh" and MT_flie == True) or (library_type == "cytassist" and config['database']['reference'] == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A"):
        min_mito = format(min(summ['mean_percent.mito_QC'].astype(float) * 100), '.2f')
        max_mito = format(max(summ['mean_percent.mito_QC'].astype(float) * 100), '.2f')

    Clustering_folder = list(glob(r"[0-9]*.Clustering"))[0]
    Marker_folder = list(glob(r"[0-9]*.Marker"))[0]
    Marker_file = list(glob('%s/top10_markers_for_each_cluster_anno.xls' % (Marker_folder)))[0]
    Marker_results = pd.read_csv(Marker_file, sep='\t')
    cluster_num = max(Marker_results['cluster'])  # 获取聚类群数
    Diffexp_folder = list(glob(r"[0-9]*.Diffexp"))

    # Define header info
    header_info = {key: config["report"][key] for key in config["report"] if
                   key in ["Project_Num", "Customer", "Species", "Sample_Num","Task_Num"]}

    header_info["订单编号"] = header_info.pop("Project_Num")
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["物种信息"] = header_info.pop("Species")
    header_info["样本数量"] = header_info.pop("Sample_Num")
    if library_type == "fresh":
        header_info["文库类型"] = "Fresh"
    elif library_type == "ffpe":
        header_info["文库类型"] = "FFPE"
    elif library_type == "cytassist":
        header_info["文库类型"] = "CytAssist"
    header_info["任务单号"] = header_info.pop("Task_Num")
    title = "10x Visium空间转录组项目结题报告"
    report = Report(title, title=title, header_info=header_info, oe_welcome='''
    ####  感谢您选择欧易生物！
    
    **关于网页报告使用有以下几点提示：**
    
    1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节，便于您快速查看相应内容。
    2. **报告正文**中的**图片**均可以点击后进行**放大查看**，且图片点击后可以包含更多的细节信息，比如左上角会显示具体的**图片数量**，右上角有相关的**图片工具**，**键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
    3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张，更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
    4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**，当表格多于20行的时候，表格会嵌入到网页之中，可以利用滚动栏进行调整设置，同时会在表格上方增加搜索设置，便于快速查询表格信息。
    5. 在报告中浏览，需要返回顶部的时候，可以使用网页**右下角的白色箭头标识**快速返回顶部。
    6. 本提示可以点击右上方的不再显示，则后续打开网页均不会显示本提示。
    ''')
    if library_type == "fresh":
        report.add_yaml_config(os.path.join(scriptdir, 'visium_report.yaml'))
    elif library_type == "ffpe":
        report.add_yaml_config(os.path.join(scriptdir, 'visium_report_FFPE.yaml'))
    elif library_type == "cytassist":
        report.add_yaml_config(os.path.join(scriptdir, 'visium_report_cytassist.yaml'))

    os.environ['OEBIO'] = os.path.join(scriptdir, "pic")

    ##
    reduct1 = config['params']['reduct1_method']
    reduct2 = config['params']['reduct2_method']
    ref_genome=config['database']['reference']

    if reduct2 == "tsne": REDUCT2 = "t-SNE"
    elif reduct2 == "umap": REDUCT2="UMAP" # "UMAP 主要基于流形理论和拓扑算法的理论，对高维数据进行降维，从而能够保留更多数据的全局结构，并且具有优越的运行性能。"
    if  reduct1 == "mnn":
        cluster_seurat_info = f"首先，我们采用 PCA 进行降维，并采用 {REDUCT2} 进行聚类展示。由于 {REDUCT2}  聚类可视化结果显示各个" \
                              f"样本切片中Spot点群样本间区分度较大，我们推测各个样本切片可能存在批次效应，因此我们使用改进版的MNN算法" \
                              f"fastMNN[^fastMNN] 算法进行批次效应剔除，并采用 {REDUCT2} 进行聚类展示，后续下游分析内容中聚类结果为去除批次效应后分析结果。"
    if  reduct1 == "pca,harmony":
        cluster_seurat_info = f"首先，我们采用 PCA 进行降维，并采用 {REDUCT2} 进行聚类展示。由于 {REDUCT2}  聚类可视化结果显示各个" \
                              f"样本切片中Spot点群样本间区分度较大，我们推测各个样本切片可能存在批次效应，因此我们改用 Harmony[^Harmony] " \
                              f"算法进行批次效应剔除，并采用 {REDUCT2} 进行聚类展示，后续下游分析内容中聚类结果为去除批次效应后分析结果。"
    if  reduct1 == "pca":
        cluster_seurat_info = f"常规基础分析报告中，我们默认采用 {reduct1.upper()} 进行降维，并采用 {REDUCT2} 进行聚类展示。"
    ################################################ 项目概况#################################################################
    summary = report.add_section('项目概况')
    report_video='https://www.bilibili.com/video/BV18u411T7UG/?spm_id_from=333.999.0.0'
    ###如果存在差异
    if (MT_flie == True and library_type == "fresh") or (library_type == "cytassist" and ref_genome == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A"):
            if sample_num >1:
                mitoinfo=f"，每个Spot中平均线粒体基因比例分布在 {min_mito}~{max_mito}"
            else:
                mitoinfo=f"，每个Spot中平均线粒体基因比例为 {min_mito}"
    else :
        mitoinfo=""
    if len(Diffexp_folder) == 0:
        if library_type == "fresh" or (library_type == "cytassist" and ref_genome == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A"):
            if sample_num > 1:
                summary_1 = summary.add_section('项目摘要_无差异',
                                            sample_num=sample_num,
                                            slices_num=slices_num,
                                            min_spot=min_spot,
                                            max_spot=max_spot,
                                            min_umi=min_umi,
                                            max_umi=max_umi,
                                            min_gene=min_gene,
                                            max_gene=max_gene,
                                            mitoinfo=mitoinfo,
                                            cluster_num=cluster_num)
            else:
                summary_1 = summary.add_section('项目摘要_无差异_样本数目为1',
                                            sample_num=sample_num,
                                            slices_num=slices_num,
                                            min_spot=min_spot,
                                            min_umi=min_umi,
                                            min_gene=min_gene,
                                            mitoinfo=mitoinfo,
                                            cluster_num=cluster_num)
        else:
            if sample_num > 1:
                summary_1 = summary.add_section('项目摘要_无差异',
                                            sample_num=sample_num,
                                            slices_num=slices_num,
                                            min_spot=min_spot,
                                            max_spot=max_spot,
                                            min_umi=min_umi,
                                            max_umi=max_umi,
                                            min_gene=min_gene,
                                            max_gene=max_gene,
                                            #min_mito='%s%%' % (min_mito),
                                            #max_mito='%s%%' % (max_mito),
                                            mitoinfo=mitoinfo,
                                            cluster_num=cluster_num)
            else:
                summary_1 = summary.add_section('项目摘要_无差异_样本数目为1',
                                            sample_num=sample_num,
                                            slices_num=slices_num,
                                            min_spot=min_spot,
                                            min_umi=min_umi,
                                            min_gene=min_gene,
                                            #min_mito='%s%%' % (min_mito),
                                            #max_mito='%s%%' % (max_mito),
                                            mitoinfo=mitoinfo,
                                            cluster_num=cluster_num)    
        summary_2 = summary.add_section('项目快捷链接_无差异',report_video=report_video) 
    else:
        diff = pd.read_csv('%s/diffexp_results_stat.xls' % (Diffexp_folder[0]), sep='\t')
        num_diff = len(diff)
        diff_num_tmp = [str(i) for i in diff.iloc[:, 4]]
        diff_num = ', '.join(diff_num_tmp)
        if library_type == "fresh" or (library_type == "cytassist" and ref_genome == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A"):
            summary_1 = summary.add_section('项目摘要_有差异',
                                            sample_num=sample_num,
                                            slices_num=slices_num,
                                            min_spot=min_spot,
                                            max_spot=max_spot,
                                            min_umi=min_umi,
                                            max_umi=max_umi,
                                            min_gene=min_gene,
                                            max_gene=max_gene,
                                            mitoinfo=mitoinfo,
                                            cluster_num=cluster_num,
                                            num_diff=num_diff,
                                            diff_num=diff_num)
        else:
            summary_1 = summary.add_section('项目摘要_有差异',
                                            sample_num=sample_num,
                                            slices_num=slices_num,
                                            min_spot=min_spot,
                                            max_spot=max_spot,
                                            min_umi=min_umi,
                                            max_umi=max_umi,
                                            min_gene=min_gene,
                                            max_gene=max_gene,
                                            #min_mito='%s%%' % (min_mito),
                                            #max_mito='%s%%' % (max_mito),
                                            mitoinfo=mitoinfo,
                                            cluster_num=cluster_num,
                                            num_diff=num_diff,
                                            diff_num=diff_num)
        summary_2 = summary.add_section('项目快捷链接_有差异', report_video=report_video)

    description=f'''[10x_visium_空间转录组报告树状图.pdf](10x_visium_空间转录组报告树状图.pdf)'''
    summary_3 = summary.add_section('项目目录结构',description=description)

    ################################################ 分析流程 #################################################################
    experi_module = report.add_section('技术简介')
    experi1 = experi_module.add_section('背景简介')
    experi2 = experi_module.add_section('技术原理')
    experi2 = experi_module.add_section('实验流程')
    workflow_module = report.add_section('生信分析流程')

    ################################################ 项目分析结果 #############################################################
    result_module = report.add_section('项目分析结果')

    ################################################ 基因定量质控 #############################################################
    samplefile = pd.read_csv('../../../config/samples.csv', index_col=0, sep=',')
    sample_name = samplefile.index[0]
    if  sample_num == 1 :
        cloupe_file_link  =  f"样本的cloupe文件：[1.SpaceRanger/{sample_name}/cloupe.cloupe](1.SpaceRanger/{sample_name}/cloupe.cloupe)"
    else:
        cloupe_file_link  =  f"aggr多样本整合分析后的cloupe文件：[1.SpaceRanger/aggr/cloupe.cloupe](1.SpaceRanger/aggr/cloupe.cloupe)"

    HE_file = list(glob(f"../../spaceranger/{sample_name}/outs/spatial/*{sample_name}*"))[0]
    SpaceRanger_folder = list(glob(r"[0-9]*.SpaceRanger"))[0]
    section1 = result_module.add_section('基因定量质控',HE_format = os.path.splitext(HE_file)[1], cloupe_file_link = cloupe_file_link)
    #section1_1 = section1.add_section('SpaceRanger分析')
    
    section2 = result_module.add_section('定量后质控')
    section2.add_table('*.Count_QC/statitics_for_QC.xls',
                         caption = '样本质控情况统计表',
                         headerdict={'sample':'样本名称',
                                    'mean_nFeature_Spatial_QC':'样本中所有Spots的UMI数平均值',
                                    'median_nFeature_Spatial_QC':'样本中所有Spots的UMI数中位数',
                                    'mean_nCount_Spatial_QC':'样本中所有Spots的基因数平均值',
                                    'median_nCount_Spatial_QC':'样本中所有Spots的UMI数中位数',
                                    'Total_Spots_QC':'样本中Spots数目'
                                    } )
    section1_2_2 = section2.add_section('表达基因比例')
    section1_2_3 = section2.add_section('UMI数目')
    if (library_type == "fresh" and MT_flie == True) or (library_type == "cytassist" and ref_genome == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A") :
        section1_2_1 = section2.add_section('线粒体基因比例')
    ################################################ 降维与聚类分析 ###########################################################
    Clustering_folder = list(glob(r"[0-9]*.Clustering"))[0]
    section3 = result_module.add_section('降维与聚类分析',cluster_seurat_info=cluster_seurat_info,Clustering_folder=Clustering_folder)
    section3_1 = section3.add_section(f'降维聚类结果')
    
    ################################################ 分组可视化
    vis_module = result_module.add_section(f'聚类结果可视化展示')
    vis_module_1 = vis_module.add_section(f'{REDUCT2}分组展示图')
    vis_module_2 = vis_module.add_section(f'{REDUCT2}分面展示图')
    vis_module_3 = vis_module.add_section(f'柱状统计图')


    ################################################ 空间特征基因鉴定#########################################################
    Marker_module = result_module.add_section('Marker基因鉴定')
    # Marker_folder = list(glob(r"[0-9]*.Marker"))[0]
    # link = list(glob(f'{Marker_folder}/all_markers_for_each_cluster_anno.xls'))[0]
    # Marker_module.add_comment('''每个Spot群中所有 Marker 基因结果表格：[{link}]({link})''',link=link)
    # link = list(glob(f'{Marker_folder}/top10_markers_for_each_cluster_anno.xls'))[0]
    # Marker_module.add_comment('''每个Spot群中 Top10 Marker 基因结果表格：[{link}]({link})''',link=link)
    Marker_module.add_table('*.Marker/top10_markers_for_each_cluster_anno.xls',
                            caption='每个 Spot 群 Top10 Marker 基因列表',
                            show_rows=10,
                            headerdict={'gene': '基因名',
                                        'p-value': 'p值',
                                        'avg_log2FC': '对数转化后的平均Foldchange值',
                                        'pct.1': '表达Marker基因的Spot在当前群中的占比',
                                        'pct.2': '表达Marker基因的Spot在其余群中的占比',
                                        'q-value': '校正后的p值',
                                        'cluster': 'Spot 群编号',
                                        'gene_diff': 'pct.1与pct.2的比值。Top10 Marker以gene_diff列为筛选依据。',
                                        'ensembl_id': 'ensembl 基因ID',
                                        'Dbxref_GeneID': ' NCBI基因ID索引号',
                                        'gene_type': '基因的类型',
                                        'gene_description': '基因描述',
                                        'TFs_Family': '转录因子家族',
                                        'GO_id': 'Gene Ontology登录号',
                                        'GO_term': 'GO 条目描述',
                                        'KEGG_id': 'KEGG 通路号',
                                        'KEGG_description': 'KEGG 通路描述'})

    section5_1 = Marker_module.add_section("Marker基因表达热图")
    section5_2 = Marker_module.add_section("Marker基因表达小提琴图")
    section5_3 = Marker_module.add_section("Marker基因空间可视化展示")
    ################################################ 空间特征基因鉴定#########################################################
    num=4
    if config['module']['celltype']:
        celltype_folder = list(glob(r"[0-9]*.Reference_CellType"))[0]
        num=num+1
        if config['params']['celltype']=='spotlight':
            section_add1 = result_module.add_section('细胞类型鉴定_spotlight', celltype_folder=celltype_folder)
        else:
            section_add1 = result_module.add_section('细胞类型鉴定_rctd', celltype_folder=celltype_folder)

        section_add1_1 = section_add1.add_section("Spot位点细胞类型组成饼图展示",celltype_folder=celltype_folder)
        section_add1_2 = section_add1.add_section("Spot位点细胞类型共定位情况",celltype_folder=celltype_folder)
        section_add1_3 = section_add1.add_section("Spot位点细胞类型组成热图展示",celltype_folder=celltype_folder)
        section_add1_4 = section_add1.add_section("Spot位点占比较高的细胞类型分布",celltype_folder=celltype_folder)
    else:
        section_add1 = result_module.add_section('细胞类型鉴定_no', tissue=header_info["物种信息"])

    ################################################## 基因差异表达分析 #########################################################

    if config['module']['diff_gene']:
        module_diff_gene_num=num+1
        if config['params']['qvalue']:
            difffiles = list(glob(f'*.Diffexp/*diff*qvalue*xls'))
            pvalue_or_padj = '矫正后p值，qvalue'
        else:
            difffiles = list(glob(f'*.Diffexp/*diff*pvalue*xls'))
            pvalue_or_padj = 'pvalue'
        diffexp_folder = list(glob(r"[0-9]*.Diffexp/"))[0]
        # diff_gene_GO3_plots = list(glob('%s/GO_enrichment/*/GO.top.*.png' % (enrichment_folder)))
        # diff_gene_GO5_plots = list(glob('%s/GO_enrichment/*/ALL_vs_DEG.GO.level2.stat.png' % (enrichment_folder[0])))
        # diff_gene_KEGG_top20_plots = list(glob('%s/KEGG_enrichment/*/KEGG.top.*.png' % (enrichment_folder[0])))
        # diff_gene_KEGG_level2_plots = list(
        #     glob('%s/KEGG_enrichment/*/ALL_vs_DEG.KEGG_Classification.png' % (enrichment_folder[0])))

        #diffexp_module = result_module.add_section('差异表达基因分析',diffexp_folder=diffexp_folder)
        diffexp_module= result_module.add_section('差异表达基因筛选',diffexp_folder=diffexp_folder,pvalue_or_padj=pvalue_or_padj)
        diffexp_module.add_table(f'*.Diffexp/diffexp_results_stat.xls',
                              caption='差异表达基因统计表',
                              headerdict={'Case': '实验组名',
                                          'Control': '对照组名',
                                          'Up_diff': '相比 Control，在Case中显著上调的差异基因数量',
                                          'Down_diff': '相比Control，在Case中显著下调的差异基因数量',
                                          'Total_diff(pvalue<0.05&FoldChange>1.5)': '显著性差异基因总数量'})
        diffexp_module.add_table('header_stat.txt', caption = '各分组差异基因数据统计结果各列说明')
        diffexp_module.add_comment("差异显著基因结果示例：")
        diffexp_module.add_table(difffiles[0],
                              caption='差异显著基因结果表格',
                              show_rows=10,
                              headerdict={'GeneID': '基因名',
                                          'p-value': 'p值',
                                          'pct.1': '表达该基因的Spot在比较组Spots群中的占比',
                                          'pct.2': '表达该基因的Spot在对照组Spots群中的占比',
                                          'q-value': '校正后的p值',
                                          'FoldChange': '差异倍数',
                                          'baseMean': '所有细胞中该基因的平均表达值',
                                          'log2FoldChange': 'log2转化后的差异倍数',
                                          'Regulation': '基因上下调描述',
                                          'ensembl_id': 'ensembl 基因ID',
                                          'Dbxref_GeneID': ' NCBI基因ID索引号',
                                          'gene_type': '基因的类型',
                                          'gene_description': '基因描述',
                                          'TFs_Family': '转录因子家族',
                                          'GO_id': 'Gene Ontology登录号',
                                          'GO_term': 'GO 条目描述',
                                          'KEGG_id': 'KEGG 通路号',
                                          'KEGG_description': 'KEGG 通路描述'})
        diffexp_module.add_table('header_diffexp.txt', caption = '差异显著基因结果各列说明')
        # section4_2 = section4.add_section('差异表达基因筛选统计结果')
        # section4_3 = section4.add_section('差异表达基因的MA图及火山图')
        enrichment_folder = list(glob(r"[0-9]*.Enrichment/"))[0]
        background_files_info=f'本次富集分析使用的背景文件见 [{enrichment_folder}/background_files]({enrichment_folder}/background_files)， 用此背景文件可以导入[云平台富集分析小工具](https://cloud.oebiotech.cn/task/detail/enrichment-oehw/)，进行灵活调整分析。'
        enrich_module = result_module.add_section('差异表达基因富集分析', background_files_info = background_files_info )
        diff_gene_GO = enrich_module.add_section('差异表达基因 GO 功能富集分析',enrichment_folder=enrichment_folder)
        link = f'{enrichment_folder}/GO_enrichment' 
        diff_gene_GO.add_comment('''差异表达基因 GO 功能富集分析结果：[{link}]({link})''',link=link)
        diff_gene_GO.add_plot('enrich.png',
			caption='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式',
			content='',
			description='其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。')
        #定义函数获取比较组名称
        def get_diffgroup_name(target,tag):
            diff_group = []
            for f in Path(target).glob(tag):
                t = str(f).split('/')[-1]
                t = t.split('-diff')[0]
                diff_group.append(t)
            return diff_group

        diff_group = get_diffgroup_name(Diffexp_folder[0], '*diff*pvalue*xls')
        diff_gene_GO.add_comment("GO 富集分析结果示例：")
        diff_gene_GO.add_table(f'{enrichment_folder}/GO_enrichment/{diff_group[0]}/enrichment-go-*-Total.xls',
			caption='富集分析结果表格',
			 show_rows=10)
			# headerdict={'id':'条目在Gene Ontology的登录号','term':'该条目的描述','category':'GO分类','ListHits':'该GO条目中差异基因数','ListTotal':'注释到GO的总差异基因数','PopHits':'注释到该条目中的所有基因数目','PopTotal':'注释到GO的总基因数','pval':'富集显著性p值','padj':'校正后的p值','Enrichment_score':'富集打分','Gene':'属于该条目的差异gene'})
        diff_gene_GO.add_table('header_enrich_go.txt', caption = 'GO 富集分析结果各列说明')
        diff_gene_GO3_plots = list(glob(f'{enrichment_folder}/GO_enrichment/*/GO.top.*.png'))
        if len(diff_gene_GO3_plots) >0:
            diff_gene_GO.add_comment('''GO 富集分析 Top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序的各 10 条）条形图展示如下：''')
            diff_gene_GO3_contents = ['']*len(diff_gene_GO3_plots)
            diff_gene_GO.add_plot(list(zip(diff_gene_GO3_plots, diff_gene_GO3_contents)),
				caption='GO 富集分析结果展示',
				description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1wS4y1z75L/)\n\n图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>Pvalue。')
            diff_gene_cohord_c3_plots = list(glob(f'{enrichment_folder}/GO_enrichment/*/GO.*chord*.png'))
        if len(diff_gene_cohord_c3_plots) >0:
            diff_gene_GO.add_comment('''富集分析和弦图，显示 p-value 最小的 10 个分类：''')
            diff_gene_cohord_c3_contents = ['']*len(diff_gene_cohord_c3_plots)
            diff_gene_GO.add_plot(list(zip(diff_gene_cohord_c3_plots,diff_gene_cohord_c3_contents)),
				caption='GO 富集 Top10 和弦图',
				description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1PF411j7Hq/)\n\n图片说明：左侧为每个分类中 |log2FC| 最大的 10 个基因，右侧反映分类组成情况，中间线条表示分类、基因对应关系。外侧热图表示对应基因的 log2FC 值。')
            diff_gene_GO_circle_plots = list(glob(f'{enrichment_folder}/GO_enrichment/*/enrichment-go-*-Total.circos.png'))
        if len(diff_gene_GO_circle_plots) >0:
            diff_gene_GO.add_comment('''GO 富集分析圈图，显示 p-value 最小的 20 个分类，从外到内共四圈：''')
            diff_gene_GO_circle_plots_contents = ['']*len(diff_gene_GO_circle_plots)
            diff_gene_GO.add_plot(list(zip(diff_gene_GO_circle_plots, diff_gene_GO_circle_plots_contents)),
				caption='GO 富集分析圈图',
				description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV15a411E7Le?spm_id_from=333.999.0.0)\n\n图片说明：第一圈：富集的分类，圈外为基因数目的坐标尺。不同的颜色代表不同的分类；第二圈：背景基因中该分类的数目以及 p-value。基因越多条形越长，值越小颜色越红，越大越蓝；第三圈：上下调基因比例条形图，浅红色代表上调基因比例，浅蓝色代表下调基因比例；下方显示具体的数值；第四圈：各分类的 RichFactor 值(该分类中前景基因数量除以背景基因数量)，背景辅助线每个小格表示 0.2。')
		
		## ========================= 2.4.9 差异基因KEGG富集分析 ========================= 
        link = f'{enrichment_folder}/KEGG_enrichment' 
        diff_gene_KEGG = enrich_module.add_section('差异表达基因 KEGG 功能富集分析')
        diff_gene_KEGG.add_comment('''差异表达基因 KEGG 功能富集分析结果：[{link}]({link})''',link=link)
		
        diff_gene_KEGG_top20_plots = list(glob(f'{enrichment_folder}/KEGG_enrichment/*/KEGG.top.*.png'))
        if len(diff_gene_KEGG_top20_plots)>0:
            diff_gene_KEGG.add_comment("KEGG 富集分析 Top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序）气泡图如下：")
            diff_gene_KEGG_top20_contents = ['']*len(diff_gene_KEGG_top20_plots)
            diff_gene_KEGG_top20_plots = list(zip(diff_gene_KEGG_top20_plots,diff_gene_KEGG_top20_contents))
            diff_gene_KEGG.add_plot(diff_gene_KEGG_top20_plots,
                caption='KEGG 富集 Top20 气泡图',
                description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV17a411E7nR?spm_id_from=333.999.0.0)\n\n图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由蓝-白-黄-红变化，其富集 p-value 值越小，显著程度越大。')
            diff_gene_cohord_c3_plots = list(glob(f'{enrichment_folder}/KEGG_enrichment/*/KEGG.*chord*.png'))
        if len(diff_gene_cohord_c3_plots)>0:
            diff_gene_KEGG.add_comment("富集分析和弦图，显示 p-value 最小的 10 个分类，图形分为左右两侧：左侧为每个分类中 |log2FC| 最大的 10 个基因，外侧热图表示对应基因的 log2FC 值：")
            diff_gene_cohord_c3_contents = ['']*len(diff_gene_cohord_c3_plots)
            diff_gene_cohord_c3_plots = list(zip(diff_gene_cohord_c3_plots, diff_gene_cohord_c3_contents))
            diff_gene_KEGG.add_plot(diff_gene_cohord_c3_plots,
                caption='KEGG 富集 Top10 和弦图',
                description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1jY4y167D7?spm_id_from=333.999.0.0)\n\n图片说明：右侧反映分类组成情况，中间线条表示分类、基因对应关系。')

            diff_gene_KEGG_circle_plots = list(glob(f'{enrichment_folder}/KEGG_enrichment/*/enrichment-kegg-*-Total.circos.png'))
        if len(diff_gene_KEGG_circle_plots)>0:
            diff_gene_KEGG.add_comment("KEGG 富集分析圈图，显示 p-value 最小的 20 个分类，从外到内共四圈：")
            diff_gene_KEGG_circle_contents = ['']*len(diff_gene_KEGG_circle_plots)
            diff_gene_KEGG_circle_plots = list(zip(diff_gene_KEGG_circle_plots,diff_gene_KEGG_circle_contents))
            diff_gene_KEGG.add_plot(diff_gene_KEGG_circle_plots,
            	caption='KEGG 富集分析圈图',
                description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1RB4y197nQ?spm_id_from=333.999.0.0)\n\n图片说明：第一圈：富集的分类，圈外为基因数目的坐标尺。不同的颜色代表不同的分类；第二圈：背景基因中该分类的数目以及 p-value。基因越多条形越长，值越小颜色越红，越大越蓝；第三圈：上下调基因比例条形图，浅红色代表上调基因比例，浅蓝色代表下调基因比例；下方显示具体的数值；第四圈：各分类的 RichFactor 值(该分类中前景基因数量除以背景基因数量)，背景辅助线每个小格表示 0.2。')

		
		## ========================= 2.4.10 差异基因PPI富集分析 =========================
        ppi_folder = list(glob(r"[0-9]*PPI"))
        if len(ppi_folder)>0: 
            ppi_module = result_module.add_section('差异表达基因蛋白网络互作分析')
            link = ppi_folder[0]
            ppi_module.add_comment('''结果目录：[{link}]({link})''',link=link)
            ppi_module.add_comment('''差异基因互作关系结果文件：''')
            ppi_table = list(glob(f'{ppi_folder[0]}/*_protein-protein-interaction.tsv' ))[0]
            ppi_module.add_table(ppi_table, 
                caption='差异基因互作关系表',show_rows=10)
            ppi_module.add_table('header_ppi.txt', caption = '差异基因互作关系结果各列说明')

			# ppi_module.add_table(ppi_table, 
				# caption='差异基因互作关系表',show_rows=10, 
				# headerdict={'stringId_A':'基因A的stringId','stringId_B':'基因B的stringId','preferredName_A':'A基因名','preferredName_B':'B基因名','ncbiTaxonId':'NCBI的物种分类标识符','score':'结合不同渠道证据的概率矫正后的综合得分','nscore':'邻近得分（从基因间核苷酸计数计算）','fscore':'融合得分（来自其他物种的融合蛋白质）','pscore':'系统谱系的共存得分（源自相似的基因缺失/存在模式）','ascore':'共表达得分（源自由DNA阵列和类似技术测量相似的mRNA表达模式）','escore':'实验得分（来自实验数据，如亲和色谱）','dscore':'数据库得分（派生自各种数据库的精选数据）','tscore':'文本挖掘得分（源自摘要中同时存在的基因/蛋白质名称）'})
			
            ppi_module.add_comment("上下调 Top25 差异表达基因互作网络图展示如下：")
            ppi_plots = list(glob(f'{ppi_folder[0]}/*string_protein-protein-interaction.new_colors.png'))
            ppi_plots = [ str(i) for i in ppi_plots ]
            ppi_contents = ['']*len(ppi_plots)
            ppi_plots = list(zip(ppi_plots,ppi_contents))
            ppi_module.add_plot(ppi_plots, caption = '差异基因互作网络图', 
                description = '图片说明：红色圆圈表示上调表达基因，绿色圆圈表示下调表达基因，节点之间的连线（或称为边）表示两蛋白之间具有相互作用，线的粗细表示相互作用关系的可靠性。')
            ppi_module.add_comment("上下调 Top25 差异表达基因互作圆形图展示如下：")
            ppi_plots = list(glob(f'{ppi_folder[0]}/*top_25_ppi_network.png'))
            ppi_plots = [ str(i) for i in ppi_plots ]
            ppi_contents = ['']*len(ppi_plots)
            ppi_plots = list(zip(ppi_plots,ppi_contents))
            ppi_module.add_plot(ppi_plots, caption = '差异基因互作圆形图', 
                description = '[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1X54y1Z7XW?spm_id_from=333.999.0.0)\n\n图片说明：红色表示上调差异表达基因，蓝色表示下调差异表达基因；关联的基因越多，基因点越大。')

        # if len(diff_gene_GO3_plots) > 0:
        #     section5_1 = section5.add_section('GO富集条形图')
        # if len(list(glob('%s/GO_enrichment/*/topGO_*.png' % (enrichment_folder[0])))) > 0:
        #     section5_2 = section5.add_section('topGO 有向无环图')
        # if len(diff_gene_GO5_plots) > 0:
        #     section5_3 = section5.add_section('差异基因和所有基因GO level2水平分布比较图')
        #     section5_4 = section5.add_section('上下调基因GO Level2水平分布比较图')

        # section6 = result_module.add_section('差异表达基因KEGG功能富集分析',module_diff_gene_num=module_diff_gene_num)
        # if len(diff_gene_KEGG_top20_plots) > 0:
        #     section6_1 = section6.add_section('KEGG富集气泡图')
        # if len(diff_gene_KEGG_level2_plots) > 0:
        #     section6_2 = section6.add_section('差异基因和所有基因在KEGG Level2水平分布比较图')
        #     section6_3 = section6.add_section('上下调基因在KEGG Level2水平分布比较图')

        # ######################################## 差异基因蛋白网络互作分析 ########################################################################
        # if config['module']['diff_gene_ppi']:
        #     section7 = result_module.add_section('样本间差异表达基因蛋白网络互作分析',module_diff_gene_num=module_diff_gene_num)
        #     section7.add_comment('''差异基因互作关系结果文件：''')
        #     ppi_table = list(glob(f'*Diffexp/PPI/*_protein-protein-interaction.tsv'))[0]
        #     section7.add_table(ppi_table,
        #                          caption='差异基因互作关系表', show_rows=10,
        #                          headerdict={'stringId_A': '基因A的stringId',
        #                                      'stringId_B': '基因B的stringId',
        #                                      'preferredName_A': 'A基因名',
        #                                      'preferredName_B': 'B基因名',
        #                                      'ncbiTaxonId': 'NCBI的物种分类标识符',
        #                                      'score': '结合不同渠道证据的概率矫正后的综合得分',
        #                                      'nscore': '邻近得分（从基因间核苷酸计数计算）',
        #                                      'fscore': '融合得分（来自其他物种的融合蛋白质）',
        #                                      'pscore': '系统谱系的共存得分（源自相似的基因缺失/存在模式）',
        #                                      'ascore': '共表达得分（源自由DNA阵列和类似技术测量相似的mRNA表达模式）',
        #                                      'escore': '实验得分（来自实验数据，如亲和色谱）',
        #                                      'dscore': '数据库得分（派生自各种数据库的精选数据）',
        #                                      'tscore': '文本挖掘得分（源自摘要中同时存在的基因/蛋白质名称）'})

        #     section7_1 = section7.add_section('差异基因互作网络图')
        #     section7_2 = section7.add_section('上下调差异基因互作网络图')

    ######################################### 附录 #############################################################################
    section8 = report.add_section('附录')
    folder = list(glob(r"[0-9]*Supplemental_material"))[0]
    section8.add_comment('''附录目录：[{link}]({link})''',link=folder)
    section8_0 = section8.add_section("报告解读视频", report_video=report_video)
    section8_0 = section8.add_section("下一步分析建议", folder=folder)
    section8_0 = section8.add_section("空间转录组常见问题（FAQ）", folder=folder)
    section8_1 = section8.add_section("Loupe Browser使用教程", folder=folder)
    section8_2 = section8.add_section("实验技术方法说明", folder=folder)
    section8_3 = section8.add_section("生信分析方法说明", folder=folder)
    section8_4 = section8.add_section('AI修图使用说明', folder=folder)
    section8_5 = section8.add_section("参考数据库链接", folder=folder)
    section8_6 = section8.add_section("分析软件列表", folder=folder)

    section9 = report.add_section('申明')
    report.write_to('分析说明.html', zip_report_name=f'{outdir}.zip')

if __name__ == "__main__":
    convert_html()
