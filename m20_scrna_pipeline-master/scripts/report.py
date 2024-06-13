#!/usr/bin/env python37
# encoding:utf-8

# Report
# 1.STARsolo
# 2.Count_QC
# 3.Clustering
# 4.Annotation
# 5.Marker
# 6.MarkerEnrich
# 7.Diffgene
# 8.DiffgeneEnrich


import click
import os
import pandas as pd
import re
import shutil
import subprocess
import sys
import time
import yaml
from glob import glob
# from oebio.app import *
from oebio.report import Report, oeweb_register


@click.command()
@click.option('-i', '--inputdir', prompt='the program path',
              help='the directory of your program. ')
@click.option('-c', '--configfile', prompt='config file', default='config/config.yaml',
              help='the config file which contain your program information.')
@click.option('-d', '--dirnum', default=0)


def multimodal_report(inputdir, configfile, dirnum):
    ##=============================================== 1. 读取项目配置文件信息 ===========================================
    """ python $0 -i result -c config/config.yaml"""
    cfg = open(configfile, 'r', encoding='utf-8').read()
    config = yaml.full_load(cfg)
    program_path = os.path.abspath(inputdir)
    # paragram from config
    program_num = dict(config['report'])['Project_Num']
    outdir = f"report/{program_num}_Report"
    log_file = open("%s/report_log.txt" % (program_path), "w")
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    # Report
    os.chdir(program_path)
    if os.path.exists(outdir):
        shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    if not os.path.exists("%s/src" % (outdir)):
        os.makedirs("%s/src" % (outdir))
        subprocess.call(f'cp {scriptdir}/src/* {outdir}/src/', shell=True)

    ##================================================2. 生成结果文件 =================================================

    # 1.STARsolo
    if os.path.exists("%s/1.STARsolo" % (program_path)):
        os.makedirs("%s/1.STARsolo" % (outdir))
        subprocess.call('cp -r %s/1.STARsolo/*/*_Solo.out %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/BarcodeRank_plot %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/Read_Distribution %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/Genetype_plot %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/STARsolo.statistic.tsv %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        os.makedirs("%s/1.STARsolo/Qualimap_result" % (outdir)),
        subprocess.call('cp -r %s/1.STARsolo/Qualimap_result/qualimap_plot.* %s/1.STARsolo/Qualimap_result/' % (program_path, outdir),
                        shell=True)
    else:
        print("Can not find 1.STARsolo results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 1.STARsolo results!!!!!!!!!!!!!!!!" + "\n")

    # 2.Count_QC
    if os.path.exists("%s/2.Count_QC" % (program_path)):
        os.makedirs("%s/2.Count_QC" % (outdir))
        subprocess.call('cp %s/2.Count_QC/*.png %s/2.Count_QC' % (program_path, outdir), shell=True)
        subprocess.call('cp %s/2.Count_QC/*.pdf %s/2.Count_QC' % (program_path, outdir), shell=True)
        subprocess.call('cp %s/2.Count_QC/*.xls %s/2.Count_QC' % (program_path, outdir), shell=True)
    else:
        print("Can not find 2.Count_QC results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 2.Count_QC results!!!!!!!!!!!!!!!!!" + "\n")

    # 3.Clustering
    if os.path.exists("%s/3.Clustering" % (program_path)):
        os.makedirs("%s/3.Clustering" % (outdir))
        subprocess.call('cp -r %s/3.Clustering/* %s/3.Clustering' % (program_path, outdir), shell=True)
    else:
        print("Can not find 3.Clustering results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 3.Clustering results!!!!!!!!!!!!!!!!!" + "\n")

    # 4.All_Marker
    if os.path.exists("%s/4.mRNA_marker" % (program_path)):
        os.makedirs("%s/4.mRNA_marker" % (outdir))
        subprocess.call('cp -r %s/4.mRNA_marker/* %s/4.mRNA_marker' % (program_path, outdir), shell=True)
    else:
        print("Can not find 4.mRNA_marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 4.mRNA_marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 4.LncRNA_Marker
    if os.path.exists("%s/5.LncRNA_marker" % (program_path)):
        os.makedirs("%s/5.LncRNA_marker" % (outdir))
        subprocess.call('cp -r %s/5.LncRNA_marker/* %s/5.LncRNA_marker' % (program_path, outdir), shell=True)
    else:
        print("Can not find 5.LncRNA_marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 5.LncRNA_marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 5.Celltyping
    if os.path.exists("%s/6.Reference_celltype" % (program_path)):
        os.makedirs("%s/6.Reference_celltype" % (outdir))
        subprocess.call('cp -r %s/6.Reference_celltype/* %s/6.Reference_celltype' % (program_path, outdir), shell=True)
    else:
        print("Can not find 6.Reference_celltype results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 6.Reference_celltype results!!!!!!!!!!!!!!!!!" + "\n")

    # 6.Diffgene
    if os.path.exists("%s/7.Diffexp" % (program_path)):
        os.makedirs("%s/7.Diffexp" % (outdir))
        subprocess.call('cp -r %s/7.Diffexp/* %s/7.Diffexp' % (program_path, outdir), shell=True)
    else:
        print("Can not find 7.Diffexp results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 7.Diffexp results!!!!!!!!!!!!!!!!!" + "\n")


    ##===================================================== 3. 生成报告================================================
    header_info = {key: config["report"][key] for key in config["report"] if
                   key in ["Project_Num", "Customer", "Species","Sample", "Executor", "Sample_Num", "Task_Num"]}
    header_info["项目编号"] = header_info.pop("Project_Num")
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["实验物种"] = header_info.pop("Species").capitalize()
    header_info["实验样本"] = header_info.pop("Sample")
    header_info["执行编号"] = header_info.pop("Executor")
    header_info["样本数目"] = header_info.pop("Sample_Num")
    header_info["任务单号"] = header_info.pop("Task_Num")

    report = Report('M20 单细胞转录组项目结题报告', title='M20 单细胞转录组项目结题报告',
                    header_info=dict(header_info))
    report.add_yaml_config(os.path.abspath(f'{scriptdir}/report.yaml'))
    os.environ['oeweb_register_token'] = 'oearray'
    Project_module = report.add_module('项目概况')
    register_info = oeweb_register(project_id=program_num, target_url='https://cloud.oebiotech.cn/task/category/scrna',
                                   note='本报告包含项目基本分析内容，如需快速做细胞表达可视化等分析(数据默认保留半年)')
    if register_info:
        Project_module.add_comment(register_info)
    report_summary_list = []
    report_hyperlink_list = []
    os.chdir(outdir)
    ## ======================================================= 实验及分析概述============================================
    Background_module = report.add_module('技术简介')
    Background_section1 = Background_module.add_section('实验流程')  ## 实验流程
    workflow_module = Background_module.add_section('生信分析流程')  ## 分析流程
    # ====================================================== 4. 分析结果 ================================================
    result_module = report.add_module('项目分析结果')
    ##====================================================== 4.1 质控分析 =================================================
    if os.path.exists("1.STARsolo"):
        section1 = result_module.add_section('基因定量质控')
        report_hyperlink_list.append('[STARsolo分析](#STAR)\n\n')
        section1_1 = section1.add_section('STARsolo分析')
        section1_1.add_table(f'1.STARsolo/STARsolo.statistic.tsv',
                             caption='样本STARsolo结果汇总统计表',
                             headerdict={'Sampleid': '样本名称',
                                         'Number of Reads': '测序Reads数目',
                                         'Reads With Valid Barcodes': '包含有效barcode的Reads占比',
                                         'Sequencing Saturation': '测序饱和度',
                                         'Q30 Bases in CB+UMI': 'Barcode+UMI序列的Q30质量',
                                         'Q30 Bases in RNA read': 'Read2的Q30质量',
                                         'Reads Mapped to Genome: Unique+Multiple': '比对到参考基因组的Reads比例',
                                         'Reads Mapped to Genome: Unique': '唯一比对到参考基因组的Reads比例',
                                         'Reads Mapped to GeneFull: Unique+Multiple GeneFull': '未采用该模式，基因比对统计基于唯一比对Reads',
                                         'Reads Mapped to GeneFull: Unique GeneFull': '唯一比对到基因区域的Reads比例',
                                         'Estimated Number of Cells': '有效细胞数目',
                                         'Unique Reads in Cells Mapped to GeneFull': '唯一比对到基因的Reads数目',
                                         'Fraction of Unique Reads in Cells': '唯一比对Reads在细胞中的占比',
                                         'Mean Reads per Cell': '细胞中Reads均值',
                                         'Median Reads per Cell': '细胞中Reads中位数',
                                         'UMIs in Cells': '细胞中UMIs数目',
                                         'Mean UMI per Cell': '细胞中UMI均值',
                                         'Median UMI per Cell': '细胞中UMI中位值',
                                         'Mean GeneFull per Cell': '细胞中基因均值',
                                         'Median GeneFull per Cell': '细胞中基因中位值',
                                         'Total GeneFull Detected': '总的基因检出数目'
                                         })
        section1_2 = section1.add_section('有效液滴指标图')
        section1_3 = section1.add_section('基因元件组成图')
        section1_4 = section1.add_section('基因类型分布图')
        section1_5 = section1.add_section('基因Reads覆盖分布图')
    ##====================================================== 4.2 基因及样本过滤 ==============================================
    if os.path.exists("2.Count_QC"):
        report_hyperlink_list.append('[基因及细胞过滤结果](#filter)\n\n')
        filter_section = result_module.add_section('基因及细胞过滤')
        QC_folder = "2.Count_QC"
        summ = pd.read_csv('%s/statitics_before_after_QC.xls' % (QC_folder), sep='\t')
        sample_num = len(summ)
        min_beforeQC_cell = min(summ.Total_cells_beforeQC)
        max_beforeQC_cell = max(summ.Total_cells_beforeQC)
        min_afterQC_cell = min(summ.Total_cells_afterQC)
        max_afterQC_cell = max(summ.Total_cells_afterQC)
        min_nFeature_afterQC = round(min(summ.mean_nFeature_RNA_afterQC.astype(float)))
        max_nFeature_afterQC = round(max(summ.mean_nFeature_RNA_afterQC.astype(float)))
        min_nCount_RNA_afterQC = round(min(summ.mean_nCount_RNA_afterQC.astype(float)))
        max_nCount_RNA_afterQC = round(max(summ.mean_nCount_RNA_afterQC.astype(float)))
        if len(summ) == 1:
            report_summary_list.append(
                f"本次分析共完成 {sample_num} 个样本的M20 单细胞转录组测序，数据经STAR[^STAR]质控过滤后的高质量细胞数为 {min_beforeQC_cell}"
                f" 个，经剔除双细胞、多细胞、异常基因表达以及低基因表达的细胞后，最终获得的细胞数目为 {min_afterQC_cell} 个，每个细胞中的平均 UMI 数为"
                f" {min_nCount_RNA_afterQC}，每个细胞中的平均基因数为 {min_nFeature_afterQC}")
        else:
            report_summary_list.append(
                f"本次分析共完成 {sample_num} 个样本的M20 单细胞转录组测序，各样本经STAR质控过滤后的高质量细胞数分布在"
                f" {min_beforeQC_cell}~{max_beforeQC_cell} 个，经剔除双细胞、多细胞、异常基因表达以及低基因表达的细胞后，最终获得的细胞数目分布在"
                f" {min_afterQC_cell}~{max_afterQC_cell} 个，每个细胞中的平均 UMI 数分布在 {min_nCount_RNA_afterQC}~{max_nCount_RNA_afterQC}，"
                f"每个细胞中的平均基因数分布在 {min_nFeature_afterQC}~{max_nFeature_afterQC}")

        report_summary_list.append(f"。")

        ##
        QC_description = f"在 STAR 初步质控的基础上进一步利用Seurat[^seurat]包对实验数据进行质控，剔除多细胞、双细胞或者未结合上细胞的数据，" \
                         f"进行下游分析。理论上大部分细胞中表达的基因数量、" \
                         f"UMI 数量、lncRNA基因以及rRNA基因数量等会集中分布在某一区域内，根据它们的分布特征可以拟合分布模型，使用该模型找" \
                         f"到其中的离域值，剔除异常数据。本项目中的质控标准为：保留基因数 、UMI数在上下两倍标准差范围内，log10GenesPerUMI大于0.7、" \
                         f"线粒体占比低于 30%、血红细胞基因占比低于5%的细胞作为高质量细胞，进行下游分析，并对lncRNA基因以及rRNA基因的占比进行统计。" \
                         f"然后使用DoubletFinder[^doubletfinder]软件进行双细胞去除，进行下游分析。\n [^doubletfinder]: Mcginnis C S , " \
                         f"Murrow L M , Gartner Z J . DoubletFinder: Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors. 2018.\n" \
                         f" [^seurat]: Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data[J]. Cell, 2021, 184(13):3573-3587.\n" \
                         f"\n 质控前后每个细胞的基因数量（nFeature）、UMI数量（nCount）、单位UMI内基因数目的比例（log10GenesPerUMI）、线粒体UMI所占比例（percent.mito）、" \
                         f"血红细胞基因所占比例（percent.HB）、lncRNA基因所占比例（percent.lncRNA）、rRNA基因所占比例（percent.rRNA）的小提琴图展示如下："
        ##
        filter_section2 = filter_section.add_section(name='基因及样本过滤结果', description=QC_description)
        filter_section2.add_comment("详细结果见目录：[2.Count_QC](./2.Count_QC)")
        filter_section2.add_plot('2.Count_QC/QC_metrics_beforeQC.png', caption='过滤低质量细胞前各质控指标的小提琴图')
        filter_section2.add_plot('2.Count_QC/QC_metrics_afterQC.png', caption='过滤低质量细胞后各质控指标的小提琴图')
        filter_section2.add_comment("定量质控前后的细胞数统计情况及过滤标准如下表所示：")
        filter_section2.add_table(summ.iloc[:, :31],
                                  caption='质控前后细胞数目统计表',
                                  headerdict={'sample': '样本名',
                                              'mean_nCount_RNA_beforeQC': '质控前RNA平均Counts数',
                                              'median_nCount_RNA_beforeQC': '质控前RNA Counts的中位数',
                                              'mean_nFeature_RNA_beforeQC': '质控前RNA平均Feature数',
                                              'median_nFeature_RNA_beforeQC': '质控前RNA feature的中位数',
                                              'mean_log10GenesPerUMI_beforeQC': '质控前单位UMI内基因数目的平均比例',
                                              'median_log10GenesPerUMI_beforeQC':'质控前单位UMI内基因数目的比例中位数',
                                              'mean_percent.mito_beforeQC':'质控前线粒体基因平均比例',
                                              'median_percent.mito_beforeQC':'质控前线粒体基因比例中位数',
                                              'mean_percent.HB_beforeQC': '质控前红细胞基因平均比例',
                                              'median_percent.HB_beforeQC': '质控前红细胞基因比例中位数',
                                              'mean_percent.lncRNA_beforeQC': '质控前lncRNA基因比例平均数',
                                              'median_percent.lncRNA_beforeQC': '质控前lncRNA基因比例中位数',
                                              'mean_percent.rRNA_beforeQC': '质控前rRNA基因比例平均数',
                                              'median_percent.rRNA_beforeQC': '质控前rRNA基因比例中位数',
                                              'Total_cells_beforeQC': '质控前细胞总数',
                                              'mean_nCount_RNA_afterQC': '质控后RNA平均counts数',
                                              'median_nCount_RNA_afterQC': '质控后RNA Counts的中位数',
                                              'mean_nFeature_RNA_afterQC': '质控后RNA平均Feature数',
                                              'median_nFeature_RNA_afterQC': '质控后RNA feature的中位数',
                                              'mean_log10GenesPerUMI_afterQC':'质控后单位UMI内基因数目的平均比例',
                                              'median_log10GenesPerUMI_afterQC':'质控后单位UMI内基因数目的比例中位数',
                                              'mean_percent.mito_afterQC':'质控后线粒体基因平均比例',
                                              'median_percent.mito_afterQC':'质控后线粒体基因比例中位数',
                                              'mean_percent.HB_afterQC': '质控后红细胞基因平均比例',
                                              'median_percent.HB_afterQC': '质控后红细胞基因比例中位数',
                                              'mean_percent.lncRNA_afterQC': '质控后lncRNA基因比例平均数',
                                              'median_percent.lncRNA_afterQC': '质控后lncRNA基因比例中位数',
                                              'mean_percent.rRNA_afterQC': '质控后rRNA基因比例平均数',
                                              'median_percent.rRNA_afterQC': '质控后rRNA基因比例中位数',
                                              'Total_cells_afterQC': '质控后细胞总数'
                                              })

    ##=============================================== 4.3 降维聚类 ==============================================
    ## 降维与聚类
    if os.path.exists("3.Clustering"):
        cluster_section = result_module.add_section('数据降维与聚类分析')
        report_hyperlink_list.append('[降维与聚类分析结果](#rna_cluster)\n\n')
        if len(summ) == 1:
            sample_des = f"仅有1个样本，并不存在批次效应"
        else:
            sample_des = f"有不同样本，但并不存在批次效应"
        if config['params']['Clustering']["RNA_reduct1"] == "pca":
            reduct1_des = f"PCA (principal component analysis, 主成分分析)[^pca] 的降维算法进行降维 "
            reduct1_name = "PCA"
        else:
            sample_des = f"有不同样本，且存在批次效应"
            reduct1_des = f"MNN (mutual nearest neighbors, 互享最近邻)[^mnn] 的降维算法剔除批次效应"
            reduct1_name = "MNN"
        if config['params']['Clustering']["RNA_reduct2"] == "tsne":
            reduct2_des = f"t-SNE (t分布式随机邻近嵌入)[^tsne] 算法进行二级降维。t-SNE创建了一个缩小的特征空间，相似的样本由附近的点建模，不相似的样本由高概率的远点建模，更好的展示了样本之间的局部关系"
            reduct2_name = "t-SNE"
        else:
            reduct2_des = f"UMAP[^umap] (统一流形逼近与投影) 算法进行二级降维。UMAP 主要基于流形理论和拓扑算法的理论，对高维数据进行降维，从而能够保留更多数据的全局结构，并且具有优越的运行性能"
            reduct2_name = "UMAP"

        cluster_section2 = cluster_section.add_section("降维与聚类分析结果")
        cluster_section2.add_comment(
            f"由于本项目{sample_des}，因此我们首先采用了{reduct1_des}，再通过 {reduct2_des}。基于{reduct1_name}降维结果通过{reduct2_name}对单细胞群聚类进行可视化，聚类算法采用SNN，最终获得最优细胞分群。\n [^pca]: Ian T. Jolliffe. PRINCIPAL COMPONENT ANALYSIS: A BEGINNER'S GUIDE — I. Introduction and application[J]. Weather, 1990, 45.\n [^tsne]: van der Maaten, Laurens. Accelerating t-SNE using Tree-Based Algorithms[J]. Journal of Machine Learning Research, 2014, 15(1):3221-3245. \n [^umap]: Mcinnes L,Healy J,Melville J.UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction[J].The Journal of Open Source Software,2018. \n [^mnn]: Haghverdi, L., Lun, A., Morgan, M. et al. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat Biotechnol 36, 421–427 (2018). https://doi.org/10.1038/nbt.4091. \n")
        cluster_section2.add_comment("详细结果见目录：[3.Clustering](./3.Clustering)")

        cluster_section2.add_plot('3.Clustering/*/*_groupby_cluster_resolution*_plot.png',
                                  caption='基因表达水平降维与聚类图',
                                  description='图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')

        ## 样本分组展示
        if (sample_num > 1):
            cluster_section3 = cluster_section.add_section('样本间降维聚类分组展示')
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/groupby-sampleid_contrast_plot.png',
                                      caption='多样本降维聚类分组展示图',
                                      description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。')
            ################
            cluster_section3.add_comment(" ")
            cluster_section3.add_comment("每个细胞群中样本占比的柱状统计图如下：")
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/groupby-clusters_summary_plot.png',
                                      caption='每个细胞群中样本占比的柱状统计图',
                                      description='图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。')

            ################
            cluster_section3.add_comment(" ")
            cluster_section3.add_comment("多个样本的降维聚类分面展示图如下：")
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.png',
                                      caption='多样本降维聚类分面展示图',
                                      description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同群的细胞以不同颜色区分。')

            ################
            cluster_section3.add_comment(" ")
            cluster_section3.add_comment("每个样本中细胞群占比的柱状统计图如下：")
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/groupby-sampleid_summary_plot.png',
                                      caption='样本间细胞群占比柱状图',
                                      description='图片说明：横坐标表示不同样本，纵坐标表示不同组别中细胞数目所占的百分比。')
            ###############
            cluster_section3.add_comment("各样本在不同细胞群中的细胞数目统计表如下：")
            cluster_section3.add_table("3.Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls",
                                       caption='样本间细胞群数目统计表')

        ## 细胞群间相似性分析
        cluster_section4 = cluster_section.add_section("细胞群间相似性分析")
        cluster_section4.add_plot('3.Clustering/clusters_correlation/coefficient_heatmap.png',
                                  caption='细胞群间相关性热图',
                                  description='图片说明：横纵坐标为不同细胞群，图中的数字为皮尔森相关性系数。该值越大，热图颜色越红，表示细胞群之间的相关性程度越高。')

    ##=============================================== 4.4 mRNA_Marker基因鉴定 ================================================
    if os.path.exists("4.mRNA_marker"):
        test_method = config["params"]['Marker']["RNA_test_method"]
        Marker_description = f"Marker 基因的定义是该基因在指定细胞群的绝大多数细胞中有较高的表达，而在其余细胞类群中只有少部分表达，且该基因在此细胞群相对于其他细胞群中是显著上调表达。\n 本研究采用  {test_method}  检验方法对指定细胞群与其余所有细胞群进行差异检验，先设置log2FC > 0 ,min.pct > 0.05导出所有可能的 marker ,然后按照P值，FoldChange 值以及 gene_diff 值筛选 top10 marker 作为每个细胞群的特异性 Marker 基因。"
        marker_section = result_module.add_section(name='Marker基因鉴定', description=Marker_description)
        marker_section1 = marker_section.add_section(name="mRNA编码基因鉴定")
        marker_section1.add_comment("mRNA编码基因在基因表达和蛋白质翻译过程中起着至关重要的作用，筛选每个细胞群的特异性mRNA编码基因是十分必要的。本研究挑选出所有mRNA编码基因，根据P值，FoldChange 值以及 gene_diff 值筛选出 top10 mRNA编码基因作为每个细胞群的特异性 Marker 基因。")
        report_hyperlink_list.append('[mRNA Marker基因鉴定结果](#mrna_marker_gene)\n\n')
        marker_section1.add_comment("详细结果见目录：[4.mRNA_marker](./4.mRNA_marker)")
        #####
        marker_result = pd.read_csv('4.mRNA_marker/top10_markers_for_each_cluster.xls', sep='\t')
        cluster_num = max(marker_result['cluster'])
        report_summary_list.append(f"通过基因表达数据进行降维聚类后共分为 {cluster_num} 群细胞，")
        ####
        marker_section1.add_table("4.mRNA_marker/top10_mRNA_genes.xls",
                                  caption='每个细胞群 Top10 mRNA Marker 基因列表',
                                  show_rows=5)
        ####
        marker_section1.add_plot("4.mRNA_marker/topmarker_gene_heatmap.png",
                                 caption='Top10 Marker 基因表达热图',
                                 description='图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')
        ####
        marker_section1.add_plot("4.mRNA_marker/*/marker_gene_featureplot.png",
                                 caption='Top10 Marker 基因在降维聚类结果中的可视化图',
                                 description='图片说明：红色越深表示该细胞中对应基因的表达量越高。')
        ####
        marker_section1.add_plot("4.mRNA_marker/*/marker_gene_violin_plot.png",
                                 caption='Top10 Marker 基因表达量小提琴图',
                                 description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')

    ##=============================================== 4.5 LncRNA_Marker基因鉴定 ================================================
    if os.path.exists("5.LncRNA_marker"):
        marker_section2 = marker_section.add_section(name='LncRNA编码基因鉴定')
        report_hyperlink_list.append('[LncRNA Marker基因鉴定结果](#lncrna_marker_gene)\n\n')
        marker_section2.add_comment("LncRNA作为蛋白编码基因和非编码基因的重要调节因子在生理和病理细胞过程中发挥着不同的作用，本研究选取所有lncRNA编码基因，同样根据P值，FoldChange 值以及 gene_diff 值筛选出 top10 lncRNA编码基因作为每个细胞群的特异性 Marker 基因。")
        marker_section2.add_comment("详细结果见目录：[5.LncRNA_marker](./5.LncRNA_marker)")
        #####
        marker_section2.add_table("5.LncRNA_marker/top10_lncRNA_genes.xls",
                                  caption='每个细胞群 Top10 lncRNA Marker 基因列表',
                                  show_rows=5)
        ####
        marker_section2.add_plot("5.LncRNA_marker/topmarker_gene_heatmap.png",
                                 caption='Top10 Marker 基因表达热图',
                                 description='图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')
        ####
        marker_section2.add_plot("5.LncRNA_marker/*/marker_gene_featureplot.png",
                                 caption='Top10 Marker 基因在降维聚类结果中的可视化图',
                                 description='图片说明：红色越深表示该细胞中对应基因的表达量越高。')
        ####
        marker_section2.add_plot("5.LncRNA_marker/*/marker_gene_violin_plot.png",
                                 caption='Top10 Marker 基因表达量小提琴图',
                                 description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')

    ##=============================================== 4.6 细胞类型鉴定 ===================================================
    if os.path.exists("6.Reference_celltype"):
        report_hyperlink_list.append('[细胞类型鉴定](#celltype)\n\n')
        celltype_file = list(glob('6.Reference_celltype/*simplified_celltyping_results.csv'))[0]
        cell_results = pd.read_csv(celltype_file, sep=',')
        celltype = ', '.join(set(cell_results.celltype))
        report_summary_list.append(f"通过对细胞类型进行注释，得到供参考的细胞类型有{celltype}。")
        rna_section5 = result_module.add_section('细胞类型鉴定')
        link = list(glob(f'6.Reference_celltype/*top.*celltyping_plot.pdf'))[0]
        refcelltype = str(link).split('/')[-1].split("_")[1]
        rna_section5.add_comment(name=refcelltype.lower())
        rna_section5.add_comment("详细结果见目录：[6.Reference_celltype](./6.Reference_celltype)")
        rna_section5.add_plot("6.Reference_celltype/*ref_*_top.main_celltyping_plot.png",
                              caption='细胞类型参考结果图',
                              description='图片说明：细胞类型结果在UMAP图上的展示，每种细胞类型以不同颜色进行区分。')
        rna_section5.add_plot("6.Reference_celltype/*ref_*_main_celltyping_heatmap.png",
                              caption='细胞类型鉴定相关性热图',
                              description='图片说明：每一行表示参考数据集中的细胞类型注释名称，每一列表示待鉴定的细胞。黄色越深表示相关性值越大，表明待鉴定的细胞类型与参考数据集中的该种细胞类型最为相似。')
        rna_section5.add_table("6.Reference_celltype/*simplified_celltyping_results.csv",
                               caption='细胞类型注释结果表格',
                               show_rows=5)

    ##=============================================== 4.7 差异基因筛选及功能富集 ============================================
    if os.path.exists("7.Diffexp"):
        diff_section = result_module.add_section('差异基因筛选及富集')
        report_hyperlink_list.append('[差异表达基因筛选及富集](#diffexp)\n\n')
        diffexp_folder = "7.Diffexp"
        diffexp_file = list(glob('%s/diffexp_results_stat.xls' % (diffexp_folder)))[0]
        diff = pd.read_csv(diffexp_file, sep='\t')
        num_diff = len(diff)
        diff_num = ', '.join([str(i) for i in diff.iloc[:, 4]])
        report_summary_list.append(f"共设有{num_diff}个差异基因分组，其检测到的差异基因数量为：{diff_num}。")
        test_method = config['params']['Diffexp']['test']
        if test_method == 'presto':
            diff_descript = f'差异表达基因使用presto[^presto]包的WilcoxAUC函数进行筛选，默认使用 presto 差异检验方法。\n [^presto]: Korsunsky I, Nathan A, Millard N, et al. Presto scales Wilcoxon and auROC analyses to millions of observations[J]. BioRxiv, 2019: 653253. \n'
        else:
            diff_descript = f'差异表达基因使用Seurat包的FindMarkers函数进行筛选，默认使用 {test_method} 差异检验方法。'
        ##差异基因
        diff_section1 = diff_section.add_section(name='差异表达基因筛选', description=diff_descript)
        diff_section1.add_comment('差异表达基因具体结果见: [7.Diffexp](./7.Diffexp)')
        diff_section1.add_comment("各分组差异基因数目统计表如下：")
        diff_section1.add_table("7.Diffexp/diffexp_results_stat.xls", caption='差异表达基因统计表')
        diff_section1.add_comment("差异显著基因结果示例：")
        difffiles = list(glob(f'{diffexp_folder}/*diff*xls'))
        diff_section1.add_table(difffiles[0],
                                caption='差异显著基因结果表格',
                                show_rows=5)
        diff_section1.add_comment("将差异基因的差异倍数（FoldChange）从大到小排列，上下调各选取 25 个基因绘制热图：")
        diff_section1.add_plot(f'{diffexp_folder}/Heatmap/*_heatmap.png',
                               caption='上下调Top 25差异基因热图',
                               description='图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足25个，'
                                           '则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。')
        ###########
        if os.path.exists("7.Diffexp/enrichment/GO_enrichment"):
            ##功能富集
            ## GO 富集分析
            diffenrich_section2 = diff_section.add_section("差异基因GO富集分析")
            diff_group = os.listdir("7.Diffexp/enrichment/GO_enrichment")
            # diff_group.remove("enrichment_go.xls")
            diffenrich_section2.add_plot('src/enrich.png',
                                  caption='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式',
                                  content='超几何分布检验计算 p 值的公式和 Enrichment score 计算公式。其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。',
                                  description='其中，N 为所有基因中具有 GO 注释的基因数目；n 为 N 中差异表达基因中具有 GO 注释的基因数目；M 为所有基因中注释为某特定 GO Term 的基因数目；m 为注释为某特定 GO Term 的差异表达基因数目。可以根据 GO 分析的结果结合生物学意义从而挑选用于后续研究的基因。')

            diffenrich_section2.add_comment("GO 富集分析结果示例：")
            diffenrich_section2.add_table(f'7.Diffexp/enrichment/GO_enrichment/{diff_group[0]}/enrichment-go-*-Total.xls',
                                   caption='富集分析结果表格',
                                   show_rows=10)
            #diffenrich_section2.add_table('src/header_enrich_go.txt', caption='GO 富集分析结果各列说明')
            diff_gene_GO3_plots = list(glob('7.Diffexp/enrichment/GO_enrichment/*/GO.top.*.png'))
            if len(diff_gene_GO3_plots) > 0:
                diffenrich_section2.add_comment(
                    '''GO 富集分析 Top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序的各 10 条）条形图展示如下：''')
                diff_gene_GO3_contents = ['GO 富集分析结果展示。图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>Pvalue。'] * len(
                    diff_gene_GO3_plots)
                diffenrich_section2.add_plot(list(zip(diff_gene_GO3_plots, diff_gene_GO3_contents)),
                                      caption='GO 富集分析结果展示',
                                      description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1wS4y1z75L/)\n\n图片说明：图中横轴为 GO 条目名称，纵轴为 -log<sub>10</sub>Pvalue。')
                diff_gene_cohord_c3_plots = list(glob('7.Diffexp/enrichment/GO_enrichment/*/GO.*chord*.png'))
            if len(diff_gene_cohord_c3_plots) > 0:
                diffenrich_section2.add_comment('''富集分析和弦图，显示 p-value 最小的 10 个分类：''')
                diff_gene_cohord_c3_contents = [
                                                   'GO 富集 Top10 和弦图。图片说明：左侧为每个分类中 |log2FC| 最大的 10 个基因，右侧反映分类组成情况，中间线条表示分类、基因对应关系。外侧热图表示对应基因的 log2FC 值。'] * len(
                    diff_gene_cohord_c3_plots)
                diffenrich_section2.add_plot(list(zip(diff_gene_cohord_c3_plots, diff_gene_cohord_c3_contents)),
                                      caption='GO 富集 Top10 和弦图',
                                      description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1PF411j7Hq/)\n\n图片说明：左侧为每个分类中 |log2FC| 最大的 10 个基因，右侧反映分类组成情况，中间线条表示分类、基因对应关系。外侧热图表示对应基因的 log2FC 值。')
                diff_gene_GO_circle_plots = list(
                    glob('7.Diffexp/enrichment/GO_enrichment/*/enrichment-go-*-Total.circos.png'))
            if len(diff_gene_GO_circle_plots) > 0:
                diffenrich_section2.add_comment('''GO 富集分析圈图，显示 p-value 最小的 20 个分类，从外到内共四圈：''')
                diff_gene_GO_circle_plots_contents = [
                                                         '差异表达基因及所有基因在 GO 富集分析圈图。图片说明：第一圈：富集的分类，圈外为基因数目的坐标尺。不同的颜色代表不同的分类；第二圈：背景基因中该分类的数目以及 p-value。基因越多条形越长，值越小颜色越红，越大越蓝；第三圈：上下调基因比例条形图，浅红色代表上调基因比例，浅蓝色代表下调基因比例；下方显示具体的数值；第四圈：各分类的 RichFactor 值(该分类中前景基因数量除以背景基因数量)，背景辅助线每个小格表示0.2。'] * len(
                    diff_gene_GO_circle_plots)
                diffenrich_section2.add_plot(list(zip(diff_gene_GO_circle_plots, diff_gene_GO_circle_plots_contents)),
                                      caption='差异表达基因及所有基因在 GO 富集分析圈图',
                                      description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV15a411E7Le?spm_id_from=333.999.0.0)\n\n图片说明：第一圈：富集的分类，圈外为基因数目的坐标尺。不同的颜色代表不同的分类；第二圈：背景基因中该分类的数目以及 p-value。基因越多条形越长，值越小颜色越红，越大越蓝；第三圈：上下调基因比例条形图，浅红色代表上调基因比例，浅蓝色代表下调基因比例；下方显示具体的数值；第四圈：各分类的 RichFactor 值(该分类中前景基因数量除以背景基因数量)，背景辅助线每个小格表示 0.2。')

        ###########
        if os.path.exists("7.Diffexp/enrichment/KEGG_enrichment"):
            ##KEGG富集分析
            diffenrich_section3 = diff_section.add_section("差异基因KEGG富集分析")

            diff_gene_KEGG_top20_plots = list(glob('7.Diffexp/enrichment/KEGG_enrichment/*/KEGG.top.*.png'))
            if len(diff_gene_KEGG_top20_plots) > 0:
                diffenrich_section3.add_comment(
                    "KEGG 富集分析 Top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log<sub>10</sub>Pvalue 由大到小排序）气泡图如下：")
                diff_gene_KEGG_top20_contents = [
                                                    'KEGG 富集 Top20 气泡图。图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。'] * len(
                    diff_gene_KEGG_top20_plots)
                diff_gene_KEGG_top20_plots = list(zip(diff_gene_KEGG_top20_plots, diff_gene_KEGG_top20_contents))
                diffenrich_section3.add_plot(diff_gene_KEGG_top20_plots,
                                        caption='KEGG 富集 Top20 气泡图',
                                        description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV17a411E7nR?spm_id_from=333.999.0.0)\n\n图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由蓝-白-黄-红变化，其富集 p-value 值越小，显著程度越大。')
                diff_gene_cohord_c3_plots = list(glob('7.Diffexp/enrichment/KEGG_enrichment/*/KEGG.*chord*.png'))
            if len(diff_gene_cohord_c3_plots) > 0:
                diffenrich_section3.add_comment(
                    "富集分析和弦图，显示 p-value 最小的 10 个分类，图形分为左右两侧：左侧为每个分类中 |log2FC| 最大的 10 个基因，外侧热图表示对应基因的 log2FC 值：")
                diff_gene_cohord_c3_contents = ['KEGG富集 Top10 和弦图。图片说明：右侧反映分类组成情况，中间线条表示分类、基因对应关系。'] * len(
                    diff_gene_cohord_c3_plots)
                diff_gene_cohord_c3_plots = list(zip(diff_gene_cohord_c3_plots, diff_gene_cohord_c3_contents))
                diffenrich_section3.add_plot(diff_gene_cohord_c3_plots,
                                        caption='KEGG 富集 Top10 和弦图',
                                        description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1jY4y167D7?spm_id_from=333.999.0.0)\n\n图片说明：右侧反映分类组成情况，中间线条表示分类、基因对应关系。')

                diff_gene_KEGG_circle_plots = list(
                    glob('7.Diffexp/enrichment/KEGG_enrichment/*/enrichment-kegg-*-Total.circos.png'))
            if len(diff_gene_KEGG_circle_plots) > 0:
                diffenrich_section3.add_comment("KEGG 富集分析圈图，显示 p-value 最小的 20 个分类，从外到内共四圈：")
                diff_gene_KEGG_circle_contents = ['差异表达基因及所有基因在 KEGG 富集分析圈图图片说明：右侧反映分类组成情况，中间线条表示分类、基因对应关系。'] * len(
                    diff_gene_KEGG_circle_plots)
                diff_gene_KEGG_circle_plots = list(zip(diff_gene_KEGG_circle_plots, diff_gene_KEGG_circle_contents))
                diffenrich_section3.add_plot(diff_gene_KEGG_circle_plots,
                                        caption='差异表达基因及所有基因在 KEGG 富集分析圈图',
                                        description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1RB4y197nQ?spm_id_from=333.999.0.0)\n\n图片说明：第一圈：富集的分类，圈外为基因数目的坐标尺。不同的颜色代表不同的分类；第二圈：背景基因中该分类的数目以及 p-value。基因越多条形越长，值越小颜色越红，越大越蓝；第三圈：上下调基因比例条形图，浅红色代表上调基因比例，浅蓝色代表下调基因比例；下方显示具体的数值；第四圈：各分类的 RichFactor 值(该分类中前景基因数量除以背景基因数量)，背景辅助线每个小格表示 0.2。')

        #############
        if os.path.exists("7.Diffexp/ppi"):
            ppi_section = result_module.add_section('差异基因蛋白网络互作分析')
            report_hyperlink_list.append('[差异基因蛋白网络互作分析](#diffexp_ppi)\n\n')
            ppi_section.add_comment('''差异基因互作关系结果文件：''')
            ppi_table = list(glob('7.Diffexp/ppi/*_protein-protein-interaction.tsv'))[0]
            ppi_section.add_table(ppi_table,
                                 caption='差异基因互作关系表', show_rows=10)
            ppi_section.add_comment("上下调 Top25 差异表达基因互作网络图展示如下：")
            ppi_plots = list(glob('7.Diffexp/ppi/*string_protein-protein-interaction.new_colors.png'))
            ppi_plots = [str(i) for i in ppi_plots]
            ppi_contents = [
                               '差异基因互作网络图。图片说明：红色圆圈表示上调表达基因，绿色圆圈表示下调表达基因，节点之间的连线（或称为边）表示两蛋白之间具有相互作用，线的粗细表示相互作用关系的可靠性。'] * len(
                ppi_plots)
            ppi_plots = list(zip(ppi_plots, ppi_contents))
            ppi_section.add_plot(ppi_plots, caption='差异基因互作网络图',
                                description='图片说明：红色圆圈表示上调表达基因，绿色圆圈表示下调表达基因，节点之间的连线（或称为边）表示两蛋白之间具有相互作用，线的粗细表示相互作用关系的可靠性。')
            ppi_section.add_comment("上下调 Top25 差异表达基因互作圆形图展示如下：")
            ppi_plots = list(glob('7.Diffexp/ppi/*top_25_ppi_network.png'))
            ppi_plots = [str(i) for i in ppi_plots]
            ppi_contents = ['差异基因互作圆形图。图片说明：红色表示上调差异表达基因，蓝色表示下调差异表达基因；关联的基因越多，基因点越大。'] * len(ppi_plots)
            ppi_plots = list(zip(ppi_plots, ppi_contents))
            ppi_section.add_plot(ppi_plots, caption='差异基因互作圆形图',
                                description='[点击此处查看完整图片解读视频](https://www.bilibili.com/video/BV1X54y1Z7XW?spm_id_from=333.999.0.0)\n\n图片说明：红色表示上调差异表达基因，蓝色表示下调差异表达基因；关联的基因越多，基因点越大。')

    ##========================================================= 4.项目摘要、链接========================================
    Project_module.add_section(name="项目摘要", description=''.join(report_summary_list))
    report_hyperlink = ''.join(report_hyperlink_list)
    Project_module.add_section(name="项目快捷链接", description=report_hyperlink)
    ##============================================================= 5.附录 ============================================
    description_module = report.add_module('附录')
    refgenome = config['database_url']['Ref_genome']['Linkage']
    refgenome_annotation = config['database_url']['Ref_genome_annotation']['Linkage']
    database = pd.DataFrame(
        [['GO Database', 'http://geneontology.org/'],
         ['KEGG Database', 'http://www.genome.jp/kegg/'],
         ['Reference genome', refgenome],
         ['Reference genome annotation', refgenome_annotation]
         ])

    database.columns = ['使用数据库', '网页链接']
    description_section5 = description_module.add_section('数据库信息', description='')
    supply_section_database_table = description_section5.add_table(database, caption='数据库信息')

    database2 = pd.DataFrame(
        [['STAR', '2.7.10b'],
         ['Seurat', '4.0.1'],
         ['SeuratDisk', '0.0.0.9019'],
         ['SingleR', '1.4.1'],
         ['presto', '1.0.0'],
         ['Biobase', '2.50.0']
         ])
    database2.columns = ['软件', '版本']
    description_section6 = description_module.add_section('数据分析软件')
    software_table = description_section6.add_table(database2, caption='数据分析软件')

    ##============================================================== 6.申明 ===========================================
    affirming_module = report.add_module('申明')

    ############################################ Generate Report HTML ###################################################
    report.write_to('Report.html')


if __name__ == "__main__":
    multimodal_report()

