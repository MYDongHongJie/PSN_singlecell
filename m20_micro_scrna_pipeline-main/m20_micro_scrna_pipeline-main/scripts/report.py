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
from oebio.app import *
from oebio.report import Report, oeweb_register


@click.command()
@click.option('-i', '--inputdir', prompt='the program path',
              help='the directory of your program. ')
@click.option('-c', '--configfile', prompt='config file', default='config/config.yaml',
              help='the config file which contain your program information.')
def multimodal_report(inputdir, configfile):
    ##=============================================== 1. 读取项目配置文件信息 ===========================================
    """ python multimodal_report.py  -i ./ -c ./config.yaml"""
    cfg = open(configfile, 'r', encoding='utf-8').read()
    config = yaml.full_load(cfg)
    program_path = os.path.abspath(inputdir)
    # paragram from config
    program_num = dict(config['report'])['Project_Num']
    report_time = time.strftime("%Y_%m_%d")
    outdir = f"report/{program_num}_Report_{report_time}"
    log_file = open("%s/report_log.txt" % (program_path), "w")
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    # Report
    os.chdir(program_path)
    if os.path.exists(outdir):
        shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    ##================================================2. 生成结果文件 =================================================

    # 1.STARsolo
    if os.path.exists("%s/1.STARsolo" % (program_path)):
        os.makedirs("%s/1.STARsolo" % (outdir))
        # name = [str(i).split('/')[-1] for i in glob("1.STARsolo/*/*_Solo.out")]
        # for j in name:
        subprocess.call('cp -r %s/1.STARsolo/*/*_Solo.out           %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/BarcodeRank_plot       %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/Read_Distribution      %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/STARsolo.statistic.tsv %s/1.STARsolo/' % (program_path, outdir),
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

    # 4.Annotation
    if os.path.exists("%s/4.Annotation" % (program_path)):
        os.makedirs("%s/4.Annotation" % (outdir))
        subprocess.call('cp -r %s/4.Annotation/* %s/4.Annotation' % (program_path, outdir), shell=True)
        subprocess.call('rm %s/4.Annotation/*/*.xml' % (outdir), shell=True)
        subprocess.call('rm %s/4.Annotation/*/*.R' % (outdir), shell=True)
        subprocess.call('rm %s/4.Annotation/*/*_tmp' % (outdir), shell=True)
    else:
        print("Can not find 4.Annotation results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 4.Annotation results!!!!!!!!!!!!!!!!!" + "\n")

    # 5.Marker
    if os.path.exists("%s/5.Marker" % (program_path)):
        os.makedirs("%s/5.Marker" % (outdir))
        subprocess.call('cp -r %s/5.Marker/* %s/5.Marker' % (program_path, outdir), shell=True)
    else:
        print("Can not find 5.Marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 5.Marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 6.MarkerEnrich
    if os.path.exists("%s/6.MarkerEnrich" % (program_path)):
        os.makedirs("%s/6.MarkerEnrich" % (outdir))
        subprocess.call('cp  -r %s/6.MarkerEnrich/* %s/6.MarkerEnrich' % (program_path, outdir), shell=True)
    else:
        print("Can not find 6.MarkerEnrich results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 6.MarkerEnrich results!!!!!!!!!!!!!!!!!" + "\n")

    # 7.Diffgene
    if os.path.exists("%s/7.Diffgene" % (program_path)):
        os.makedirs("%s/7.Diffgene" % (outdir))
        subprocess.call('cp -r %s/7.Diffgene/* %s/7.Diffgene' % (program_path, outdir), shell=True)
    else:
        print("Can not find 7.Diffgene results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 7.Diffgene results!!!!!!!!!!!!!!!!!" + "\n")

    # 8.DiffgeneEnrich
    if os.path.exists("%s/8.DiffgeneEnrich" % (program_path)):
        os.makedirs("%s/8.DiffgeneEnrich" % (outdir))
        subprocess.call('cp -r %s/8.DiffgeneEnrich/* %s/8.DiffgeneEnrich/' % (program_path, outdir), shell=True)
        subprocess.call('cp -r %s/8.DiffgeneEnrich/* %s/8.DiffgeneEnrich/' % (program_path, outdir), shell=True)
    else:
        print("Can not find 8.DiffgeneEnrich results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 8.DiffgeneEnrich results!!!!!!!!!!!!!!!!!" + "\n")
    # report
    if not os.path.exists("%s/src" % (outdir)):
        os.makedirs("%s/src" % (outdir))
        subprocess.call(f'cp {scriptdir}/src/* {outdir}/src/', shell=True)

    ##===================================================== 3. 生成报告================================================
    header_info = {key: config["report"][key] for key in config["report"] if
                   key in ["Project_Num", "Customer", "Species", "Executor", "Sample_Num", "Task_Num"]}
    header_info["项目编号"] = header_info.pop("Project_Num")
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["实验物种"] = header_info.pop("Species")
    header_info["执行编号"] = header_info.pop("Executor")
    header_info["样本数目"] = header_info.pop("Sample_Num")
    header_info["任务单号"] = header_info.pop("Task_Num")

    report = Report('M20 单细胞微生物转录组结题报告', title='M20 单细胞微生物转录组结题报告',
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
    ## ======================================================= 1. 实验及分析概述============================================
    Background_module = report.add_module('技术简介')
    Background_section1 = Background_module.add_section('实验流程')  ## 实验流程
    workflow_module = Background_module.add_section('生信分析流程')  ## 分析流程
    # ====================================================== 3. 分析结果 ================================================
    result_module = report.add_module('项目分析结果')
    ##====================================================== 3.1 质控分析 =================================================
    if os.path.exists("1.STARsolo"):
        section1 = result_module.add_section('基因定量质控')
        report_hyperlink_list.append('[STARsolo分析](#STAR)\n\n')
        section1_1 = section1.add_section('STARsolo分析')
        section1_1.add_table(f'1.STARsolo/STARsolo.statistic.tsv',
                             caption='样本STARsolo结果汇总统计表',
                             headerdict={'sample': '样本名称',
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

    ##====================================================== 3.3 基因及样本过滤 ==============================================
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
                f"本次分析共完成 {sample_num} 个样本的M20 单细胞微生物转录组测序，数据经STAR质控过滤后的高质量细胞数为 {min_beforeQC_cell}"
                f" 个，经剔除异常基因表达以及低基因表达的细胞后，最终获得的细胞数目为 {min_afterQC_cell} 个，每个细胞中的平均 UMI 数为"
                f" {min_nCount_RNA_afterQC}，每个细胞中的平均基因数为 {min_nFeature_afterQC}")
        else:
            report_summary_list.append(
                f"本次分析共完成 {sample_num} 个样本的M20 单细胞微生物转录组测序，各样本经STAR质控过滤后的高质量细胞数分布在"
                f" {min_beforeQC_cell}~{max_beforeQC_cell} 个，经剔除异常基因表达以及低基因表达的细胞后，最终获得的细胞数目分布在"
                f" {min_afterQC_cell}~{max_afterQC_cell} 个，每个细胞中的平均 UMI 数分布在 {min_nCount_RNA_afterQC}~{max_nCount_RNA_afterQC}，"
                f"每个细胞中的平均基因数分布在 {min_nFeature_afterQC}~{max_nFeature_afterQC}")

        if "rRNA.mito_higher" in summ.columns:
            min_percent_mito_afterQC = format(min(summ["mean_rRNA.mito_afterQC"].astype(float)), '.4f')
            max_percent_mito_afterQC = format(max(summ["mean_rRNA.mito_afterQC"].astype(float)), '.4f')
            rRNA_txt = f"、rRNA基因所占比例（rRNA.mito）"
            rRNA_descript = f"、rRNA基因比例过高"
            if len(summ) == 1:
                report_summary_list.append(f"，每个细胞中平均rRNA基因比例为 {min_percent_mito_afterQC}")
            else:
                report_summary_list.append(
                    f"，每个细胞中平均rRNA基因比例分布在 {min_percent_mito_afterQC}~{max_percent_mito_afterQC}")
        else:
            rRNA_descript = ""
            rRNA_txt = ""
        if "tRNA.mito_higher" in summ.columns:
            min_percent_mito_afterQC = format(min(summ["mean_tRNA.mito_afterQC"].astype(float)), '.4f')
            max_percent_mito_afterQC = format(max(summ["mean_tRNA.mito_afterQC"].astype(float)), '.4f')
            tRNA_txt = f"、tRNA基因所占比例（tRNA.mito）"
            tRNA_descript = f"、tRNA基因比例过高"
            if len(summ) == 1:
                report_summary_list.append(f"，每个细胞中平均tRNA基因比例为 {min_percent_mito_afterQC}")
            else:
                report_summary_list.append(
                    f"，每个细胞中平均tRNA基因比例分布在 {min_percent_mito_afterQC}~{max_percent_mito_afterQC}")
        else:
            tRNA_txt = ""
            tRNA_descript = ""
        report_summary_list.append(f"。")


        ##
        QC_description = f"在 STAR 初步质控的基础上进一步利用Seurat[^seurat]包对实验数据进行质控。理论上大部分细胞中表达的基因数量、" \
                         f"UMI 数量、rRNA基因以及tRNA基因数量等会集中分布在某一区域内，根据它们的分布特征可以拟合分布模型，使用该模型找" \
                         f"到其中的离域值，剔除异常数据。本项目中的质控标准为：1).去除在所有细胞中均未表达的基因，2).过滤细胞中少于10个" \
                         f"基因表达{rRNA_descript}{tRNA_descript}的细胞，以同时满足以上标准的细胞作为高质量细胞，进行下游分析。\n" \
                         f" [^seurat]: Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data[J]. Cell, 2021, 184(13):3573-3587.\n" \
                         f"\n 质控前后每个细胞的基因数量（nFeature_RNA）、UMI 数量（nCount_RNA）{rRNA_txt}{tRNA_txt}的小提琴图展示如下："
        ##
        filter_section2 = filter_section.add_section(name='基因及样本过滤结果', description=QC_description)
        filter_section2.add_comment("详细结果见目录：[2.Count_QC](./2.Count_QC)")
        filter_section2.add_plot('2.Count_QC/QC_metrics_beforeQC.png', caption='过滤低质量细胞前各质控指标的小提琴图')
        filter_section2.add_comment(
            "过滤低质量细胞前各质控指标的小提琴图：[2.Count_QC/QC_metrics_beforeQC.pdf](./2.Count_QC/QC_metrics_beforeQC.pdf)")
        filter_section2.add_plot('2.Count_QC/QC_metrics_afterQC.png', caption='过滤低质量细胞后各质控指标的小提琴图')
        filter_section2.add_comment(
            "过滤低质量细胞后各质控指标的小提琴图：[2.Count_QC/QC_metrics_afterQC.pdf](./2.Count_QC/QC_metrics_afterQC.pdf)")
        filter_section2.add_comment("定量质控前后的细胞数统计情况及过滤标准如下表所示：")
        filter_section2.add_table(summ.iloc[:, :19],
                                  caption='质控前后细胞数目统计表',
                                  headerdict={'sample': '样本名',
                                              'mean_nCount_RNA_beforeQC': '质控前RNA平均Counts数',
                                              'median_nCount_RNA_beforeQC': '质控前RNA Counts的中位数',
                                              'mean_nFeature_RNA_beforeQC': '质控前RNA平均Feature数',
                                              'median_nFeature_RNA_beforeQC': '质控前RNA feature的中位数',
                                              'mean_rRNA.mito_beforeQC': '质控前rRNA基因比例平均数',
                                              'median_rRNA.mito_beforeQC': '质控前rRNA基因比例中位数',
                                              'mean_tRNA.mito_beforeQC': '质控前tRNA基因比例平均数',
                                              'median_tRNA.mito_beforeQC': '质控前tRNA基因比例中位数',
                                              'Total_cells_beforeQC': '质控前细胞总数',
                                              'mean_nCount_RNA_afterQC': '质控后RNA平均counts数',
                                              'median_nCount_RNA_afterQC': '质控后RNA Counts的中位数',
                                              'mean_nFeature_RNA_afterQC': '质控后RNA平均Feature数',
                                              'median_nFeature_RNA_afterQC': '质控后RNA feature的中位数',
                                              'mean_rRNA.mito_afterQC': '质控后rRNA基因比例平均数',
                                              'median_rRNA.mito_afterQC': '质控后rRNA基因比例中位数',
                                              'mean_tRNA.mito_afterQC': '质控后tRNA基因比例平均数',
                                              'median_tRNA.mito_afterQC': '质控后tRNA基因比例中位数',
                                              'Total_cells_afterQC': '质控后细胞总数'
                                              })

    ##=============================================== 3.4 降维聚类 ==============================================
    ## 降维与聚类
    if os.path.exists("3.Clustering"):
        cluster_section = result_module.add_section('数据降维与聚类分析')
        report_hyperlink_list.append('[降维与聚类分析结果](#rna_cluster)\n\n')
        if len(summ) == 1:
            sample_des = f"仅有1个样本，并不存在批次效应"
        else:
            sample_des = f"有不同样本，但并不存在批次效应"

        if config['params']['Clustering']["RNA_reduct1"] == "pca":
            reduct1_des = f"PCA(principal component analysis, 主成分分析)[^pca]的降维算法进行降维 "
            reduct1_name = "PCA"
        else:
            sample_des = f"有不同样本，且存在批次效应"
            reduct1_des = f"MNN(mutual nearest neighbors, 互享最近邻)的降维算法剔除批次效应"
            reduct1_name = "MNN"
        if config['params']['Clustering']["RNA_reduct2"] == "tsne":
            reduct2_des = f"t-SNE(t分布式随机邻近嵌入)[^tsne]算法进行二级降维。t-SNE创建了一个缩小的特征空间，相似的样本由附近的点建模，不相似的样本由高概率的远点建模，更好的展示了样本之间的局部关系"
            reduct2_name = "t-SNE"
        else:
            reduct2_des = f"UMAP[^umap](统一流形逼近与投影)算法进行二级降维。UMAP 主要基于流形理论和拓扑算法的理论，对高维数据进行降维，从而能够保留更多数据的全局结构，并且具有优越的运行性能"
            reduct2_name = "UMAP"

        cluster_section2 = cluster_section.add_section("降维与聚类分析结果")
        cluster_section2.add_comment(
            f"由于本项目{sample_des}，因此我们首先采用了{reduct1_des}，再通过 {reduct2_des}。基于{reduct1_name}降维结果通过{reduct2_name}对单细胞群聚类进行可视化，聚类算法采用SNN，最终获得最优细胞分群。\n [^pca]: Ian T. Jolliffe. PRINCIPAL COMPONENT ANALYSIS: A BEGINNER'S GUIDE — I. Introduction and application[J]. Weather, 1990, 45.\n [^tsne]: van der Maaten, Laurens. Accelerating t-SNE using Tree-Based Algorithms[J]. Journal of Machine Learning Research, 2014, 15(1):3221-3245. \n [^umap]: Mcinnes L,Healy J,Melville J.UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction[J].The Journal of Open Source Software,2018. \n")
        cluster_section2.add_comment("详细结果见目录：[3.Clustering](./3.Clustering)")

        cluster_section2.add_plot('3.Clustering/*/*_groupby_cluster_resolution*_plot.png',
                                  caption='基因表达水平降维与聚类图',
                                  description='图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')
        cluster_section2.add_comment(
            "基因表达水平降维与聚类图：[3.Clustering/\*/\*_groupby_clusters_resolution*_plot.png](./3.Clustering/)")

        ## 样本分组展示
        if (sample_num > 1):
            cluster_section3 = cluster_section.add_section('样本间降维聚类分组展示')
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/groupby-sampleid_contrast_plot.png',
                                      caption='多样本降维聚类分组展示图',
                                      description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。')
            cluster_section3.add_comment(
                "基因表达水平多个样本的降维聚类分组展示图：[3.Clustering/visualize_cluster_by_\*/\*groupby-sampleid_contrast_plot.pdf](./3.Clustering/)")
            ################
            cluster_section3.add_comment(" ")
            cluster_section3.add_comment("每个细胞群中样本占比的柱状统计图如下：")
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/groupby-*.res.*_summary_plot.png',
                                      caption='每个细胞群中样本占比的柱状统计图',
                                      description='图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。')

            ################
            cluster_section3.add_comment(" ")
            cluster_section3.add_comment("多个样本的降维聚类分面展示图如下：")
            cluster_section3.add_plot('3.Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.png',
                                      caption='多样本降维聚类分面展示图',
                                      description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同群的细胞以不同颜色区分。')
            cluster_section3.add_comment(
                "多个样本的降维聚类分面展示图如下：[3.Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.pdf](./3.Clustering/)")

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
            cluster_section3.add_comment(
                "各样本在不同细胞群中的细胞数目统计表：[3.Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls](./3.Clustering/)")

    ##=============================================== 3.5 基因功能注释 ==============================================
    ##基因注释
    if os.path.exists("4.Annotation"):
        annotation_section = result_module.add_section('基因功能注释')
        report_hyperlink_list.append('[基因功能注释](#rna_annotation)\n\n')
        ##KEGG注释结果
        annotation_section2 = annotation_section.add_section('KEGG注释结果')
        annotation_section2.add_comment("KEGG注释详细结果文件见目录：[4.Annotation/KEGG](./4.Annotation/KEGG/)")
        annotation_section2.add_comment(
            "使用diamond[^diamond]与KEGG数据库比对之后，对KEGG注释结果汇总，统计结果如下表所示：\n [^diamond]: Buchfink B, Reuter K, Drost HG, Sensitive protein alignments at tree-of-life scale using DIAMOND[J], Nature Methods 18, 366–368 (2021). \n")
        annotation_section2.add_table("4.Annotation/KEGG/*.KEGG.gene.anno.xls",
                                      caption='KEGG数据库注释结果',
                                      show_rows=5)
        annotation_section2.add_comment(
            "根据功能分级，通常将 KEGG 分为三个层级，level1 包含六个分类：Metabolism、Genetic Information Processing、Environmental Information Processing、Cellular Processes、Organismal Systems 和 Human Diseases（具体物种注释可能有删减）。level2 包含 Cell growth and death、Transcription 和 Development 等 44 个分类（具体物种注释可能有删减），level3 即为常规富集使用的数百个 Pathway，从 level1 到 level3 功能更具体，反之，更概括。")
        annotation_section2.add_comment("所有基因在 KEGG Level2 水平分布图如下：")
        annotation_section2.add_plot('4.Annotation/KEGG/*.classification.png',
                                     caption='所有基因在 KEGG Level2 水平分布图',
                                     description='图片说明：横轴表示gene数量，纵轴表示Level2 pathway通路名，柱子右边数字代表注释到该Level2 pathway下的gene数量。')
        annotation_section2.add_table("4.Annotation/KEGG/*.KEGG.classification.xls",
                                      caption='KEGG Level2 水平注释结果',
                                      show_rows=5)
        ##GO注释结果
        annotation_section3 = annotation_section.add_section('GO注释结果')
        annotation_section3.add_comment("GO注释详细结果见目录：[4.Annotation/GO](./4.Annotation/GO/)")
        annotation_section3.add_comment("使用diamond与GO数据库比对之后，对GO注释结果汇总，统计结果如下表所示：")
        annotation_section3.add_table("4.Annotation/GO/*.GO.gene.anno.xls",
                                      caption='GO数据库注释结果',
                                      show_rows=5)
        annotation_section3.add_comment("GO功能分类图如下：")
        annotation_section3.add_plot('4.Annotation/GO/*.GO.classification.stat.png',
                                     caption='GO功能分类图',
                                     description='图片说明：横轴表示GO功能分类，左边纵轴表示注释到该类gene数量占比，右边纵轴表示注释到该类gene数量。')
        annotation_section3.add_table("4.Annotation/GO/*.GO.classification.stat.xls",
                                      caption='GO功能分类结果统计表',
                                      show_rows=5)
        ##CARD注释结果
        annotation_section4 = annotation_section.add_section('CARD注释结果')
        annotation_section4.add_comment("CARD注释详细结果见目录：[4.Annotation/CARD](./4.Annotation/CARD/)")
        annotation_section4.add_table("4.Annotation/CARD/*.anno.xls",
                                      caption='CARD数据库注释结果统计表',
                                      show_rows=5)
        annotation_section4.add_comment("CARD数据库注释top10 term展示图如下：")
        annotation_section4.add_plot('4.Annotation/CARD/*.png',
                                     caption='CARD数据库注释top10 term展示图',
                                     description='图片说明：不同颜色的扇形大小代表注释到该term的基因数量占比，左上角图例后面的数字是注释到该term的gene数量。')
        ##CAZy注释结果
        annotation_section5 = annotation_section.add_section('CAZy注释结果')
        annotation_section5.add_comment("CAZy注释详细结果见目录：[4.Annotation/CAZy](./4.Annotation/CAZy/)")
        annotation_section5.add_comment("CAZy注释分类分布图如下：")
        annotation_section5.add_plot('4.Annotation/CAZy/CAZy_class_plot.png',
                                     caption='CAZy注释class分布图',
                                     description='图片说明：横坐标是CAZy分类，纵坐标是注释到对应分类的基因数量。')
    ##=============================================== 3.6 Marker基因鉴定及功能富集 ============================================
    if os.path.exists("5.Marker"):
        marker_section = result_module.add_section('Marker基因鉴定及功能富集')
        report_hyperlink_list.append('[Marker 基因鉴定结果](#rna_marker_gene)\n\n')
        test_method = config["params"]['Marker']["RNA_test_method"]
        Marker_description = f"Marker 基因的定义是该基因在指定细胞群的绝大多数细胞中有较高的表达，而在其余细胞类群中只有少部分表达，且该基因在此细胞群相对于其他细胞群中是显著上调表达。\n 本研究采用 {test_method} 检验方法对指定细胞群与其余所有细胞群进行差异检验，先设置log2FC > 0 ,min.pct > 0.05导出所有可能的marker,然后按照P值，FoldChange值以及gene_diff值筛选top10 marker作为每个细胞群的特异性 Marker 基因。"
        marker_section1 = marker_section.add_section(name='Marker 基因鉴定', description=Marker_description)
        #####
        marker_result = pd.read_csv('5.Marker/top10_markers_for_each_cluster.xls', sep='\t')
        cluster_num = max(marker_result['cluster'])
        report_summary_list.append(f"通过基因表达数据进行降维聚类后共分为 {cluster_num} 群细胞，")
        ####
        marker_section1.add_comment(
            "每个细胞群中 Top10 Marker 基因结果表格：[5.Marker/top10_markers_for_each_cluster.xls](./5.Marker/top10_markers_for_each_cluster.xls)")
        marker_section1.add_table("5.Marker/top10_markers_for_each_cluster.xls",
                                  caption='每个细胞群 Top10 Marker 基因列表',
                                  show_rows=5)
        ####
        marker_section1.add_comment(
            "Top10 Marker 基因表达热图：[5.Marker/topmarker_gene_heatmap.pdf](./5.Marker/topmarker_gene_heatmap.pdf)")
        marker_section1.add_plot("5.Marker/topmarker_gene_heatmap.png",
                                 caption='Top10 Marker 基因表达热图',
                                 description='图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')
        ####
        marker_section1.add_comment(
            "Top10 Marker 基因在降维聚类结果中的可视化图：[5.Marker/*/marker_gene_featureplot.pdf](./5.Marker/)")
        marker_section1.add_plot("5.Marker/*/marker_gene_featureplot.png",
                                 caption='Top10 Marker 基因在降维聚类结果中的可视化图',
                                 description='图片说明：红色越深表示该细胞中对应基因的表达量越高。')
        ####
        marker_section1.add_comment(
            "Top10 Marker 基因表达量小提琴图：[5.Marker/*/marker_gene_violin_plot.pdf](./5.Marker/)")
        marker_section1.add_plot("5.Marker/*/marker_gene_violin_plot.png",
                                 caption='Top10 Marker 基因表达量小提琴图',
                                 description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')
        #### Mark基因富集分析
        if os.path.exists("6.MarkerEnrich"):
            marker_section2 = marker_section.add_section('Marker 基因富集')
            report_hyperlink_list.append('[Marker 基因富集结果](#rna_marker_enrich)\n\n')
            marker_section2.add_comment("Marker基因富集结果文件见：[6.MarkerEnrich](./6.MarkerEnrich/)")
            ### Mark基因KEGG富集分析
            marker_section3 = marker_section2.add_section('KEGG富集分析')
            marker_section3.add_comment(
                "KEGG 富集分析 top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log10Pvalue 由大到小排序）气泡图如下：")
            marker_section3.add_plot("6.MarkerEnrich/KEGG_enrichment/*/KEGG.top.Total.png",
                                     caption='KEGG富集 top20 气泡图',
                                     description='图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。')
            ### Mark基因GO富集分析
            marker_section4 = marker_section2.add_section('GO富集分析')
            marker_section4.add_comment(
                "GO 富集分析 top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log10pValue 由大到小排序的各 10 条）条形图展示如下：")
            marker_section4.add_plot("6.MarkerEnrich/GO_enrichment/*/GO.top.Total.png",
                                     caption='GO富集分析结果展示',
                                     description='图片说明：图中横轴为 GO 条目名称，纵轴为 -log10pValue。')
    ##=============================================== 3.7 差异基因筛选及功能富集 ============================================
    if os.path.exists("7.Diffgene"):
        diff_section = result_module.add_section('差异基因筛选及富集')
        report_hyperlink_list.append('[差异表达基因筛选及富集](#diffexp)\n\n')
        diffexp_folder = "7.Diffgene"
        diffexp_file = list(glob('%s/diffexp_results_stat.xls' % (diffexp_folder)))[0]
        diff = pd.read_csv(diffexp_file, sep='\t')
        num_diff = len(diff)
        diff_num = ', '.join([str(i) for i in diff.iloc[:, 4]])
        report_summary_list.append(f"共设有{num_diff}个差异基因分组，其检测到的差异基因数量为：{diff_num}。")
        test_method = config['params']['Diffexp']['test']
        if test_method == 'presto':
            diff_descript = f'差异表达基因使用presto包的WilcoxAUC函数进行筛选，默认使用 presto 差异检验方法。'
        else:
            diff_descript = f'差异表达基因使用Seurat包的FindMarkers函数进行筛选，默认使用 {test_method} 差异检验方法。'
        ##差异基因
        diff_section1 = diff_section.add_section(name='差异表达基因筛选', description=diff_descript)
        diff_section1.add_comment('差异表达基因具体结果见: [7.Diffgene](./7.Diffgene)')
        diff_section1.add_comment("各分组差异基因数目统计表如下：")
        diff_section1.add_table("7.Diffgene/diffexp_results_stat.xls", caption='差异表达基因统计表')
        diff_section1.add_comment("差异显著基因结果示例：")
        difffiles = list(glob(f'{diffexp_folder}/*diff*pval*xls'))
        diff_section1.add_table(difffiles[0],
                                caption='差异显著基因结果表格',
                                show_rows=5)
        diff_section1.add_comment("将差异基因的差异倍数（FoldChange）从大到小排列，上下调各选取 25 个基因绘制热图：")
        # topdiff_heatmap_plots = list(glob(f'{diffexp_folder}/*_heatmap.png'))
        # topdiff_heatmap_plots = [str(i) for i in topdiff_heatmap_plots]
        # topdiff_heatmap_prefix = [re.findall("*top(.*)_heatmap.png", i)[0] for i in topdiff_heatmap_plots]
        # topdiff_heatmap_contents = [
        #     j + ' 上下调 Top25 差异基因热图。图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足20个，则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。'
        #     for j in topdiff_heatmap_prefix]
        # list(zip(topdiff_heatmap_plots, topdiff_heatmap_contents)),
        diff_section1.add_plot( f'{diffexp_folder}/*_heatmap.png',
                               caption='上下调Top 25差异基因热图',
                               description='图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足25个，'
                                           '则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。')
        if os.path.exists("8.DiffgeneEnrich"):
            ##功能富集
            ## GO 富集分析
            diffenrich_section2 = diff_section.add_section("差异基因GO富集分析")
            diffenrich_section2.add_comment(
                "GO富集分析具体结果文件为：[8.DiffgeneEnrich/GO_enrichment/](./8.DiffgeneEnrich/GO_enrichment/)")
            diffenrich_section2.add_comment("GO富集分析汇总表如下：")
            diffenrich_section2.add_table("8.DiffgeneEnrich/GO_enrichment/enrichment_go.xls",
                                          caption='GO富集分析汇总表')
            diffenrich_section2.add_comment("GO富集分析结果如下：")
            diffenrich_table = list(glob("8.DiffgeneEnrich/GO_enrichment/*/enrichment-go-*_*-vs-*-Total.xls"))[0]
            diffenrich_section2.add_table(diffenrich_table,
                                          caption='GO富集分析结果',
                                          show_rows=5)
            diffenrich_section2.add_comment(
                "GO 富集分析 top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log10pValue 由大到小排序的各 10 条）条形图展示如下")
            diffenrich_section2.add_plot("8.DiffgeneEnrich/GO_enrichment/*/GO.top.*.png",
                                         caption='GO富集分析结果展示',
                                         description='图片说明：图中横轴为 GO 条目名称，纵轴为 -log10pValue。')
            diffenrich_section2.add_comment(
                "根据功能分级，一般将 GO 分为三个层级，level1 包含三个条目：biological process、cellular component和molecular function，level2 包含 biological adhesion、cell 和 binding 等 64 个条目，level3 即为常规富集使用的数万个条目。从 level1到 level3 功能更具体，反之，更概括。")
            diffenrich_section2.add_comment("上调差异基因和下调差异基因在 GO Level2 水平分布比较图如下：")
            diffenrich_section2.add_plot("8.DiffgeneEnrich/GO_enrichment/*/GO.level2.stat.Total.png",
                                         caption='上调差异基因和下调差异基因在 GO Level2 水平分布比较图',
                                         description='图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')
            ##KEGG富集分析
            diffenrich_section3 = diff_section.add_section("差异基因KEGG富集分析")
            diffenrich_section3.add_comment(
                "KEGG是有关 Pathway 的主要公共数据库，利用KEGG数据库对差异蛋白编码基因进行 Pathway 分析（结合 KEGG 注释结果），并用超几何分布检验的方法计算每个 Pathway 条目中差异基因富集的显著性。")
            diffenrich_section3.add_comment(
                "差异基因KEGG富集分析结果文件如下：[8.DiffgeneEnrich/KEGG_enrichment](./8.DiffgeneEnrich/KEGG_enrichment)")
            diffenrich_section3.add_comment(
                "KEGG 富集分析 top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log10Pvalue 由大到小排序）气泡图如下：")
            diffenrich_section3.add_plot("8.DiffgeneEnrich/KEGG_enrichment/*/KEGG.top.*.png",
                                         caption='KEGG富集 top20 气泡图',
                                         description='图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。')
            diffenrich_section3.add_comment(
                "根据功能分级，通常将 KEGG 分为三个层级，level1 包含六个分类：Metabolism、Genetic Information Processing、Environmental Information Processing、Cellular Processes、Organismal Systems 和 Human Diseases（具体物种注释可能有删减）。level2 包含 Cell growth and death、Transcription 和 Development 等 44 个分类（具体物种注释可能有删减），level3 即为常规富集使用的数百个 Pathway，从 level1 到 level3 功能更具体，反之，更概括。")
            diffenrich_section3.add_comment("差异表达基因及所有基因在 KEGG Level2 水平分布比较图如下：")
            diffenrich_section3.add_plot("8.DiffgeneEnrich/KEGG_enrichment/*/ALL_vs_DEG.KEGG_Classification.png",
                                         caption='差异表达基因及所有基因在 KEGG Level2 水平分布比较图',
                                         description='图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。')
            diffenrich_section3.add_comment("上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图如下：")
            diffenrich_section3.add_plot("8.DiffgeneEnrich/KEGG_enrichment/*/Up_vs_Down.KEGG_Classification.png",
                                         caption='上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图',
                                         description='图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。')

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
         ['CARD', 'https://card.mcmaster.ca/'],
         ['CAZy', 'http://www.cazy.org/'],
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
         ['diamond', '2.0.15'],
         ['HMMER', '3.2.1']
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
