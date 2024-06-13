#!/usr/bin/env python37
# encoding:utf-8

# Report
# 1.Xenium
# 2.Count_QC
# 3.Clustering
# 4.Visualize_cluster_by_clusters
# 5.Marker
# 6.Reference_Celltype


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
    outdir = f"report/{program_num}_QCReport"
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

    # 1.Xenium
    if os.path.exists("%s/1.Xenium" % (program_path)):
        os.makedirs("%s/1.Xenium" % (outdir))
        name = [str(i).split('/')[-2] for i in glob("1.Xenium/*/outs")]
        for j in name:
            os.makedirs("%s/1.Xenium/%s" % (outdir, j))
            subprocess.call('cp %s/1.Xenium/%s/*.png %s/1.Xenium/' % (program_path,j,outdir),shell=True)
            #subprocess.call('ln -s %s/1.Xenium/%s/outs/* %s/1.Xenium/%s ' % (program_path, j, outdir, j),shell=True)
        subprocess.call('cp %s/1.Xenium/summary_all.csv %s/1.Xenium/' % (program_path, outdir), shell=True)
    else:
        print("Can not find 1.Xenium results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 1.Xenium results!!!!!!!!!!!!!!!!" + "\n")

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

    report = Report('10x Xenium 组织原位分析项目质控报告', title='10x Xenium 组织原位分析项目质控报告',
                    header_info=dict(header_info))
    report.add_yaml_config(os.path.abspath(f'{scriptdir}/report.yaml'))
    print(f'{scriptdir}/report.yaml')
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
    # ====================================================== 4. 分析结果 =================================================
    result_module = report.add_module('项目分析结果')
    ##====================================================== 4.1 质控分析 ================================================
    if os.path.exists("1.Xenium"):
        section1 = result_module.add_section('基因定量质控')
        report_hyperlink_list.append('[Xenium Ranger分析](#xenium)\n\n')
        section1_1 = section1.add_section('Xenium Ranger分析')
        summ_all = pd.read_table("1.Xenium/summary_all.csv", sep=",")
        summ_all = summ_all[["region_name",
                             "num_cells_detected",
                             "cells_per_100um2",
                             "median_genes_per_cell",
                             "median_transcripts_per_cell",
                             "total_high_quality_decoded_transcripts",
                             "decoded_transcripts_per_100um2",
                             "fraction_transcripts_decoded_q20",
                             "adjusted_negative_control_probe_rate"]]
        section1_1.add_plot(f'1.Xenium/*.png',
                            caption='Xenium Ranger分析结果概述', content='Xenium Ranger分析结果概述')
        section1_1.add_table(summ_all,
                             caption='样本Xenium Ranger结果汇总统计表')
        
    # ##====================================================== 4.2 基因及样本过滤 ===========================================
    # if os.path.exists("2.Count_QC"):
    #     report_hyperlink_list.append('[定量后质控](#filter)\n\n')
    #     filter_section = section1.add_section('定量后质控')
    #     filter_section.add_comment("详细结果见目录：[2.Count_QC](./2.Count_QC)")
    #     QC_folder = "2.Count_QC"
    #     summ = pd.read_table('%s/statitics_for_QC.xls' % (QC_folder), sep='\t')
    #     sample_num = len(summ)
    #     min_cell = min(summ.Total_cells_QC)
    #     max_cell = max(summ.Total_cells_QC)
    #     min_nFeature = round(min(summ.mean_nFeature_xenium_QC.astype(float)))
    #     max_nFeature = round(max(summ.mean_nFeature_xenium_QC.astype(float)))
    #     min_nCount = round(min(summ.mean_nCount_xenium_QC.astype(float)))
    #     max_nCount = round(max(summ.mean_nCount_xenium_QC.astype(float)))
    #     if len(summ) == 1:
    #         report_summary_list.append(
    #             f"本次分析共完成 {sample_num} 个样本的10x Xenium 组织原位分析测序，数据经Xenium Ranger重新运行，用10x最新算法分割细胞，最终获得{max_cell}个细胞，每个细胞中的平均UMI数为{max_nCount}，平均基因数为{max_nFeature}。")
    #     else:
    #         report_summary_list.append(
    #             f"本次分析共完成 {sample_num} 个样本的10x Xenium 组织原位分析测序，数据经Xenium Ranger重新运行，用10x最新算法分割细胞，最终获得的细胞数目分布在{min_cell}~{max_cell}个，每个细胞中的平均UMI数分布在{min_nCount}~{max_nCount}，平均基因数分布在{min_nFeature}~{max_nFeature}。")
    #     ###
    #     filter_section.add_table(summ,caption='样本质控情况统计表')
    #     filter_section.add_section('基因表达数目统计')
    #     filter_section.add_section('UMI数目统计')
    # ##====================================================== 4.3 降维聚类 =================================================
    # if os.path.exists("3.Clustering"):
    #     report_hyperlink_list.append('[降维与聚类分析结果](#rna_cluster)\n\n')
    #     if len(summ) == 1:
    #         sample_des = f"仅有1个样本，并不存在批次效应"
    #     else:
    #         sample_des = f"有不同样本，但并不存在批次效应"
    #     if config['params']['Clustering']["RNA_reduct1"] == "pca":
    #         reduct1_des = f"PCA 的降维算法进行降维 "
    #         reduct1_name = "PCA"
    #     else:
    #         sample_des = f"有不同样本，且存在批次效应"
    #         reduct1_des = f"MNN (mutual nearest neighbors, 互享最近邻)[^fastMNN] 的降维算法剔除批次效应"
    #         reduct1_name = "MNN"
    #     if config['params']['Clustering']["RNA_reduct2"] == "tsne":
    #         reduct2_des = f"t-SNE (t分布式随机邻近嵌入) 算法进行二级降维。t-SNE创建了一个缩小的特征空间，相似的样本由附近的点建模，不相似的样本由高概率的远点建模，更好的展示了样本之间的局部关系"
    #         reduct2_name = "t-SNE"
    #     else:
    #         reduct2_des = f"UMAP (统一流形逼近与投影) 算法进行二级降维。UMAP 主要基于流形理论和拓扑算法的理论，对高维数据进行降维，从而能够保留更多数据的全局结构，并且具有优越的运行性能"
    #         reduct2_name = "UMAP"
    #
    #     cluster_seurat_info = f"由于本项目{sample_des}，因此我们首先采用了{reduct1_des}，再通过 {reduct2_des}。基于{reduct1_name}降维结果通过{reduct2_name}对单细胞群聚类进行可视化，聚类算法采用SNN，最终获得最优细胞分群。"
    #     cluster_section = result_module.add_section('数据降维与聚类分析',cluster_seurat_info=cluster_seurat_info)
    #     cluster_section.add_comment("详细结果见目录：[3.Clustering](./3.Clustering)")
    #
    #     cluster_section2 = cluster_section.add_section("降维与聚类分析结果")
    #
    #     cluster_section2.add_plot('3.Clustering/*/*_groupby_cluster_resolution*_plot.png',
    #                               caption='基因表达水平降维与聚类图',
    #                               description='图片说明：左图横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。右图反映的是不同群在样本空间切片位置上的具体分布，颜色与上图对应。')
    #
    #     ## 样本分组展示
    # ##===================================================== 4.4 降维聚类可视化 =================================================
    # if os.path.exists("4.Visualize_cluster_by_clusters"):
    #     if (sample_num > 1):
    #         cluster_section3 = cluster_section.add_section('聚类结果分组展示')
    #         report_hyperlink_list.append('[聚类结果分组展示](#cluster_vis)\n\n')
    #         cluster_section3.add_section('样本分组展示图')
    #         cluster_section3.add_section('样本分面展示图')
    #         cluster_section3.add_section('样本占比统计图')
    #
    # ##===================================================== 4.5 Marker基因鉴定 ===========================================
    # if os.path.exists("5.Marker"):
    #     test_method = config["params"]['Marker']["RNA_test_method"]
    #     Marker_description = f"Marker 基因的定义是该基因在指定细胞群的绝大多数细胞中有较高的表达，而在其余细胞类群中只有少部分表达，且该基因在此细胞群相对于其他细胞群中是显著上调表达。\n 本研究采用  {test_method}  检验方法对指定细胞群与其余所有细胞群进行差异检验，先设置log2FC > 0 ,min.pct > 0.05导出所有可能的 marker ,然后按照P值，FoldChange 值以及 gene_diff 值筛选 top10 marker 作为每个细胞群的特异性 Marker 基因。"
    #     marker_section = result_module.add_section(name='Marker基因鉴定', description=Marker_description)
    #     marker_section.add_comment("详细结果见目录：[5.Marker](./5.Marker)")
    #     #####
    #     marker_result = pd.read_csv('5.Marker/top10_markers_for_each_cluster.xls', sep='\t')
    #     cluster_num = max(marker_result['cluster'])
    #     report_summary_list.append(f"通过基因表达数据进行降维聚类后共分为 {cluster_num} 群细胞，")
    #     ####
    #     marker_section.add_table("5.Marker/top10_markers_for_each_cluster.xls",
    #                               caption='每个细胞群 Top10 Marker 基因列表',
    #                               show_rows=5)
    #     ####
    #     marker_section.add_plot("5.Marker/topmarker_gene_heatmap.png",
    #                              caption='Top10 Marker 基因表达热图',
    #                              description='图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中红色表示高表达，蓝色表示低表达。')
    #     ####
    #     marker_section.add_plot("5.Marker/mark*/*Spatialfeatureplot.png",
    #                              caption='Top10 Marker 基因在降维聚类结果中的可视化图',
    #                              description='图片说明：红色越深表示该细胞中对应基因的表达量越高。')
    #     ####
    #     marker_section.add_plot("5.Marker/mark*/topmarker_gene_violin_plot.png",
    #                              caption='Top10 Marker 基因表达量小提琴图',
    #                              description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')
    #
    # ##===================================================== 4.6 细胞类型鉴定 ================================================
    # if os.path.exists("6.Reference_Celltype"):
    #     report_hyperlink_list.append('[细胞类型鉴定](#celltype)\n\n')
    #     celltype_file = list(glob('6.Reference_Celltype/*_celltyping_results.csv'))[0]
    #     cell_results = pd.read_csv(celltype_file, sep=',')
    #     celltype = cell_results['celltype'].unique()
    #     celltype_num = len(cell_results['celltype'].unique())
    #     celltype = ', '.join(celltype)
    #     report_summary_list.append(f"通过对细胞类型进行注释，得到供参考的细胞类型有{celltype_num}种，分别为{celltype}。")
    #     rna_section5 = result_module.add_section('细胞类型鉴定')
    #     rna_section5.add_comment("详细结果见目录：[6.Reference_Celltype](./6.Reference_Celltype)")
    #
    #     rna_section5.add_plot("6.Reference_Celltype/*celltyping_plot.png",
    #                           caption='细胞类型参考结果图',
    #                           description='图片说明：细胞类型结果在UMAP图上的展示，每种细胞类型以不同颜色进行区分。')
    #     rna_section5.add_plot("6.Reference_Celltype/*celltyping_spatical_plot.png",
    #                           caption='细胞类型空间分布图',
    #                           description='图片说明：不同颜色代表不同细胞类型，该图片反应不同样本的细胞类型分布情况。')
    #     rna_section5.add_table("6.Reference_Celltype/*celltyping_results.csv",
    #                            caption='细胞类型注释结果表格',
    #                            show_rows=5)

    ##==================================================== 4.项目摘要、链接==================================================
    Project_module.add_section(name="项目摘要", description=''.join(report_summary_list))
    report_hyperlink = ''.join(report_hyperlink_list)
    Project_module.add_section(name="项目快捷链接", description=report_hyperlink)
    ##============================================================= 5.附录 ============================================
    description_module = report.add_module('附录')

    refgenome_annotation = config['database_url']['Ref_genome_annotation']['Linkage']
    database = pd.DataFrame(
        [['Pannel Database', 'https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/panel-design/pre-designed-xenium-gene-expression-panels'],
         ['Reference genome annotation', refgenome_annotation]
         ])

    database.columns = ['使用数据库', '网页链接']
    description_section5 = description_module.add_section('数据库信息', description='')
    supply_section_database_table = description_section5.add_table(database, caption='数据库信息')

    database2 = pd.DataFrame(
        [['Xenium Ranger', '1.7.0'],
         ['Seurat', '5.0.0'],
         ['SeuratObject','5.0.0'],
         ['sp', '1.6.0'],
         ['ScType','1.0']
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

