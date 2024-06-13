#!/usr/bin/env python37
# encoding:utf-8
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
def report(inputdir, configfile):
    ##=============================================== 1. 读取项目配置文件信息 ===========================================
    """ python multimodal_report.py  -i ./ -c ./config.yaml"""
    cfg = open(configfile, 'r', encoding='utf-8').read()
    config = yaml.full_load(cfg)
    program_path = os.path.abspath(inputdir)
    # paragram from config
    program_num = dict(config['report'])['Project_Num']
    report_time = time.strftime("%Y_%m_%d")
    outdir = f"report/{program_num}_QC_Report_{report_time}"
    log_file = open("%s/report_single_species_log.txt" % (program_path), "w")
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
        subprocess.call('cp -r %s/1.STARsolo/Read_Distribution     %s/1.STARsolo/' % (program_path, outdir),
                        shell=True)
        subprocess.call('cp -r %s/1.STARsolo/STARsolo.statistic.tsv  %s/1.STARsolo/' % (program_path, outdir),
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

    report = Report('M20 单细胞微生物转录组质控报告', title='M20 单细胞微生物转录组质控报告',
                    header_info=dict(header_info))
    report.add_yaml_config(os.path.abspath(f'{scriptdir}/report.yaml'))
    report_summary_list = []
    report_hyperlink_list = []
    os.chdir(outdir)
    ## ======================================================= 1. 实验及分析概述============================================
    Background_module = report.add_module('技术简介')
    Background_section1 = Background_module.add_section('实验流程')  ## 实验流程
    #workflow_module = Background_module.add_section('生信分析流程')  ## 分析流程
    # ====================================================== 3. 分析结果 ================================================
    result_module = report.add_module('质控分析结果')
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
                f"本次分析共完成 {sample_num} 个样本的M20 单细胞微生物转录组测序，数据经STAR质控过滤后的高"
                f"质量细胞数为 {min_beforeQC_cell} 个，经剔除异常基因表达以及低基因表达的细胞后，最终获得的"
                f"细胞数目为 {min_afterQC_cell} 个，每个细胞中的平均 UMI 数为 {min_nCount_RNA_afterQC}，"
                f"每个细胞中的平均基因数为 {min_nFeature_afterQC}")
        else:
            report_summary_list.append(
                f"本次分析共完成 {sample_num} 个样本的M20 单细胞微生物转录组测序，各样本经STAR质控过滤后的高质量细胞数分布在 "
                f"{min_beforeQC_cell}~{max_beforeQC_cell} 个，经剔除异常基因表达以及低基因表达的细胞后，最终获得的细胞数目分布在 "
                f"{min_afterQC_cell}~{max_afterQC_cell} 个，每个细胞中的平均 UMI 数分布在 {min_nCount_RNA_afterQC}~{max_nCount_RNA_afterQC}，"
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
        report_summary_list.append(f"。")
        ##
        QC_description = f"在 STAR 初步质控的基础上进一步利用Seurat[^seurat]包对实验数据进行质控。理论上大部分细胞中表达的基因数量、" \
                         f"UMI 数量、rRNA基因以及tRNA基因数量等会集中分布在某一区域内，根据它们的分布特征可以拟合分布模型，使用该模型找到" \
                         f"其中的离域值，剔除异常数据。本项目中的质控标准为：1.去除在所有细胞中均未表达的基因，2.过滤细胞中少于10个基因" \
                         f"表达{rRNA_descript}{tRNA_descript}的细胞，以同时满足以上标准的细胞作为高质量细胞，进行下游分析。" \
                         f"\n [^seurat]: Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data[J]. Cell, 2021, 184(13):3573-3587.\n" \
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
    ##============================================================= 5.附录 ============================================
    description_module = report.add_module('附录')
    ####get database url information
    refgenome = config['database_url']['Ref_genome']['Linkage']
    refgenome_annotation = config['database_url']['Ref_genome_annotation']['Linkage']
    database = pd.DataFrame(
        [
         ['Reference genome', refgenome],
         ['Reference genome annotation', refgenome_annotation]
         ])
    database.columns = ['使用数据库', '网页链接']
    description_section5 = description_module.add_section('数据库信息', description='')
    supply_section_database_table = description_section5.add_table(database, caption='数据库信息')

    ##============================================================== 6.申明 ===========================================
    affirming_module = report.add_module('申明')

    ############################################ Generate Report HTML ###################################################
    report.write_to('Report.html')


if __name__ == "__main__":
    report()
