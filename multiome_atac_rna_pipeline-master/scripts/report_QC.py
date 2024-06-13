#!/usr/bin/env python37
# encoding:utf-8
import os,time,sys,subprocess,click,shutil,yaml,re
from oebio.report import Report,oeweb_register
from oebio.app import *
from glob import glob
import pandas as pd


@click.command()
@click.option('-i','--inputdir', prompt='the program path',
                         help='the directory of your program. ')
@click.option('-c','--configfile', prompt='config file', default='config/config.yaml',
                         help='the config file which contain your program information.')


def multimodal_report(inputdir,configfile):
    ##=============================================== 1. 读取项目配置文件信息 ===========================================
    """ python multimodal_report.py  -i ./ -c ./config.yaml"""
    cfg = open(configfile,'r',encoding='utf-8').read()
    config = yaml.full_load(cfg)
    program_path = os.path.abspath(inputdir)
    #paragram from config
    program_num = dict(config['report'])['Project_Num']
    report_time = time.strftime("%Y_%m_%d")  # report_time=config['report']['Project_End_Time']
    outdir = f"report/{program_num}_Report_{report_time}"
    log_file = open("%s/report_log.txt" % (program_path),"w")
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    # Report
    os.chdir(program_path)
    if os.path.exists(outdir):
        shutil.rmtree(outdir,ignore_errors=True)
    os.makedirs(outdir)
    ##================================================2. 生成结果文件 =================================================

    # 1.CellRanger
    if os.path.exists("%s/cellranger" % (program_path)):
        os.makedirs("%s/1.Cellranger" % (outdir))
        name = [str(i).split('/')[-2] for i in glob("cellranger/*/outs")]
        subprocess.call('cp %s/cellranger/summary.csv %s/1.Cellranger/' % (program_path, outdir), shell=True)

        for j in name:
            os.makedirs("%s/1.Cellranger/%s" % (outdir, j))
            subprocess.call('cp %s/cellranger/%s/*.png %s/1.Cellranger/' % (program_path,j,outdir),shell=True)
            subprocess.call('ln -s %s/cellranger/%s/outs/* %s/1.Cellranger/%s ' % (program_path, j, outdir, j),shell=True)

    else:
        print("Can not find 1.Cellranger results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 1.Cellranger results!!!!!!!!!!!!!!!!" + "\n")

    # report
    if not os.path.exists("%s/src" % (outdir)):
        os.makedirs("%s/src" % (outdir))
        subprocess.call(f'cp {scriptdir}/src/report/* {outdir}/src/',shell=True)

    ##===================================================== 3. 生成报告================================================
    header_info = {key:config["report"][key] for key in config["report"] if
                   key in ["Project_Num","Customer","Species","Executor","Sample_Num","Task_Num"]}
    header_info["项目编号"] = header_info.pop("Project_Num")
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["实验物种"] = header_info.pop("Species")
    header_info["执行编号"] = header_info.pop("Executor")
    header_info["样本数目"] = header_info.pop("Sample_Num")
    header_info["任务单号"] = header_info.pop("Task_Num")

    report = Report('10x 单细胞多组学（ATAC + 基因表达）质控报告', title = '10x 单细胞多组学（ATAC + 基因表达）质控报告',header_info = dict(header_info))
    report.add_yaml_config(os.path.abspath(f'{scriptdir}/src/report.ymal'))
    os.environ['oeweb_register_token'] = 'oearray'
    Project_module = report.add_module('项目概况')
    register_info = oeweb_register(project_id=program_num, target_url='https://cloud.oebiotech.cn/task/category/scrna',
                                   note='本报告包含项目基本分析内容，如需快速做细胞表达可视化等分析(数据默认保留半年)')
    if register_info:
        Project_module.add_comment(register_info)
    report_summary_list = []

    ## ======================================================= 3.1 背景信息 ============================================
    Background_module = report.add_module('10x 单细胞多组学（ATAC + 基因表达）简介')
    ## 实验流程
    Background_section1 = Background_module.add_section('实验技术流程')
    Background_section1.add_plot( "%s/src/pipeline1.png"%(outdir),
                                  caption='10x 单细胞多组学（ATAC + 基因表达）文库结构示意图',
                                  content='10x 单细胞多组学（ATAC + 基因表达）文库结构示意图')
    Background_section1.add_plot( "%s/src/pipeline2.png"%(outdir),
                                  caption='10x 单细胞多组学（ATAC + 基因表达）文库构建流程示意图',
                                  content="10x 单细胞多组学（ATAC + 基因表达）文库构建流程示意图")

    ## 分析流程
    workflow_module = Background_module.add_module('10x 单细胞多组学（ATAC + 基因表达）生物信息质控分析流程')
    workflow_module.add_plot( "%s/src/qc_bio.png"%(outdir),
                              caption='10x 单细胞多组学（ATAC + 基因表达）生物信息质控分析流程图',
                              content="10x 单细胞多组学（ATAC + 基因表达）生物信息质控分析流程图")

    ##====================================================== 3.2 质控分析 =============================================
    result_module = report.add_module('质控分析结果')

    #################### Cell Ranger 标准分析 ####################
    qc_section2 = result_module.add_section('Cellranger ARC 标准质控')
    cellranger_out_plot = qc_section2.add_plot("%s/1.Cellranger/*.png"%(outdir),
                                                           caption='样本 Cell Ranger ARC 质控结果概述',
                                                           content='样本 Cell Ranger ARC 质控结果概述')
    cellranger_summary = qc_section2.add_table("%s/1.Cellranger/summary.csv"%(outdir),caption="样本细胞质量统计结果")
    #==============================
    summ = pd.read_csv('%s/1.Cellranger/summary.csv' % (outdir), sep=',')
    sample_num = len(summ)
    max_number_cells = max(summ["Estimated number of cells"])
    min_number_cells = min(summ["Estimated number of cells"])
    max_media_gene = max(summ["GEX Median genes per cell"])
    min_media_gene = min(summ["GEX Median genes per cell"])
    max_media_peak = max(summ["ATAC Median high-quality fragments per cell"])
    min_media_peak = min(summ["ATAC Median high-quality fragments per cell"])
    max_tss = max(summ["ATAC TSS enrichment score"])
    min_tss = min(summ["ATAC TSS enrichment score"])
    if len(summ) == 1:
        report_summary_list.append(f"本次分析共完成 {sample_num} 个样本的10x 单细胞多组学（ATAC + 基因表达）测序，Cell Ranger ARC 定量质控的高质量细胞数为 {min_number_cells} 个,获得的基因中位数为{min_media_gene}，获得的高质量peaks的中位数为{min_media_peak}，TSS的富集分数为{min_tss},其他指标详见本报告。")
    else:
        report_summary_list.append(f"本次分析共完成 {sample_num} 个样本的10x 单细胞多组学（ATAC + 基因表达）测序，Cell Ranger ARC 定量质控的高质量细胞数分布在 {min_number_cells}~{max_number_cells} 个,获得基因的中位数分布在{min_media_gene}~{max_media_gene}，获得的高质量peaks的中位数分布在{min_media_peak}~{max_media_peak}，TSS的富集分数分布在{min_tss}~{max_tss},其他指标详见本报告。")

    cellranger_out_table = qc_section2.add_table("%s/src/multimodal_qc_stat.txt" % (outdir), caption='样本细胞质量统计结果说明')


    ##========================================================= 4.项目摘要、链接========================================
    Project_module.add_section(name="项目摘要",description=''.join(report_summary_list))

    ##============================================================= 5.附录 ============================================
    description_module = report.add_module('附录')
    description_section1 = description_module.add_section('Loupe Browser 使用教程')
    description_section3 = description_module.add_section('10x 单细胞多组学（ATAC + 基因表达）生信分析方法说明')
    description_section4 = description_module.add_section('10x 单细胞多组学（ATAC + 基因表达）常见问题（FAQ）')
    
    refgenome = config['database_url']['Ref_genome']['Linkage']
    refgenome_annotation = config['database_url']['Ref_genome_annotation']['Linkage']
    database = pd.DataFrame(
        [['GO Database', 'http://geneontology.org/'],
         ['KEGG Database','http://www.genome.jp/kegg/'],
         ['Reference genome',refgenome],
         ['Reference genome annotation',refgenome_annotation]
         ])    


    database.columns = ['使用数据库', '网页链接']
    description_section5 = description_module.add_section('数据库信息', description='')
    supply_section_database_table = description_section5.add_table(database, caption='数据库信息')

    database2 = pd.DataFrame(
        [['Cell Ranger ARC', '2.0.0']])
    database2.columns = ['软件', '版本']
    description_section6 = description_module.add_section('数据分析软件')
    software_table = description_section6.add_table(database2, caption='数据分析软件')

    ##============================================================== 6.申明 ===========================================
    affirming_module = report.add_module('申明')

    ############################################ Generate Report HTML ###################################################
    report.write_to('%s/Report.html'%(outdir))
if __name__ == "__main__":
    multimodal_report()
