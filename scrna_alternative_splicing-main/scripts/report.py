#!/usr/bin/env python37
# encoding:utf-8

from email import header, message
import click
import os
import pandas as pd
import numpy as np
import re
import shutil
import subprocess
import sys
import time
import yaml
from glob import glob
import shutil

# from oebio.app import *
from oebio.report import Report, oeweb_register

@click.command()
@click.option('-i', '--inputdir', prompt='the program path', help='the directory of your program. ')
@click.option('-c', '--configfile', prompt='config file', default='config/config.yaml', help='the config file which contain your program information.')
@click.option('-d', '--dirnum', default=0)



def multimodal_report(inputdir, configfile, dirnum):

    #######################  1.读取项目配置文件信息   #######################

    """ python $0 -i result -c config/config.yaml"""
    cfg = open(configfile, 'r', encoding='utf-8').read()
    config = yaml.full_load(cfg)
    program_path = os.path.abspath(inputdir)

    # samples
    metadatafile = config['metadatafile']
    if os.path.exists(metadatafile):
        samples = pd.read_csv(metadatafile, dtype=str).set_index("sampleid",drop=False)
        samples.dropna(how='all', inplace=True)
        samples = samples[ ~np.array([s.startswith("#") for s in samples.sampleid.to_list()])]
        samples = samples.set_index("sampleid",drop=False)
        samples.index.names = ["sample_id"]
    else:
        print(f"您尚未准备好样本文件{metadatafile}")
        os._exit(0)

    # diff groups
    diff_DEG_group = config['diff_DEG_group']
    if os.path.exists(diff_DEG_group):
        diff_groups = pd.read_csv(diff_DEG_group, dtype=str, delimiter=",")
        diff_groups.dropna(how='all', inplace=True)
        diff_groups = diff_groups[ ~np.array([s.startswith("#") for s in diff_groups.treatment.to_list()]) ]
        diff_groups.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        diff_groups.index = diff_groups["type"] + "_" + diff_groups["treatment_name"].replace(" ", "_", regex=True) + "-vs-" + diff_groups["control_name"].replace(" ", "_", regex=True)
        diff_groups['groupname'] = diff_groups.index
    else:
        print(f"您尚未准备好差异分组设置文件 {diff_DEG_group}")
        os._exit(0)

    # project
    if os.path.exists('result/cellranger'):
        project='10X'
    if os.path.exists('result/1.STARsolo'):
        project='M20'

    # 从config文件读取差异筛选参数
    program_num = dict(config['report'])['Project_Num']
    outdir = f"report/{program_num}_Report"
    scriptdir = os.path.dirname(os.path.realpath(__file__))

    log2fc_gene = config["params"]["splicing"]["plot"]["log2fc_gene"]
    min_gene_expression = config["params"]["splicing"]["analysis"]["min_gene_expression"]
    min_gene_expression_name = "min.log2EP-" + str(min_gene_expression)
    dimension_type = config["params"]["splicing"]["prepare"]["dimension_type"]
    top_DE_junctions = config["params"]["splicing"]["plot"]["top_DE_junctions"]

    filter_hard = config["params"]["splicing"]["analysis"]["filter_hard"]
    if filter_hard:
        filter_gene_in_percent_cell = round(config["params"]["splicing"]["analysis"]["filter_gene_in_percent_cell"] * 100,0)
        filter_sj_in_percent_cell = round(config["params"]["splicing"]["analysis"]["filter_sj_in_percent_cell"] * 100,0)
    # else:


    pvalue_sj = config["params"]["splicing"]["plot"]["pvalue_sj"]
    # 基因差异显著性可以采用 pvalue，也可以采用 qvalue
    if "pvalue_gene" in config["params"]["splicing"]["plot"].keys():
        pval_or_padj = "pvalue.gene-" + str(config["params"]["splicing"]["plot"]["pvalue_gene"])
        pval_or_padj_stand = "p-value"
        pval_or_padj_value = config["params"]["splicing"]["plot"]["pvalue_gene"]
    if "qvalue_gene" in config["params"]["splicing"]["plot"].keys():
        pval_or_padj = "qvalue.gene-" + str(config["params"]["splicing"]["plot"]["qvalue_gene"])
        pval_or_padj_stand = "q-value"
        pval_or_padj_value = config["params"]["splicing"]["plot"]["qvalue_gene"]
    if "pvalue_gene" in config["params"]["splicing"]["plot"].keys() & "qvalue_gene" in config["params"]["splicing"]["plot"].keys():
        print("Warning: 参数 pvalue_gene 和 qvalue_gene 二选一")
        os._exit(0)
    # junction差异总量可以使用 log2fc，也可以使用PSI
    if "delta_sj" in config["params"]["splicing"]["plot"].keys():
        delat_or_log2fc_stand = "PSI"
        delat_or_log2fc = "delta.sj-" + str(config["params"]["splicing"]["plot"]["delta_sj"])
        delat_or_log2fc_value = config["params"]["splicing"]["plot"]["delta_sj"]
    if "log2fc_sj" in config["params"]["splicing"]["plot"].keys():
        delat_or_log2fc = "log2fc.sj-" + str(config["params"]["splicing"]["plot"]["log2fc_sj"])
        delat_or_log2fc_stand = "log2fc(PSI)"
        delat_or_log2fc_value = config["params"]["splicing"]["plot"]["log2fc_sj"]
    
    diff_long_name = "pvalue.sj-" + str(pvalue_sj) + "-" + delat_or_log2fc + "-" + pval_or_padj + "-" + "log2fc.gene-" + str(log2fc_gene) + "-" + min_gene_expression_name
    volcan_long_name = "pvalue.sj-" + str(pvalue_sj) + "-" + delat_or_log2fc + "-" + min_gene_expression_name
    
    number_diff_group = len(diff_groups.index)
    content_diff_group = ', '.join(diff_groups.index)


    # Report
    # if os.path.exists(f'')
    if os.path.exists(outdir):
        print(f'删除已经存在的 {outdir} 文件夹，并重新创建该文件夹')
        shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    if not os.path.exists("%s/src" % (outdir)):
        os.makedirs("%s/src" % (outdir))
        subprocess.call(f'cp {scriptdir}/../Documents/MARVEL-droplet-DSJ-pipeline.png {outdir}/src/', shell=True)
        subprocess.call(f'cp {scriptdir}/../Documents/marvel-pipeline-{project}.png {outdir}/src/', shell=True)
        subprocess.call(f'cp {scriptdir}/../Documents/差异分析表格说明.csv {outdir}/src/', shell=True)
        subprocess.call(f"cp {scriptdir}/../Documents/结构图表格说明.csv {outdir}/src/", shell=True)

    ######################  2.生成结果文件   #######################
 
    # 报告默认不把junction矩阵发给老师
    # # 1.STARsolo
    # if os.path.exists("%s/Splicing/STARSolo" % (program_path)):
    #     print(f"正在拷贝STARsolo比对结果到 {outdir}")
    #     for sample in samples.index:
    #         copypath = outdir + "/1.STARsolo/" + sample + "/SJ/"
    #         os.makedirs(copypath)
    #         rawfiles = glob("result/Splicing/STARSolo/" + sample + "/" + sample + "_Solo.out/SJ/raw/*.gz")
    #         for file in rawfiles:
    #             copyfile = copypath + "/" + file.split("/")[-1]
    #             shutil.copy2(file, copyfile)
    # else:
    #     print("Can not find result/Splicing/STARSolo results!!!!!!!!!!!!!!!!")

    # Splicing
    print(f"正在拷贝剪切结果到 {outdir}")
    for diff_group in diff_groups.index:
        rawpath = "result/Splicing/" + diff_group
        copypath = outdir + "/Splicing" 
        if not os.path.exists(copypath):
            os.makedirs(copypath)
        if os.path.exists(rawpath):
            subprocess.call(f'cp -r {rawpath} {copypath}', shell=True)
            if len(glob(f"{copypath}/{diff_group}/scatter/*/scatter.done")) > 0:
                subprocess.call(f'rm {copypath}/{diff_group}/scatter/*/scatter.done', shell=True)
            if len(glob(f"{copypath}/{diff_group}/MARVEL.RData")) > 0:
                subprocess.call(f"rm {copypath}/{diff_group}/MARVEL.RData", shell=True)
        else:
            print(f"Can not find {rawpath} in your path")

    #######################  3.生成网页版报告  #######################

    print("开始生成网页版报告")
    header_info = {}
    # header_info["原项目编号"] = config['report']['Raw_Task_Num']
    header_info["原项目编号"] = ','.join(set(samples['raw_task_number'].to_list()))
    header_info["项目编号"] = config['report']["Project_Num"]
    header_info["客户姓名"] = config['report']["Customer"]
    header_info["实验物种"] = config['report']["Species"].capitalize()
    header_info["实验样本"] = config['report']["Sample"]
    header_info["执行编号"] = config['report']["Executor"]
    header_info["样本数目"] = len(samples.index)
    header_info["任务单号"] = config['report']['Task_Num']
    header_info["比较组数量"] = len(diff_groups.index)

    report = Report('单细胞可变剪切分析结题报告', title='单细胞可变剪切分析结题报告',
                    header_info=dict(header_info))
    report.add_yaml_config(os.path.abspath(f'{scriptdir}/report.yaml'))
    os.environ['oeweb_register_token'] = 'oearray'
    Project_module = report.add_module('项目概况')

    samplesize = len(samples.index)
    sample_names = ", ".join(samples.index)
    diffgroup_size = len(diff_groups.index)
    diffgroup_names = ", ".join(diff_groups.index)
    summary = Project_module.add_section(f'项目摘要', samplesize=samplesize,
                               sample_names=sample_names,
                               diffgroup_size=diffgroup_size,
                               diffgroup_names=diffgroup_names,                        
                               )
    for diff_group in diff_groups.index:
        candidate_data = pd.read_csv(f'{outdir}/Splicing/{diff_group}/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.csv')
        total_data = pd.read_csv(f'{outdir}/Splicing/{diff_group}/diff_SJ_exp/splicing_junction_differect_expression.total.csv')
        total_gene_number = len(set(total_data.gene_name.to_list()))
        number_candidate_data = len(set(candidate_data.gene_name.to_list()))
        summary.add_comment(f'{diff_group} 共有 {total_gene_number} 个基因参与分析，其中共检测到 {number_candidate_data} 个基因含有差异表达的junction。')
    # register_info = oeweb_register(project_id=program_num, target_url='https://cloud.oebiotech.cn/task/category/scrna',
    #                                note='本报告包含项目基本分析内容，如需快速做细胞表达可视化等分析(数据默认保留半年)')
    # if register_info:
    #     Project_module.add_comment(register_info)
    report_summary_list = []
    report_hyperlink_list = []

    #############################
    # 报告具体内容
    #############################

    ## ============ 项目信息

    os.chdir(outdir)

    ## ============ 基本原理
    iterations = config["params"]["splicing"]["analysis"]["iterations"]
    report.add_section('MARVEL检测可变剪切基本原理')

    ## ============ 生信流程
    dimension_type = config["params"]["splicing"]["prepare"]["dimension_type"]
    section_pipeline = report.add_section('生信分析' )
    analysis_pipeline = section_pipeline.add_section('分析流程', dimension_type=dimension_type, iterations=iterations)
    analysis_pipeline.add_plot(f'src/marvel-pipeline-{project}.png', caption='生信分析流程图')
    
    #  ============ 基因筛选标准
    section_pipeline.add_section('基因筛选的标准', dimension_type=dimension_type, 
                                                 filter_gene_in_percent_cell=filter_gene_in_percent_cell,
                                                 filter_sj_in_percent_cell=filter_sj_in_percent_cell,
                                                 pvalue_sj=pvalue_sj,
                                                 delat_or_log2fc_stand=delat_or_log2fc_stand,
                                                 delat_or_log2fc_value=delat_or_log2fc_value,
                                                 pval_or_padj_stand=pval_or_padj_stand,
                                                 pval_or_padj_value=pval_or_padj_value,
                                                 log2fc_gene=log2fc_gene )
    
    ##  ============ 基因分类标准
    section_pipeline.add_section('基因分类标准')

    ## ============ 分析结果
    result_module = report.add_module('项目分析结果')

    ### ============ 添加 1.STARsolo
    result_module.add_section('STARsolo分析')

    if os.path.exists('Splicing'):
        
        ### 可变剪切分析
        if len(diff_groups.index) > 1:
            resp="分别"
        else:
            resp=""
        section_diff = result_module.add_section('可变剪切分析', number_diff_group=number_diff_group, content_diff_group=content_diff_group, resp=resp)

        ## 添加差异分析总表
        section_total_table = section_diff.add_section("junction差异表达分析总表")  
        for diff_group in diff_groups.index:
            section_total_table.add_comment(f'详细结果见目录：[Splicing/{diff_group}/diff_SJ_exp](./Splicing/{diff_group}/diff_SJ_exp)')
        table_1 = diff_groups.index[0]
        table_diff_total = f"Splicing/{table_1}/diff_SJ_exp/splicing_junction_differect_expression.total.csv"
        section_total_table.add_table(table_diff_total, caption="差异可变剪切分析总表", show_rows=11)
        section_total_table.add_table('src/差异分析表格说明.csv', caption="junction差异表达分析总表表格说明")

         ## 添加火山图
        section_volcan_png = section_diff.add_section('junction差异表达火山图')
        for diff_group in diff_groups.index:
            section_volcan_png.add_comment(f"详细结果见目录：[Splicing/{diff_group}/volcano](./Splicing/{diff_group}/volcano)")
        section_volcan_png.add_comment('junction表达火山图见下图：')
        section_volcan_png.add_plot(f"Splicing/*/volcano/diff_splicing_volcano.{volcan_long_name}.png", 
                                    caption='横坐标为junction在两种细胞类型间的平均表达量，纵坐标为junction差异水平')

        ## 添加筛选后候选基因表格
        section_candidate_table = section_diff.add_section('候选可变剪切基因列表')
        for diff_group in diff_groups.index:
            section_candidate_table.add_comment(f'详细结果见目录：[Splicing/{diff_group}/diff_SJ_exp](./Splicing/{diff_group}/diff_SJ_exp)')
        table_candidate = f"Splicing/{table_1}/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.csv"
        section_candidate_table.add_table(table_candidate, caption="候选可变剪切基因列表", show_rows=11)

        ## 添加候选基因分类统计图
        section_catalog_png = section_diff.add_section('候选基因分类统计')
        for diff_group in diff_groups.index:
            section_catalog_png.add_comment(f"详细结果见目录：[Splicing/{diff_group}/diff_SJ_exp](./Splicing/{diff_group}/diff_SJ_exp)")
        section_catalog_png.add_comment('候选基因分类统计见下图：')
        section_catalog_png.add_plot(f"Splicing/*/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.png", caption="四种分类基因占比饼图")

        # 从散点图开始，图片以候选基因的类型再分别展示
        ## 添加散点图
        
        section_scatter_png = section_diff.add_section('候选基因PSI散点图', top_DE_junctions=top_DE_junctions)
        for diff_group in diff_groups.index:
            catalogs_gene = glob(f'Splicing/{diff_group}/scatter/{diff_long_name}/*')
            catalogs_gene = [catalog.split("/")[-1] for catalog in catalogs_gene]
            for catalog in catalogs_gene:
                section_scatter_png.add_comment('<br>')
                section_scatter_png.add_comment(f'{diff_group} 类型 {catalog} 候选基因PSI散点图展示如下：')
                section_scatter_png.add_comment(f'详细结果见目录：[Splicing/{diff_group}/scatter/{diff_long_name}/{catalog}](./Splicing/{diff_group}/scatter/{diff_long_name}/{catalog})')
                section_scatter_png.add_plot(f"Splicing/{diff_group}/scatter/{diff_long_name}/{catalog}/*.png", 
                                             caption='（A）细胞分类示意图；（B）基因在细胞类型间的表达量对比；（C）差异表达 junction 在细胞类型间的 PSI 值')


        ## 添加表达量面板图
        section_tabulate_png = section_diff.add_section('候选基因表达量和junction表达量在细胞群水平对比图', top_DE_junctions=top_DE_junctions)
        for diff_group in diff_groups.index:
            catalogs_gene = glob(f'Splicing/{diff_group}/scatter/{diff_long_name}/*')
            catalogs_gene = [catalog.split("/")[-1] for catalog in catalogs_gene]
            for catalog in catalogs_gene:
                section_tabulate_png.add_comment('<br>')
                section_tabulate_png.add_comment(f'{diff_group} 类型 {catalog} 候选基因PSI散点图展示如下：')
                section_tabulate_png.add_comment(f'详细结果见目录：[Splicing/{diff_group}/tabulate/{diff_long_name}/{catalog}](./Splicing/{diff_group}/tabulate/{diff_long_name}/{catalog})')
                section_tabulate_png.add_plot(f"Splicing/{diff_group}/tabulate/{diff_long_name}/{catalog}/*.png", 
                                             caption='（A）细胞分类示意图；（B）基因在细胞类型间的表达量对比；（C）差异表达 junction 在细胞类型间的 PSI 值')
        section_tabulate_png.add_comment('<br>')
        section_tabulate_png.add_comment("以细胞群为单位统计基因和junction表达量的详细信息总结在下表中：")
        for diff_group in diff_groups.index:
            section_tabulate_png.add_comment(f"详细结果见目录：[Splicing/{diff_group}/tabulate/{diff_long_name}](./Splicing/{diff_group}/tabulate/{diff_long_name})")
        section_tabulate_png.add_table(f'Splicing/{table_1}/tabulate/{diff_long_name}/gene_SJ_expression.csv', caption="基因和junction在细胞群水平的表达量统计", show_rows=11)


        ## 添加候选基因结构图
        section_structure_png = section_diff.add_section('候选基因转录本和差异junction位置结构图', top_DE_junctions=top_DE_junctions)
        for diff_group in diff_groups.index:
            catalogs_gene = glob(f'Splicing/{diff_group}/scatter/{diff_long_name}/*')
            catalogs_gene = [catalog.split("/")[-1] for catalog in catalogs_gene]
            for catalog in catalogs_gene:
                section_structure_png.add_comment('<br>')
                section_structure_png.add_comment(f'{diff_group} 类型 {catalog} 候选基因PSI散点图展示如下：')
                section_structure_png.add_comment(f'详细结果见目录：[Splicing/{diff_group}/structure/{diff_long_name}/{catalog}](./Splicing/{diff_group}/structure/{diff_long_name}/{catalog})')
                section_structure_png.add_plot(f"Splicing/{diff_group}/structure/{diff_long_name}/{catalog}/*.png", 
                                             caption='黑色和灰色方框表示外显子区域，其中灰色表示CDS区')

        section_structure_png.add_comment("<br>")
        section_structure_png.add_comment("转录本结构注释文件：")
        section_structure_png.add_comment("**注意**：正链转录本序号从左往右数，负链转录本序号从右往左数")
        section_structure_png.add_comment(f"详细结果见目录：[Splicing/{diff_group}/structure/{diff_long_name}](./Splicing/{diff_group}/structure/{diff_long_name})")
        section_structure_png.add_table(f'Splicing/{diff_group}/structure/{diff_long_name}/structure_transtripts.csv', caption="转录本结构注释", show_rows=11 )
        section_structure_png.add_comment("<br>")
        section_structure_png.add_comment("junction结构注释文件：")
        section_structure_png.add_table(f'Splicing/{diff_group}/structure/{diff_long_name}/structure_junctions.csv', caption="junction结构注释", show_rows=11)

        section_structure_png.add_comment("<br>")
        section_structure_png.add_comment('转了本和junction的结构表格说明：')
        section_structure_png.add_table(f'src/结构图表格说明.csv', caption='结构图表格说明')

        ## 添加富集分析结果
        # 由于筛选到的基因假阳性较多，暂时没有添加注释结果到报告中

        ## 添加附录

        ## 添加申明
        affirming_module = report.add_module('申明')

    report.write_to('Report.html')


if __name__ == "__main__":
    multimodal_report()
