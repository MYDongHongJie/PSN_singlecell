#!/usr/bin/env python37
# encoding:utf-8

# Report
# 1.Cellranger
# 2.Count_QC
# 3.RNA_Clustering
# 4.RNA_Marker
# 5.RNA_Reference_celltype
# 6.ATAC_Clustering
# 7.ATAC_Marker
# 8.ATAC_Reference_celltype
# 9.WNN_Clustering
# 10.WNN_RNA_Marker
# 11.WNN_ATAC_Marker
# 12.WNN_Reference_celltype
# 13.Diffexp_RNA
# 14.Diffexp_ATAC
# 15.enrichment

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

        for j in name:
            os.makedirs("%s/1.Cellranger/%s" % (outdir, j))
            subprocess.call('cp %s/cellranger/%s/*.png %s/1.Cellranger/' % (program_path,j,outdir),shell=True)
            subprocess.call('ln -s %s/cellranger/%s/outs/* %s/1.Cellranger/%s ' % (program_path, j, outdir, j),shell=True)

    else:
        print("Can not find 1.Cellranger results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 1.Cellranger results!!!!!!!!!!!!!!!!" + "\n")

    # 2.Count_QC
    if os.path.exists("%s/2.Count_QC" % (program_path)):
        os.makedirs("%s/2.Count_QC" % (outdir))
        subprocess.call('cp %s/2.Count_QC/*.png %s/2.Count_QC' % (program_path, outdir),shell=True)
        subprocess.call('cp %s/2.Count_QC/*.pdf %s/2.Count_QC' % (program_path, outdir),shell=True)
        subprocess.call('cp %s/2.Count_QC/*.xls %s/2.Count_QC' % (program_path,outdir),shell=True)
    else:
        print("Can not find 2.Count_QC results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 2.Count_QC results!!!!!!!!!!!!!!!!!" + "\n")

    # 3.RNA_Cluster
    if os.path.exists("%s/3.RNA_Clustering" % (program_path)):
        os.makedirs("%s/3.RNA_Clustering" % (outdir))
        subprocess.call('cp -r %s/3.RNA_Clustering/* %s/3.RNA_Clustering' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/3.RNA_Clustering/*.log;rm %s/3.RNA_Clustering/*/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 3.RNA_Clustering results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 3.RNA_Clustering results!!!!!!!!!!!!!!!!!" + "\n")

    # 4.RNA_Marker
    if os.path.exists("%s/4.RNA_Marker" % (program_path)):
        os.makedirs("%s/4.RNA_Marker" % (outdir))
        subprocess.call('cp -r %s/4.RNA_Marker/* %s/4.RNA_Marker' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/4.RNA_Marker/*.log %s/4.RNA_Marker/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 4.RNA_Marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 4.RNA_Marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 5.RNA_Reference_celltype
    if os.path.exists("%s/5.RNA_Reference_celltype" % (program_path)):
        os.makedirs("%s/5.RNA_Reference_celltype" % (outdir))
        subprocess.call('cp %s/5.RNA_Reference_celltype/* %s/5.RNA_Reference_celltype' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/5.RNA_Reference_celltype/celltyping*' % (outdir),shell=True)
    else:
        print("Can not find 5.RNA_Reference_celltype results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 5.RNA_Reference_celltype results!!!!!!!!!!!!!!!!!" + "\n")

    # 6.ATAC_Clustering
    if os.path.exists("%s/6.ATAC_Clustering" % (program_path)):
        os.makedirs("%s/6.ATAC_Clustering" % (outdir))
        subprocess.call('cp  -r %s/6.ATAC_Clustering/* %s/6.ATAC_Clustering' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/6.ATAC_Clustering/*.log;rm %s/6.ATAC_Clustering/*/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 6.ATAC_Clustering results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 6.ATAC_Clustering results!!!!!!!!!!!!!!!!!" + "\n")

    # 7.ATAC_Marker
    if os.path.exists("%s/7.ATAC_Marker" % (program_path)):
        os.makedirs("%s/7.ATAC_Marker" % (outdir))
        subprocess.call('cp -r %s/7.ATAC_Marker/* %s/7.ATAC_Marker' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/7.ATAC_Marker/*.log %s/7.ATAC_Marker/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 7.ATAC_Marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 7.ATAC_Marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 8.ATAC_Reference_celltype
    if os.path.exists("%s/8.ATAC_Reference_celltype" % (program_path)):
        os.makedirs("%s/8.ATAC_Reference_celltype" % (outdir))
        subprocess.call('cp %s/8.ATAC_Reference_celltype/visualize_cluster_by_celltype/groupby-celltype_resolution0.4_contrast_plot.pdf* %s/8.ATAC_Reference_celltype/humanref_hpca_top.main_celltyping_plot.pdf' % (program_path, outdir),shell=True)
        subprocess.call('cp %s/8.ATAC_Reference_celltype/visualize_cluster_by_celltype/groupby-celltype_resolution0.4_contrast_plot.png* %s/8.ATAC_Reference_celltype/humanref_hpca_top.main_celltyping_plot.png' % (program_path,outdir),shell=True)
    else:
        print("Can not find 8.ATAC_Reference_celltype results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 8.ATAC_Reference_celltype results!!!!!!!!!!!!!!!!!" + "\n")
    
    # 9.WNN_Clustering
    if os.path.exists("%s/9.WNN_Clustering" % (program_path)):
        os.makedirs("%s/9.WNN_Clustering" % (outdir))
        subprocess.call('cp -r %s/9.WNN_Clustering/* %s/9.WNN_Clustering' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/9.WNN_Clustering/*.log;rm %s/9.WNN_Clustering/*/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 9.WNN_Clustering results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 9.WNN_Clustering results!!!!!!!!!!!!!!!!!" + "\n")

    # 10.WNN_RNA_Marker
    if os.path.exists("%s/10.WNN_RNA_Marker" % (program_path)):
        os.makedirs("%s/10.WNN_RNA_Marker" % (outdir))
        subprocess.call('cp  -r %s/10.WNN_RNA_Marker/* %s/10.WNN_RNA_Marker' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/10.WNN_RNA_Marker/*.log %s/10.WNN_RNA_Marker/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 10.WNN_RNA_Marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 10.WNN_RNA_Marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 11.WNN_ATAC_Marker
    if os.path.exists("%s/11.WNN_ATAC_Marker" % (program_path)):
        os.makedirs("%s/11.WNN_ATAC_Marker" % (outdir))
        subprocess.call('cp  -r %s/11.WNN_ATAC_Marker/* %s/11.WNN_ATAC_Marker' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/11.WNN_ATAC_Marker/*.log %s/11.WNN_ATAC_Marker/*.sessioninfo' % (outdir,outdir),shell=True)
    else:
        print("Can not find 11.WNN_ATAC_Marker results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 11.WNN_ATAC_Marker results!!!!!!!!!!!!!!!!!" + "\n")

    # 12.WNN_Reference_celltype
    if os.path.exists("%s/12.WNN_Reference_celltype" % (program_path)):
        os.makedirs("%s/12.WNN_Reference_celltype" % (outdir))
        subprocess.call('cp %s/12.WNN_Reference_celltype/* %s/12.WNN_Reference_celltype' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/12.WNN_Reference_celltype/celltyping*' % (outdir),shell=True)
    else:
        print("Can not find 12.WNN_Reference_celltype results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 12.WNN_Reference_celltype results!!!!!!!!!!!!!!!!!" + "\n")
    
    # 13.Diffexp_RNA
    if os.path.exists("%s/13.Diffexp_RNA" % (program_path)):
        os.makedirs("%s/13.Diffexp_RNA" % (outdir))
        subprocess.call('cp %s/13.Diffexp_RNA/* %s/13.Diffexp_RNA' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/13.Diffexp_RNA/diffexp*' % (outdir),shell=True)
    else:
        print("Can not find 13.Diffexp_RNA results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 13.Diffexp_RNA results!!!!!!!!!!!!!!!!!" + "\n")
    
    # 14.Diffexp_ATAC
    if os.path.exists("%s/14.Diffexp_ATAC" % (program_path)):
        os.makedirs("%s/14.Diffexp_ATAC" % (outdir))
        subprocess.call('cp %s/14.Diffexp_ATAC/* %s/14.Diffexp_ATAC' % (program_path, outdir),shell=True)
        #subprocess.call('rm %s/14.Diffexp_ATAC/diffexp*' % (outdir),shell=True)
    else:
        print("Can not find 14.Diffexp_ATAC results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 14.Diffexp_ATAC results!!!!!!!!!!!!!!!!!" + "\n")
    
    # 15.enrichment
    if os.path.exists("%s/15.enrichment" % (program_path)):
        os.makedirs("%s/15.enrichment" % (outdir))
        subprocess.call('cp -r %s/15.enrichment/GO_enrichment %s/15.enrichment/' % (program_path, outdir),shell=True)
        subprocess.call('cp -r %s/15.enrichment/KEGG_enrichment %s/15.enrichment/' % (program_path,outdir),shell=True)
    else:
        print("Can not find 15.enrichment results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find 15.enrichment results!!!!!!!!!!!!!!!!!" + "\n")

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

    report = Report('10x 单细胞多组学（ATAC + 基因表达）结题报告', title = '10x 单细胞多组学（ATAC + 基因表达）结题报告',header_info = dict(header_info))
    report.add_yaml_config(os.path.abspath(f'{scriptdir}/src/report.ymal'))
    os.environ['oeweb_register_token'] = 'oearray'
    Project_module = report.add_module('项目概况')
    register_info = oeweb_register(project_id=program_num, target_url='https://cloud.oebiotech.cn/task/category/scrna',
                                   note='本报告包含项目基本分析内容，如需快速做细胞表达可视化等分析(数据默认保留半年)')
    if register_info:
        Project_module.add_comment(register_info)
    report_summary_list = []
    report_hyperlink_list = []

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
    workflow_module = report.add_module('10x 单细胞多组学（ATAC + 基因表达）生物信息标准分析流程')
    workflow_module.add_plot( "%s/src/pipeline.png"%(outdir),
                              caption='10x 单细胞多组学（ATAC + 基因表达）生物信息标准分析流程图',
                              content="10x 单细胞多组学（ATAC + 基因表达）生物信息标准分析流程图")

    ##====================================================== 3.2 质控分析 =============================================
    result_module = report.add_module('项目分析结果')
    qc_section1 = result_module.add_section('质控分析')

    #################### Cell Ranger 标准分析 ####################
    qc_section2 = qc_section1.add_section('Cellranger ARC 标准质控')
    report_hyperlink_list.append('[样本质控结果](#cellranger)\n\n')
    cellranger_out_plot = qc_section2.add_plot("%s/1.Cellranger/*.png"%(outdir),
                                                           caption='样本 Cell Ranger ARC 质控结果概述',
                                                           content='样本 Cell Ranger ARC 质控结果概述。 # 图片说明：Estimated number of cells：有效细胞数；ATAC Median high-quality fragments per cell：每个细胞中 ATAC 的测序片段数量中位值；GEX Median genes per cell：每个细胞中检测到的基因数中位值')

    cellranger_out_table = qc_section2.add_table("%s/src/multimodal_qc_stat.txt" % (outdir), caption='样本细胞质量统计结果说明')

    ##====================================================== 3.3 QC 分析 ==============================================
    if os.path.exists("%s/2.Count_QC" % (outdir)):
        QC_folder = "%s/2.Count_QC" % (outdir)
        summ = pd.read_csv('%s/statitics_before_after_QC.xls' % (QC_folder),sep='\t')
        sample_num = len(summ)
        min_beforeQC_cell = min(summ.Total_cells_beforeQC)
        max_beforeQC_cell = max(summ.Total_cells_beforeQC)
        min_afterQC_cell = min(summ.Total_cells_afterQC)
        max_afterQC_cell = max(summ.Total_cells_afterQC)
        min_nFeature_afterQC = round(min(summ.mean_nFeature_RNA_afterQC.astype(float)))
        max_nFeature_afterQC = round(max(summ.mean_nFeature_RNA_afterQC.astype(float)))
        min_nCount_RNA_afterQC = round(min(summ.mean_nCount_RNA_afterQC.astype(float)))
        max_nCount_RNA_afterQC = round(max(summ.mean_nCount_RNA_afterQC.astype(float)))
        min_percent_fragment_peaks_afterQC = format(min(summ.mean_atac_fragments_afterQC.astype(float)),'.4f')
        max_percent_fragment_peaks_afterQC = format(max(summ.mean_atac_fragments_afterQC.astype(float)),'.4f')
        if len(summ)==1:
            report_summary_list.append(f"本次分析共完成 {sample_num} 个样本的10x 单细胞多组学（ATAC + 基因表达）测序，Cell Ranger ARC 定量质控的高质量细胞数为 {min_beforeQC_cell} 个，经剔除双细胞、多细胞和凋亡细胞等质控后，最终获得的细胞数目为 {min_afterQC_cell} 个，每个细胞中的平均 UMI 数为 {min_nCount_RNA_afterQC}，每个细胞中的平均基因数为 {min_nFeature_afterQC}，每个细胞中 Peak 片段数量为 {min_percent_fragment_peaks_afterQC}。")
        else:
            report_summary_list.append(f"本次分析共完成 {sample_num} 个样本的10x 单细胞多组学（ATAC + 基因表达）测序，各样本 Cell Ranger ARC 定量质控的高质量细胞数分布在 {min_beforeQC_cell}~{max_beforeQC_cell} 个，经剔除双细胞、多细胞和凋亡细胞等质控后，最终获得的细胞数目分布在 {min_afterQC_cell}~{max_afterQC_cell} 个，每个细胞中的平均 UMI 数分布在 {min_nCount_RNA_afterQC}~{max_nCount_RNA_afterQC}，每个细胞中的平均基因数分布在 {min_nFeature_afterQC}~{max_nFeature_afterQC}，每个细胞中 Peak 片段数量分布在 {min_percent_fragment_peaks_afterQC}~{max_percent_fragment_peaks_afterQC}")
        if  "percent.mito_higher" in summ.columns:
            min_percent_mito_afterQC = format(min(summ["mean_percent.mito_afterQC"].astype(float)),'.4f')
            max_percent_mito_afterQC = format(max(summ["mean_percent.mito_afterQC"].astype(float)),'.4f')
            percent_txt = f"线粒体基因所占比例（percent.mito）、"
            percent_descript = f"；5.线粒体基因比例低于 15%"
            if len(summ)==1:
                report_summary_list.append(f"，每个细胞中平均线粒体基因比例为 {min_percent_mito_afterQC}。")
            else:
                report_summary_list.append(f"，每个细胞中平均线粒体基因比例分布在 {min_percent_mito_afterQC}~{max_percent_mito_afterQC}。")
        else:
            report_summary_list.append(f"。")

        ##
        QC_description=f"在 Cell Ranger ARC 初步质控的基础上进一步对实验数据进行质控，剔除多细胞、双细胞或者未结合上的低质量细胞的数据。理论上大部分细胞中表达的基因数量、UMI 数量、fragments的数目、线粒体基因数量等会集中分布在某一区域内，根据它们的分布特征可以拟合分布模型，使用该模型找到其中的离域值，剔除异常数据。本项目中的质控标准为：1.保留细胞基因数大于200和 UMI 数大于1000的细胞；2.转录起始位点富集分数（TSS.enrichment）大于2；3.fragments数在1000~50000范围内；4.缠绕核小体的片段与非核小体片段（nucleosome_signal）的比例小于2{percent_descript}，以同时满足以上标准的细胞作为高质量细胞，进行下游分析。" \
            f"\n 质控前后每个细胞的基因数量（nFeature_RNA）、UMI 数量（nCount_RNA）、{percent_txt}染色体fragments的数目（atac_fragments）、转录起始位点富集分数（TSS.enrichment）、缠绕核小体的片段与非核小体片段的比例（nucleosome_signal）的小提琴图展示如下："
        ##
        qc_section3 = qc_section1.add_section(name='过滤低质量细胞',description=QC_description)
        qc_section3.add_comment("详细结果见目录：[2.Count_QC](./2.Count_QC)")
        qc_section3.add_plot( '%s/2.Count_QC/QC_metrics_beforeQC.png'%(outdir),caption='过滤低质量细胞前各质控指标的小提琴图')
        qc_section3.add_comment(
            "过滤低质量细胞前各质控指标的小提琴图：[2.Count_QC/QC_metrics_beforeQC.pdf](./2.Count_QC/QC_metrics_beforeQC.pdf)")
        qc_section3.add_plot('%s/2.Count_QC/QC_metrics_afterQC.png'%(outdir),caption='过滤低质量细胞后各质控指标的小提琴图')
        qc_section3.add_comment(
            "过滤低质量细胞后各质控指标的小提琴图：[2.Count_QC/QC_metrics_afterQC.pdf](./2.Count_QC/QC_metrics_afterQC.pdf)")
        qc_section3.add_comment("定量质控前后的细胞数统计情况及过滤标准如下表所示：")
        qc_section3.add_table(summ.iloc[:,:30],
                              caption = '质控前后细胞数目统计表',
                              headerdict={'sample':'样本名',
                                          'mean_nFeature_RNA_beforeQC':'质控前RNA平均Feature数',
                                          'median_nFeature_RNA_beforeQC':'质控前RNA feature的中位数',
                                          'mean_nCount_RNA_beforeQC':'质控前RNA平均Counts数',
                                          'median_nCount_RNA_beforeQC':'质控前RNA Counts的中位数',
                                          'mean_percent.mito_beforeQC':'质控前平均线粒体比例',
                                          'median_percent.mito_beforeQC':'质控前线粒体比例中位数',
                                          'mean_atac_fragments_beforeQC':'质控前ATAC的平均fragments数',
                                          'median_atac_fragments_beforeQC':'质控前ATAC fragments的中位数',
                                          'mean_TSS.enrichment_beforeQC':'质控前TSS平均富集分数',
                                          'median_TSS.enrichment_beforeQC':'质控前TSS富集分数的中位数',
                                          'mean_nucleosome_signal_beforeQC':'质控前核小体片段的平均比例',
                                          'median_nucleosome_signal_beforeQC':'质控前核小体片段比例中位数',
                                          'mean_nCount_ATAC_beforeQC':'质控前ATAC的平均Counts数',
                                          'median_nCount_ATAC_beforeQC':'质控前ATAC的Counts数的中位数',
                                          'Total_cells_beforeQC':'质控前总的细胞数'})

   ##=============================================== 3.4基因表达水平分析 ==============================================
    rna_section1 = result_module.add_section('基因表达水平分析')

    ## 降维与聚类
    if os.path.exists("%s/3.RNA_Clustering" % (outdir)):
        report_hyperlink_list.append('[基因表达水平降维与聚类分析结果](#rna_cluster)\n\n')
        if len(summ)==1:
            sample_des=f"仅有1个样本，并不存在批次效应"
        else:
            sample_des = f"有不同样本，但并不存在批次效应"
        if config["cluster_params"]["RNA_reduct1"]=="pca":
            reduct1_des=f"PCA(principal component analysis, 主成分分析)的降维算法进行降维"
            reduct1_name="PCA"
        else:
            sample_des = f"有不同样本，且存在批次效应"
            reduct1_des=f"MNN(mutual nearest neighbors, 互享最近邻) 的降维算法剔除批次效应"
            reduct1_name="MNN"
        if config["cluster_params"]["RNA_reduct2"]=="tsne":
            reduct2_des=f"t-SNE(t分布式随机邻近嵌入)算法进行二级降维。t-SNE创建了一个缩小的特征空间，相似的样本由附近的点建模，不相似的样本由高概率的远点建模，更好的展示了样本之间的局部关系"
            reduct2_name="t-SNE"
        else:
            reduct2_des=f"UMAP(统一流形逼近与投影)算法进行二级降维。UMAP 主要基于流形理论和拓扑算法的理论，对高维数据进行降维，从而能够保留更多数据的全局结构，并且具有优越的运行性能"
            reduct2_name="UMAP"

        rna_section2 = rna_section1.add_section("基因表达数据降维与聚类分析")
        rna_section2.add_comment(f"由于本项目{sample_des}，因此我们首先采用了{reduct1_des}，再通过 {reduct2_des}。基于{reduct1_name}降维结果通过{reduct2_name}对单细胞群聚类进行可视化，聚类算法采用SNN，最终获得最优细胞分群。")
        rna_section2.add_comment("详细结果见目录：[3.RNA_Clustering](./3.RNA_Clustering)")

        rna_section2.add_plot('%s/3.RNA_Clustering/*/*_groupby_cluster_resolution*_plot.png'%(outdir),
                              caption='基因表达水平降维与聚类图',
                              description='图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')
        rna_section2.add_comment("基因表达水平降维与聚类图：[3.RNA_Clustering/\*/\*_groupby_clusters_resolution*_plot.png](./3.RNA_Clustering/)")

    ## 样本分组展示
        if (sample_num > 1):
            rna_section3 = rna_section1.add_section('样本间降维聚类分组展示')
            rna_section3.add_plot('%s/3.RNA_Clustering/visualize_cluster_by_*/groupby-sampleid_contrast_plot.png' % (outdir),
                                  caption='多样本降维聚类分组展示图',
                                  description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。')
            rna_section3.add_comment("基因表达水平多个样本的降维聚类分组展示图：[3.RNA_Clustering/visualize_cluster_by_\*/\*groupby-sampleid_contrast_plot.pdf](./3.RNA_Clustering/)")
        ################
            rna_section3.add_comment(" ")
            rna_section3.add_comment("每个细胞群中样本占比的柱状统计图如下：")
            rna_section3.add_plot('%s/3.RNA_Clustering/visualize_cluster_by_*/groupby-*.res.*_summary_plot.png' % (outdir),
                                  caption='每个细胞群中样本占比的柱状统计图',
                                  description='图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。')
       
        ################
            rna_section3.add_comment(" ")
            rna_section3.add_comment("多个样本的降维聚类分面展示图如下：")
            rna_section3.add_plot('%s/3.RNA_Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.png' % (outdir),
                                  caption='多样本降维聚类分面展示图',
                                  description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同群的细胞以不同颜色区分。')
            rna_section3.add_comment("多个样本的降维聚类分面展示图如下：[3.RNA_Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.pdf](./3.RNA_Clustering/)")
        
        ################
            rna_section3.add_comment(" ")
            rna_section3.add_comment("每个样本中细胞群占比的柱状统计图如下：")
            rna_section3.add_plot('%s/3.RNA_Clustering/visualize_cluster_by_*/groupby-sampleid_summary_plot.png' % (outdir),
                                  caption='样本间细胞群占比柱状图',
                                  description='图片说明：横坐标表示不同样本，纵坐标表示不同组别中细胞数目所占的百分比。')
        ###############
            rna_section3.add_comment("各样本在不同细胞群中的细胞数目统计表如下：")
            rna_section3.add_table("%s/3.RNA_Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls" % (outdir),
                                   caption='样本间细胞群数目统计表')
            rna_section3.add_comment("各样本在不同细胞群中的细胞数目统计表：[3.RNA_Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls](./3.RNA_Clustering/)")

    ## Marker 基因鉴定
    if os.path.exists("%s/4.RNA_Marker" % (outdir)):
        report_hyperlink_list.append('[基因表达水平 Marker 基因鉴定结果](#rna_marker_gene)')
        test_method=config["marker_params"]["RNA_test_method"]
        Marker_description=f"Marker 基因的定义是该基因在指定细胞群的绝大多数细胞中有较高的表达，而在其余细胞类群中只有少部分表达，且该基因在此细胞群相对于其他细胞群中是显著上调表达。\n 本研究采用 {test_method} 检验方法对指定细胞群与其余所有细胞群进行差异检验，从而筛选得到每个细胞群的特异性 Marker 基因"
        rna_section4 = rna_section1.add_section(name='Marker 基因鉴定',description=Marker_description)
        #####
        marker_result=pd.read_csv('%s/4.RNA_Marker/top10_markers_for_each_cluster.xls' % (outdir), sep='\t')
        cluster_num = max(marker_result['cluster'])
        report_summary_list.append(f"通过基因表达数据进行降维聚类后共分为 {cluster_num} 群细胞，")
        ####
        rna_section4.add_comment("每个细胞群中 Top10 Marker 基因结果表格：[4.RNA_Marker/top10_markers_for_each_cluster.xls](./4.RNA_Marker/top10_markers_for_each_cluster.xls)")
        rna_section4.add_table("%s/4.RNA_Marker/top10_markers_for_each_cluster.xls"%(outdir),caption='每个细胞群 Top10 Marker 基因列表', show_rows=5)
        ####
        rna_section4.add_comment("Top10 Marker 基因表达热图：[4.RNA_Marker/topmarker_gene_heatmap.pdf](./4.RNA_Marker/topmarker_gene_heatmap.pdf)")
        rna_section4.add_plot("%s/4.RNA_Marker/topmarker_gene_heatmap.png"%(outdir),
                              caption='Top10 Marker 基因表达热图',
                              description='图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')
        ####
        rna_section4.add_comment("Top10 Marker 基因在降维聚类结果中的可视化图：[4.RNA_Marker/marker_gene_featureplot.pdf](./4.RNA_Marker/)")
        rna_section4.add_plot("%s/4.RNA_Marker/*/marker_gene_featureplot.png"%(outdir),
                              caption='Top10 Marker 基因在降维聚类结果中的可视化图',
                              description='图片说明：红色越深表示该细胞中对应基因的表达量越高。')
        ####
        rna_section4.add_comment("Top10 Marker 基因表达量小提琴图：[4.RNA_Marker/marker_gene_violin_plot.pdf](./4.RNA_Marker/)")
        rna_section4.add_plot("%s/4.RNA_Marker/*/marker_gene_violin_plot.png"%(outdir),
                              caption='Top10 Marker 基因表达量小提琴图',
                              description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')
    
    ## 细胞类型鉴定
    if os.path.exists("%s/5.RNA_Reference_celltype" % (outdir)):
        report_hyperlink_list.append('[基因表达水平细胞类型鉴定结果](#rna_celltype)\n\n')
        celltype_folder = "%s/5.RNA_Reference_celltype" % (outdir)
        celltype_file = list(glob('%s/*simplified_celltyping_results.csv' % (celltype_folder)))[0]
        cell_results = pd.read_csv(celltype_file,sep=',')
        celltype = ', '.join(set(cell_results.celltype))
        report_summary_list.append(f"通过对细胞类型进行注释，得到供参考的细胞类型有{celltype}。")

        rna_section5 = rna_section1.add_section('细胞类型鉴定')
        link = list(glob(f'report/{program_num}_Report_{report_time}/5.RNA_Reference_celltype/*top.*celltyping_plot.pdf'))[0]
        refcelltype = str(link).split('/')[-1].split("_")[1]
        rna_section5.add_comment(name=refcelltype.lower())

        rna_section5.add_comment("数据集方法细胞类型鉴定参考结果见：[5.RNA_Reference_celltype/\*ref_\*_top.main_celltyping_plot.pdf](./5.RNA_Reference_celltype/*ref_*_top.main_celltyping_plot.pdf)")
        rna_section5.add_plot("%s/5.RNA_Reference_celltype/*ref_*_top.main_celltyping_plot.png"%(outdir),
                              caption='细胞类型参考结果图',
                              description='图片说明：细胞类型结果在UMAP图上的展示，每种细胞类型以不同颜色进行区分。')
        rna_section5.add_comment("细胞类型鉴定相关性热图结果见：[5.RNA_Reference_celltype/\*ref_\*_main_celltyping_heatmap.pdf](./5.RNA_Reference_celltype/*ref_*_main_celltyping_heatmap.pdf)")
        rna_section5.add_plot("%s/5.RNA_Reference_celltype/*ref_*_main_celltyping_heatmap.png"%(outdir),
                              caption='细胞类型鉴定相关性热图',
                              description='图片说明：每一行表示参考数据集中的细胞类型注释名称，每一列表示待鉴定的细胞。黄色越深表示相关性值越大，表明待鉴定的细胞类型与参考数据集中的该种细胞类型最为相似。')
        rna_section5.add_comment("细胞类型注释表格见：[5.RNA_Reference_celltype/*ref_*_main_simplified_celltyping_results.csv](./5.RNA_Reference_celltype/*ref_*_main_simplified_celltyping_results.csv)")
        rna_section5.add_table("%s/5.RNA_Reference_celltype/*ref_*_main_simplified_celltyping_results.csv"%(outdir),
                               caption='细胞类型注释结果表格',
                               show_rows=5)

    ##==========================================  3.5 染色质可及性水平分析 ===============================================
    atac_section1 = result_module.add_section('染色质可及性水平分析')

    ## 降维与聚类分析
    if os.path.exists("%s/6.ATAC_Clustering" % (outdir)):
        report_hyperlink_list.append('[染色质可及性水平降维与聚类分析结果](#atac_cluster)\n\n')
        atac_section2 = atac_section1.add_section('染色质可及性水平降维与聚类分析')
        atac_section2.add_plot('%s/6.ATAC_Clustering/umap_Dimension_Reduction/*_groupby_cluster_resolution*_plot.png'%(outdir),
                               caption='染色质可及性水平降维与聚类图',
                               description='图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')
        atac_section2.add_comment("染色质可及性水平降维与聚类图：[6.ATAC_Clustering/umap_Dimension_Reduction/\*_groupby_clusters_resolution\*_plot.png](./6.ATAC_Clustering/)")
    
    ## 染色质可及性水平样本分组展示
        if (sample_num > 1):
            atac_section3 = atac_section1.add_section('样本间降维聚类分组展示')
            atac_section3.add_plot('%s/6.ATAC_Clustering/visualize_cluster_by_*/groupby-sampleid_contrast_plot.png' % (outdir),
                                   caption='多样本降维聚类分组展示图',
                                   description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。')
            atac_section3.add_comment("染色质可及性水平多个样本的降维聚类分组展示图：[6.ATAC_Clustering/visualize_cluster_by_\*/\*groupby-sampleid_contrast_plot.pdf](./6.ATAC_Clustering/)")
        ####
            atac_section3.add_comment(" ")
            atac_section3.add_comment("每个细胞群中样本占比的柱状统计图如下：")
            atac_section3.add_plot('%s/6.ATAC_Clustering/visualize_cluster_by_*/groupby-*.res.*_summary_plot.png' % (outdir),
                                   caption='每个细胞群中样本占比的柱状统计图',
                                   description='图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。')
        ####
            atac_section3.add_comment(" ")
            atac_section3.add_comment("多个样本的降维聚类分面展示图如下：")
            atac_section3.add_plot('%s/6.ATAC_Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.png' % (outdir),
                                   caption='多个样本的降维聚类分面展示图',
                                   description='图片说明：横纵坐标分别代表降维第一和第二主成分，不同细胞群的细胞以不同颜色区分。')
            atac_section3.add_comment("多个样本的降维聚类分面展示图：[6.ATAC_Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.pdf](./6.ATAC_Clustering/)")
        ####
            atac_section3.add_comment(" ")
            atac_section3.add_comment("每个样本中细胞群占比的柱状统计图如下：")
            atac_section3.add_plot('%s/6.ATAC_Clustering/visualize_cluster_by_*/groupby-sampleid_summary_plot.png' % (outdir),
                                   caption='每个样本中细胞群占比的柱状统计图',
                                   description='图片说明：横坐标表示不同样本，纵坐标表示不同组别中细胞数目所占的百分比。')
        ####
            atac_section3.add_comment(" ")
            atac_section3.add_comment("各样本在不同细胞群中的细胞数目统计表如下：")
            atac_section3.add_table("%s/6.ATAC_Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls" %(outdir),caption='样本间细胞群数目统计表')
            atac_section3.add_comment("各样本在不同细胞群中的细胞数目统计表：[6.ATAC_Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls](./6.ATAC_Clustering/)")

    ## 染色质可及性水平 Marker Peak 鉴定
    if os.path.exists("%s/7.ATAC_Marker" % (outdir)):
        report_hyperlink_list.append('[染色质可及性水平 Marker Peak 鉴定结果](#atac_marker_peak)\n\n')
        test_ATAC_method = config["marker_params"]["ATAC_test_method"]
        atac_section4 = atac_section1.add_section('染色质可及性水平 Marker Peak 鉴定',test_ATAC_method=test_ATAC_method)
        ####
        marker_result = pd.read_csv('%s/7.ATAC_Marker/top10_markers_for_each_cluster.xls' % (outdir),sep='\t')
        cluster_num = max(marker_result['cluster'])
        report_summary_list.append(f"通过染色质可及性数据降维聚类后共分为 {cluster_num} 群细胞，")
        ####
        atac_section4.add_table("%s/7.ATAC_Marker/top10_markers_for_each_cluster.xls"%(outdir),
                                caption='每个细胞群 Top10 Marker Peak 列表',
                                show_rows=5)
        atac_section4.add_comment("每个细胞群中所有 Marker Peak 结果表格：[7.ATAC_Marker/all_markers_for_each_cluster.xls](./7.ATAC_Marker/all_markers_for_each_cluster.xls)")
        atac_section4.add_comment("每个细胞群中 Top10 Marker Peak 结果表格：[7.ATAC_Marker/top10_markers_for_each_cluster.xls](./7.ATAC_Marker/top10_markers_for_each_cluster.xls)")
        atac_section4.add_plot("%s/7.ATAC_Marker/topmarker_gene_heatmap.png"%(outdir),
                               caption='Top10 Marker Peak 开放程度热图',
                               description='图片说明：横坐标为细胞群，纵坐标为 Marker Peak，图中黄色表示高开放，紫色表示低开放。')
        atac_section4.add_comment("Top10 Marker Peak 开放程度热图：[7.ATAC_Marker/topmarker_gene_heatmap.pdf](./7.ATAC_Marker/topmarker_gene_heatmap.pdf)")
        atac_section4.add_plot("%s/7.ATAC_Marker/*/marker_gene_featureplot.png"%(outdir),
                               caption='Top10 Marker Peak 在降维聚类结果中的可视化图',
                               description='图片说明：红色越深表示该细胞中对应基因的开放程度越高。')
        atac_section4.add_comment("Top10 Marker Peak 在降维聚类结果中的可视化图：[7.ATAC_Marker/*/marker_gene_featureplot.pdf](./7.ATAC_Marker/)")
        atac_section4.add_plot("%s/7.ATAC_Marker/*/marker_gene_violin_plot.png"%(outdir),
                               caption='Top10 Marker Peak 开放程度小提琴图',
                               description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的 Peak 开放程度。')
        atac_section4.add_comment("Top10 Marker Peak 开放程度小提琴图：[7.ATAC_Marker/*/marker_gene_violin_plot.pdf](./7.ATAC_Marker/)")
    
    ## 染色质可及性水平细胞类型鉴定
    if os.path.exists("%s/8.ATAC_Reference_celltype" % (outdir)):
        report_hyperlink_list.append('[染色质可及性水平细胞类型映射展示](#atac_celltype)\n\n')
        atac_section5 = atac_section1.add_section('染色质可及性水平细胞类型映射展示')
        atac_section5.add_comment("染色质可及性水平细胞类型鉴定结果见：[8.ATAC_Reference_celltype/\*ref_\*_top.main_celltyping_plot.pdf](./8.ATAC_Reference_celltype/*ref_*_top.main_celltyping_plot.png)")
        atac_section5.add_plot("%s/8.ATAC_Reference_celltype/*ref_*_top.main_celltyping_plot.png"%(outdir),
                               caption='细胞类型参考结果图',
                               description='图片说明：细胞类型结果在UMAP图上的展示，每种细胞类型以不同颜色进行区分。')

    ##=================================================== 3.6 WNN 整合分析 =============================================
    wnn_section1 = result_module.add_section('WNN 整合分析')

    ## WNN降维聚类
    if os.path.exists("%s/9.WNN_Clustering" % (outdir)):
        report_hyperlink_list.append('[WNN 整合分析结果](#wnn)\n\n')
        wnn_section2 = wnn_section1.add_section('降维和聚类分析')
        wnn_section2.add_plot('%s/9.WNN_Clustering/wnn_groupby_cluster_*_plot.png'%(outdir),
                          caption='WNN 整合分析降维与聚类图',
                          description='图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。',
                          content='WNN 整合分析降维与聚类图。图片说明：横纵坐标分别代表降维的第一和第二主成分，图中的每个点代表一个细胞，不同群的细胞以不同颜色区分。')
        wnn_section2.add_comment('WNN 整合分析降维与聚类图:[9.WNN_Clustering/wnn_groupby_cluster_*_plot.pdf](./9.WNN_Clustering/)')
    ####
        if (sample_num > 1):
            wnn_section2.add_comment('WNN多样本降维聚类分组展示图：')
            wnn_section2.add_plot('%s/9.WNN_Clustering/visualize_cluster_by_*/groupby-sampleid_contrast_plot.png' % (outdir),
                                  caption='多样本降维聚类分组展示图',
                                  description='图片说明：横纵坐标分别代表降维的第一和第二主成分，不同样本来源的细胞以不同颜色区分。')
            wnn_section2.add_comment("多个样本的降维聚类分组展示图：[9.WNN_Clustering/visualize_cluster_by_\*/\*groupby-sampleid_contrast_plot.pdf](./9.WNN_Clustering/)")
        ####
            wnn_section2.add_comment(" ")
            wnn_section2.add_comment("每个细胞群中样本占比的柱状统计图如下：")
            wnn_section2.add_plot('%s/9.WNN_Clustering/visualize_cluster_by_*/groupby-*.res.*_summary_plot.png' % (outdir),
                                  caption='每个细胞群中样本占比的柱状统计图',
                                  description='图片说明：横坐标表示不同细胞群，纵坐标表示不同组别中细胞数目所占的百分比。')
        ####
            wnn_section2.add_comment(" ")
            wnn_section2.add_comment("多个样本的降维聚类分面展示图如下：")
            wnn_section2.add_plot('%s/9.WNN_Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.png' % (outdir),
                                  caption='多个样本的降维聚类分面展示图',
                                  description='图片说明：横纵坐标分别代表降维第一和第二主成分，不同细胞群的细胞以不同颜色区分。')
            wnn_section2.add_comment("多个样本的降维聚类分面展示图：[9.WNN_Clustering/visualize_cluster_by_*/splitby-sampleid_split_plot.pdf](./9.WNN_Clustering/)")
        ####
            wnn_section2.add_comment(" ")
            wnn_section2.add_comment("每个样本中细胞群占比的柱状统计图如下：")
            wnn_section2.add_plot('%s/9.WNN_Clustering/visualize_cluster_by_*/groupby-sampleid_summary_plot.png' % (outdir),
                                  caption='每个样本中细胞群占比的柱状统计图',
                                  description='图片说明：横坐标表示不同样本，纵坐标表示不同组别中细胞数目所占的百分比。')
        ####
            wnn_section2.add_comment(" ")
            wnn_section2.add_comment("各样本在不同细胞群中的细胞数目统计表如下：")
            wnn_section2.add_table("%s/9.WNN_Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls" % (outdir),caption='样本间细胞群数目统计表')
            wnn_section2.add_comment("各样本在不同细胞群中的细胞数目统计表：[9.WNN_Clustering/visualize_cluster_by_*/clust_cond_freq_info.xls](./9.WNN_Clustering/)")

    ## Marker RNA
    if os.path.exists("%s/10.WNN_RNA_Marker" % (outdir)):
        wnn_cluster_file = list(glob('%s/10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls' % (outdir)))[0]
        wnn_cluster_num = max(pd.read_csv(wnn_cluster_file,sep='\t')['cluster'])
        report_summary_list.append(f"通过WNN 整合分析后共分为 {wnn_cluster_num} 群细胞。")
        wnn_section3 = wnn_section1.add_section('Marker基因鉴定')
        wnn_section3.add_comment('WNN 整合后 Marker 基因鉴定：')
        wnn_section3.add_table("%s/10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls"%(outdir),
                               caption='每个细胞群 Top10 Marker 基因列表',
                               show_rows=5)
        wnn_section3.add_comment("每个细胞群中所有 Marker 基因结果表格：[10.WNN_RNA_Marker/all_markers_for_each_cluster.xls](./10.WNN_RNA_Marker/all_markers_for_each_cluster.xls)")
        wnn_section3.add_comment("每个细胞群中 Top10 Marker 基因结果表格：[10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls](./10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls)")
        wnn_section3.add_plot("%s/10.WNN_RNA_Marker/topmarker_gene_heatmap.png"%(outdir),
                              caption='Top10 Marker 基因表达热图',
                              description='图片说明：横坐标为细胞群，纵坐标为 Marker 基因，图中黄色表示高表达，紫色表示低表达。')
        wnn_section3.add_comment("Top10 Marker 基因表达热图：[10.WNN_RNA_Marker/topmarker_gene_heatmap.pdf](./10.WNN_RNA_Marker/topmarker_gene_heatmap.pdf)")
        wnn_section3.add_plot("%s/10.WNN_RNA_Marker/*/marker_gene_featureplot.png"%(outdir),
                              caption='Top10 Marker 基因在降维聚类结果中的可视化图',
                              description='图片说明：红色越深表示该细胞中对应基因的表达量越高。')
        wnn_section3.add_comment("Top10 Marker 基因在降维聚类结果中的可视化图：[10.WNN_RNA_Marker/marker_gene_featureplot.pdf](./10.WNN_RNA_Marker/)")
        wnn_section3.add_plot("%s/10.WNN_RNA_Marker/*/marker_gene_violin_plot.png"%(outdir),
                              caption='Top10 Marker 基因表达量小提琴图',
                              description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的基因表达值。')
        wnn_section3.add_comment("Top10 Marker 基因表达量小提琴图：[10.WNN_RNA_Marker/*/marker_gene_violin_plot.pdf](./10.WNN_RNA_Marker/)")
    
    ### Mark ATAC
    if os.path.exists("%s/11.WNN_ATAC_Marker" % (outdir)):
        wnn_section4 = wnn_section1.add_section('Marker peak鉴定')
        wnn_section4.add_comment('WNN 整合后 Marker Peak 鉴定：')
        wnn_section4.add_table("%s/11.WNN_ATAC_Marker/top10_markers_for_each_cluster.xls"%(outdir),
                               caption='每个细胞群 Top10 Marker Peak 列表',
                               show_rows=5)
        wnn_section4.add_comment("每个细胞群中所有 Marker Peak 结果表格：[11.WNN_ATAC_Marker/all_markers_for_each_cluster.xls](./11.WNN_ATAC_Marker/all_markers_for_each_cluster.xls)")
        wnn_section4.add_comment("每个细胞群中 Top10 Marker Peak 结果表格：[11.WNN_ATAC_Marker/top10_markers_for_each_cluster.xls](./11.WNN_ATAC_Marker/top10_markers_for_each_cluster.xls)")
        wnn_section4.add_plot("%s/11.WNN_ATAC_Marker/topmarker_gene_heatmap.png"%(outdir),
                              caption='Top10 Marker Peak 开放程度热图',
                              description='图片说明：横坐标为细胞群，纵坐标为 Marker Peak，图中黄色表示高开放，紫色表示低开放。')
        wnn_section4.add_comment("Top10 Marker Peak 开放程度热图：[11.WNN_ATAC_Marker/topmarker_gene_heatmap.pdf](./11.WNN_ATAC_Marker/topmarker_gene_heatmap.pdf)")
        wnn_section4.add_plot("%s/11.WNN_ATAC_Marker/*/marker_gene_featureplot.png"%(outdir),
                              caption='Top10 Marker Peak 在降维聚类结果中的可视化图',
                              description='图片说明：红色越深表示该细胞中对应基因的开放程度越高。')
        wnn_section4.add_comment("Top10 Marker Peak 在降维聚类结果中的可视化图：[11.WNN_ATAC_Marker/*/marker_gene_featureplot.pdf](./11.WNN_ATAC_Marker/)")
        wnn_section4.add_plot("%s/11.WNN_ATAC_Marker/*/marker_gene_violin_plot.png"%(outdir),
                              caption='Top10 Marker Peak 开放程度小提琴图',
                              description='图片说明：横坐标为细胞群编号，纵坐标为标准化后的 Peak 开放程度。')
        wnn_section4.add_comment("Top10 Marker Peak 开放程度小提琴图：[11.WNN_ATAC_Marker/*/marker_gene_violin_plot.pdf](./11.WNN_ATAC_Marker/)")
    ### WNN celltype
    if os.path.exists("%s/12.WNN_Reference_celltype" % (outdir)):
        wnn_section5 = wnn_section1.add_section("WNN 整合后细胞类型映射展示")
        wnn_section5.add_comment("WNN 整合后细胞类型映射展示结果见：[12.WNN_Reference_celltype/*ref_\*_top.main_celltyping_plot.pdf](./12.WNN_Reference_celltype/*ref_*_top.main_celltyping_plot.pdf)")
        wnn_section1.add_plot('%s/12.WNN_Reference_celltype/*ref_*_top.main_celltyping_plot.png'%(outdir),
                              caption='细胞类型参考结果图',
                              description = '图片说明：细胞类型注释结果在降维聚类结果图上的展示，每种细胞类型以不同颜色区分。')

    ##=========================================== 3.7 差异表达基因/peak筛选 ===============================================
    if os.path.exists("%s/13.Diffexp_RNA" % (outdir)):
        report_hyperlink_list.append('[差异表达基因/peak筛选](#diffexp)\n\n')
        diffexp_folder = "%s/13.Diffexp_RNA" % (outdir)
        diffexp_file = list(glob('%s/diffexp_results_stat.xls' % (diffexp_folder)))[0]
        diff = pd.read_csv(diffexp_file,sep='\t')
        num_diff = len(diff)
        diff_num = ', '.join([str(i) for i in diff.iloc[:,4]])
        report_summary_list.append(f"共设有{num_diff}个差异基因分组，其检测到的差异基因数量为：{diff_num}。")
        test_wnn_method = config['diffexp_params']['test']
        if test_wnn_method=='presto':
            diff_descript=f'差异表达基因使用presto包的WilcoxAUC函数进行筛选，默认使用 presto 差异检验方法；差异表达peak同样采用presto包的WilcoxAUC函数进行筛选，采用presto进行差异检验。'
        else:
            diff_descript=f'差异表达基因使用Seurat包的FindMarkers函数进行筛选，默认使用 {test_wnn_method} 差异检验方法；差异表达peak同样采用Seurat包的FindMarkers函数进行筛选，采用{test_wnn_method}进行差异检验。'
        diff_section1 = result_module.add_section(name = '差异表达基因/peak筛选',description=diff_descript)
        diff_section1.add_comment('差异表达基因具体结果见: [13.Diffexp_RNA](./13.Diffexp_RNA)')
        diff_section1.add_comment('差异peak具体结果见: [14.Diffexp_ATAC](./14.Diffexp_ATAC)')
    ##差异表达基因
        diff_section2 = diff_section1.add_section("差异表达基因筛选")
        diff_section2.add_comment("所有差异表达基因结果见：[13.Diffexp_RNA/\*_\*-vs-*-all_diffexp_genes_anno.xls](./13.Diffexp_RNA/*_*-vs-*-all_diffexp_genes_anno.xls)")
        diff_section2.add_comment("基于以上差异结果，进一步根据差异倍数（FoldChange）及显著性检验（pvalue）筛选显著的差异表达基因，结果见: [13.Diffexp_RNA/\*_\*-vs-\*-diff-pval-\*-FC-\*_anno.xls](./13.Diffexp_RNA/*_*-vs-*-diff-pval-*-FC-*_anno.xls)")
        diff_section2.add_comment("各分组差异基因数目统计表如下：")
        diff_section2.add_table("%s/13.Diffexp_RNA/diffexp_results_stat.xls"%(outdir),caption='差异表达基因统计表')
        diff_section2.add_comment("差异显著基因结果示例：")
        difffiles = list(glob(f'{diffexp_folder}/*diff*pval*xls'))
        diff_section2.add_table(difffiles[0],
                                caption='差异显著基因结果表格',
                                show_rows=5)
        diff_section2.add_comment("将差异基因的差异倍数（FoldChange）从大到小排列，上下调各选取 25 个基因绘制热图：")
        topdiff_heatmap_plots = list(glob(f'{diffexp_folder}/*_heatmap.png'))
        topdiff_heatmap_plots = [ str(i) for i in topdiff_heatmap_plots]
        topdiff_heatmap_prefix = [ re.findall("top(.*)_heatmap.png",i)[0] for i in topdiff_heatmap_plots]
        topdiff_heatmap_contents = [j + ' 上下调 Top25 差异基因热图。图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足20个，则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。'
            for j in topdiff_heatmap_prefix]
        diff_section2.add_plot(list(zip(topdiff_heatmap_plots, topdiff_heatmap_contents)),
                               caption='上下调Top 25差异基因热图',
                               description='图片说明：横坐标为差异分组信息，纵坐标为上下调 Top25 基因（如果上下调差异基因不足25个，则绘制全部基因；线粒体基因和核糖体基因默认不进行绘图）。图中黄色表示高表达，紫色表示低表达。')
    ## 差异peak
    if os.path.exists("%s/14.Diffexp_ATAC" % (outdir)):
        diff_section3 = diff_section1.add_section("差异peak筛选")
        diff_section3.add_comment("所有差异peak结果见：[14.Diffexp_ATAC/\*_\*-vs-\*-all_diffexp_genes_anno.xls](./14.Diffexp_ATAC/*_*-vs-*-all_diffexp_genes_anno.xls)")
        diff_section3.add_comment("基于以上差异结果，进一步根据差异倍数（FoldChange）及显著性检验（pvalue）筛选显著的差异表达基因，结果见: [14.Diffexp_ATAC/\*_\*-vs-\*-diff-pval-\*-FC-\*_anno.xls](./14.Diffexp_ATAC/*_*-vs-*-diff-pval-*-FC-*_anno.xls)")
        diff_section3.add_comment("各分组差异peak数目统计表如下：")
        diff_section3.add_table("%s/14.Diffexp_ATAC/diffexp_results_stat.xls"%(outdir),caption='差异表达基因统计表')
    
    ##==================================================== 3.8差异表达基因富集分析=======================================
    if os.path.exists("%s/15.enrichment/GO_enrichment/" % (outdir)):
        report_hyperlink_list.append('[差异表达基因富集分析](#diffexp_GO_KEGG)\n\n')
        diffenrich_section1 = result_module.add_section('差异表达基因富集分析')

    ## GO 富集分析
        diffenrich_section2 = diffenrich_section1.add_section("差异基因GO富集分析")
        diffenrich_section2.add_comment("GO富集分析具体结果文件为：[15.enrichment/GO_enrichment/](./15.enrichment/GO_enrichment/)")
        diffenrich_section2.add_comment("GO富集分析汇总表如下：")
        diffenrich_section2.add_table("%s/15.enrichment/GO_enrichment/enrichment_go.xls"%(outdir),caption='GO富集分析汇总表')
        diffenrich_section2.add_comment("GO富集分析结果如下：")
        diffenrich_table = list(glob("%s/15.enrichment/GO_enrichment/*/enrichment-go-*_*-vs-*-Total.xls"%(outdir)))[0]
        diffenrich_section2.add_table(diffenrich_table,
                                      caption='GO富集分析结果',
                                      show_rows=5)
        diffenrich_section2.add_comment("GO 富集分析 top30 （筛选三种分类中对应差异基因数目大于 2 的 GO 条目，按照每个条目对应的 -log10pValue 由大到小排序的各 10 条）条形图展示如下")
        diffenrich_section2.add_plot("%s/15.enrichment/GO_enrichment/*/GO.top.*.png"%(outdir),
                                     caption='GO富集分析结果展示',
                                     description='图片说明：图中横轴为 GO 条目名称，纵轴为 -log10pValue。')
        diffenrich_section2.add_comment("使用 fisher 算法分别对样本间差异基因进行 CC，BP，MF 富集分析，并使用 topGO对富集到的 Term 绘制有向无环图。topGO 有向无环图能直观展示差异表达基因富集的 GO 节点（Term）及其层级关系，是差异表达基因 GO 富集分析的结果图形化展示，分支代表的包含关系，从上至下所定义的功能描述范围越来越具体。")
        diffenrich_section2.add_plot("%s/15.enrichment/GO_enrichment/*/topGO_*_*.png"%(outdir),
                                     caption='差异基因topGO有向无环示例图展示',
                                     description='图片说明：对每个 GO Term 进行富集，最显著的 10 个节点用矩形表示。矩形的颜色代表富集显著性，从黄色到红色显著性越来越高。每个节点的基本信息显示在相应的图形中，为 GO ID 和 GO Term。')
        diffenrich_section2.add_comment("根据功能分级，一般将 GO 分为三个层级，level1 包含三个条目：biological process、cellular component和molecular function，level2 包含 biological adhesion、cell 和 binding 等 64 个条目，level3 即为常规富集使用的数万个条目。从 level1到 level3 功能更具体，反之，更概括。")
        #diffenrich_section2.add_comment("差异基因和所有基因在 GO Level2 水平分布比较图如下：")
        #diffenrich_section2.add_plot("%s/15.enrichment/GO_enrichment/*/ALL_vs_DEG.GO.level2.stat.png"%(outdir),
        #                         caption='差异表达基因及所有基因在 GO Level2 水平分布比较图',
        #                         description='图片说明：蓝色表示所有基因富集的 GO Level2 条目，红色表示差异基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')
        diffenrich_section2.add_comment("上调差异基因和下调差异基因在 GO Level2 水平分布比较图如下：")
        diffenrich_section2.add_plot("%s/15.enrichment/GO_enrichment/*/GO.level2.stat.Total.png"%(outdir),
                                     caption='上调差异基因和下调差异基因在 GO Level2 水平分布比较图',
                                     description='图片说明：红色表示上调差异表达基因富集的 GO Level2 条目，绿色表示下调差异表达基因富集的 GO Level2 条目，横轴为条目名称，纵轴表示对应条目的基因数量和其百分比。')

    ## KEGG 富集分析
    if os.path.exists("%s/15.enrichment/KEGG_enrichment" % (outdir)):
        diffenrich_section3 = diffenrich_section1.add_section("差异基因KEGG富集分析")
        diffenrich_section3.add_comment("KEGG是有关 Pathway 的主要公共数据库，利用KEGG数据库对差异蛋白编码基因进行 Pathway 分析（结合 KEGG 注释结果），并用超几何分布检验的方法计算每个 Pathway 条目中差异基因富集的显著性。")
        diffenrich_section3.add_comment("差异基因KEGG富集分析结果文件如下：[15.enrichment/KEGG_enrichment](./15.enrichment/KEGG_enrichment)")
        diffenrich_section3.add_comment("KEGG 富集分析 top20（筛选对应差异基因数目大于 2 的 Pathway 条目，按照每个条目对应的 -log10Pvalue 由大到小排序）气泡图如下：")
        diffenrich_section3.add_plot("%s/15.enrichment/KEGG_enrichment/*/KEGG.top.*.png"%(outdir),
                                     caption='KEGG富集 top20 气泡图',
                                     description='图片说明：图中横轴 Enrichment Score 为富集分值，气泡越大的条目包含的差异蛋白编码基因数目越多，气泡颜色由紫-蓝-绿-红变化，其富集 pValue 值越小，显著程度越大。')
        diffenrich_section3.add_comment("根据功能分级，通常将 KEGG 分为三个层级，level1 包含六个分类：Metabolism、Genetic Information Processing、Environmental Information Processing、Cellular Processes、Organismal Systems 和 Human Diseases（具体物种注释可能有删减）。level2 包含 Cell growth and death、Transcription 和 Development 等 44 个分类（具体物种注释可能有删减），level3 即为常规富集使用的数百个 Pathway，从 level1 到 level3 功能更具体，反之，更概括。")
        diffenrich_section3.add_comment("差异表达基因及所有基因在 KEGG Level2 水平分布比较图如下：")
        diffenrich_section3.add_plot("%s/15.enrichment/KEGG_enrichment/*/ALL_vs_DEG.KEGG_Classification.png"%(outdir),
                                     caption='差异表达基因及所有基因在 KEGG Level2 水平分布比较图',
                                     description='图片说明：横轴是注释到各 Level2 通路的基因（差异表达基因）和所有注释到 KEGG 通路的基因（差异表达基因）总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该Level2 Pathway下的差异表达基因数量。')
        diffenrich_section3.add_comment("上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图如下：")
        diffenrich_section3.add_plot("%s/15.enrichment/KEGG_enrichment/*/Up_vs_Down.KEGG_Classification.png"%(outdir),
                                     caption='上调差异表达基因及下调差异表达基因在 KEGG Level2 水平分布图',
                                     description='图片说明：横轴是注释到各 Level2 通路的上调（下调）差异表达基因和所有注释到 KEGG 通路的上调（下调）差异表达基因总数的比值（%），纵轴表示 Level2 Pathway 的名称，柱子右边数字代表注释到该 Level2 Pathway 的上调（下调）差异表达基因数量。')

    ##========================================================= 4.项目摘要、链接========================================
    Project_module.add_section(name="项目摘要",description=''.join(report_summary_list))
    report_hyperlink_list.append('[10x 单细胞多组学（ATAC + 基因表达）常见问题（FAQ）](#faq)\n\n[Loupe-Cell-Browser使用教程](supplemental_material/Loupe-Cell-Browser使用教程.html)')
    report_hyperlink = ''.join(report_hyperlink_list)
    Project_module.add_section(name="项目快捷链接",description=report_hyperlink)
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
        [['Cell Ranger ARC', '2.0.0'],
         ['Seurat', '4.0.1'],
         ['SeuratDisk', '0.0.0.9019'],
         ['Signac','1.7.0'],
         ['harmony','1.0'],
         ['batchelor','1.6.3']
          ])
    database2.columns = ['软件', '版本']
    description_section6 = description_module.add_section('数据分析软件')
    software_table = description_section6.add_table(database2, caption='数据分析软件')

    ##============================================================== 6.申明 ===========================================
    affirming_module = report.add_module('申明')

    ############################################ Generate Report HTML ###################################################
    report.write_to('%s/Report.html'%(outdir))
if __name__ == "__main__":
    multimodal_report()
