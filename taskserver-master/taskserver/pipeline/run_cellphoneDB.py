import pandas as pd 
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger
logger = getLogger('oe.cloud.sc.qsub')


from taskserver.tools.module_cmd import module_cmd

def task_cellphoneDB(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    wkdir =workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/cellphoneDB")
    # #获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/*.rds"
    #seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/bclust/*.rds"

    celltype = f"--celltype {input.loc['parameters', 'celltype']} " if input.loc['parameters', 'celltype'] !="" else ""
    var2use = f"--var2use {input.loc['parameters', 'var2use']} " if input.loc['parameters', 'var2use'] !="" else ""
    levels4var = f"--levels4var {input.loc['parameters', 'levels4var']} " if input.loc['parameters', 'levels4var'] !="" else ""
    module = input.loc['task', 'type']
    outdir = f"{wkdir}/{module}"
    environment = input.loc["project","environment"]
    assay = input.loc["parameters","assay"]
    species = input.loc["parameters","species"]
    column4cell = input.loc["parameters","column4cell"]
    plot = f"--plot {input.loc['parameters','plot']} " if input.loc['parameters', 'plot'] !="" else f" --plot network,circos,dotplot,chorddiagram,heatmap,barplot"  
    topn = f"--topn {input.loc['parameters','topn']} " if input.loc['parameters', 'topn'] !="" else ""  


    cmd=f"module purge && module load OESingleCell/2.0.0  && "\
        f"Rscript /public/scRNA_works/works/liuhongyan/Test/oecloude/RunCellphonedb.R "\
        f"  -i {seurat_ob} "\
        f"  -f seurat " \
        f"  -s {species} "\
        f"  -o {outdir} "\
        f"  --assay {assay} "\
        f"  -p 10 "\
        f"  -c {column4cell}" \
        f" {celltype} {var2use} {levels4var}"

    logger.info(cmd)
    logger.info("step1: 开始 cellphoneDB 分析")
    ### env ####
    ret = subprocess.run(cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        encoding="utf-8")
    logger.info(ret.stdout)  
    logger.info(ret.stderr)  


    seurat_ob=f"{wkdir}/{module}/cellphonedb_results.rds"
    cmd=f"module purge && module load OESingleCell/2.0.0  && "\
            f"Rscript /public/scRNA_works/works/liuhongyan/Test/oecloude/visualize_cellcomm.R "\
            f"  -i {seurat_ob} "\
            f"  -f rds " \
            f"  -o {outdir} "\
            f" {plot} {topn}" 

    logger.info(cmd)
    logger.info("step2: 开始 cellphoneDB visualization 分析")
    ### env ####
    ret = subprocess.run(cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        encoding="utf-8")
    logger.info(ret.stdout)  
    logger.info(ret.stderr)  

# ==================================================================================================================


    logger.info("step2: 构建output/download/分析结果")
    ## change xls to tsv
    os.rename(f"{workdir}/{module}/cell_comm_annotation.xls",
                f"{workdir}/{module}/cell_comm_annotation.tsv")
    os.rename(f"{workdir}/{module}/cell_comm_summary.xls",
                f"{workdir}/{module}/cell_comm_summary.tsv")

    if not os.path.exists(f'{workdir}/output/download/'):
        os.makedirs(f'{workdir}/output/download/')
    cmd = f"ln -s {workdir}/{module}/*.png {workdir}/{module}/*.pdf {workdir}/{module}/*.tsv {workdir}/{module}/cellphoneDB* {workdir}/output/download/"
    os.system(cmd)
    
    ## 图形 tsv ### 
   # ==================================================================================================================

    ## =================================================
    ## 2. amending output.json
    ## =================================================
    logger.info("step2 生成output.json.tsv") 

    df = pd.DataFrame(columns=['task_type', 'result_module', 'input', 'type', 'file', 'title','downloadName', 'downloadPath'])
    df.loc[len(df.index)] = ["cellphoneDB","data","cell_comm_annotation.tsv","tsv","cell_comm_annotation.tsv","communication","communication","download/cell_comm_annotation.tsv"]
    df.loc[len(df.index)] = ["cellphoneDB","summary","cell_comm_summary.tsv","tsv","cell_comm_summary.tsv","communication","communication","download/cell_comm_annotation.tsv"]
    if "dotplot" in plot.split(",") :
        df.loc[len(df.index)] = ["cellphoneDB","diagram","cell_comm_dotplot.pljson","oe","cell_comm_dotplot.pljson","cell_comm_dotplot","",""] 

    if "barplot" in plot.split(",") :
        df.loc[len(df.index)] = ["cellphoneDB","diagram","cell_comm_histogram_plot.pljson","oe","cell_comm_histogram_plot.pljson","cell_comm_histogram_plot","",""]
   
    df.to_csv(f"{workdir}/{module}/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')

