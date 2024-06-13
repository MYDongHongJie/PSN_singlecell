import os
import sys
import re
import json
import subprocess
import pandas as pd
import requests
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')


def task_Scorpius(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # ==================================================================================================================
    logger.info("step1:解析命令行，执行分析")
    ## 1.获取seurat_ob路径及参数信息
    env_module = input.loc['project', 'environment']
    seurat_ob = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/bclust/data_ob_v3.rds"
    colorby = input.loc['parameters', 'colorby']
    reduct = input.loc['parameters', 'reduct']
    assay = input.loc['parameters', 'assay']
    qvalue = input.loc['parameters', 'qvalue']
    extraGene = input.loc['parameters', 'extraGene']
    reverse = input.loc['parameters', 'reverse']
    #groupby = input.loc['parameters', 'groupby']
    #which_cells = input.loc['parameters', 'which_cells']
    #gmt = input.loc['parameters', 'gmt']
    cellmeta = input.loc['base', 'tasks'][0]['filename']

    add_params = ""
    if cellmeta != "":
        add_params = f" --cloudmeta  {cellmeta} "
    #add_gmt = ""
    #if gmt != "":
    #    add_gmt = f" --gmt  {gmt} "


    ## 2.获取命令行
    cmd1 = f"module purge && source /home/lipeng/miniconda3/bin/activate base && Rscript /public/scRNA_works/works/lipeng/script/SCORPIUS/cloud/Scorpius_v1.4.R " \
        f"-i {seurat_ob} " \
        f"-o {workdir}/{module} " \
        f"-q {qvalue} " \
        f"-a {assay} " \
        f"-c {colorby} " \
        f"-x {extraGene} " \
        f"-r {reduct} " \
        f"--reverse {reverse} " \
        #f"-u {which_cells} " \
        #f"-g {groupby} "
    ## 3.执行分析
    with module_cmd(env_module) as p:
        status = p(cmd1, projectid, taskid)

    # =================================================================================================================
    logger.info("step2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd = f"ln -s {workdir}/{module}/*.tsv {workdir}/output/ && ln -s {workdir}/{module}/Scorpius_*xls {workdir}/output/download && ln -s {workdir}/{module}/Scorpius_*png {workdir}/output/download && ln -s {workdir}/{module}/Scorpius_*pdf {workdir}/output/download "
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    # logger.info("step3:生成output.json.xls")
    #output_cfg = f"{workdir}/scMetabolism/output.json.tsv"
    #output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    return status
