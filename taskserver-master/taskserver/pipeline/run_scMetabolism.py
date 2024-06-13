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


def task_scMetabolism(input, output_cfg, projectid, taskid, workdir):
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
    seurat_ob = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    contrast = input.loc['parameters', 'diff']
    collapseby = input.loc['parameters', 'collapseby']
    topn = input.loc['parameters', 'topn']
    pvalue = input.loc['parameters', 'pvalue']
    type = input.loc['parameters', 'type']
    method = input.loc['parameters', 'method']
    #gmt = input.loc['parameters', 'gmt']
    cellmeta = input.loc['base', 'tasks'][0]['filename']

    add_params = ""
    if cellmeta != "":
        add_params = f" --cloudmeta  {cellmeta} "
    #add_gmt = ""
    #if gmt != "":
    #    add_gmt = f" --gmt  {gmt} "


    ## 2.获取命令行
    cmd1 = f"/home/lipeng/miniconda3/envs/Seurat3.1/bin/Rscript /public/scRNA_works/works/jhyu/test/boxplot/cloud_scMET/bak/scmetabolism_scFEA.v1.7_test.R " \
        f"-v {seurat_ob} " \
        f"-o {workdir}/{module} " \
        f"-q {collapseby} " \
        f"-n {topn} " \
        f"-t {type} " \
        f"-p {pvalue} " \
        f"--cpu 6 " \
        f"-d {contrast} " \
        f"-m {method} "
    ## 3.执行分析
    with module_cmd(env_module) as p:
        status = p(cmd1, projectid, taskid)

    # =================================================================================================================
    logger.info("step2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
        os.makedirs(f'{workdir}/output/download/1.Diff_Term')
    cmd = f"ln -s {workdir}/{module}/{{average_*.tsv,metabolism_*.tsv}} {workdir}/output/download && ln -s {workdir}/{module}/*.png {workdir}/output/download && ln -s {workdir}/{module}/1.Diff_Term/{{*.tsv,*.png}} {workdir}/output/download/1.Diff_Term"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    # logger.info("step3:生成output.json.xls")
    #output_cfg = f"{workdir}/output/download/output.json.tsv"
    #output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    return status
