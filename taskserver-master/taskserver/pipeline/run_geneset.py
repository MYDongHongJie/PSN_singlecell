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





def task_geneset(input,projectid="项目ID", taskid="任务ID", workdir="运行目录"):
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
    logger.info("step1:解析命令行")
    ## 1.获取filterd.h5seurat路径及参数信息
    rds = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    environment="OESingleCell/2.0.0"
    genelist=input.loc["parameters","genelist"]
    groupby=input.loc["parameters","groupby"]
    pointsize=input.loc["parameters","pointsize"]
    reduct=input.loc["parameters","reduct"]
    pvalue=input.loc["parameters","pvalue"]
    assay=input.loc["parameters","assay"]
    method=input.loc["parameters","method"] 
    if method == "geneset":
    
        ## 2.获取命令行
        cmd=f"Rscript /public/scRNA_works/works/jhyu/test/boxplot/cloud_addmodulescore/cloud/visualize_markers.R"\
            f" -v {rds} "\
            f" -x {workdir}/input/{genelist} " \
            f" -m {method} " \
            f" -s {pointsize} " \
            f" --reduct {reduct} " \
            f" -e {assay}" \
            f" -o {workdir}/{module} " \
            f" -g {groupby} " \
            f" -p {pvalue}"
        logger.info(cmd)
        logger.info("step2:geneset 分析")
        with module_cmd(environment) as p:
            status=p(cmd, projectid, taskid)
        # =================================================================================================================
        logger.info("step3:构建output/download/分析结果")
        if not os.path.exists(f'{workdir}/output/download'):
            os.makedirs(f'{workdir}/output/download')
        cmd1= f"ln -s {workdir}/{module}/geneset_visualization/{{*png,*tsv}}  {workdir}/output/download/"
        with module_cmd(environment) as p:
            status = p(cmd1, projectid, taskid)
        return status
    if method == "ggstatsplot":
        ## 2.获取命令行
        cmd=f"Rscript /public/scRNA_works/works/jhyu/test/boxplot/cloud_addmodulescore/cloud/visualize_markers.R"\
            f" -v {rds} "\
            f" -x {workdir}/input/{genelist} " \
            f" -m {method} " \
            f" -e {assay}" \
            f" -o {workdir}/{module} " \
            f" -g {groupby} " \
            f" -p {pvalue}"
        
        logger.info(cmd)
        logger.info("step2:ggstatsplot 分析")
        with module_cmd(environment) as p:
            status=p(cmd, projectid, taskid)
        # =================================================================================================================
        logger.info("step3:构建output/download/分析结果")
        if not os.path.exists(f'{workdir}/output/download'):
            os.makedirs(f'{workdir}/output/download')
        cmd1= f"ln -s {workdir}/{module}/ggstatsplot_visualization/{{*png,ggstatsplot_plot.tsv}}  {workdir}/output/download/"
        with module_cmd(environment) as p:
            status = p(cmd1, projectid, taskid)
        return status
