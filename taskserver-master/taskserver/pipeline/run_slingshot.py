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


def task_slingshot(input, projectid, taskid, workdir):
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
    logger.info("step3:解析命令行，执行分析")
    ## 1.获取filterd.h5seurat路径及参数信息
    # env_module = input.loc['project', 'environment']
    env_module = "OESingleCell/3.0.d"
    rds = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    groupby = input.loc['parameters', 'groupby']
    reduct = input.loc['parameters', 'reduct']
    # seurat = input.loc['parameters', 'seurat']
    # start = str(input.loc['parameters', 'start'])
    # end = str(input.loc['parameters', 'end'])
    cellmeta = input.loc['base', 'tasks'][0]['filename']

    add_params = ""
    # if var2use != "" and levels4var != "":
    #     add_params = f" --predicate {var2use} %in%  {levels4var}"
    if input.loc['parameters', 'seurat']:
        add_params = add_params + f" -b TRUE"
    if input.loc['parameters', 'start'] != "":
        add_params = add_params + f" -a {str(input.loc['parameters', 'start'])}"
    if input.loc['parameters', 'end'] != "":
        add_params = add_params + f" -z {str(input.loc['parameters', 'end'])}"

    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "

    logger.info("step3.1:开始分析")
    ## 2.获取命令行
    cmd1 = f"Rscript /public/scRNA_works/works/guokaiqi/test/sc_cloud/slingshot.R " \
        f"-i {rds} " \
        f"-o {workdir}/{module} " \
        f"-g {groupby} " \
        f"--reduct {reduct} " \
        f"{add_params} " 
        ## 3.执行分析
    with module_cmd(env_module) as p:
        status = p(cmd1, projectid, taskid)

    # =================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd= f"ln -s {workdir}/{module}/{{*.pdf,*png,*docx}}  {workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    # logger.info("step3.3:生成output.json.xls")
    return status
