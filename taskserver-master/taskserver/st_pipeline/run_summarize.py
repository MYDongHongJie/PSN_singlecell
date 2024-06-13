import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_summarize(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
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
    # env_module = "OESingleCell/v_3.0.0_visium"
    env_module = "OESingleCell/v_3.0.0_visium_produce"
    input_rds = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    assay = input.loc['parameters', 'assay']
    reduct = input.loc['parameters', 'reduct']
    library_size = input.loc['parameters','library_size']
    resolution = input.loc['parameters', 'resolution']
    ptsize = input.loc['parameters', 'ptsize']
    spointsize = input.loc['parameters', 'spointsize']
    groupby = input.loc['parameters','groupby']
    propby = input.loc['parameters','propby']
    splitby = input.loc['parameters','splitby']
    groups = input.loc['parameters','groups']
    if library_size == "cytassist":
        crop = "--crop TRUE"
    else:
        crop = " "
    #print(propby)
    ## 2.获取命令行
    cmd = f"/public/scRNA_works/works/guochy/ST_taskserver/exec/scVis " \
          f"-i {input_rds} " \
          f"-f rds " \
          f"--reduct {assay}_{reduct} " \
          f"-o {workdir}/{module} " \
          f"--assay {assay} " \
          f"--image TRUE " \
          f"dimplot " \
          f"      --resolution {resolution}" \
          f"      --ptsize {ptsize} " \
          f"      --spointsize {spointsize} " \
          f"      --groupby {groupby} " \
          f"      {crop}" \
          f"      --splitby {splitby} " \
          f"      --groups {groups} " \
          f"      --common.legend FALSE " \
          f"      --propby {propby} "

    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    ## 3.执行分析
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd= f"ln -s {workdir}/{module}/visualize_cluster_by_*  {workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    return status
