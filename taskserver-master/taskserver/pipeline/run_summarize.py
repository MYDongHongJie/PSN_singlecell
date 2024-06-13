import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd
from glob import glob
logger = getLogger('oe.cloud.sc.qsub')

def task_summarize(input,  projectid, taskid, workdir):
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
    env_module=input.loc['project', 'environment']
    #bclust_h5seurat = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/bclust/bclust.h5seurat"
    ####进行输入对象rds或者h5seurat格式的识别和判断
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files

    rds_input = search_files(f'{base_wd}/{base_module_name}', "rds")
    h5seurat_input = search_files(f'{base_wd}/{base_module_name}', "h5seurat")
    seurat_ob = ""
    format_ob = ""
    if (len(rds_input) == 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    elif (len(rds_input) != 0 and len(h5seurat_input) == 0):
        seurat_ob = rds_input[0]
        format_ob = "rds"
    elif (len(rds_input) != 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    else:
        print(f"输入seurat对象为空，请检查！！！")
    #seurat_ob_realpath=os.readlink(seurat_ob)
    ####以上完成输入对象的识别及赋值

    # base_taskid = input.loc["base", "tasks"][0]["taskId"]
    # base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    # base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    # reduct=base_json.loc['parameters','reduct2']

    groupby = input.loc['parameters','summ_groupby']
    #propby = input.loc['parameters','propby'].replace("[","").replace("]","").replace("\"","")
    #propby =','.join([str(elem) for elem in input.loc['parameters','propby']])
    propby = input.loc['parameters','propby']
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta {workdir}/input/{cellmeta} "
    #print(propby)
    ## 2.获取命令行
    cmd = f"sctool " \
          f"-i {seurat_ob} " \
          f"-f h5seurat " \
          f"-o {workdir}/{module} " \
          f"--assay RNA " \
          f"{add_params} " \
          f"summarize " \
          f"      --method barplot " \
          f"      --groupby {groupby} " \
          f"      --propby {propby} "
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    ## 3.执行分析
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd= f"ln -s {workdir}/{module}/{{*.pdf,*png,*xls}}  {workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    return status
