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


def task_visualize_marker(input,projectid="项目ID", taskid="任务ID", workdir="运行目录"):
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
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    if input.loc['parameters','var2use'] != "" and input.loc['parameters','levels4var'] != "":
        var2use = input.loc["parameters","var2use"]
        levels4var = input.loc["parameters","levels4var"]
        add_params = add_params + f" --var2use  {var2use} --levels4var {levels4var}"

    environment="OESingleCell/2.0.0"
    ##解析传入基因列表
    if input.loc['parameters','plot_genes_file'] == "":
        topn=input.loc["parameters","topn"]
        base_taskid = input.loc["base", "tasks"][0]["taskId"]
        base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
        base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
        base_module_name = base_json.loc['task', 'type']
        if base_module_name == "findallmarkers":
            rds = f"{workdir}/../../{base_json.loc['base', 'tasks'][0]['projectId']}/{base_json.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
            base_file=f"{base_wd}/{base_module_name}/all_markers_for_each_cluster_anno.tsv"
            file_param = f"-l {base_file}  -n  {topn}  -c  gene_diff "
        # elif base_module_name == "diffexp":
            # base_file=f"{base_wd}/{base_module_name}/group_*-vs-*-diff-*-FC-*_anno.tsv"
            # file_param = f"-d {base_file}"
        else:
            print("未找到合适的基因列表，请核对上游依赖任务或上传基因列表。")
    else:
        rds = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
        file_param = f" -x {workdir}/input/{input.loc['parameters','plot_genes_file']}"

    groupby=input.loc["parameters","groupby"]
    reduct=input.loc["parameters","reduct"]
    assay=input.loc["parameters","assay"]

    method=",".join(input.loc["parameters", "method"])
    ## 2.获取命令行
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " /opt/softwares/envs/oesinglecell_v3/bin/Rscript /visualize_markers_script/visualize_markers.R "
    else:
        tmp_envs="module purge && module load OESingleCell/2.0.0 && Rscript /public/scRNA_works/works/guokaiqi/test/sc_cloud/visualize_markers/visualize_markers.R"
 
    cmd=f"{tmp_envs}"\
        f" -v {rds} "\
        f" {file_param} " \
        f" -m {method} " \
        f" --reduct {reduct} " \
        f" -e RNA" \
        f" -o {workdir}/{module} " \
        f" -g {groupby} " \
        f" {add_params} " 
        
    logger.info(cmd)
    logger.info("step2:visualize_marker 分析")
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)
    # =================================================================================================================
    logger.info("step3:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd1= f"ls {workdir}/{module}/ | grep -v tsv$ | xargs -i ln -s {workdir}/{module}/{{}} {workdir}/output/download/"
    with module_cmd(environment) as p:
        status = p(cmd1, projectid, taskid)
    return status




