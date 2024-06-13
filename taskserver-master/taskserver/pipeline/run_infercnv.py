import pandas as pd
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from glob import glob
from oebio.utils.log import getLogger
logger = getLogger('oe.cloud.sc.qsub')


from taskserver.tools.module_cmd import module_cmd


def task_infercnv(input,projectid="项目ID", taskid="任务ID", workdir="运行目录"):
    wkdir = workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/infercnv")
    cellmeta=input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {wkdir}/input/{cellmeta} "
    else:
        pass

    ####获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files
    seurat_ob = search_files(f'{base_wd}/{base_module_name}', "rds")[0]
    
    outdir=f"{wkdir}/infercnv"

    #infercnv分析参数
    gtf = f"{input.loc['parameters', 'database']}/genes/genes.gtf"
    celltype=input.loc["parameters","celltype"]
    ref_group_name = input.loc["parameters", "ref_group_name"]
    assay=input.loc["parameters","assay"]
    cutoff=input.loc["parameters","cutoff"]
    mode=input.loc["parameters","analysis_mode"]
    cpu=input.loc["parameters","ncores"]
    hclust_method=input.loc["parameters","hclust_method"]
    Pval=input.loc["parameters","Pval"]

    if input.loc["parameters", "DoHMM"] :
        add_params = add_params + f"--doHMM TRUE"
    else:
        pass

    if input.loc["parameters", "denoise"] :
        pass
    else:
        add_params = add_params + f"--denoise FALSE"

    #script_dir="."
    environment="OESingleCell/2.0.0"

    cmd=f"Rscript /home/xuejuan/Script/infercnv/infercnv.R "\
        f" -i {seurat_ob} "\
        f" -f seurat " \
        f" --assay {assay} " \
        f" -l {celltype} " \
        f" -r {ref_group_name} " \
        f" -z {gtf} " \
        f" -o {outdir} " \
        f" -u {mode}" \
        f" --clusting2use {hclust_method} " \
        f" --cutoff {cutoff}" \
        f" --pval {Pval}" \
        f" --ncores {cpu}" \
        f" {add_params} "

    logger.info(cmd)
    logger.info("infercnv:执行分析")
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    ## 2.链接rds
    logger.info("链接rds")
    if not os.path.exists(f'{wkdir}/infercnv/data_ob_v3.rds'):
        cmd2= f"ln -s {seurat_ob}  {wkdir}/infercnv/data_ob_v3.rds"
        with module_cmd(environment) as p:
            status = p(cmd2, projectid, taskid)
        return status
    else:
        logger.info("存在data_ob_v3.rds")

    ## 3.执行分析
    logger.info("构建output/download/分析结果")
    if not os.path.exists(f'{wkdir}/output/download'):
        os.makedirs(f'{wkdir}/output/download')
    cmd3= f"ln -s {wkdir}/infercnv/*.pdf  {wkdir}/output/download/;ln -s {wkdir}/infercnv/*png {wkdir}/output/;cp {wkdir}/infercnv/clustering_results.tsv {wkdir}/output/"
    with module_cmd(environment) as p:
        status = p(cmd3, projectid, taskid)
    return status
