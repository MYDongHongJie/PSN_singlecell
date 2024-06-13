# -*-coding:utf-8-*-
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

#===================change all the input into standard str form===========================================================================
def params_setting(d):
    params = []
    for key in d:
        param = f'{key} {d[key]}'
        params.append(param)
    return params
#===================input DiffusionMap input and run the cmds===================================================================================
def task_DiffusionMap(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    if not os.path.exists( f"{wkdir}/output/download"):
        os.makedirs( f"{wkdir}/output/download")
    if not os.path.exists( f"{wkdir}/DiffusionMap"):
        os.makedirs( f"{wkdir}/DiffusionMap")
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    root1 = input.loc["parameters","root1"]
    root2 = str(input.loc["parameters","root2"])
    add_params=""
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"
    outdir=f"{wkdir}/DiffusionMap"
    environment=input.loc["project","environment"]
    # environment = "OESingleCell/2.0.0"
    groupby=input.loc["parameters","groupby"]
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    if root1 != "" and root2 != "":
        add_params = add_params + f" -r {root1}:{root2} "

    cmd=f"module purge && source /home/lipeng/miniconda3/bin/activate ggstatsplot && "\
        f"Rscript /public/scRNA_works/works/chenhaoruo/script/scrna_cloud/DiffusionMap/diffusion_map_v0.3.R "\
        f"  -i {seurat_ob} "\
        f"  -o {outdir} "\
        f"  --assay RNA "\
        f"  -c {groupby} "\
        f"  {add_params} "\

    logger.info(cmd)
    logger.info("开始运行")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    d = {"task_type":["DiffusionMap"],"result_module":["data"], 
    "input":["Diffusion_reduct_metadata.tsv"], 
    "type":["tsv"], 
    "file":["Diffusion_reduct_metadata.tsv"],
    "title":["DiffusionMap"], 
    "downloadName":["Diffusion_reduct_metadata"], 
    "downloadPath":["download/Diffusion_reduct_metadata.tsv"]}
    df = pd.DataFrame(d)
    df.to_csv(f"{wkdir}/DiffusionMap/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")

    cmd_ln = f"ln -s {wkdir}/DiffusionMap/* {wkdir}/output/download/"  +"\n"+\
        f"rm -f {wkdir}/output/download/*.rds"+"\n"+\
        f"mv {wkdir}/DiffusionMap/*.rds {wkdir}/DiffusionMap/data_ob_v3.rds"+"\n"

    logger.info("完成")
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status ###run_cmd use the default function from xiufeng.yang
    
