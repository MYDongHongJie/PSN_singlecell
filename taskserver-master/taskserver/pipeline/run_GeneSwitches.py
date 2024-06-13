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

def task_GeneSwitches(input, projectid, taskid, workdir):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/GeneSwitches")

    #拟时序的下游分析，所以文件读取路径在拟时序结果文件夹中
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/*seurat_withpseudotime.rds"
    monocle_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/pseudotime_results.rds"
    outdir=f"{wkdir}/GeneSwitches"

    environment=input.loc["project","environment"]
    bn_cutoff=input.loc["parameters","bn_cutoff"]
    reference=input.loc["parameters","reference"]
    topn=int(input.loc["parameters","topn"])

    add_params = ""
    if input.loc["parameters","states"] == "" :
        pass
    else:
        add_params = add_params + f" --states {states} "

    if input.loc["parameters", "showtf"] :
        add_params = add_params + f"--showtf TRUE"
    else:
        add_params = add_params + f"--showtf FALSE"

    environment = "OESingleCell/2.0.0"
    cmd = f"Rscript /home/xuejuan/Script/GeneSwitches/GeneSwitches.R" \
          f" -i {seurat_ob}"  \
          f" -m {monocle_ob}"  \
          f" --reference {reference}" \
          f" --bn_cutoff {bn_cutoff}" \
          f" --topn {topn}" \
          f" -o {outdir}" \
          f" {add_params} "

    logger.info(cmd)
    logger.info("step1: 开始GeneSwitches分析")
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)

#################################################################################
    logger.info("step2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd2 = f"ln -s {workdir}/GeneSwitches/{{*.pdf,*png,*xls,*txt}} {workdir}/output/download/ ;" \
           f"ln -s {workdir}/GeneSwitches/enrichment/ {workdir}/output/download/"
    with module_cmd(environment) as p:
        status = p(cmd2, projectid, taskid)
    return status
