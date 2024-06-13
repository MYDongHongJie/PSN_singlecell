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


def task_transrds(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    if not os.path.exists( f"{wkdir}/output/download"):
        os.makedirs( f"{wkdir}/output/download")
    if not os.path.exists( f"{wkdir}/transrds"):
        os.makedirs( f"{wkdir}/transrds")
    #input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{wkdir}/input/cell.tsv"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"
    outdir=f"{wkdir}/transrds"
    # environment=input.loc["project","environment"]
    environment="OESingleCell/2.0.0"
    inTaxid=int(input.loc["parameters","inTaxid"])
    outTaxid=int(input.loc["parameters","outTaxid"])
    cmd=f"Rscript /opt/softwares/oe_packages/script/homologene_transformed.R "\
        f"  -v {seurat_ob} "\
        f"  --assay RNA "\
        f"  -o {outdir} "\
        f"  -i {inTaxid} "\
        f"  -t {outTaxid}  "\
        f"  {add_params} "\
        # f"  --topn {topn}  "\
        # f"  -topby {topby} "\

    logger.info(cmd)
    logger.info("开始运行")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    cmd_ln = f"ln -s {wkdir}/transrds/homologene.*.xls {wkdir}/output/download/" +"\n"+\
        f"mv {wkdir}/transrds/*.rds {wkdir}/transrds/data_ob_v3.rds"+"\n"


    logger.info("完成")
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status ###run_cmd use the default function from xiufeng.yang
