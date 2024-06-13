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


def task_scheatmap(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    if not os.path.exists( f"{wkdir}/output/download"):
        os.makedirs( f"{wkdir}/output/download")
    if not os.path.exists( f"{wkdir}/scheatmap"):
        os.makedirs( f"{wkdir}/scheatmap")
    #input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{wkdir}/input/cell.tsv"
    genelist = f"{wkdir}/input/{input.loc['parameters', 'genelist']}"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"
    outdir=f"{wkdir}/scheatmap"
    environment=input.loc["project","environment"]
    groupby=input.loc["parameters","groupby"]
    # topn=int(input.loc["parameters","topn"])
    rowcluster=input.loc["parameters","rowcluster"]
    colcluster=input.loc["parameters","colcluster"]
    gaps_row=str(input.loc["parameters","gaps_row"])
    gaps_col=str(input.loc["parameters","gaps_col"])
    showname=input.loc["parameters","showname"]
    # topby=input.loc["parameters","topby"]

    cmd=f"Rscript /opt/softwares/oe_packages/script/sc_heatmap.r "\
        f"  -i {seurat_ob} "\
        f"  -g {genelist} "\
        f"  -o {outdir} "\
        f"  -c {groupby} "\
        f"  -r {rowcluster}  "\
        f"  -l {colcluster}  "\
        f"  -w {gaps_row} "\
        f"  -L {gaps_col} "\
        f"  --showname {showname} "\
        f"  {add_params} "\
        # f"  --topn {topn}  "\
        # f"  -topby {topby} "\

    logger.info(cmd)
    logger.info("开始绘制均值热图")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    cmd_ln = f"ln -s {wkdir}/scheatmap/heatmap.png {wkdir}/output/download/heatmap.png" +"\n"+\
        f"ln -s {wkdir}/scheatmap/heatmap.pdf {wkdir}/output/download/heatmap.pdf" +"\n"+\
        f"ln -s {wkdir}/scheatmap/heatmap_count.xls {wkdir}/output/download/heatmap_count.xls" +"\n"+\
        f"ln -s {wkdir}/scheatmap/heatmap_count_scaled.xls {wkdir}/output/download/heatmap_count_scaled.xls" +"\n"+\
        f"if [ -e {wkdir}/scheatmap/filtered_gene.xls ];then ln -s {wkdir}/scheatmap/filtered_gene.xls {wkdir}/output/download/filtered_gene.xls;fi"+"\n"

    logger.info("完成")
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status ###run_cmd use the default function from xiufeng.yang
