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


def task_infercnv_vis(input,projectid="项目ID", taskid="任务ID", workdir="运行目录"):
    wkdir = workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/infercnv_vis")

    #获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"

    outdir=f"{wkdir}/infercnv_vis"

    #infercnv绘图参数
    reduct=input.loc["parameters","reduct"]
    vismethod=','.join(str(x) for x in list(set(input.loc["parameters","vismethod"])))
    groupby=input.loc["parameters","groupby"]
    splitby=input.loc["parameters","splitby"]
    colormapping=input.loc["parameters","colormapping"]

    #infercnv_vis rds路径
    if input.loc["parameters","infercnv_output_path"] == "" :
        infercnv_output_path=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/infercnv/"
    else:
        infercnv_output_path=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{int(input.loc['parameters','infercnv_output_path'])}/infercnv/"

    #script_dir="."
    environment="OESingleCell/2.0.0"
    # cmd=f"source /opt/softwares/oe_packages/script/bashrc.sh && "\
    cmd=f"Rscript /opt/softwares/oe_packages/script/infercnv_vis.R" \
        f" -i {seurat_ob}" \
        f" -f seurat " \
        f" -l {infercnv_output_path} " \
        f" -g {groupby} " \
        f" -m {colormapping}" \
        f" -p {vismethod} " \
        f" --reduct {reduct} " \
        f" -y {splitby}"  \
        f" -o {outdir}"
    logger.info(cmd)
    logger.info("infercnv:绘图分析")
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    #执行分析
    logger.info("构建output/download/分析结果")
    if not os.path.exists(f'{wkdir}/output/download'):
        os.makedirs(f'{wkdir}/output/download')
    cmd2= f"ln -s {wkdir}/infercnv_vis/{{*.pdf,*png,*xls}}  {wkdir}/output/download/; ln -s {wkdir}/infercnv_vis/*png {wkdir}/output/ ; cp {wkdir}/infercnv_vis/cnv_result.xls {wkdir}/output/cnv_result.tsv"
    with module_cmd(environment) as p:
        status = p(cmd2, projectid, taskid)
    return status
