# -*-coding:utf-8-*-
import glob
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

def find_loom(wkdir,input):
    logger.info("寻找loom文件")
    # os.makedirs(f"{wkdir}/velocity/loom")
    loom = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/loom/loom"
    if not os.path.exists(loom):
        logger.info("上级任务不是loom，在项目目录寻找loom")
        loom = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/loom"
    check_loom = glob.glob(os.path.join(f'{loom}', '*.loom'))
    if check_loom:
        logger.info("loom文件存在，开始处理")
    else:
        logger.info("loom文件不存在，请检查loom任务")
        sys.exit(1)
    return loom

def task_velocity(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir  = workdir
    environment = "OESingleCell/2.0.0" #input.loc['project', 'environment']
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/velocity")
    cellmeta=input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" -m  {wkdir}/input/{cellmeta} "
  ### velocity
    #读取rds路径
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"
    if not os.path.exists(seurat_ob):
        logger.info("RDS不存在，请从降维聚类开始进行loom和velocity分析")
    if input.loc['parameters', 'pointsize'] != "" :
        pointsize = input.loc["parameters", "pointsize"]
    else :
         pointsize = 0.5
    if input.loc["parameters", "ident2use"] != "" and  input.loc["parameters", "which_cells"] != "" :
        ident2use = input.loc["parameters", "ident2use"]
        which_cells = input.loc["parameters", "which_cells"]
        add_params = add_params + f"-q  {which_cells}  -u  {ident2use} "     
    else:
        logger.info("ident2use and which_cells are supporting parameters, which must be provided at the same time")
    reduct_vis = input.loc["parameters","reduct_vis"]
    reduct1 =  input.loc["parameters","reduct1"]
    groupby = input.loc["parameters","groupby"]
    #运行cmd主命令
    if not os.path.exists(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/seurat_Velocity.rds"):
        logger.info("未检测到seurat_Velocity.rds,即将生成seurat_Velocity.rds,较耗时！")
        loom = find_loom(wkdir,input)
        loom_files = glob.glob(os.path.join(loom, '*.loom'))
        loom_result_string = ','.join(loom_files)
        cmd = f" Rscript  /opt/softwares/oe_packages/script/velecto_seurat.R" \
          f"  -i  {seurat_ob} " \
          f"  -s  {pointsize} " \
          f"  -l  {loom_result_string} " \
          f"  -o  {wkdir}/velocity " \
          f"  -r {reduct_vis}" \
          f"  --reduct1 {reduct1} " \
          f"  -g {groupby} " \
          f"{add_params}"
    else:
        logger.info("检测到adata_with_scvelo.h5ad,将基于此次的结果进行分析")
        h5ad=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/seurat_Velocity.rds"  
        os.symlink(h5ad,f"{wkdir}/velocity/seurat_Velocity.rds")
        cmd=f"Rscript  /opt/softwares/oe_packages/script/velecto_seurat.R " \
            f" -i {h5ad} " \
            f" -o {wkdir}/velocity " \
            f" -g {groupby}  "\
            f"  -r {reduct_vis}" \
            f" --reduct1 {reduct1}  "\
            f" {add_params}" \
            f"--replot T " \
            f"  -s  {pointsize} " 
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)
    # cmd=f"python  /opt/softwares/oe_packages/script/velocity_plot.py {groupby}" 
    # with module_cmd(environment) as p:
    #     status = p(cmd, projectid, taskid)

    logger.info("The drawing has been completed, and the indexing of the results is started")
    cmd_ln = f"ln -s {wkdir}/velocity/*{{pdf,png}} {wkdir}/output/download/" + "\n"+\
             f"ln -s {wkdir}/velocity/*分析说明.docx  {wkdir}/output/download/" + "\n"+\
             f"ln -s {wkdir}/velocity/*png {wkdir}/output/"
    with module_cmd(environment) as p:
        status=p(cmd_ln, projectid, taskid)
    return status
