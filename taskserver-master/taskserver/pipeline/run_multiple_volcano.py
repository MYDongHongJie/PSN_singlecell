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

#===================input multiple_volcano input and run the cmds===================================================================================
def task_multiple_volcano(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir  = workdir
    environment = input.loc['project', 'environment']
    if not os.path.exists( f"{wkdir}/output/download"):
        os.makedirs( f"{wkdir}/output/download")
    if not os.path.exists( f"{wkdir}/multiple_volcano"):
        os.makedirs( f"{wkdir}/multiple_volcano")

    # genelist = input.loc["parameters","genelist"] if input.loc["parameters","genelist"] != "" else ""
    # if not genelist:
    #     logger.info("There is no file input")
    #
    # topn = f"--topn {int(input.loc['parameters','topn'])} " if input.loc['parameters','topn' ] != "" else ""
    # padj = f"--padj {float(input.loc['parameters','padj'])}" if input.loc['parameters','padj'] !="" else ""
    # pointsize = f"--pointsize {float(input.loc['parameters','pointsize'])} " if input.loc["parameters","pointsize"] != "" else ""
    # cluster = f"--cluster {input.loc['parameters','cluster']}" if input.loc['parameters','cluster'] !="" else ""
    # topby = input.loc["parameters","topby"]
    # logger.info(f"{genelist}")
    if input.loc['parameters','genelist'] =="":
        base_taskid = input.loc["base", "tasks"][0]["taskId"]
        base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
        base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
        base_module_name = base_json.loc['task', 'type']
        base_file=f"{base_wd}/{base_module_name}/all_markers_for_each_cluster_anno.tsv"
        add_param = f"-i {base_file}"
    else:
        add_param = f" --genelist {workdir}/input/{input.loc['parameters','genelist']}"
    topn=input.loc['parameters','topn']
    pvalue=input.loc['parameters','pvalue']
    fc=input.loc['parameters','fc']
    ###使用方式,数据输入有两种，1.上传的输入数据列名为gene,p_val,FC,cluster列的tab分隔文件；2.上游依赖任务模块为findallmarkers的结果文件all_markers_for_each_cluster_anno.tsv
    # Rscirpt multiple_volcano.R -i xxx.file -o ./ -n NULL -p 0.05 -f 1.5
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " /opt/softwares/envs/oesinglecell_v4/bin/Rscript /opt/softwares/oe_packages/cloud_report_script/scripts_use/run_multiple_volcano/multiple_volcano.r "
    else:
        tmp_envs="module purge && module load OESingleCell/2.0.0 && Rscript /public/scRNA_works/works/xujingmei/cloud/multiple_volcano/multiple_volcano.r"
    topn_param=f"-n {int(topn)}" if topn!="" else ""
    cmd = f" {tmp_envs}" \
            f"  {add_param} " \
            f"  -o  {wkdir}/multiple_volcano/  "\
            f" {topn_param}"\
            f" -p {pvalue}"\
            f" -f {fc}"\

    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    logger.info("The drawing has been completed, and the indexing of the results is started")

    cmd_ln = f"ln -s {wkdir}/multiple_volcano/*.pdf  {wkdir}/output/download/"+"\n"+\
             f"ln -s {wkdir}/multiple_volcano/*.png  {wkdir}/output/download/"+"\n"+\
             f"ln -s {wkdir}/multiple_volcano/*.png  {wkdir}/output/"+"\n"+\
             f"ln -s {wkdir}/multiple_volcano/*.tsv  {wkdir}/output/download/"+"\n"+\
             f"ln -s {wkdir}/multiple_volcano/*.tsv  {wkdir}/output/"

    d = {
            "task_type":["multiple_volcano"]*2,
            "result_module":["diagram","data"],
            "input":["multiple_volcano_plot.png","multiple_volcano_plot_diffexp_gene_result.tsv"],
            "type":["image","gene"],
            "file":["multiple_volcano_plot.png","multiple_volcano_plot_diffexp_gene_result.tsv"],
            "title":["multiple_vlocano","multiple_volcano_plot_diffexp_gene_result"],
            "downloadName":["",""] ,
            "downloadPath":["",""]
        }

    df = pd.DataFrame(d)
    df.to_csv(f"{workdir}/multiple_volcano/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status


