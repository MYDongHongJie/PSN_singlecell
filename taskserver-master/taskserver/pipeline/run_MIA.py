# -*-coding:utf-8-*-
from glob import glob

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
def task_MIA(input,projectid="projectid", taskid="taskid", workdir="workdir"):
    wkdir  = workdir
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    environment = "/home/xfyang/modulefiles/OESingleCell/beta"
    if input.loc["parameters", "queryexp"] != "":
        queryexp = f"{wkdir}/input/{input.loc['parameters','queryexp']}"
    else:
        # if user didn't specify the rds path we'll try to find it in the default path, this default path can be introduced from the parameters stored in input.json, so the structure can be discussed later
        queryexp_input_json = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json"
        input_info = pd.read_json(f'{queryexp_input_json}', orient="index")
        queryexp = f"{wkdir}/../../{input_info.loc['base', 'tasks'][0]['projectId']}/{input_info.loc['base', 'tasks'][0]['taskId']}/bclust/data_ob_v3.rds"
    if input.loc["parameters", "refexp"] != "":
        refexp = f"{wkdir}/input/{input.loc['parameters', 'refexp']}"
    else:
        refexp = queryexp
    if input.loc["parameters", "input"] != "":
        plot_input = f"{wkdir}/input/{input.loc['parameters', 'input']}"
    else:
        plot_input = f"{wkdir}/MIA/celltype_enrichment.tsv"
    queryby = input.loc["parameters","queryby"]
    queryfilter = input.loc["parameters","queryfilter"]
    queryassay = input.loc["parameters", "queryassay"]
    queryprefix = input.loc["parameters", "queryprefix"] #shared by mia analysis and results visualize
    if input.loc["parameters", "querymarker"] != "":
        querymarker = f"{wkdir}/input/{input.loc['parameters', 'querymarker']}"
    else:
        querymarker = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/output/all_markers_for_each_cluster_anno.tsv"
    refby = input.loc["parameters", "refby"]
    refmarker = input.loc["parameters", "refmarker"]
    reffilter = input.loc["parameters", "reffilter"]
    refassay = input.loc["parameters", "refassay"]
    refprefix = input.loc["parameters", "refprefix"] #shared by mia analysis and results visualize
    min_pct1 = int(input.loc["parameters", "min_pct1"])
    max_pct2 = int(input.loc["parameters", "max_pct2"])
    pct_fold = int(input.loc["parameters", "pct_fold"])
    avg_log2FC = float(input.loc["parameters", "avg_log2FC"])
    pvalue = float(input.loc["parameters", "pvalue"])
    test = input.loc["parameters", "test"]
    # below params for results visualization
    reorder = input.loc["parameters", "reorder"]
    plot = input.loc["parameters", "plot"]
    # analyze steps: 1. MIA 2. plot, use plot param to determine whether to draw plots
    logger.info("MIA analysis start!")
    step1 = f"module purge && module load {environment} && Rscript  /home/guochenyue/pipeline/run_mia.R" \
          f" -o {wkdir}/MIA" \
          f" mia" \
          f" --queryexp {queryexp}" \
          f" --queryby {queryby}" \
          f" --queryfilter {queryfilter}" \
          f" --queryassay {queryassay}" \
          f" --queryprefix {queryprefix}" \
          f" --querymarker {querymarker}" \
          f" --refexp {refexp}" \
          f" --refby {refby}" \
          f" --refmarker {wkdir}/input/{refmarker}" \
          f" --reffilter {reffilter}" \
          f" --refassay {refassay}" \
          f" --refprefix {refprefix}" \
          f" --min_pct1 {min_pct1}" \
          f" --max_pct2 {max_pct2}" \
          f" --pct_fold {pct_fold}" \
          f" --avg_log2FC {avg_log2FC}" \
          f" --pvalue {pvalue}" \
          f" --test {test}"
    with module_cmd(input.loc['project', 'environment']) as p:
        status = p(step1, projectid, taskid)
    if plot == "yes":
        logger.info("MIA plot start!")
        step2 = f"module purge && module load {environment} && Rscript  /home/guochenyue/pipeline/run_mia.R" \
                f" -o {wkdir}/MIA" \
                f" plot" \
                f" --input {plot_input}" \
                f" --reorder {reorder}" \
                f" --refprefix {refprefix}" \
                f" --queryprefix {queryprefix}"
        with module_cmd(input.loc['project', 'environment']) as p:
            status = p(step2, projectid, taskid)


    logger.info("MIA analysis finished!")
    # generate output.json.tsv so that the diagrams can be uploaded to obs
    d = {"task_type":["MIA"]*3,"result_module":["diagram", "diagram","data"],
    "input":["celltype_heatmap.png", "celltype_sankeyNetwork.html", "celltype_enrichment.tsv"],
    "type":["image","html","tsv"],
    "file":["celltype_heatmap.png", "celltype_sankeyNetwork.html", "celltype_enrichment.tsv"],
    "title":["celltype_heatmap", "celltype_sankeyNetwork", "celltype_enrichment"],
    "downloadName":["celltype_heatmap", "celltype_sankeyNetwork", "celltype_enrichment"],
    "downloadPath":["download/celltype_heatmap.png", "download/celltype_sankeyNetwork.html", "download/celltype_enrichment.tsv"]}
    df = pd.DataFrame(d)
    df.to_csv(f"{wkdir}/MIA/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd_ln = f"ln -s {workdir}/{module}/" \
          f"{{celltype_enrichment*.tsv," \
          f"celltype_heatmap.*," \
          f"celltype_sankeyNetwork.*}}  "  \
          f"{workdir}/output/download/"
    with module_cmd(input.loc['project', 'environment']) as p:
        status = p(cmd_ln, projectid, taskid)
    cmd_cp = f"cp -rf {workdir}/{module}/celltype_heatmap.png {workdir}/output/celltype_heatmap.png"
    with module_cmd(input.loc['project', 'environment']) as p:
        status = p(cmd_cp, projectid, taskid)
    return status
