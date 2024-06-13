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
def task_Nichenet(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/Nichenet")
    environment = "/home/xfyang/modulefiles/OESingleCell/dev"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    seurat_ob = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    outdir=f"{wkdir}/Nichenet"
    assay=input.loc["parameters","assay"]
    celltype=input.loc["parameters","celltype"]
    species=input.loc["parameters","species"]
    if input.loc["parameters", "database"] != "":
        database =  input.loc["parameters", "database"]
    else:
        database = "/data/database/NicheNet"
    receiver=input.loc["parameters","receiver"]
    sender=input.loc["parameters","sender"]
    genelist = input.loc["parameters", "genelist"]
    topn=int(input.loc["parameters","topn"])
    if input.loc["parameters", "pct"] != "":
        pct =  float(input.loc["parameters", "pct"])
    else:
        pct = 0.1
    if input.loc["parameters", "ligand_visiual_cutoff"] != "":
        ligand_visiual_cutoff =  float(input.loc["parameters", "ligand_visiual_cutoff"])
    else:
        ligand_visiual_cutoff = 0.25


    cmd=f"/public/scRNA_works/works/ziqingzhen/test/Nichenet_cloud/exec/sctool "\
        f"  -i {seurat_ob} "\
        f"  -f rds "\
        f"  -o {outdir} "\
        f"  --assay {assay} "\
        f"{add_params} " \
        f" nichenet "\
        f"  --celltype {celltype}  "\
        f"  --species {species} "\
        f"  --database {database}  "\
        f"  --receiver {receiver} "\
        f"  --sender {sender} "\
        f"  --genelist {wkdir}/input/{genelist} "\
        f"  --topn {topn} "\
        f"  --pct {pct} "\
        f"  --ligand_visiual_cutoff {ligand_visiual_cutoff} "

    logger.info(cmd)
    logger.info(" Nicenet analysis start!")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    # generate output.json.tsv so that the diagrams can be uploaded to obs
    d = {"task_type":["Nichenet"]*5,
    "result_module":["diagram","diagram","data","data","data"], 
    "input":["combined_plot_receptor.png","combined_plot_target.png","vis_ligand_receptor_network.tsv", "ligand_activities.tsv", "vis_ligand_target.tsv"], 
    "type":["image","image","gene","gene","gene"], 
    "file":["combined_plot_receptor.png","combined_plot_target.png","vis_ligand_receptor_network.tsv", "ligand_activities.tsv", "vis_ligand_target.tsv"],
    "title":["combined_plot_receptor","combined_plot_target","vis_ligand_receptor_network","ligand_activities","vis_ligand_target"], 
    "downloadName":["combined_plot_receptor","combined_plot_target","ligand_receptor", "ligand_activities", "ligand_target"], 
    "downloadPath":["combined_plot_receptor.png", "combined_plot_target.png","ligand_activities.tsv", "vis_ligand_receptor_network.tsv","vis_ligand_target.tsv"]}
    df = pd.DataFrame(d)
    #df = pd.DataFrame.from_dict(d, orient='index')
    df.to_csv(f"{wkdir}/Nichenet/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")
    cmd_ln = f"ln -s {wkdir}/Nichenet/{{[vl]*.tsv,*.pdf,*.png,*.docx}}  {wkdir}/output/download/;cp {wkdir}/Nichenet/*.png  {wkdir}/output/"
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)

    return status
