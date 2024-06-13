#corscatter
import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from glob import glob
from taskserver.tools.module_cmd import module_cmd

def task_corscatter(input, projectid="项目ID", taskid="任务ID", workdir="分析目录"):

    logger = getLogger('oe.cloud.sc.qsub')
    wkdir=workdir
    environment="OESingleCell/2.0.0"
    if not os.path.exists(f"{wkdir}/output/download"):
        os.makedirs(f"{wkdir}/output/download")
    if not os.path.exists(f"{wkdir}/corscatter"):
        os.makedirs(f"{wkdir}/corscatter")


    ##################### Input files #####################

    #Input rds file
    #input_rds=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    ####获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files
    input_rds = search_files(f'{base_wd}/{base_module_name}', "rds")[0]
    #Input genelist
    if input.loc["parameters","genelist"] != "":
        genelist=input.loc["parameters","genelist"]
    else:logger.info("There must be genelist input")


    ##################### Get parameters  ##################

    #Point size, the default is 0.5
    if input.loc["parameters","ptsize"] != "":
        ptsize=float(input.loc["parameters","ptsize"])
    else:ptsize=0.5

    #Groupby, the default is clusters
    if input.loc["parameters","groupby"] != "" :
        groupby=input.loc["parameters","groupby"]
    else:
        groupby="clusters"

    #
    if input.loc["parameters", "var2use"] != "" and  input.loc["parameters", "levels4var"] != "" :
        var2use= input.loc["parameters", "var2use"]
        levels4var = input.loc["parameters", "levels4var"]
        subset = f" -q  {var2use}  -u  {levels4var} "
    elif input.loc["parameters", "var2use"] == "" and  input.loc["parameters", "levels4var"] == "" :
        subset=""
    else:
        logger.info("var2use and levels4var are supporting parameters, which must be provided at the same time")
        subset = ""

    #Order, the order of the clusters can be set
    if input.loc["parameters","order"] != "" :
        order=f" -w {input.loc['parameters','order']}"
    else: order=""

    #The method of reduction, the default value is umap
    if input.loc["parameters","reduct"] != "":
        reduct=input.loc["parameters","reduct"]
    else:
         reduct="umap"

    #Assay category, the default is RNA
    if input.loc["parameters","assay"] != "":
        assay=input.loc["parameters","assay"]
    else:assay="RNA"

    ###################### Run code  ######################
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " /opt/softwares/envs/oesinglecell_v3/bin/Rscript  "
    else:
        tmp_envs="module purge && module load OESingleCell/2.0.0 && Rscript"
    cmd=f"{tmp_envs}   /opt/softwares/oe_packages/cloud_report_script/scripts_use/run_corscatter/visualize_markers.R "\
        f" -v {input_rds}"\
        f" -x {wkdir}/input/{genelist}"\
        f" -o {wkdir}/corscatter"\
        f" -m corscatter"\
        f" -s {ptsize}"\
        f" -e {assay}"\
        f" -g {groupby}"\
        f"{subset}"\
        f" --reduct {reduct}"\
        f"{order}"

    with module_cmd(environment) as p:
            status = p(cmd, projectid, taskid)

    ###################### Link files to download directory  ######################
    cmd_ln = f"ln -s {wkdir}/corscatter/geneset_corscatter/{{*.png,*.pdf}} {wkdir}/output/download/"
    with module_cmd(environment) as p:
            status = p(cmd_ln, projectid, taskid)

    ####################### cp png to output ##############
    cmd_png=f"ln -s {wkdir}/corscatter/geneset_corscatter/geneset_corscatter.png {wkdir}/output/"
    with module_cmd(environment) as p:
            status = p(cmd_png, projectid, taskid)

    ####################### cp genelist to output ##############
    filename=genelist.split('.')[0]
    cmd_genelist=f"cp  {wkdir}/input/{genelist} {wkdir}/output/{filename}.tsv"
    with module_cmd(environment) as p:
            status = p(cmd_genelist, projectid, taskid)

    ####################### generate output.json.tsv  ######################
    d={"task_type":["corscatter","corscatter"],"result_module":["diagram","data"],"input":["geneset_corscatter.png",f"{filename}.tsv"],"type":["image","gene"],"file":[f"geneset_corscatter.png",f"{filename}.tsv"],"title":["corscatter","genelist"],"downloadName":["",""],"downloadPath":["",""]}
    df=pd.DataFrame(d)
    df.to_csv( f"{wkdir}/corscatter/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")

    return status

