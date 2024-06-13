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

def task_cellchat(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    wkdir =workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/cellchat")
    #获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/*.rds"

    celltype = f"--celltype {input.loc['parameters', 'celltype']} " if input.loc['parameters', 'celltype'] !="" else ""
    groupby = f"--groupby {input.loc['parameters', 'groupby']} " if input.loc['parameters', 'groupby'] !="" else ""
    contrast = f"--contrast {input.loc['parameters', 'contrast']} " if input.loc['parameters', 'contrast'] !="" else ""
    module = input.loc['task', 'type']
    outdir=f"{wkdir}/{module}"
    environment=input.loc["project","environment"]
    assay=input.loc["parameters","assay"]
    species=input.loc["parameters","species"]
    column4cell=input.loc["parameters","column4cell"]

    cmd=f"module purge && source /home/luyao/miniconda3/bin/activate  scrna && "\
        f"Rscript /public/scRNA_works/works/liuhongyan/Test/oecloude/CellChat.R "\
        f"  -i {seurat_ob} "\
        f"  -s {species} "\
        f"  -o {outdir} "\
        f"  --assay {assay} "\
        f"  -p 10 "\
        f"  -c {column4cell}" \
        f" {celltype} {groupby} {contrast}"

    logger.info(cmd)
    logger.info("step1: 开始 cellchat 分析")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    # with module_cmd(environment) as p:
    #     status=p(cmd, projectid, taskid)
    ## 其中 projectid, taskid 是冗余的形参。
    ### env ####
    ret = subprocess.run(cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        encoding="utf-8")
    logger.info(ret.stdout)  
    logger.info(ret.stderr)  

# ==================================================================================================================
    logger.info("step2: 构建output/download/分析结果")
    ## change xls to tsv
    os.rename(f"{workdir}/{module}/communication.xls",
                f"{workdir}/{module}/communication.tsv")

    if not os.path.exists(f'{workdir}/output/download/'):
        os.makedirs(f'{workdir}/output/download/')
    cmd = f"ln -s {workdir}/{module}/*.png {workdir}/{module}/*.pdf {workdir}/{module}/communication.tsv {workdir}/{module}/CellChat细胞通讯分析说明* {workdir}/{module}/signaling_pathway_visualize {workdir}/output/download/"
    os.system(cmd)

    
    ## 图形 tsv ### 
   # ==================================================================================================================

    ## =================================================
    ## 2. amending output.json
    ## =================================================
    logger.info("step2 生成output.json.tsv") 

    df = pd.DataFrame(columns=['task_type', 'result_module', 'input', 'type', 'file', 'title','downloadName', 'downloadPath'])
    df.loc[len(df.index)] = ["cellchat","data","communication.tsv","ligand_receptor_pair","communication.tsv","communication","communication","download/communication.tsv"]
    df.loc[len(df.index)] = ["cellchat","diagram","significant_interactions_bubble_plot.pljson","oe","significant_interactions_bubble_plot.pljson","significant_interactions_bubble_plot","",""]  
   
    if input.loc['parameters', 'groupby'] != "":
        for contrast in input.loc['parameters', 'contrast'].split("+"):
            case = contrast.split(":")[0]
            con = contrast.split(":")[1]
            df.loc[len(df.index)] = ["cellchat","diagram",f"diff_{case}-vs-{con}_bubble_plot.pljson","oe",f"diff_{case}-vs-{con}_bubble_plot.pljson",f"diff_{case}-vs-{con}_bubble_plot","",""]

    df.to_csv(f"{workdir}/{module}/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')

