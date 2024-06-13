import pandas as pd
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger
from glob import glob
logger = getLogger('oe.cloud.sc.qsub')


from taskserver.tools.module_cmd import module_cmd


def task_celltyping(input,projectid="项目ID", taskid="任务ID", workdir="运行目录"):
    wkdir = workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/celltyping")
    #input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{wkdir}/input/cell.tsv"
    cellmeta=input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    else:
        pass
    #seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/bclust/bclust.h5seurat"
    ####进行输入对象rds或者h5seurat格式的识别和判断
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files

    rds_input = search_files(f'{base_wd}/{base_module_name}', "rds")
    h5seurat_input = search_files(f'{base_wd}/{base_module_name}', "h5seurat")
    seurat_ob = ""
    format_ob = ""
    if (len(rds_input) == 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    elif (len(rds_input) != 0 and len(h5seurat_input) == 0):
        seurat_ob = rds_input[0]
        format_ob = "rds"
    elif (len(rds_input) != 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    else:
        print(f"输入seurat对象为空，请检查！！！")
    #seurat_ob_realpath=os.readlink(seurat_ob)
    ####以上完成输入对象的识别及赋值

    outdir=f"{wkdir}/celltyping"
    # celltypingdb=input.loc["parameters","ref_dataset"]
    # celltypingdb_ttt=celltypingdb.split(",")
    # celltypingdb_list=",".join(["/data/database/celltype_refdata/logNorm_rds/"+x+".rds" for x in celltypingdb_ttt])
    celltypingdb_list=input.loc["parameters","ref_dataset"]
    celltypingdb_ttt =celltypingdb_list.split("/")[-1].split(".")[0]
    pointsize=1
    reduct2=input.loc["parameters","reduct2"]
    species=input.loc["parameters","species_celltype"]
    assay=input.loc["parameters","assay"]
    level=input.loc["parameters","LEVEL"]
    script_dir="."
    environment=input.loc["project","environment"]
    ## 针对弹性云，更改默认环境以运行singlR
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " source /opt/softwares/oe_packages/.bashrc  &&  "
    else:
        tmp_envs=" "
    cmd=f"{tmp_envs}  sctool "\
        f" -i {seurat_ob} "\
        f" -f h5seurat "\
        f" -j  5 " \
        f" -o {outdir} "\
        f" -d h5seurat "\
        f" --assay {assay} "\
        f" --dataslot counts,data,scale.data "\
        f"{add_params} " \
        f" celltyping "\
        f" -r {celltypingdb_list} "\
        f" --usecluster F "\
        f" --demethod classic "\
        f" --annolevel {level} "\
        f" -v {pointsize} "\
        f" -n 25 "\
        f" --reduct {reduct2} "\
        f" --species {species}"

    logger.info(cmd)
    logger.info("celltyping:执行分析")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    seurat_ob_new=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['task', 'id']}/celltyping/celltyping.h5seurat"
    cwd=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['task', 'id']}/celltyping/"
    projectname=input.loc["project","name"]
    #v4tov3_script="/data/software/conda_envs/OESingleCell_3.0.0_scRNAcloud/lib/R/library/OESingleCell/exec/v4tov3.r"

    cmd2=f"  v4tov3  -i {seurat_ob_new} -f h5seurat -o {cwd} "
    logger.info(cmd2)
    logger.info("celltyping: v4tov3_step1")
    with module_cmd(environment) as p:
        status=p(cmd2, projectid, taskid)
    ## 针对弹性云，更改seurat_3.1.2环境以运行
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " source /opt/softwares/oe_packages/.bashrc_seurat_3.1.2  &&  "
    else:
        tmp_envs=" "
    cmd3= f" {tmp_envs}  v4tov3 -o {cwd}"
    logger.info(cmd3)
    logger.info("celltyping: v4tov3_step2")
    with module_cmd("OESingleCell/2.0.0") as p:
        status=p(cmd3, projectid, taskid)

    # cmd4=f"ssh scrna@192.168.10.24  '[ -d /oecloud/userfiles/{projectname}  ] && echo Directory Exists ||   mkdir /oecloud/userfiles/{projectname}' "
    # logger.info(cmd4)
    # logger.info("celltyping:mkdir in cu04")
    # with module_cmd(environment) as p:
    #     status=p(cmd4, projectid, taskid)
    #
    # cmd5=f"scp {cwd}/data_ob_v3.rds scrna@192.168.10.24:/oecloud/userfiles/{projectname}/"
    # logger.info(cmd5)
    # logger.info("celltyping: scp data_ob_v3.rds to cu04")
    # with module_cmd(environment) as p:
    #     status=p(cmd5, projectid, taskid)


    #cmd_ln = f"ln -s {wkdir}/celltyping/metadata.tsv {wkdir}/output/"+"\n"
    cmd_ln =""
    #print(cmd_ln)
    for i in [celltypingdb_ttt]:
        cmd_ln = cmd_ln + "\n" + \
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_celltyping_plot.pdf {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_celltyping_plot.png {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_celltyping_heatmap.pdf {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_celltyping_heatmap.png {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_celltyping_results.xls {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_simplified_celltyping_results.csv {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_top.{level}_celltyping_plot.pdf {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_top.{level}_celltyping_plot.png {wkdir}/output/download/" + "\n"+\
            f"ln -s {wkdir}/celltyping/{species}ref_*_{level}_celltyping_statistics.tsv {wkdir}/output/download/" + "\n"
			#f"ln -s {wkdir}/celltyping/{species}ref_{i}_top.{level}_celltyping_data.tsv {wkdir}/output/{species}ref_{i}_top.{level}_celltyping_data.tsv"+"\n"+\
            #f"ln -s {wkdir}/celltyping/{species}ref_{i}_{level}_celltyping_data.tsv {wkdir}/output/{species}ref_{i}_{level}_celltyping_data.tsv"
    #print(cmd_ln)
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status ###run_cmd use the default function from xiufeng.yang
