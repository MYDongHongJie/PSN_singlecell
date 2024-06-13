import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd
from glob import glob
logger = getLogger('oe.cloud.sc.qsub')
from taskserver.tools.module_cmd import module_cmd


def task_cellcycle(input, projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/cellcycle")
    # input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    # input_meta=f"{wkdir}/input/cell.tsv"
    #cellmeta = f"{wkdir}/input/{input.loc['base', 'tasks'][0]['filename']}"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    env_module=input.loc['project', 'environment']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {wkdir}/input/{input.loc['base', 'tasks'][0]['filename']} "
    #seurat_ob = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/bclust/bclust.h5seurat"
    ####获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files
    seurat_ob = search_files(f'{base_wd}/{base_module_name}', "rds")[0]

    outdir = f"{wkdir}/cellcycle"
    species = input.loc["parameters", "species"]
    method = input.loc["parameters", "method"]
    splitby = ",".join(input.loc["parameters", "splitby"])
    bartype = input.loc["parameters", "bartype"]
    reduct = input.loc["parameters", "reduct"]
    module = input.loc['task', 'type']
    environment = input.loc['project', 'environment']
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " /opt/softwares/envs/oesinglecell_v3/bin/Rscript  "
    else:
        tmp_envs="module purge && module load OESingleCell/3.0.d && Rscript"
    cmd = f"{tmp_envs} /opt/softwares/oe_packages/cloud_report_script/scripts_use/run_cellcycle/cellcycle.R " \
          f"  -i {seurat_ob} " \
          f"  -o {outdir} " \
          f"  --species {species}" \
          f"  --informat rds"\
          f"  --method {method}" \
          f"  --splitby   {splitby}" \
          f"  --bartype {bartype}" \
          f"  --reduct {reduct}" \
          f"    {add_params} "
    logger.info(cmd)
    logger.info("step1: 开始cellcycle分析")
    # with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    logger.info("step2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    filelist = f"*.tsv,*.pdf,*.png"
    cmd = f"ln -s {workdir}/{module}/{{{filelist}}} {workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    logger.info("step3:convert seurat object to v3")
    projectname = input.loc["project", "name"]
    ## 针对弹性云，不必运行下面代码
    # if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
    #     tmp_envs = " /opt/softwares/envs/oesinglecell_v3/bin/Rscript /opt/softwares/oe_packages/cloud_report_script/scripts_use/v4tov3.R "
    # else:
    #     tmp_envs=" v4tov3.R"
    # cmd = f"{tmp_envs}   -i {workdir}/{module}/cellcycle_seurat.h5seurat  -f h5seurat -o {workdir}/{module}/"
    # logger.info("step3.1: v4tov3_step1")
    # with module_cmd(env_module) as p:
    #     status = p(cmd, projectid, taskid)
    # 
    # cmd = f" {tmp_envs}     -o {workdir}/{module} "
    # logger.info("step3.2:  v4tov3_step2")
    # with module_cmd("OESingleCell/2.0.0") as p:
    #     status = p(cmd, projectid, taskid)
