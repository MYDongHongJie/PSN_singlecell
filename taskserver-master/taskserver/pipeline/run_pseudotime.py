import os, sys
import json
import subprocess
import pandas as pd
import requests
import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')
from taskserver.tools.module_cmd import module_cmd


def task_pseudotime(input, projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/pseudotime")
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {wkdir}/input/{cellmeta} "
    #seurat_ob = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/bclust/data_ob_v3.rds"
    types = ('*.rds', '*.h5seurat')
    seurat_file = []
    for files in types:
        seurat_file.extend(glob.glob(f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*" + files))
    rds_input = [element for element in seurat_file if "rds" in element]
    h5seurat_input = [element for element in seurat_file if "h5seurat" in element]
    if (len(rds_input) == 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
    elif (len(rds_input) != 0 and len(h5seurat_input) == 0):
        seurat_ob = rds_input[0]
    elif (len(rds_input) != 0 and len(h5seurat_input) != 0):
        seurat_ob = rds_input[0]
    else:
        print(f"输入seurat对象为空，请检查！！！")
    outdir = f"{wkdir}/pseudotime"

    #colorby = ",".join(input.loc["parameters", "colorby"])
    colorby = input.loc["parameters", "colorby"]
    design = input.loc["parameters", "design"]
    # if input.loc["parameters", "groupby"] != "":
        # add_params = add_params + f" --groupby  {input.loc['parameters', 'groupby']} "
    # if input.loc["parameters", "which_cells"] != "":
        # add_params = add_params + f" --which_cells  {input.loc['parameters', 'which_cells']} "
    if input.loc["parameters", "downsample"] != "":
        add_params = add_params + f" --downsample  {int(input.loc['parameters', 'downsample'])} "
    #if input.loc["parameters", "root"] != "":
    #    add_params = add_params + f" --root  {input.loc['parameters', 'root']} "
    if input.loc["parameters", "rds"] :
        add_params = add_params + f" --rds {wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/pseudotime/pseudotime_results.rds  "
    assay = input.loc["parameters", "assay"]
    module = input.loc['task', 'type']
    env_module = input.loc['project', 'environment']

    environment = "OESingleCell/2.0.0"
    cmd = f"source deactivate  && source activate /opt/softwares/mamba-forge/envs/monocle2 &&  /opt/softwares/mamba-forge/envs/monocle2/bin/Rscript /mnt/monocle2/monocle.R " \
          f"  -i {seurat_ob} " \
          f"  -f seurat"\
          f"  -o {outdir} " \
          f"  --colorby {colorby}" \
          f"  --design {design}" \
          f"    {add_params} "
    logger.info(cmd)
    logger.info("step1: 开始pseudotime分析")
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    logger.info("step2:生成绘图json文件")
    cmd2=f"cd {outdir} \n"
    colorby2=colorby.split(",")
    colorby2.append("Pseudotime")
    #for i in colorby2:
    # for i in colorby.split(","):
        # cmd2=cmd2+f"  source deactivate  && /opt/softwares/conda_envs/oesinglecell/bin/python  /mnt/monocle2/pseudotime.py {i} \n"
    # with module_cmd(environment) as p:
        # status = p(cmd2, projectid, taskid)
#################################################################################
    logger.info("step3:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    #filelist = f"*.tsv,*.pdf,*.png,*.docx"
    filelist = f"Pseudotime_State_count_info.tsv,*.pdf,*.png,*.docx"
    cmd3 = f"ln -s {workdir}/{module}/{{{filelist}}} {workdir}/output/download/ \n"
    # cmd3=cmd3+f"ln -s {workdir}/{module}/*json {workdir}/output/"
    #修改pljson 为 tsv 进行可视化，type 为 monocle2 ：20240123
    cmd3=cmd3+f"ln -s {workdir}/{module}/pseudotime*tsv {workdir}/output/"
    with module_cmd(env_module) as p:
        status = p(cmd3, projectid, taskid)
