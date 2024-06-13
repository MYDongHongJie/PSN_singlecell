import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_bclust(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.tsv文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作目录,pojectid层级
    :return: status
    """
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    ## 1.获取filterd.h5seurat路径及参数信息
    # env_module = input.loc['project', 'environment']
    # env_module = "OESingleCell/v_3.0.0_visium"
    env_module = "OESingleCell/v_3.0.0_visium_produce"
    filtered_rds = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    reduct1 = input.loc['parameters', 'reduct1']
    reduct2 = input.loc['parameters', 'reduct2']
    resolution = input.loc['parameters', 'resolution']
    clusteringuse = input.loc['parameters', 'clusteringuse']
    assay = input.loc['parameters', 'assay']
    perplexity = input.loc['parameters', 'perplexity']
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    var2use = input.loc['parameters', 'var2use']
    levels4var = input.loc['parameters', 'levels4var']
    ##
    component = input.loc['parameters', 'component']
    if component == "":
        reduct1_component = {"mnn": 10 ,"pca": None }
        component = reduct1_component.get(input.loc['parameters', 'reduct1'])

    component_parmas = f" --component {component} " if not component is None else f"  "
    ##
    if reduct1 == "mnn" or reduct1 == "harmony":
        reduct1_params = f" --batchid {input.loc['parameters', 'batchid']} "
    else:
        reduct1_params = ""

    ##
    add_params = ""
    add_params2 = ""
    if var2use != "" and levels4var != "":
        levels4var = '"{0}"'.format('","'.join(levels4var.split(",")))
        add_params = f" --predicate '{var2use} %in% c({levels4var})'"
        add_params2 = f" --dosubset TRUE --assay4SCT Spatial "
        
    if cellmeta != "":
        add_params2 = f" --dosubset TRUE --assay4SCT Spatial  --cloudmeta {workdir}/input/{cellmeta} "

    ## 2.获取命令行
    cmd1 = f"/public/scRNA_works/works/guokaiqi/test/st_cloud/exec-dev/sctool " \
           f"-i {filtered_rds} " \
           f"-f rds " \
           f"-o {workdir}/{module} " \
           f"-d rds " \
           f"--assay SCT " \
           f"--image TRUE " \
           f"--prefix singlecell_object.clustering_resolution{resolution}  " \
           f"--update F " \
           f"{add_params} " \
           f"bclust " \
           f"      --reduct1 {reduct1} {reduct1_params} " \
           f"      --reduct2 {reduct2} " \
           f"      {component_parmas} " \
           f"      --clusteringuse {clusteringuse} " \
           f"      --resolution {resolution}" \
           f"      --rerun T " \
           f"      --cloudcfg {output_cfg} " \
           f"{add_params2} " 

    # cmd2 = f"sctool " \
           # f"-i {workdir}/{module}/bclust.h5seurat " \
           # f"-f h5seurat " \
           # f"-o {workdir}/{module} " \
           # f"--assay RNA " \
           # f"summarize " \
           # f"      --reduct {reduct2} " \
           # f"      -c clusters " \
           # f"      --palette customecol2 " \
           # f"      -b sampleid,group" \
           # f"      --dosummary F  "

    ## 3.执行分析
    logger.info("step3.1:执行分析")
    with module_cmd(env_module) as p:
        status = p(cmd1, projectid, taskid)

    # with module_cmd(env_module) as p:
        # status = p(cmd2, projectid, taskid)

    # ==================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    filelist = f"{reduct1}_Dimension_Reduction,{reduct2}_Dimension_Reduction"
    if reduct1 == "mnn" or reduct1 == "harmony":
        filelist = filelist + "," + str(input.loc['parameters', 'batchid']) + "-batchid.xls"
    cmd = f"ln -s {workdir}/{module}/{{{filelist}}} {workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
   # ==================================================================================================================
    # logger.info("step4.convert seurat object to v3")
    # projectname=input.loc["project","name"]
    # cmd=f"v4tov3.r -i {workdir}/{module}/bclust.h5seurat  -f h5seurat -o {workdir}/{module}/"
    # logger.info("bclust: v4tov3_step1")
    # with module_cmd(env_module) as p:
        # status=p(cmd, projectid, taskid)

    # cmd= f"v4tov3.r  -o {workdir}/{module} "
    # logger.info("bclust:  v4tov3_step2")
    # with module_cmd(env_module,"OESingleCell/2.0.0") as p:
        # status=p(cmd, projectid, taskid)

    # cmd=f"ssh scrna@192.168.10.24  '[ -d /oecloud/userfiles/{projectname}  ] && echo Directory Exists || mkdir /oecloud/userfiles/{projectname}' "
    # logger.info("bclust:check or mkdir in cu04")
    # with module_cmd(env_module) as p:
        # status=p(cmd, projectid, taskid)

    # cmd=f"scp {workdir}/{module}/data_ob_v3.rds scrna@192.168.10.24:/oecloud/userfiles/{projectname}/"
    # logger.info("bclust: scp data_ob_v3.rds to cu04")
    # with module_cmd(env_module) as p:
        # status=p(cmd, projectid, taskid)


    return status
