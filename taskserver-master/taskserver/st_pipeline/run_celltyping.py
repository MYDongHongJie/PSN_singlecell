import pandas as pd 
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_celltyping(input, output_cfg, projectid="项目ID", taskid="任务ID", workdir="运行目录"):
    # workdir = workdir
    # os.makedirs( f"{workdir}/output/download")
    # os.makedirs( f"{workdir}/celltyping")
    logger.info("解析命令行")
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # env_module = input.loc['project', 'environment']
    # env_module = "OESingleCell/v_3.0.0_visium"
    env_module = "OESingleCell/v_3.0.0_visium_produce"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    else:
        pass

    seurat_ob = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    ref_rds = input.loc['parameters', 'ref_singlecell']
    refcelltype = input.loc['parameters', 'refcelltype']
    refassay = input.loc['parameters', 'refassay']
    refmarker = "auto" if input.loc['parameters', 'refmarker'] == "" else input.loc['parameters', 'refmarker']
    # celltype_visual 
    groupby = input.loc["parameters","groupby"]
    piescale = input.loc["parameters","piescale"]
    # miscname= "spotlight_results" if input.loc['parameters', 'celltyping_method'] == "spotlight" else "RCTD_results"
    miscname= f"{input.loc['parameters', 'celltyping_method']}_results"
    library_size = input.loc["parameters","library_size"]
    if library_size == "cytassist":
        crop_param = "--crop TRUE"
    else:
        crop_param = " "
    if input.loc['parameters', 'celltyping_method'] == "spotlight":
        cmd = f"sctool "\
              f" -i {seurat_ob} "\
              f" -f rds "\
              f" -o {workdir}/{module} "\
              f" -d rds "\
              f" --assay Spatial "\
              f" --dataslot counts "\
              f" --prefix 'spatial_seurat_celltype' "\
              f"{add_params} " \
              f" st_deconv "\
              f" --refexp {ref_rds}"\
              f" --refcelltype {refcelltype} "\
              f" --refassay {refassay} "\
              f" --refmarker {refmarker} "
    else:
        cmd = f"sctool "\
              f" -i {seurat_ob} "\
              f" -f rds "\
              f" -o {workdir}/{module} "\
              f" -d rds "\
              f" --assay Spatial "\
              f" --dataslot counts "\
              f" --prefix 'spatial_seurat_celltype' "\
              f"{add_params} " \
              f" RCTD "\
              f" --refexp {ref_rds}"\
              f" --refcelltype {refcelltype} "\
              f" --refassay {refassay} "\
              f" --doublet_mode FALSE "

    ##  spotlight/RCTD  celltype visualization
    cmd_Vis = f"/public/scRNA_works/works/guochy/ST_taskserver/exec/scVis "\
              f" -i {workdir}/{module}/spatial_seurat_celltype.rds "\
              f" -f rds "\
              f" -o {workdir}/{module} "\
              f" --assay Spatial "\
              f" spatialpieplot "\
              f" --group.by {groupby} "\
              f" --piescale {piescale} "\
              f" --misclist {miscname} "\
              f"{crop_param} "\
              f" --pt.alpha TRUE "

    logger.info(cmd)
    logger.info("celltyping:执行分析")
    with module_cmd(env_module) as p:
        status=p(cmd, projectid, taskid)
    logger.info(cmd_Vis)
    logger.info("celltyping:可视化")
    with module_cmd(env_module) as p:
        status=p(cmd_Vis, projectid, taskid)

    # ==================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd_ln = f"""shopt -s extglob  
             ln -s {workdir}/{module}/!(*rds|*tsv) {workdir}/output/download/"""

    with module_cmd(env_module) as p:
        status = p(cmd_ln, projectid, taskid)
#    print([os.path.abspath(os.path.join(f"{workdir}/{module}", file)) for file in os.listdir(f"{workdir}/{module}") if file.endswith(".png")])
    png_files = " ".join([os.path.abspath(os.path.join(f"{workdir}/{module}", file)) for file in os.listdir(f"{workdir}/{module}") if file.endswith(".png")])

    cmd_image_ln = f"ln -s {png_files} {workdir}/output/"
    with module_cmd(env_module) as p:
        status = p(cmd_image_ln, projectid, taskid)
    # ==================================================================================================================
    # #cmd_ln = f"ln -s {workdir}/celltyping/metadata.tsv {workdir}/output/"+"\n"
    # cmd_ln =""
    # #print(cmd_ln)
    # for i in celltypingdb_ttt:
        # cmd_ln = cmd_ln + "\n" + \
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_plot.pdf {workdir}/output/download/{species}ref_{i}_{level}_celltyping_plot.pdf"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_plot.png {workdir}/output/download/{species}ref_{i}_{level}_celltyping_plot.png"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_heatmap.pdf {workdir}/output/download/{species}ref_{i}_{level}_celltyping_heatmap.pdf"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_heatmap.png {workdir}/output/download/{species}ref_{i}_{level}_celltyping_heatmap.png"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_results.xls {workdir}/output/download/{species}ref_{i}_{level}_celltyping_results.xls"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_simplified_celltyping_results.csv {workdir}/output/download/{species}ref_{i}_{level}_simplified_celltyping_results.csv"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_top.{level}_celltyping_plot.pdf {workdir}/output/download/{species}ref_{i}_top.{level}_celltyping_plot.pdf"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_top.{level}_celltyping_plot.png {workdir}/output/download/{species}ref_{i}_top.{level}_celltyping_plot.png"+"\n"+\
            # f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_statistics.tsv {workdir}/output/download/{species}ref_{i}_{level}_celltyping_statistics.xls"+"\n"
			# #f"ln -s {workdir}/celltyping/{species}ref_{i}_top.{level}_celltyping_data.tsv {workdir}/output/{species}ref_{i}_top.{level}_celltyping_data.tsv"+"\n"+\
            # #f"ln -s {workdir}/celltyping/{species}ref_{i}_{level}_celltyping_data.tsv {workdir}/output/{species}ref_{i}_{level}_celltyping_data.tsv"
    # #print(cmd_ln)        
    # with module_cmd(environment) as p:
        # status = p(cmd_ln, projectid, taskid)
    # ==================================================================================================================
    return status ###run_cmd use the default function from xiufeng.yang
