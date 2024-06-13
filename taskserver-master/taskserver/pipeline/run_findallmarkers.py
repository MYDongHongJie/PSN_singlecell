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


def task_findallmarkers(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/findallmarkers")
    #input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{wkdir}/input/cell.tsv"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
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
    outdir=f"{wkdir}/findallmarkers"
    reduct2=input.loc["parameters","reduct2"]
    environment=input.loc["project","environment"]
    database=input.loc["parameters","database"]
    assay=input.loc["parameters","assay"]
    col=input.loc["parameters","col"]
    Significance=input.loc["parameters","marker_Significance"]
    Significance_threshold=input.loc["parameters","marker_Significance_threshold"]
    add_params_sign=""
    if Significance=="pvalue":
        add_params_sign=add_params_sign+f"  -p {Significance_threshold} "
    elif Significance=="FDR":
        add_params_sign=add_params_sign+f"  -q {Significance_threshold} "
    groupby=input.loc["parameters","marker_groupby"]
    method=input.loc["parameters","marker_method"]
    topnum=int(input.loc["parameters","topn_marker"])

    cmd=f"sctool "\
        f"  -i {seurat_ob} "\
        f"  -f {format_ob} "\
        f"  -o {outdir} "\
        f"  --assay {assay} "\
        f"  --dataslot data,counts,scale.data  "\
        f"  -j 10 "\
        f"{add_params} " \
        f" findallmarkers "\
        f"  -c 2 "\
        f"  -N {topnum} "\
        f"  -k 1 "\
        f"{add_params_sign} "\
        f"  --strict F "\
        f"  -p {Significance_threshold} "\
        f"  -e {method} "\
        f"  -n {groupby} "\
        f"  --reduct2 {reduct2}"\

    logger.info(cmd)
    logger.info("step1: 开始计算marker")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    ##############################################################################################################################
    ##############################################################################################################################      
    topn_marker=f"{wkdir}/findallmarkers/top{topnum}_markers_for_each_cluster.tsv"
    cmd2=f"scVis "\
        f" -i {seurat_ob} "\
        f" -f {format_ob} "\
        f" -o {outdir} "\
        f" -t 10 "\
        f" --assay RNA "\
        f" --slot data,scale.data "\
        f"{add_params} " \
        f" --reduct {reduct2} "\
        f" heatmap "\
        f" -l {topn_marker} "\
        f" -c gene_diff "\
        f" -n {topnum} "\
        f" -g {groupby} "\
        f" --group_colors {col} "\
        f" --sample_ratio 0.8 "\
        f" --style seurat "\


    logger.info(cmd2)
    logger.info("step2: top10 marker 画热图")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd2, projectid, taskid)

    ##############################################################################################################################
    ##############################################################################################################################  

    cmd3=f"sctool" \
        f"  -i {seurat_ob} "\
        f"  -f {format_ob} "\
        f"  -o {outdir} "\
        f"  -j 10 "\
        f"  --assay RNA "\
        f"  --dataslot data "\
        f"{add_params} " \
        f"  visualize "\
        f"  -l {topn_marker} "\
        f"  -g {groupby} "\
        f"  --reduct {reduct2} "\
        f"  --topn  {topnum}  "\
        f"  --topby gene_diff "\
        f"  -m vlnplot,featureplot "\
        f"  --vcolors {col} "\
        f"  --ccolors spectral "\
        f"  --dodge F"

    logger.info(cmd3)
    logger.info("step3: top10 marker 画小提琴图和散点图")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd3, projectid, taskid)

    ##############################################################################################################################
    ##############################################################################################################################           
    all_DEG=f"{wkdir}/findallmarkers/all_markers_for_each_cluster.tsv"

    cmd4=f"sctool annotation "\
        f"-g {all_DEG} "\
        f"--anno {database}/annotation/gene_annotation.xls"

    logger.info(cmd4)
    logger.info("step4:  allmarkers 注释")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd4, projectid, taskid)

    ##############################################################################################################################
    ##############################################################################################################################  
    cmd5=f"sctool annotation "\
        f"-g {topn_marker} "\
        f"--anno {database}/annotation/gene_annotation.xls"
        
    logger.info(cmd5)
    logger.info("step5:  topN markers 注释")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd5, projectid, taskid)

    cluster_number_file=pd.read_csv(f"{wkdir}/findallmarkers/clusters_number.tsv",sep="/t")
    cluster_num=np.shape(cluster_number_file)[0]

    cmd_ln = f"ln -s {wkdir}/findallmarkers/all_markers_for_each_cluster_anno.tsv {wkdir}/output/download/all_markers_for_each_cluster_anno.xls" +"\n"+\
        f"ln -s {wkdir}/findallmarkers/top{topnum}_markers_for_each_cluster_anno.tsv {wkdir}/output/download/top{topnum}_markers_for_each_cluster_anno.xls" +"\n"+\
        f"ln -s {wkdir}/findallmarkers/{{topmarker_gene_heatmap.png,topmarker_gene_heatmap.pdf}} {wkdir}/output/download/"+"\n"
        #f"ln -s {wkdir}/findallmarkers/topmarker_gene_heatmap.tsv {wkdir}/output/"+"\n"
    #print(cmd_ln)
    for i in cluster_number_file.cluster:
        cmd_ln = cmd_ln + "\n" + \
            f"ln -s {wkdir}/findallmarkers/markers_vis4cluster{i} {wkdir}/output/download/"+"\n"
            # f"ln -s {wkdir}/findallmarkers/markers_vis4cluster{i}/marker_gene_featureplot.png {wkdir}/output/download/featureplot_cluster{i}.png"+"\n"+\
            # f"ln -s {wkdir}/findallmarkers/markers_vis4cluster{i}/marker_gene_violin_plot.pdf {wkdir}/output/download/violin_cluster{i}.pdf"+"\n"+\
            # f"ln -s {wkdir}/findallmarkers/markers_vis4cluster{i}/marker_gene_violin_plot.png {wkdir}/output/download/violin_cluster{i}.png"+"\n"
            #f"ln -s {wkdir}/findallmarkers/featureplot_cluster{i}.tsv {wkdir}/output/featureplot_cluster{i}.tsv"+"\n"+\
            #f"ln -s {wkdir}/findallmarkers/vlnplot_cluster{i}.tsv {wkdir}/output/vlnplot_cluster{i}.tsv"+"\n"
    #print(cmd_ln)        
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status ###run_cmd use the default function from xiufeng.yang
