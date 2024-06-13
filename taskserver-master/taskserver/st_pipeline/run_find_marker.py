import os, sys
import json
from glob import glob
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_find_marker(input, output_cfg, projectid, taskid, workdir):
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
    #input = pd.read_json(f'{workdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{workdir}/input/cell.tsv"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    seurat_ob = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    #outdir=f"{workdir}/{module}"
    environment="OESingleCell/v_3.0.0_visium_produce"
    database=input.loc["parameters","database"]
    assay=input.loc["parameters","assay"]
    library_size=input.loc["parameters","library_size"] 
    #Significance=input.loc["parameters","marker_Significance"]
    if input.loc["parameters", "Significance_threshold"] != "":
        Significance_threshold =  input.loc["parameters", "Significance_threshold"]
    else:
        Significance_threshold = 0.05
    if library_size == "cytassist":
        crop = "--crop TRUE"
    else:
        crop = " "
    topnum =int(input.loc["parameters", "topnum"])
    #strict=input.loc["parameters","strict"]
    method=input.loc["parameters","method"]
    groupby=input.loc["parameters","groupby"]
    #add_params_sign=""
     # if Significance=="pvalue":
         # add_params_sign=add_params_sign+f"  -p {Significance_threshold} "
     # elif Significance=="FDR":
         # add_params_sign=add_params_sign+f"  -q {Significance_threshold} "
    # if input.loc["parameters", "min_pct1"] != "":
        # min_pct1 =  float(input.loc["parameters", "min_pct1"])
    # else:
        # min_pct1 = 0.5
    # if input.loc["parameters", "max_pct2"] != "":
        # max_pct2 =  float(input.loc["parameters", "max_pct2"])
    # else:
        # max_pct2 = 0.5
    #pct_fold =int(input.loc["parameters", "pct_fold"])
    #avg_log2FC =int(input.loc["parameters", "avg_log2FC"])

    cmd1=f"sctool "\
        f"  -i {seurat_ob} "\
        f"  -f rds "\
        f"  -o {workdir}/{module} "\
        f"  --image TRUE "\
        f"  --assay {assay} "\
        f"{add_params} " \
        f" findallmarkers "\
        f"  -e {method} "\
        f"  -n {groupby} "\
        f"--anno {database}annotation/gene_annotation.xls"\
        f"  -t 0.5 "\
        f"  -T 0.5 "\
        f"  -c 2 "\
        f"  -N 10 "\
        f"  -k 1 "\
        f"  -p {Significance_threshold}"\
        f"  --strict FALSE "


    logger.info(cmd1)
    logger.info("step1.1: 开始计算marker")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd1, projectid, taskid)
    marker_files = glob(
        f"{workdir}/{module}/*.xls")
    for file in marker_files:
        print(file)
        os.rename(file, file.replace('.xls', '.tsv'))
##############################################################################################################################
    reduct2=input.loc["parameters","reduct2"]
    topn_marker=f"{workdir}/{module}/top{topnum}_markers_for_each_cluster.tsv"
    topby=input.loc["parameters","topby"]
    #crop =input.loc["parameters","crop"]
    #spatialpointsize=input.loc["parameters","spatialpointsize"]
    cmd2=f"scVis " \
        f"  -i {seurat_ob} "\
        f"  -f rds "\
        f"  --assay {assay} "\
        f"  -o {workdir}/{module} "\
        f"  --image TRUE "\
        f"{add_params} " \
        f"  markervis   "\
        f"  -l {topn_marker} "\
        f"  --topn  {topnum}  "\
        f"  --topby  {topby}  "\
        f"  -m vlnplot,featureplot,heatmap "\
        f"  -s 0 "\
        f"  --reduct {reduct2} "\
        f"  -p 1.2 "\
        f"  -a 1 "\
        f"  {crop} "

    logger.info(cmd2)
    logger.info("step2: top10 marker 画小提琴图、散点图、热图")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd2, projectid, taskid)
    logger.info("step3 生成output.json.tsv。")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t",na_filter=False) #, encoding='gb2312')
    df= pd.DataFrame(columns=output_df.columns)
    df.loc[len(df.index)] = ["find_marker", "data", "all_markers_for_each_cluster_anno.tsv",
                  "gene", "all_markers_for_each_cluster_anno.tsv", "all", "", ""]
    df.loc[len(df.index)] = ["find_marker", "data", "top10_markers_for_each_cluster_anno.tsv",
                  "gene", "top10_markers_for_each_cluster_anno.tsv", "top10", "", ""]
    df.loc[len(df.index)] = ["find_marker", "diagram", "topmarker_gene_heatmap.png",
                  "image", "topmarker_gene_heatmap.png", "topmarker_gene_heatmap", "", ""]
    df.to_csv(f"{workdir}/{module}/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')
    logger.info("step1.3:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    # cmd_ln1= f"ln -s {workdir}/{module}/{{all_markers_for_each_cluster_anno.tsv,top10_markers_for_each_cluster_anno.tsv}}  {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*.tsv"+"\n"
    cmd_ln = f"ln -s {workdir}/{module}/all_markers_for_each_cluster_anno.tsv {workdir}/output/download/all_markers_for_each_cluster_anno.xls" +"\n"+\
        f"ln -s {workdir}/{module}/top{topnum}_markers_for_each_cluster_anno.tsv {workdir}/output/download/top{topnum}_markers_for_each_cluster_anno.xls" +"\n"+\
        f"ln -s {workdir}/{module}/{{topmarker_gene_heatmap.png,topmarker_gene_heatmap.pdf}} {workdir}/output/"+"\n"+\
        f"ln -s {workdir}/{module}/{{topmarker_gene_heatmap.png,topmarker_gene_heatmap.pdf}} {workdir}/output/download/"+"\n"
        #f"ln -s {workdir}/findallmarkers/topmarker_gene_heatmap.tsv {workdir}/output/"+"\n"
    #print(cmd_ln)
    cmd_ln = cmd_ln + "\n" + \
            f"ln -s {workdir}/{module}/markers_vis4cluster* {workdir}/output/download/"+"\n"

    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status
