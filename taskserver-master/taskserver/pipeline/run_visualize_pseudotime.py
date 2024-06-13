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

def task_visualize_pseudotime(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/visualize_pseudotime")
    #input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{wkdir}/input/cell.tsv"
    cellmeta = f"{wkdir}/input/{input.loc['base', 'tasks'][0]['filename']}"
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    monocle_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/pseudotime/pseudotime_results.rds"
    outdir=f"{wkdir}/visualize_pseudotime"
    environment=input.loc["project","environment"]
    species=input.loc["parameters","species"]
    root=input.loc["parameters","root"]
    show_branch=bool(input.loc["parameters","show_branch"])
    ordering=bool(input.loc["parameters","ordering"])
    genelist=input.loc["parameters","genelist"]
    ordering_gene_topn=int(input.loc["parameters","ordering_gene_topn"])
    heatmap=bool(input.loc["parameters","heatmap"])
    heatmap_clusters=int(input.loc["parameters","heatmap_clusters"])
    heatmap_branchpoint=input.loc["parameters","heatmap_branchpoint"]
    expressplot=bool(input.loc["parameters","expressplot"])
    expressplot_groupby=input.loc["parameters","expressplot_groupby"]
    expressplot_branchpoint=input.loc["parameters","expressplot_branchpoint"]
    trajectoryplot=bool(input.loc["parameters","trajectoryplot"])
    trajectory_groupby=input.loc["parameters","trajectory_groupby"]
    treeplot=bool(input.loc["parameters","treeplot"])
    treeplot_groupby=input.loc["parameters","treeplot_groupby"]
    expressplot_line=bool(input.loc["parameters","expressplot_line"])
    expressplot_line_groupby=input.loc["parameters","expressplot_line_groupby"]
    module_expressplot_line=bool(input.loc["parameters","module_expressplot_line"])
    module_expressplot_line_branchpoint=input.loc["parameters","module_expressplot_line_branchpoint"]

    add_params=""
    if genelist !="":
        ####上传文件对应input.json中的值是obs的全路径
        upload_file=genelist.split("/")[-1]
        #genelist=f"{wkdir}/{upload_file}"
        add_params = add_params + f" --genelist  {wkdir}/input/{upload_file} "
    if expressplot_branchpoint != "":
        add_params = add_params + f" --expressplot_branchpoint  {int(expressplot_branchpoint)} "
    if module_expressplot_line_branchpoint != "":
        add_params = add_params + f" --module_expressplot_line_branchpoint  {int(module_expressplot_line_branchpoint)} "
    if heatmap_branchpoint != "":
        add_params = add_params + f" --heatmap_branchpoint  {int(heatmap_branchpoint)} "
    if root != "":
        add_params = add_params + f" --root  {int(root)} "

    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " /opt/softwares/envs/oesinglecell_v3/bin/Rscript   /opt/softwares/oe_packages/cloud_report_script/scripts_use/run_visualize_pseudotime/visualize_pseudotime_20230215.R"
    else:
        tmp_envs="module purge && module load OESingleCell/2.0.0 && Rscript /public/scRNA_works/works/liuxuan/Yunpingtai/Scripts/visualize_pseudotime_20230215.R"
    cmd=f"{tmp_envs}   "\
        f"  -i {monocle_ob} "\
        f"  -o {outdir} "\
        f"  --species {species}"\
        f"  --show_branch {show_branch}"\
        f"  --ordering_gene_topn {ordering_gene_topn}"\
        f"  --heatmap {heatmap}"\
        f"  --expressplot {expressplot}"\
        f"  --expressplot_groupby {expressplot_groupby}"\
        f"  --trajectoryplot {trajectoryplot}"\
        f"  --trajectory_groupby {trajectory_groupby}"\
        f"  --treeplot {treeplot}"\
        f"  --treeplot_groupby {treeplot_groupby}"\
        f"  --expressplot_line {expressplot_line}"\
        f"  --expressplot_line_groupby {expressplot_line_groupby}"\
        f"  --module_expressplot_line {module_expressplot_line}"\
        f"  --ordering {ordering}"\
        f"  {add_params} " \
    #   f"  --expressplot_branchpoint {expressplot_branchpoint}"\

    logger.info(cmd)
    logger.info("step1: 开始拟时序可视化绘图")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    #####需要构建output文件夹内容，主要包括1) output/download/中供提供给客户的结果文件; 2) output/ 下提供供云平台前端展示的绘图数据、表格数据、统计数据
    link_rds=f"{outdir}/pseudotime_results.rds"
    os.symlink(monocle_ob, link_rds)
    if expressplot_line is True:
        os.makedirs( f"{wkdir}/output/download/expressplot_line")
        os.system('cp -s %s/visualize_pseudotime/expressplot_module* %s/output/download/expressplot_line'%(wkdir,wkdir))
    if expressplot is True:
        os.makedirs( f"{wkdir}/output/download/pseudotime_gene_express")
        os.system('cp -s %s/visualize_pseudotime/pseudotime_gene_express* %s/output/download/pseudotime_gene_express'%(wkdir,wkdir))
    if heatmap is True:
        os.makedirs( f"{wkdir}/output/download/pseudotime_heatmap")
        os.system('cp -s %s/visualize_pseudotime/pseudotime_heatmap* %s/output/download/pseudotime_heatmap'%(wkdir,wkdir))
    if treeplot is True:
        os.makedirs( f"{wkdir}/output/download/pseudotime_treeplot")
        os.system('cp -s %s/visualize_pseudotime/pseudotime_treeplot* %s/output/download/pseudotime_treeplot'%(wkdir,wkdir))
    if trajectoryplot is True:
        os.makedirs( f"{wkdir}/output/download/trajectory_gene")
        os.system('cp -s %s/visualize_pseudotime/trajectory_gene* %s/output/download/trajectory_gene'%(wkdir,wkdir))
    if f"{root}" != "":
        os.makedirs( f"{wkdir}/output/download/cell_trajectory")
        os.system('cp -s %s/visualize_pseudotime/cell_trajectory* %s/output/download/cell_trajectory'%(wkdir,wkdir))
    df = pd.DataFrame(columns = ['task_type', 'result_module', 'input','type','file','title','downloadName','downloadPath'])
    os.system('cp -s %s/visualize_pseudotime/*.png %s/output/'%(wkdir,wkdir))
    DownloadPath = f"{wkdir}/output/download/"
    Download_Dir = os.listdir(DownloadPath)
    for j in Download_Dir:
        files = os.listdir("%s/output/download/%s"%(wkdir,j))
        #print(files)
        for i in range(len(files)):
            file = files[i].split(".")[0]
            format = files[i].split(".")[1]
            if format == "png":
                df.loc[len(df)] = ['visualize_pseudotime','diagram',files[i],'image',files[i], files[i].split(".")[0],files[i].split(".")[0],"/download/%s/%s.png,/download/%s/%s.pdf"%(j,files[i].split(".")[0],j,files[i].split(".")[0])]
            elif format == "tsv":
                df.loc[len(df)] = ['visualize_pseudotime','data',files[i],'tsv',files[i], files[i].split(".")[0],files[i].split(".")[0],"/download/%s/%s"%(j,files[i])]
            #elif format == "csv":
            #    df.loc[len(df)] = ['visualize_pseudotime','data',files[i],'image',files[i], files[i].split(".")[0],'visualize_pseudotime',"/download/%s/%s"%(j,files[i])]
    df.to_csv(f"{outdir}/output.json.tsv",sep="\t",index=None)
#            elif format == "png":
#                df.loc[len(df)] = ['visualize_pseudotime','diagram','','visualize_pseudotime','', files[i].split(".")[0],'visualize_pseudotime',"/download/%s/%s"%(j,files[i])]

