import os, sys
import json
import subprocess
import pandas as pd
import requests
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_monocle3(input, projectid, taskid, workdir):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/monocle3")

    #rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    #seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"
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

    outdir = f"{wkdir}/monocle3"

    add_params = ""
    groupby = input.loc["parameters", "groupby"] #设置分类的列
    if input.loc["parameters", "dimension"] != "":
        add_params = add_params + f" --dimension  {input.loc['parameters', 'dimension']} " #使用原来的umap方法,tsne也适用
    if input.loc["parameters", "rootby"] != "":
        add_params = add_params + f" --rootby  {input.loc['parameters', 'rootby']} " #选择细胞起源的clusters或者new_celltype等
        if input.loc["parameters", "root"] != "":
            if input.loc["parameters", "rootby"] == "clusters" :
                add_params = add_params + f" --root  {int(input.loc['parameters', 'root'])} " #选择细胞起源对应的簇
            else:
                add_params = add_params + f" --root  {input.loc['parameters', 'root']}"
        else:
            logger.info("细胞起源簇和细胞起源必须同时使用")
    if input.loc["parameters", "vismethod"] == []:
        pass
    else:
        vismethod=','.join(input.loc['parameters', 'vismethod'])
        add_params = add_params + f" --vismethod  {vismethod} " #选择进行的下游分析

    if input.loc["parameters", "extraGene"] != "":
        extraGene = f"{wkdir}/input/{input.loc['parameters', 'extraGene']}"
        add_params = add_params + f" --extraGene  {extraGene} " #添加外源基因的文件
    module = input.loc['task', 'type']
    #env_module = input.loc['project', 'environment']

    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        tmp_envs = " /opt/softwares/envs/oesinglecell_v4/bin/Rscript /opt/softwares/oe_packages/cloud_report_script/scripts_use/run_monocle3/monocle3.R "
    else:
        tmp_envs="module purge && module load OESingleCell/3.0.0 && Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/monocle3.R"
    #environment = "OESingleCell/3.0.0"
    print(tmp_envs)
    env_module = input.loc['project', 'environment']
    cmd = f"{tmp_envs}  " \
          f"    -i {seurat_ob}"  \
          f"    --groupby {groupby}"  \
          f"    -o {outdir}" \
          f"    {add_params}"

    logger.info(cmd)
    logger.info("step1: 开始monocle3分析")
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)

#################################################################################
    #生成output.json.tsv
    output_file=os.listdir(f"{outdir}")
    df=pd.DataFrame(columns = ['task_type', 'result_module', 'input','type','file','title','downloadName','downloadPath'])
    for file in output_file :
        if file.split(".")[-1]=='png':
            file_name='.'.join(file.split(".")[0:-1])
            df.loc[len(df)]=["monocle3","diagram",f"{file_name}.png","image",f"{file_name}.png",f"{file_name}",f"{file_name}",f"download/{file_name}.pdf,download/{file_name}.png"]
            df=df.sort_values(by=["downloadName"],ascending=True)
        else:
            pass
    df.to_csv(f"{outdir}/output.json.tsv",sep="\t",index=None)

    logger.info("step2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    #filelist = f"*.pdf,*.png,*.docx"
    #cmd3 = f"ln -s {workdir}/{module}/{{{filelist}}} {workdir}/output/download/ "
    cmd3=f"find {workdir}/{module}/ -type f ! -name '*.rds' ! -name '*.json.tsv' -exec ln -s '{{}}' {workdir}/output/download/ \;"
    with module_cmd(env_module) as p:
        status = p(cmd3, projectid, taskid)
    cmd4=f"ln -s {workdir}/{module}/*.png {workdir}/output/"
    with module_cmd(env_module) as p:
        status = p(cmd4, projectid, taskid)
    return status
