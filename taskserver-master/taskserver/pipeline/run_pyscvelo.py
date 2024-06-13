import pandas as pd
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger
import glob
logger = getLogger('oe.cloud.sc.qsub')


from taskserver.tools.module_cmd import module_cmd

def find_loom(wkdir,input):
    logger.info("寻找loom文件")
    os.makedirs(f"{wkdir}/velocity/loom")
    loom = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/loom/loom"
    check_loom = glob.glob(os.path.join(f'{loom}', '*.loom'))
    if check_loom:
        logger.info("loom文件存在，开始处理")
        # 将loom链接到当前目录
        os.symlink(loom,f'{wkdir}/loom/loom')
    else:
        logger.info("loom文件不存在，请检查loom任务")
        sys.exit(1)
    return loom

def task_pyscvelo(input,projectid="项目ID", taskid="任务ID", workdir="运行目录"):
    wkdir = workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/pyscvelo")
    cellmeta=input.loc['base', 'tasks'][0]['filename']
    add_params = ""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {wkdir}/input/{cellmeta} "
    else:
        pass
    #获取rds路径，根据base 的taskid,去读取base的input.json，从base的input.json中获取base的任务类型
    rds_json = pd.read_json(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    seurat_ob=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"
    outdir=f"{wkdir}/pyscvelo"
    #pyscvelo分析参数
    #gtf=f"{input.loc['parameters', 'database']}/genes/genes.gtf" #loom 的参数
    reduction=input.loc["parameters","reduction"]
    if reduction != "umap" and reduction != "tsne":
        logger.info("reduction参数错误，只能是umap或者tsne")
        reduction = "umap"
    groupby=input.loc["parameters", "groupby"]
    assay=input.loc["parameters","assay"]
    if input.loc["parameters","genelist"] != "" :
        genelist=input.loc["parameters","genelist"]
        add_params = add_params + f"--genelist {genelist}"
    if input.loc["parameters", "onlygeneplot"] :
        add_params = add_params + f"--onlygeneplot TRUE"
    environment="/public/dev_scRNA/yfang/dev_fangying"
    yourloomdir = find_loom(wkdir,input)
    add_params = add_params + f" --loom_dir {yourloomdir} "
    if not os.path.exists(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/adata_with_scvelo.h5ad"):
        logger.info("未检测到adata_with_scvelo.h5ad,即将生成adata_with_scvelo.h5ad,较耗时！")
        cmd=f"source /opt/softwares/oe_packages/script/scvelo/bash.sh && " \
            f"/opt/softwares/oe_packages/script/scvelo/sctool " \
            f" -i {seurat_ob} "\
            f" -f rds "\
            f" --assay {assay} "\
            f" -o {outdir} "\
            f" pyscvelo "\
            f" --groupby {groupby}  "\
            f" --reduction {reduction}  "\
            f" {add_params}"
    else:
        logger.info("检测到adata_with_scvelo.h5ad,将基于此次的结果进行分析")
        h5ad=f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/adata_with_scvelo.h5ad"  
        cmd=f"source /opt/softwares/oe_packages/script/scvelo/bash.sh && " \
            f"/opt/softwares/oe_packages/script/scvelo/sctool " \
            f" -i {h5ad} " \
            f" -o {outdir} " \
            f" --groupby {groupby}  "\
            f" --reduction {reduction}  "\
            f" {add_params}"

    logger.info(cmd)
    logger.info("pyscvelo:执行分析")
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)

    ## 2.链接rds
    logger.info("链接rds") 
    if not os.path.exists(f'{outdir}/adata_with_scvelo.h5ad'):
        cmd2= f"ln -s {h5ad}  {outdir}/adata_with_scvelo.h5ad"
        with module_cmd(environment) as p:
            status = p(cmd2, projectid, taskid)
        return status
    else:
        logger.info("存在adata_with_scvelo.h5ad")

    ## 3.生成output.json
    logger.info("生成output.json")
    df = pd.DataFrame(columns=['task_type', 'result_module', 'input', 'type', 'file', 'title','downloadName', 'downloadPath'])
    df.loc[len(df.index)] = ["pyscvelo","data","latent_time.tsv","tsv","latent_time.tsv","latent_time","",""]
    df.loc[len(df.index)] = ["pyscvelo","diagram","paga_velocity-graph.png","image","paga_velocity-graph.png","paga_velocity-graph","",""]  
    df.loc[len(df.index)] = ["pyscvelo","diagram",f"Top-likelihood_genes_alonetime_groupby_{groupby}.png","image",f"Top-likelihood_genes_alonetime_groupby_{groupby}.png",f"Top-likelihood_genes_alonetime_groupby_{groupby}","",""]
    df.loc[len(df.index)] = ["pyscvelo","diagram",f"Top-likelihood_genes_groupby_{groupby}.png","image",f"Top-likelihood_genes_groupby_{groupby}.png",f"Top-likelihood_genes_groupby_{groupby}","",""]
    df.loc[len(df.index)] = ["pyscvelo","diagram",f"proportions_spliced_unspliced_counts_groupby_{groupby}.png","image",f"proportions_spliced_unspliced_counts_groupby_{groupby}.png",f"proportions_spliced_unspliced_counts_groupby_{groupby}","",""]
    
    df.to_csv(f"{workdir}/pyscvelo/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')

    logger.info("构建output/download/分析结果")
    if not os.path.exists(f'{wkdir}/output/download'):
        os.makedirs(f'{wkdir}/output/download')

    resultList = ["paga_velocity-graph","Speed_and_coherence",f"Top-likelihood_genes_alonetime_groupby_{groupby}",  f"Top-likelihood_genes_groupby_{groupby}",f"proportions_spliced_unspliced_counts_groupby_{groupby}"] 
    for i in resultList:
        cmd_cp = f"convert  -verbose -density 500 -trim  {wkdir}/pyscvelo/{i}.pdf  -quality 100  -flatten {wkdir}/pyscvelo/{i}.png" + "\n" + \
                 f"ln -s {wkdir}/pyscvelo/{i}.png {wkdir}/output/"
        with module_cmd(environment) as p:
            status = p(cmd_cp, projectid, taskid)
    
    cmd3= f"ln -s {wkdir}/pyscvelo/latent_time.xls {wkdir}/output/latent_time.tsv"+"\n" +\
          f"ln -s {wkdir}/pyscvelo/*{{pdf,png,svg,xls}}  {wkdir}/output/download/"+"\n" +\
          f"ln -s /public/scRNA_works/Documents/scVelo分析说明.docx  {wkdir}/output/download/"+"\n" +\
          f"ln -s {wkdir}/pyscvelo/{{clustersspecific_velocity_genes,clustersspecific_top_gene}} {wkdir}/output/download/" 
          
    with module_cmd(environment) as p:
        status = p(cmd3, projectid, taskid)
    return status
