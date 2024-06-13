#######空转出报告流程##########
import os, sys,re
import yaml
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')
from taskserver.tools.module_cmd import module_cmd


def task_report(input, projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/report")
    account = input.loc["parameters", "account"]
    password=input.loc["parameters", "password"]
    environment="OESingleCell/v_3.0.0_visium_produce"
    base_taskid=input.loc["base", "tasks"][0]["taskId"]
    ####测试/public/cloud_scRNA替换为/public/cloud_scRNA


    base_wd=f'{wkdir}/../../{base_taskid[0:10]}/{base_taskid}'
    base_json=pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_workflow=base_json.loc["workflow","taskInSeries"]
    os.makedirs(f"{wkdir}/report/spaceranger")
    os.makedirs(f"{wkdir}/report/sctools")
    os.makedirs(f"{wkdir}/report/count_qc")
    os.makedirs(f"{wkdir}/report/cluster_seurat")
    os.makedirs(f"{wkdir}/report/marker")
    os.makedirs(f"{wkdir}/report/diffexp")
    os.makedirs(f"{wkdir}/report/diffexp/ppi")
    os.makedirs(f"{wkdir}/report/celltype")
    #os.makedirs(f"{wkdir}/report/rds")
    cmd=""
    Sample_Num=""
    for item in base_workflow:
        taskseries_id = item["id"]
        # if item["type"]=="beforeQC":
        #     beforeQC_folder = f'/public/cloud_scRNA/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
        #     ###读取beforeQC任务的input.json路径
        #     beforeQC_json=pd.read_json(f'/public/cloud_scRNA/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}/input/input.json', orient="index", dtype={"id": 'str'})
        #     Sample_Num=len(beforeQC_json.loc["base","tasks"])
        #     #beforeQC_json=pd.read_json(f'input_beforeqc.json', orient="index", dtype={"id": 'str'})##本地测试
        #     ###找到beforeQC任务的base spaceranger任务，并将spaceranger结果cp到report目录
        #     for sample in beforeQC_json.loc["base","tasks"]:
        #         spaceranger_task_id=sample['taskId']
        #         spaceranger_res_temp=f'/public/cloud_scRNA/20{spaceranger_task_id[0:2]}/{spaceranger_task_id[2:4]}/{spaceranger_task_id[0:10]}/{spaceranger_task_id}/output/download'
        #         cmd=cmd + f'cp -r {spaceranger_res_temp}/*  {wkdir}/report/spaceranger/ \n'
        #     ###将beforeQC任务结果中的seurat.h5seurat cp到report/sctools文件夹
        #     cmd = cmd + f'cp -r {beforeQC_folder}/beforeQC/seurat.h5seurat {wkdir}/report/sctools/ \n'
        #     ###将beforeQC任务结果中的beforeQC文件夹内容cp到report/count_qc文件夹
        #     cmd = cmd + f'cp -r {beforeQC_folder}/output/download/beforeQC/* {wkdir}/report/count_qc/ \n'
        #     ###将beforeQC任务结果中的aggr文件夹cp到report/spaceranger文件夹
        #     if os.path.exists(f'{beforeQC_folder}/output/download/aggr'):
        #         cmd = cmd + f'cp -r {beforeQC_folder}/output/download/aggr {wkdir}/report/spaceranger/ \n'
        if item["type"] == "Count_QC":
            QC_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            ###读取QC任务的input.json路径,然后将获取base的样本数
            QC_json=pd.read_json(f'{QC_folder}/input/input.json', orient="index", dtype={"id": 'str'})
            Sample_Num=len(QC_json.loc["base","tasks"])
            ###将spaceranger文件夹需要的结果软链过去
            cmd = cmd + f'ln -s {QC_folder}/Count_QC/spaceranger/* {wkdir}/report/spaceranger/ \n'
            cmd = cmd + f'ln -s {QC_folder}/Count_QC/aggr {wkdir}/report/spaceranger/ \n'
            ###将sctools文件夹需要的create rds软链过去
            cmd = cmd + f'ln -s {QC_folder}/Count_QC/create/*.rds  {wkdir}/report/sctools/ \n'
            ###将count_qc 文件夹需要的结果软链过去
            cmd = cmd + f'ln -s {QC_folder}/Count_QC/*.rds  {wkdir}/report/count_qc/ \n'
            cmd = cmd + f'ln -s {QC_folder}/output/download/*  {wkdir}/report/count_qc/ && rm -rf {wkdir}/report/count_qc/aggr \n'
        if item["type"] == "bclust":
            bclust_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            #cmd = cmd + f'cp -r {bclust_folder}/bclust/data_ob_v3.rds {wkdir}/report/rds/ \n'
            cmd = cmd + f'ln -s {bclust_folder}/bclust/* {wkdir}/report/cluster_seurat/ &&rm -rf {wkdir}/report/cluster_seurat/*.tsv *.xls \n'
        if item["type"] == "summarize":
            summarize_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'ln -s  {summarize_folder}/summarize/visualize_cluster_by_clusters {wkdir}/report/cluster_seurat/ \n'
        # if item["type"]=="correlation":
        #     correlation_folder = f'/public/cloud_scRNA/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
        #     os.makedirs(f'{wkdir}/report/cluster_seurat/clusters_correlation')
        #     cmd = cmd + f'cp -r {correlation_folder}/output/download/* {wkdir}/report/cluster_seurat/clusters_correlation/ \n'
        if item["type"]=="find_marker":
            findallmarkers_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'ln -s {findallmarkers_folder}/output/download/* {wkdir}/report/marker/ \n'
        if item["type"]=="diffexp":
            diffexp_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'ln -s {diffexp_folder}/diffexp/* {wkdir}/report/diffexp/ && rm -rf {wkdir}/report/diffexp/output.json.tsv &&rename tsv xls {wkdir}/report/diffexp/*\n'
        if item["type"]=="enrichment":
            enrichment_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'ln -s {enrichment_folder}/output/download/* {wkdir}/report/diffexp/ \n'
        if item["type"]=="PPI":
            PPI_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'ln -s {PPI_folder}/PPI/* {wkdir}/report/diffexp/ppi/ && rm -rf {wkdir}/report/diffexp/ppi/output.json.tsv && rename tsv xls {wkdir}/report/diffexp/ppi/* \n'
        if item["type"]=="celltyping":
            celltyping_folder = f'{wkdir}/../../{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'ln -s {celltyping_folder}/celltyping/* {wkdir}/report/celltype/ && rm -rf {wkdir}/report/celltype/output.json.tsv \n'
    cmd=cmd + f'cp -r /public/scRNA_works/works/guochy/ST_taskserver/ST_snakemake_template/report/scripts /public/scRNA_works/works/guochy/ST_taskserver/ST_snakemake_template/st_pipeline/config {wkdir}/report/ \n'
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)
    with open(f'{wkdir}/report/config/config.yaml','r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
    config['report']['Project_Num'] =input.loc["parameters","Project_Num"]
    config['report']['Task_Num'] = input.loc["parameters", "Task_Num"]
    config['report']['Customer'] = input.loc["parameters", "Customer"]
    config['report']['Sales'] = input.loc["parameters", "Sales"]
    config['report']['Executor'] = input.loc["parameters", "Executor"]
    Species=input.loc['project',"species"]+"-"+input.loc['project',"sampletype"]
    config['report']['Species'] = Species
    #Sample_Num=len(input.loc['project',"samples"])
    config['report']['Sample_Num'] = Sample_Num
    config['report']['Project_Path'] = f"{workdir}/report/"
    config['module']['spaceranger_report'] = False
    # ###获取Chemistry参数
    # beforeQC_folder_id=[x for x in base_workflow if x['type']=="beforeQC"][0]["id"]
    # beforeQC_folder_json = pd.read_json(f'/public/cloud_scRNA/20{beforeQC_folder_id[0:2]}/{beforeQC_folder_id[2:4]}/{beforeQC_folder_id[0:10]}/{beforeQC_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    # spaceranger_id= beforeQC_folder_json.loc["base", "tasks"][0]['taskId']
    # spaceranger_json = pd.read_json(f'/public/cloud_scRNA/20{spaceranger_id[0:2]}/{spaceranger_id[2:4]}/{spaceranger_id[0:10]}/{spaceranger_id}/input/input.json',orient="index", dtype={"id": 'str'})
    # if spaceranger_json.loc["parameters","Chemistry"] =="":
    #     config['report']['Others'] = "3' scRNA"
    # else:
    #     config['report']['Others'] = spaceranger_json.loc["parameters","Chemistry"]
    # config['report']['Reference'] = spaceranger_json.loc["parameters", "database"]
    ###获取blust 参数
    blust_folder_id=[x for x in base_workflow if x['type']=="bclust"][0]["id"]
    blust_folder_json = pd.read_json(f'{wkdir}/../../{blust_folder_id[0:10]}/{blust_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['params']['reduct1_method'] = blust_folder_json.loc["parameters", "reduct1"]
    config['params']['reduct2_method'] = blust_folder_json.loc["parameters", "reduct2"]
    config['params']['clusteringuse'] = blust_folder_json.loc["parameters", "clusteringuse"]
    config['params']['resolution'] = float(blust_folder_json.loc["parameters", "resolution"])
    config['params']['batchid'] = blust_folder_json.loc["parameters", "batchid"]
    #config['params']['component'] = blust_folder_json.loc["parameters", "component"] if blust_folder_json.loc["parameters", "component"] else "NULL"
    ###获取QC参数
    qc_folder_id=[x for x in base_workflow if x['type']=="Count_QC"][0]["id"]
    qc_folder_json = pd.read_json(f'{wkdir}/../../{qc_folder_id[0:10]}/{qc_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['database']['reference']=qc_folder_json.loc["parameters", "database"]
    config['params']['library_type']=qc_folder_json.loc["parameters", "library_size"]
    config['params']['sct_split'] = qc_folder_json.loc["parameters", "sct_split"]
    config['params']['nvfeatures'] = float(qc_folder_json.loc["parameters", "nvfeatures"])
    config['params']['normmeth'] = qc_folder_json.loc["parameters", "normmeth"]
    config['params']['sct_split'] = qc_folder_json.loc["parameters", "sct_split"]
    ###获取marker gene差异分析任务的参数
    findallmarkers_folder_id=[x for x in base_workflow if x['type']=="find_marker"][0]["id"]
    findallmarkers_folder_json = pd.read_json(f'{wkdir}/../../{findallmarkers_folder_id[0:10]}/{findallmarkers_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    #config['params']['marker_test'] = findallmarkers_folder_json.loc["parameters", "method"]
    ####获取差异分析任务的参数
    diffexp_folder_id=[x for x in base_workflow if x['type']=="diffexp"][0]["id"]
    diffexp_folder_json = pd.read_json(f'{wkdir}/../../{diffexp_folder_id[0:10]}/{diffexp_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['params']['foldchange'] =float(diffexp_folder_json.loc["parameters", "foldchange"])
    config['params']['pvalue'] = float(diffexp_folder_json.loc["parameters", "DEG_Significance_threshold"])
    #####获取PPI参数
    PPI_folder_id=[x for x in base_workflow if x['type']=="PPI"][0]["id"]
    PPI_folder_json = pd.read_json(f'{wkdir}/../../{PPI_folder_id[0:10]}/{PPI_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['params']['speciesid'] = float(PPI_folder_json.loc["parameters", "species_ppi"])
    ###根据项目信息和参数信息更新config.yaml文件
    with open(f'{wkdir}/report/config/config.yaml','w',encoding='utf-8') as yw:
        yaml.dump(config,yw)
    ###根据项目分组比较信息更新diff_group.tsv文件
    fd=pd.read_csv(f'{wkdir}/report/config/diff_group.tsv',sep="\t")
    diffexp_folder_id  = [x for x in base_workflow if x['type'] == "diffexp"][0]["id"]
    diffexp_folder_json = pd.read_json(f'{wkdir}/../../{diffexp_folder_id[0:10]}/{diffexp_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    contrast=diffexp_folder_json.loc["parameters","contrast"]
    line_num=0
    for item in contrast:
        item_temp=item.split(":")
        fd.loc[line_num]=[item_temp[1],item_temp[1],item_temp[2],item_temp[2],item_temp[0]]
        line_num=line_num + 1
    fd.to_csv(f'{wkdir}/report/config/diff_group.tsv',sep="\t",index=False)

    # cmd2=cmd2+f'cd {wkdir}/report/ \n '\
    #     f'/public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {wkdir}/report/scripts/report/scRNA_report.py -i . -c {wkdir}/report/config/config.yaml &> report_stdout.txt \n'

###########################
    cmd2=f"module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {wkdir}/report/scripts/step1.report_files_copy.py  -i {wkdir}/report/ -c {wkdir}/report/config/config.yaml &> report_stdout.txt \n"
    with module_cmd(environment) as p:
        status = p(cmd2, projectid, taskid)
    report_name=os.listdir(f"{wkdir}/report/report/")
    cmd3=f'export CLOUD=True \n ' \
        f'export CLOUDUSER={account} \n'\
        f'export PASSWORD={password} \n'
    cmd3=cmd3+f"module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {wkdir}/report/scripts/step2.report_html_convert.py -c {wkdir}/report/config/config.yaml  -d {wkdir}/report/report/{report_name[0]}/ &> {wkdir}/report/report_stdout.txt \n"
    with module_cmd(environment) as p:
        status = p(cmd3, projectid, taskid)
    with open(f"{wkdir}/report/report_stdout.txt",'r') as f:
        txt=f.readlines()
        link=""
        code=""
        for line in txt:
            print(line)
            #line=log.iloc[line_num,]
            if "欧易生信云报告访问链接：" in line:
                #print(line)
                link=line.split("：")[1].split("，")[0]
                code=line.split("：")[2].split("(")[0]
                link = link.replace("http://192.168.20.201:9999/","https://cloud.oebiotech.com/")
            else:
                pass
    with open(f"{wkdir}/report/output.json.tsv",'w') as fo , open(f"{wkdir}/report/link.txt",'w') as fl:
        fl.write(f'欧易云报告系统的链接：{link}\n提取码：{code}\n以上为欧易云报告系统的链接和提取码，提取码在两周内有效，请及时查看并确认无误！')
        fo.write("task_type	result_module	input	type	file	title	downloadName	downloadPath\n")
        fo.write(f"report\tdiagram\t{link}\tlink\t{link}\t欧易云报告\t\t\n")
        fo.write(f"report\tdiagram\tlink.txt\ttext\tlink.txt\t链接和提取码\t\t")
        
