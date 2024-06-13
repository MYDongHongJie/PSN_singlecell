import os, sys,re
import yaml
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
ROOT_DIR = os.environ.get('SC_ROOT_DIR')

def task_report(input, projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/report")
    account = input.loc["parameters", "account"]
    password=input.loc["parameters", "password"]
    environment=input.loc['project', 'environment']
    base_taskid=input.loc["base", "tasks"][0]["taskId"]
    base_wd=f'{ROOT_DIR}/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json=pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_workflow=base_json.loc["workflow","taskInSeries"]
    os.makedirs(f"{wkdir}/report/cellranger")
    os.makedirs(f"{wkdir}/report/create")
    os.makedirs(f"{wkdir}/report/Count_QC")
    os.makedirs(f"{wkdir}/report/Clustering")
    os.makedirs(f"{wkdir}/report/Marker")
    os.makedirs(f"{wkdir}/report/Diffexp")
    os.makedirs(f"{wkdir}/report/Diffexp/ppi")
    os.makedirs(f"{wkdir}/report/Reference_celltype")
    os.makedirs(f"{wkdir}/report/rds")
    cmd=""
    Sample_Num=""
    for item in base_workflow:
        taskseries_id = item["id"]
        if item["type"]=="beforeQC" or item["type"] == "beforeQC_txy" :
            beforeQC_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            ###读取beforeQC任务的input.json路径
            beforeQC_json=pd.read_json(f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}/input/input.json', orient="index", dtype={"id": 'str'})
            Sample_Num=len(beforeQC_json.loc["base","tasks"])
            #beforeQC_json=pd.read_json(f'input_beforeqc.json', orient="index", dtype={"id": 'str'})##本地测试
            ###找到beforeQC任务的base cellranger任务，并将cellranger结果cp到report目录
            for sample in beforeQC_json.loc["base","tasks"]:
                cellranger_task_id=sample['taskId']
                cellranger_res_temp=f'{ROOT_DIR}/20{cellranger_task_id[0:2]}/{cellranger_task_id[2:4]}/{cellranger_task_id[0:10]}/{cellranger_task_id}/output/download'
                cmd=cmd + f'cp -r {cellranger_res_temp}/*  {wkdir}/report/cellranger/ \n'
            ###将beforeQC任务结果中的seurat.h5seurat cp到report/create文件夹
            cmd = cmd + f'cp -r {beforeQC_folder}/beforeQ*/seurat.h5seurat {wkdir}/report/create/ \n'
            ###将beforeQC任务结果中的beforeQC文件夹内容cp到report/Count_QC文件夹
            cmd = cmd + f'cp -r {beforeQC_folder}/output/download/beforeQ*/* {wkdir}/report/Count_QC/ \n'
            ###将beforeQC任务结果中的aggr文件夹cp到report/cellranger文件夹
            if os.path.exists(f'{beforeQC_folder}/output/download/aggr'):
                cmd = cmd + f'cp -r {beforeQC_folder}/output/download/aggr {wkdir}/report/cellranger/ \n'
        if item["type"] == "cellranger_aggr":
            cellranger_aggr_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            ###将beforeQC任务结果中的aggr文件夹cp到report/cellranger文件夹
            if os.path.exists(f'{cellranger_aggr_folder}/output/download/aggr'):
                cmd = cmd + f'cp -r {cellranger_aggr_folder}/output/download/aggr {wkdir}/report/cellranger/ \n'
        if item["type"] == "QC":
            QC_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {QC_folder}/output/download/* {wkdir}/report/Count_QC/ \n'
        if item["type"] == "bclust":
            bclust_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {bclust_folder}/bclust/data_ob_v3.rds {wkdir}/report/rds/ \n'
            cmd = cmd + f'cp -r {bclust_folder}/output/download/* {wkdir}/report/Clustering/ \n'
        if item["type"] == "summarize":
            summarize_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {summarize_folder}/output/download/* {wkdir}/report/Clustering/visualize_cluster_by_clusters/ \n'
        if item["type"]=="correlation":
            correlation_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            os.makedirs(f'{wkdir}/report/Clustering/clusters_correlation')
            cmd = cmd + f'cp -r {correlation_folder}/output/download/* {wkdir}/report/Clustering/clusters_correlation/ \n'
        if item["type"]=="findallmarkers":
            findallmarkers_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {findallmarkers_folder}/output/download/* {wkdir}/report/Marker/ \n'
        if item["type"]=="diffexp":
            diffexp_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {diffexp_folder}/output/download/* {wkdir}/report/Diffexp/ \n'
        if item["type"]=="enrichment":
            enrichment_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {enrichment_folder}/output/download/enrichment {wkdir}/report/Diffexp/ \n'
        if item["type"]=="PPI":
            PPI_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {PPI_folder}/output/download/* {wkdir}/report/Diffexp/ppi/ \n'
        if item["type"]=="celltyping":
            celltyping_folder = f'{ROOT_DIR}/20{taskseries_id[0:2]}/{taskseries_id[2:4]}/{taskseries_id[0:10]}/{taskseries_id}'
            cmd = cmd + f'cp -r {celltyping_folder}/output/download/* {wkdir}/report/Reference_celltype/ \n'
    if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        print("from /opt/softwares/oe_packages/report/scripts/")
        cmd=cmd + f'cp -r  /opt/softwares/oe_packages/report/scripts/  /opt/softwares/oe_packages/report/config/ {wkdir}/report/ \n'
    else:
        print("from /public/scRNA_works/works/xujingmei/cloud/report_cloud/test/report/scripts/")
        cmd = cmd + f'cp -r /public/scRNA_works/works/xujingmei/cloud/report_cloud/test/report/scripts/  /public/scRNA_works/works/xujingmei/cloud/report/config/ {wkdir}/report/ \n'
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)
    with open(f'{wkdir}/report/config/config.yaml','r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
    config['report']['Library'] = input.loc["parameters", "Library"]
    config['report']['Task_Num'] = input.loc["parameters", "Task_Num"]
    # config['report']['Project_Num'] =input.loc["parameters","Project_Num"] #项目号改为从任务单号截取获取
    config['report']['Project_Num'] =input.loc["parameters", "Task_Num"].split("-")[0]
    config['report']['Customer'] = input.loc["parameters", "Customer"]
    config['report']['Sales'] = input.loc["parameters", "Sales"]
    config['report']['Executor'] = input.loc["parameters", "Executor"]
    Species=input.loc['project',"species"]+"-"+input.loc['project',"sampletype"]
    config['report']['Species'] = Species
    #Sample_Num=len(input.loc['project',"samples"])
    config['report']['Sample_Num'] = Sample_Num
    ###获取Chemistry参数
    beforeQC_folder_id=[x for x in base_workflow if x['type']=="beforeQC_txy"][0]["id"]
    beforeQC_folder_json = pd.read_json(f'{ROOT_DIR}/20{beforeQC_folder_id[0:2]}/{beforeQC_folder_id[2:4]}/{beforeQC_folder_id[0:10]}/{beforeQC_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    cellranger_id= beforeQC_folder_json.loc["base", "tasks"][0]['taskId']
    cellranger_json = pd.read_json(f'{ROOT_DIR}/20{cellranger_id[0:2]}/{cellranger_id[2:4]}/{cellranger_id[0:10]}/{cellranger_id}/input/input.json',orient="index", dtype={"id": 'str'})
    if cellranger_json.loc["parameters","Chemistry"] =="":
        config['report']['Others'] = "3' scRNA"
    else:
        config['report']['Others'] = cellranger_json.loc["parameters","Chemistry"]
    config['report']['Reference'] = cellranger_json.loc["parameters", "database"]
    ###获取reduct1和reduct2参数
    blust_folder_id=[x for x in base_workflow if x['type']=="bclust"][0]["id"]
    blust_folder_json = pd.read_json(f'{ROOT_DIR}/20{blust_folder_id[0:2]}/{blust_folder_id[2:4]}/{blust_folder_id[0:10]}/{blust_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['report_params']['reduct1'] = blust_folder_json.loc["parameters", "reduct1"]
    config['report_params']['reduct2'] = blust_folder_json.loc["parameters", "reduct2"]
    ###获取rmdoublets_method参数
    qc_folder_id=[x for x in base_workflow if x['type']=="QC"][0]["id"]
    qc_folder_json = pd.read_json(f'{ROOT_DIR}/20{qc_folder_id[0:2]}/{qc_folder_id[2:4]}/{qc_folder_id[0:10]}/{qc_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['report_params']['rmdoublets_methods'] = qc_folder_json.loc["parameters", "rmdoublets_method"]
    ###获取marker gene差异分析方法参数
    findallmarkers_folder_id=[x for x in base_workflow if x['type']=="findallmarkers"][0]["id"]
    findallmarkers_folder_json = pd.read_json(f'{ROOT_DIR}/20{findallmarkers_folder_id[0:2]}/{findallmarkers_folder_id[2:4]}/{findallmarkers_folder_id[0:10]}/{findallmarkers_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['report_params']['marker_test'] = findallmarkers_folder_json.loc["parameters", "marker_method"]
    ####获取差异分析方法
    diffexp_folder_id=[x for x in base_workflow if x['type']=="diffexp"][0]["id"]
    diffexp_folder_json = pd.read_json(f'{ROOT_DIR}/20{diffexp_folder_id[0:2]}/{diffexp_folder_id[2:4]}/{diffexp_folder_id[0:10]}/{diffexp_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    config['report_params']['diffexp_test'] = diffexp_folder_json.loc["parameters", "DEGtest"]
    ###根据项目信息和参数信息更新config.yaml文件
    with open(f'{wkdir}/report/config/config.yaml','w',encoding='utf-8') as yw:
        yaml.dump(config,yw)
    ###根据项目分组比较信息更新diffexp.tsv文件
    fd=pd.read_csv(f'{wkdir}/report/config/diffexp.tsv',sep="\t")
    diffexp_folder_id  = [x for x in base_workflow if x['type'] == "diffexp"][0]["id"]
    diffexp_folder_json = pd.read_json(f'{ROOT_DIR}/20{diffexp_folder_id[0:2]}/{diffexp_folder_id[2:4]}/{diffexp_folder_id[0:10]}/{diffexp_folder_id}/input/input.json', orient="index", dtype={"id": 'str'})
    contrast=diffexp_folder_json.loc["parameters","contrast"]
    line_num=0
    for item in contrast:
        item_temp=item.split(":")
        fd.loc[line_num]=[item_temp[1],item_temp[1],item_temp[2],item_temp[2],item_temp[0]]
        line_num=line_num + 1
    fd.to_csv(f'{wkdir}/report/config/diffexp.tsv',sep="\t",index=False)
    cmd2=f'export CLOUD=True \n ' \
         f'export CLOUDUSER={account} \n'\
         f'export PASSWORD={password} \n'
    cmd2=cmd2+f'cd {wkdir}/report/ \n '\
        f'python {wkdir}/report/scripts/report/scRNA_report.py -i . -r {wkdir}/report/rds/*.rds -c {wkdir}/report/config/config.yaml &> report_stdout.txt \n'
    with module_cmd(environment) as p:
        status = p(cmd2, projectid, taskid)
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
                ###将内网前缀http://192.168.20.201:9999改为公网https://cloud.oebiotech.com/spa#/
                link = link.replace('http://10.51.0.100:9999/', 'https://cloud.oebiotech.com/')
                print(link)
                code=line.split("：")[2].split("(")[0]
            else:
                pass
    with open(f"{wkdir}/report/output.json.tsv",'w') as fo , open(f"{wkdir}/report/link.txt",'w') as fl:
        fl.write(f'欧易云报告系统的链接：{link}\n提取码：{code}\n以上为欧易云报告系统的链接和提取码，提取码在两周内有效，请及时查看并确认无误！')
        fo.write("task_type	result_module	input	type	file	title	downloadName	downloadPath\n")
        fo.write(f"report\tdiagram\t{link}\tlink\t{link}\t欧易云报告\t\t\n")
        fo.write(f"report\tdiagram\tlink.txt\ttext\tlink.txt\t链接和提取码\t\t")

