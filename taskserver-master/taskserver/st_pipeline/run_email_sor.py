import os
import sys
import re
import json
import subprocess
import pandas as pd
import yaml
import requests
import time
import csv
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')


def task_email_sor(input,  projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    ## 1.获取filterd.h5seurat路径及参数信息
    #env_module = input.loc['project', 'environment']
    env_module = "OESingleCell/v_3.0.0_visium_produce"
    projectid = input.loc['project', 'name']
    customer = input.loc['parameters', 'customer']
    species = input.loc['project', 'species']
    taskid = input.loc['task', 'id']
    raw_data_obs_address = ''   ###
    metrics_summary = []
    samples_html = []
    samples_png = []
    obs_folder = []
    sample_meta=[]
    ###构建包含所有样本的spaceranger结果目录
    if not os.path.exists(f'{workdir}/{module}/spaceranger'):
        os.makedirs(f'{workdir}/{module}/spaceranger')
    for task in input.loc['base', 'tasks']:
        cmd_ln=f"ln -s {workdir}/../../{task['projectId']}/{task['taskId']}/spaceranger/* {workdir}/{module}/spaceranger/ && rm {workdir}/{module}/spaceranger/output.json.tsv"
        logger.info("将spaceranger结果软链到当期目录")
        with module_cmd(env_module) as p:
            status = p(cmd_ln, projectid, taskid)
        spaceranger_folder = glob(f"{workdir}/../../{task['projectId']}/{task['taskId']}/spaceranger/*/")[0]
        metrics_summary_csv = f"{spaceranger_folder}/outs/metrics_summary.csv"
        samples_png.append(
            glob(f"{spaceranger_folder}/outs/*.png")[0])
        samplename = spaceranger_folder.split('/')[-2]
        copyfile(f"{spaceranger_folder}/outs/web_summary.html",
                 f'{workdir}/{module}/{samplename}_summary.html')
        samples_html.append(f"{workdir}/{module}/{samplename}_summary.html")
        cmd = f"echo -e \"sampleid\\n{samplename}\" | paste - {metrics_summary_csv} | sed 's/\\t/,/g' - > {workdir}/{module}/{samplename}_metrics_summary.csv "
        with module_cmd(env_module) as p:
            status = p(cmd, projectid, taskid)
        metrics_summary.append(
            f"{workdir}/{module}/{samplename}_metrics_summary.csv")
        #obs_folder.append(
        #    f"obs://scrna-cloud/projects/{task['projectId'][0:4]}/{task['projectId']}/tasks/{task['taskId']}/output/download/")
        ####读取每个sample.csv的metadata文件
        base_folder=f"{workdir}/../../{task['projectId']}/{task['taskId']}/config/samples.csv"
        sample_metadata_temp=pd.read_csv(base_folder,header=0)
        if sample_meta:
            sample_meta.append(sample_metadata_temp.loc[0,].tolist())
        else:
            sample_meta.append(sample_metadata_temp.columns.tolist())
            sample_meta.append(sample_metadata_temp.loc[0,].tolist())
    ###将合并后的samples.csv metadata输出文件
    with open(f'{workdir}/{module}/samples.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for row in sample_meta:
            writer.writerow(row)

    metrics_summary_input = " ".join(metrics_summary)
    samples_html_input = ",".join(samples_html)
    samples_png_input = ",".join(samples_png)
    num_of_reads = str(int(input.loc['parameters', 'num_of_reads']))
    #include_intron = True if input.loc['parameters', 'include_intron'] else False
    project_manager = input.loc['parameters', 'project_manager']
    laboratory_manager = input.loc['parameters', 'laboratory_manager']
    tolist = input.loc['parameters', 'tolist']
    password = input.loc['parameters', 'password']
    #remark = input.loc['parameters', 'remark']
    generate_QC_report = True if input.loc['parameters', 'generate_QC_report'] else False
    if projectid.count('-') > 1 or '-b' in projectid:
        projectid_nobatch = re.sub('-\d+$|-b\d+$','',projectid)
    else:
        projectid_nobatch = projectid
    date = time.strftime('%Y%m%d',time.localtime(time.time()))

    ## 2.获取命令行
    cmd = f"cat {metrics_summary_input} | awk '!a[$0]++' > {workdir}/{module}/metrics_summary_stat.csv "
    #qc_report_folder = f"obs://oe-scrna/Analysis_Report/{projectid_nobatch}/{projectid}_QC_Report{date}.tar.gz"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)

    sample_num=len(input.loc['base','tasks'])
    # py_object = {'report': {'Project_Num':projectid, 'Customer':customer, 'Species':species,'Project_Type': '单细胞空间','Sample_Num':sample_num,'Rawdata_Path':{'image_obs_address':raw_data_obs_address,'image_obs_address':""}},'email': {'email_to': tolist, 'Project_management': project_manager, 'Laboratory': laboratory_manager, 'remark': remark, 'generate_QC_report':generate_QC_report, 'qc_report_folder':qc_report_folder}, 'spaceranger_params': {'include_introns': include_intron, 'raw_data_obs_address': raw_data_obs_address, 'Number_of_Reads': num_of_reads},'params':{'email_to': tolist,'raw_reads': num_of_reads}}
    # file = open(f'{workdir}/{module}/config.yaml', 'w', encoding='utf-8')
    # yaml.dump(py_object, file)

    cmd_cp=f'cp -r /public/scRNA_works/works/guochy/ST_taskserver/ST_snakemake_template/report/scripts /public/scRNA_works/works/guochy/ST_taskserver/ST_snakemake_template/st_pipeline/config {workdir}/{module}/'
    logger.info("将标准config.yaml 复制过来")
    with module_cmd(env_module) as p:
        status = p(cmd_cp, projectid, taskid)

    with open(f'{workdir}/{module}/config/config.yaml','r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
    config['report']['Project_Num'] =projectid_nobatch
    config['report']['Task_Num'] =projectid
    config['params']['library_type']=input.loc["parameters", "library_type"]
    config['module']['spaceranger_report'] =generate_QC_report
    config['module']['report'] = False
    config['report']['Customer'] = input.loc["parameters", "customer"]
    config['report']['Sales'] = input.loc["parameters", "Sales"]
    config['report']['Executor'] = input.loc["parameters", "Executor"]
    Species=input.loc['project',"species"]+"-"+input.loc['project',"sampletype"]
    config['report']['Species'] = Species
    config['report']['Project_Type'] ='单细胞空间'
    config['report']['Sample_Num'] = sample_num
    config['report']['Rawdata_Path']['image_obs_address'] = ""###暂不上传该项内容，设为空
    config['report']['Rawdata_Path']['image_obs_address'] = ""###暂不上传该项内容，设为空
    # config['email']['email_to'] = tolist
    # config['email']['Project_management'] = project_manager
    # config['email']['Laboratory'] = laboratory_manager
    config['params']['email_to'] = tolist
    config['params']['raw_reads'] = num_of_reads
    ###根据项目信息和参数信息更新config.yaml文件
    with open(f'{workdir}/{module}/config/config.yaml','w',encoding='utf-8') as yw:
        yaml.dump(config,yw)

    cmd_summary=f"Rscript /public/scRNA_works/works/guochy/ST_taskserver/scripts/qc_summary.R -s  {workdir}/{module}/samples.csv -i  {workdir}/{module}/spaceranger -o  {workdir}/{module} -c  {workdir}/{module}/config/config.yaml  "
    logger.info("生成summary网页报告结果and发送邮件")
    with module_cmd(env_module) as p:
        status = p(cmd_summary, projectid, taskid)
    


    cmd_sor=f"module purge && /data/software/conda_envs/OESingleCell_3.0.0_visium/bin/python /public/scRNA_works/works/guochy/ST_taskserver/scripts/sor_upload_visium.py -i {workdir}/{module}/config/config.yaml -m {workdir}/{module}/samples.csv -l {workdir}/{module}/ -s 'Project,Sample,QC,Database,Software,Cpu_Mem'"
    logger.info("SOR 上传")
    with module_cmd(env_module) as p:
        status = p(cmd_sor, projectid, taskid)



    if generate_QC_report:
        #    with module_cmd("obsutil/5.2.12") as p:
        #        status = p(cmd1, projectid, taskid)
        #cmd2 = f"cd {workdir}/{module} && tar -czvhf {projectid}_QC_Report{date}.tar.gz {projectid}_QC_Report{date} && obsutil cp {projectid}_QC_Report{date}.tar.gz obs://oe-scrna/Analysis_Report/{projectid_nobatch}/  -r -f -flat -vlength -vmd5  && rm -rf {workdir}/output/download/* && ln -s {workdir}/{module}/{projectid}_QC_Report{date}.tar.gz {workdir}/output/download/  && echo -e \"质控结果已反馈至{tolist}，请注意查收！\n质控报告已上传至 {qc_report_folder}，请及时下载。\"  > {workdir}/{module}/结果.txt "
        cmd2 = f"module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {workdir}/{module}/scripts/step1.report_files_copy.py  -i {workdir}/{module}/ -c {workdir}/{module}/config/config.yaml"
        with module_cmd("obsutil/5.2.12") as p:
            status = p(cmd2, projectid, taskid)
        report_name=os.listdir(f"{workdir}/{module}/report/")
        cmd3=f'export CLOUD=True \n ' \
        f'export CLOUDUSER={tolist} \n'\
        f'export PASSWORD={password} \n'
        cmd3=cmd3+f"module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {workdir}/{module}/scripts/step2.spaceranger_html_convert.py -c {workdir}/{module}/config/config.yaml  -d {workdir}/{module}/report/{report_name[0]}/ &> {workdir}/{module}/report/report_stdout.txt \n"
        with module_cmd(env_module) as p:
            status = p(cmd3, projectid, taskid)
        with open(f"{workdir}/{module}/report/report_stdout.txt",'r') as f:
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
        with open(f"{workdir}/{module}/output.json.tsv",'w') as fo , open(f"{workdir}/{module}/link.txt",'w') as fl:
            fl.write(f'欧易云报告系统的链接：{link}\n提取码：{code}\n以上为欧易云报告系统的链接和提取码，提取码在两周内有效，请及时查看并确认无误！')
            fo.write("task_type	result_module	input	type	file	title	downloadName	downloadPath\n")
            fo.write(f"spaceranger_email\tdiagram\t{link}\tlink\t{link}\t欧易云报告\t\t\n")
            fo.write(f"spaceranger_email\tdiagram\tlink.txt\ttext\tlink.txt\t链接和提取码\t\t")
    else:
        cmd_res = f"echo -e \"质控结果已反馈至{tolist}，请注意查收！\"  > {workdir}/{module}/link.txt "
        with module_cmd(env_module) as p:
            status = p(cmd_res, projectid, taskid)
        df= pd.DataFrame(columns=['task_type','result_module','input','type','file','title','downloadName',	'downloadPath'])
        df.loc[0] = ["spaceranger_email", "diagram", "link.txt",
                    "text", "link.txt", "邮件反馈结果", "", ""]
        df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',
                index=False, header=True, encoding='utf-8')

    # =================================================================================================================
    # logger.info("step3.2:构建output/download/分析结果")
    # if not os.path.exists(f'{workdir}/output/download'):
    #     os.makedirs(f'{workdir}/output/download')
    # cmd = f"ln -s {workdir}/{module}/* {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*.tsv"
    # with module_cmd(env_module) as p:
    #     status = p(cmd, projectid, taskid)
    # logger.info("step3.3:生成output.json.xls")

    return status
