import os
import sys
import re
import json
import subprocess
import pandas as pd
import yaml
import requests
import time
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd
ROOT_DIR = os.environ.get('SC_ROOT_DIR')
logger = getLogger('oe.cloud.sc.qsub')


def task_email_sor(input, output_cfg, projectid, taskid, workdir):
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
    ROOT_DIR = os.environ.get('SC_ROOT_DIR')
    ## 1.获取filterd.h5seurat路径及参数信息
    env_module = input.loc['project', 'environment']
    projectid = input.loc['project', 'name']
    TaskNum=input.loc['parameters', 'Task_Num']
    customer = input.loc['parameters', 'customer']
    species1 = input.loc['project', 'species']
    species2 = input.loc['project', 'sampletype']
    species = str(species1) + "_" + str(species2)
    taskid = input.loc['task', 'id']
    raw_data_obs_address = ''   ###
    metrics_summary = []
    samples_html = []
    samples_png = []
    obs_folder = []
    for task in input.loc['base', 'tasks']:
        #cellranger_folder = glob(f"{workdir}/../../{task['projectId']}/{task['taskId']}/cellranger/*/")[0]
        base_taskid=task['taskId']
        cellranger_folder = glob(f'{ROOT_DIR}/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}/cellranger/*/')[0]
        metrics_summary_csv = f"{cellranger_folder}/outs/metrics_summary.csv"
        samples_png.append(
            glob(f"{cellranger_folder}/outs/*.png")[0])
        samplename = cellranger_folder.split('/')[-2]
        copyfile(f"{cellranger_folder}/outs/web_summary.html",
                 f'{workdir}/{module}/{samplename}_summary.html')
        samples_html.append(f"{workdir}/{module}/{samplename}_summary.html")
        cmd = f"echo -e \"sampleid\\n{samplename}\" | paste - {metrics_summary_csv} | sed 's/\\t/,/g' - > {workdir}/{module}/{samplename}_metrics_summary.csv "
        with module_cmd(env_module) as p:
            status = p(cmd, projectid, taskid)
        metrics_summary.append(
            f"{workdir}/{module}/{samplename}_metrics_summary.csv")
        obs_folder.append(
            f"obs://scrna-cloud/projects/{task['projectId'][0:4]}/{task['projectId']}/tasks/{task['taskId']}/output/download/")

    metrics_summary_input = " ".join(metrics_summary)
    samples_html_input = ",".join(samples_html)
    samples_png_input = ",".join(samples_png)
    num_of_reads = str(int(input.loc['parameters', 'num_of_reads']))
    include_intron = True if input.loc['parameters', 'include_intron'] else False
    project_manager = input.loc['parameters', 'project_manager']
    laboratory_manager = input.loc['parameters', 'laboratory_manager']
    tolist = input.loc['parameters', 'tolist']
    remark = input.loc['parameters', 'remark']
    generate_QC_report = True if input.loc['parameters', 'generate_QC_report'] else False
    if projectid.count('-') > 1 or '-b' in projectid:
        projectid_nobatch = re.sub('-\d+$|-b\d+$','',projectid)
    else:
        projectid_nobatch = projectid

    date = time.strftime('%Y_%m_%d',time.localtime(time.time()))


    ## 2.获取命令行
    cmd = f"cat {metrics_summary_input} | awk '!a[$0]++' > {workdir}/{module}/metrics_summary_stat.csv "
    qc_report_folder = f"obs://oe-scrna/Analysis_Report/{projectid_nobatch}/{projectid_nobatch}_QC_Report_{date}.tar.gz"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    ##
    py_object = {'report': {'Project_Num':projectid, 'Customer':customer, 'Species':species, 'Task_Num':TaskNum},
                 'email': {'emailto': tolist, 'Project_management': project_manager, 'Laboratory': laboratory_manager, 'remark': remark, 'generate_QC_report':generate_QC_report, 'qc_report_folder':qc_report_folder}, 'cellranger_params': {'include_introns': include_intron, 'raw_data_obs_address': raw_data_obs_address, 'Number_of_Reads': num_of_reads}}
    file = open(f'{workdir}/{module}/config.yaml', 'w', encoding='utf-8')
    yaml.dump(py_object, file)

    cmd_email = f"oesc_email.py --config {workdir}/{module}/config.yaml --stat {workdir}/{module}/metrics_summary_stat.csv --html {samples_html_input} --png {samples_png_input}"

    ## 3.执行分析
    with module_cmd(env_module) as p:
        status = p(cmd_email, projectid, taskid)

    if generate_QC_report:
        account = input.loc["parameters", "account"]
        password = input.loc["parameters", "password"]
        ## 更新config文件
        py_object['report']['Project_Num'] = projectid_nobatch
        py_object['report']['Executor'] = input.loc["parameters", "Executor"]
        py_object['report']['Sample_Num'] = len(obs_folder)
        task_base=input.loc['base', 'tasks'][0] ###上游某个cellranger任务
        base_taskid = task_base['taskId']
        cellranger_input = pd.read_json(f"{ROOT_DIR}/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}/input/input.json", orient="index", dtype={"id": 'str'})
        py_object['cellranger_params']['envmodules'] = cellranger_input.loc['parameters', 'cellranger_version']
        py_object['cellranger_params']['Reference'] = cellranger_input.loc["parameters", "database"]
        with open(f'{workdir}/{module}/config.yaml', 'w', encoding='utf-8') as f:
            yaml.dump(py_object, f)
        ## 创建 QC report
        if not os.path.exists(f'{workdir}/{module}/{projectid}_QC_Report_{date}/1.CellRanger'):
            os.makedirs(f'{workdir}/{module}/{projectid}_QC_Report_{date}/1.CellRanger')
        if not os.path.exists(f'{workdir}/output/download'):
            os.makedirs(f'{workdir}/output/download')
        for file in obs_folder:
            cmd_download = f"obsutil cp  {file}  {workdir}/{module}/{projectid}_QC_Report_{date}/1.CellRanger/  -r -f -flat -vlength -vmd5 \n "
            with module_cmd("obsutil/5.2.12") as p:
                status = p(cmd_download, projectid, taskid)

        cmd_export = f'export CLOUD=True \n ' \
                     f'export CLOUDUSER={account} \n'\
                     f'export PASSWORD={password} \n'
        if not os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
            cmd_export = cmd_export + f'cp -r /opt/softwares/report/scripts/  {workdir}/{module}/ \n'
            cmd_report = cmd_export + f"cd {workdir}/{module}/ \n" \
                                      f"python {workdir}/{module}/scripts/report/Cellranger_report.py "\
                                      f" -i {workdir}/{module}/  -c {workdir}/{module}/config.yaml &> report_stdout.txt \n "
        else:
            cmd_export = cmd_export + f'cp -r /public/scRNA_works/works/xujingmei/cloud/report_cloud/test/report/scripts/ {workdir}/{module}/ \n'
            cmd_report = cmd_export + f"cd {workdir}/{module}/ \n" \
                                      f"module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {workdir}/{module}/scripts/report/Cellranger_report.py "\
                                      f" -i {workdir}/{module}/  -c {workdir}/{module}/config.yaml &> report_stdout.txt \n "
        with module_cmd(env_module) as p:
            status = p(cmd_report, projectid, taskid)
        with open(f"{workdir}/{module}/report_stdout.txt",'r') as f:
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
                    ###将内网前缀http://192.168.20.201:9999改为公网https://cloud.oebiotech.com/
                    link = link.replace('http://192.168.20.201:9999/', 'https://cloud.oebiotech.com/')
                    print(link)
                else:
                    pass
        cmd2 = f"cd {workdir}/{module} && tar -czvhf {projectid}_QC_Report_{date}.tar.gz {projectid}_QC_Report_{date} && rm -rf {workdir}/output/download/* && ln -s {workdir}/{module}/{projectid}_QC_Report_{date}.tar.gz {workdir}/output/download/  && echo -e \"质控结果已反馈至{tolist}，请注意查收！该项目为只质控不分析项目，完整QC report已上传云报告。\"  > {workdir}/{module}/结果.txt "
        with module_cmd("obsutil/5.2.12") as p:
            status = p(cmd2, projectid, taskid)
    else:
        cmd_res = f"echo -e \"质控结果已反馈至{tolist}，请注意查收！\"  > {workdir}/{module}/结果.txt "
        with module_cmd(env_module) as p:
            status = p(cmd_res, projectid, taskid)

    cmd_sor = f"oesc_sor.py  -i {workdir}/{module}/config.yaml  -l {workdir}/{module}/logs/upload_cellranger_sor  -s 'Project, QC, Cpu_Mem' "
    ## 3.执行分析
    with module_cmd(env_module) as p:
        status = p(cmd_sor, projectid, taskid)

    # =================================================================================================================
    # logger.info("step3.2:构建output/download/分析结果")
    # if not os.path.exists(f'{workdir}/output/download'):
    #     os.makedirs(f'{workdir}/output/download')
    # cmd = f"ln -s {workdir}/{module}/* {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*.tsv"
    # with module_cmd(env_module) as p:
    #     status = p(cmd, projectid, taskid)
    # logger.info("step3.3:生成output.json.xls")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    df= pd.DataFrame(columns=output_df.columns)
    df.loc[0] = ["cellranger_email", "diagram", "结果.txt",
                 "text", "结果.txt", "邮件反馈结果", "", ""]
    if generate_QC_report:
        with open(f"{workdir}/{module}/link.txt",'w') as fl:
            fl.write(f'欧易云报告系统的链接：{link}\n提取码：{code}\n以上为欧易云报告系统的链接和提取码，提取码在两周内有效，请及时查看并确认无误！')
            # fo.write("task_type	result_module	input	type	file	title	downloadName	downloadPath\n")
        df.loc[1] = ["cellranger_email","diagram",link,"link",link,"欧易云报告","",""]
        df.loc[2] = ["cellranger_email","diagram","link.txt","text","link.txt","链接和提取码","",""]

    df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',
              index=False, header=True, encoding='utf-8')
    return status
