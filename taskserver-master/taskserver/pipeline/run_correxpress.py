import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_correxpress(input, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """

    wkdir  = workdir
    environment = input.loc['project', 'environment']
    #os.makedirs(f"{wkdir}/output/download")
    if not os.path.exists(f'{workdir}/correxpress'):
        os.makedirs(f'{workdir}/correxpress')
    #os.makedirs(f"{wkdir}/correxpress")

    if input.loc["parameters", "file1.txt"] != "":
        file1 = f"{wkdir}/input/{input.loc['parameters', 'file1.txt']}"
    else :
        logger.info("There is no file input.")

    if input.loc["parameters", "file2.txt"] != "":
        file2 = f"{wkdir}/input/{input.loc['parameters', 'file2.txt']}"

    cor = float(input.loc["parameters", "cor"])
    pvalue = float(input.loc["parameters", "pvalue"])
    method = input.loc["parameters", "method"]
    
    if input.loc["parameters", "file2.txt"] == "":
        cmd = f"oebio correxpress -cor {cor} -p {pvalue} -m {method} {file1} -o {workdir}/correxpress/correxpress.tsv && rename {workdir}/correxpress/correxpress_cor*p*.tsv {workdir}/correxpress/correxpress.tsv {workdir}/correxpress/correxpress_cor*p*.tsv"
    else :
        cmd = f"oebio correxpress -cor {cor} -p {pvalue} -m {method} {file1} {file2} -o {workdir}/correxpress/correxpress.tsv && rename {workdir}/correxpress/correxpress_cor*p*.tsv {workdir}/correxpress/correxpress.tsv {workdir}/correxpress/correxpress_cor*p*.tsv"
        
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)
        

    ## 3.执行分析
    logger.info("The analysis has been completed, and the indexing of the results is started")
    #logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    # cmd_ln = f"ln -s {wkdir}/correxpress/*.xlsx {wkdir}/output/download/ && rename *.xlsx correxpress.tsv {workdir}/output/download/*.xlsx"

    # with module_cmd(environment) as p:
        # status = p(cmd_ln, projectid, taskid)
        
    # =================================================================================================================
    
    # cmd = f"mv {workdir}/{module}/cor_file1_cor_*_cor*p*.xlsx {workdir}/output/ && ln -s {workdir}/{module}/* {workdir}/output/download/ && rename xlsx tsv {workdir}/output/download/*.xlsx"
    # with module_cmd(env_module) as p:
        # status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    # logger.info("step3.3:生成output.json.xls")
              
    # generate output.json.tsv so that the diagrams can be uploaded to obs
    d = {"task_type":["correxpress"],"result_module":["data"], 
    "input":["correxpress.tsv"], 
    "type":["tsv"], 
    "file":["correxpress.tsv"],
    "title":["correxpress"], 
    "downloadName":["correxpress"], 
    "downloadPath":["download/correxpress.tsv"]}
    df = pd.DataFrame(d)
    df.to_csv(f"{wkdir}/correxpress/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")
    cmd_ln = f"ln -s {wkdir}/correxpress/correxpress.tsv {wkdir}/output/download/ "
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status

