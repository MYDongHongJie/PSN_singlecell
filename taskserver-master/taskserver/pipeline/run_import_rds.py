import os, sys
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')
from taskserver.tools.module_cmd import module_cmd


def task_import_rds(input, projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir = workdir
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/import_rds")
    rds_path=input.loc["parameters", "rds_path"]
    environment = input.loc['project', 'environment']
    print(rds_path)
    cmd=f"ln -s {rds_path} {wkdir}/import_rds/"
    print(cmd)
    logger.info(cmd)
    with module_cmd(environment) as p:
        status = p(cmd, projectid, taskid)

    # cmd2=f"Rscript /public/scRNA_works/works/xujingmei/cloud/import_rds/extract_rds_structure.r -i {rds_path} -o {wkdir}/import_rds/"
    # logger.info(cmd2)
    # with module_cmd(environment) as p:
    #     status = p(cmd2, projectid, taskid)

    with open(f"{wkdir}/import_rds/rds_path.tsv", 'w') as fo:
        fo.write("rdspath\n")
        fo.write(f"{rds_path}")
    ###构建一个只有表头的空的output.json.tsv
    with open(f"{wkdir}/import_rds/output.json.tsv", 'w') as fo:
        fo.write("task_type	result_module	input	type	file	title	downloadName	downloadPath\n")
        fo.write("import_rds	data	rds_path.tsv	cell	rds_path.tsv	导入的rds路径		")
    logger.info("step2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
