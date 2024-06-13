# -*-coding:utf-8-*-
from glob import glob
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

def task_volcano(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir  = workdir
    environment = "OESingleCell/2.0.0" #input.loc['project', 'environment']
    os.makedirs(f"{wkdir}/output/download")
    os.makedirs(f"{wkdir}/volcano")

    if input.loc["parameters", "genelist"] != "":
        genemetas = input.loc["parameters", "genelist"]
    else  :
        if input.loc['base', 'tasks'][0]['filename'] != "":
            genemetas =  f"{wkdir}/input/{input.loc['base', 'tasks'][0]['filename']}"
        else :
             genemetas = glob(f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/output/download/*-vs-*all_diffexp_genes*.xls")
    if genemetas == "":
        logger.info("There is no file input, and there are also files under the difference file of the project，无差异文件输入，对应项目的差异文件夹下也无文件")

    pvalue = int(input.loc["parameters", "pvalue"])
    foldchange = int(input.loc["parameters", "foldchange"])
    symbol_topn = f"--symbol_topn {int(input.loc['parameters', 'symbol_topn'])} " if input.loc['parameters', 'symbol_topn' ] != "" else ""
    symbol_gene = f"--symbol_gene {input.loc['parameters', 'symbol_gene']} " if input.loc['parameters', 'symbol_gene'] != "" else ""
    symbol_fc = f"--symbol_fc {int(input.loc['parameters', 'symbol_fc'])} " if input.loc['parameters', 'symbol_fc'] != "" else ""
    if symbol_gene != "" and symbol_topn != "" or symbol_fc !="":
        logger.info("存在多个参数，优先使用symbol_gene")
    for genemeta in genemetas:
        cmd = f" Rscript   /public/scRNA_works/works/donghongjie/project/test/sc_Volcano.r" \
              f"  -i  {genemeta} " \
              f"  -f  {foldchange} " \
              f"  -p  {pvalue} " \
              f"  -o  {wkdir}/volcano "\
              f"{symbol_gene}  {symbol_topn}  {symbol_fc}"
        with module_cmd(environment) as p:
            status = p(cmd, projectid, taskid)


    logger.info("The drawing has been completed, and the indexing of the results is started")

    cmd_ln = f"ln -s {wkdir}/volcano/*.pdf {wkdir}/output/download/" + "\n"+\
             f"ln -s {wkdir}/volcano/*.png  {wkdir}/output/download/"
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status

