import os
import sys
import re
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')


def set_gsea_params(input, projectid, taskid, workdir):
    '''
    从 input.json 文件获取 R 脚本的命令行参数
    '''

    #module = input.loc['task', 'type']
    global environment
    environment = "OESingleCell/2.0.0"

    if not os.path.exists(f'{workdir}/gsea'):
        os.makedirs(f'{workdir}/gsea')

    if not os.path.exists(f"{workdir}/output/download"):
        os.makedirs(f"{workdir}/output/download")

    # gsea脚本参数字典
    gsea_param = {}

    # 前端 获取物种，人、小鼠 单选框（暂定）
    if not input.loc["parameters", "species"] in ("Human" or "Mouse"):
        sys.exit()
    else:
        species = input.loc["parameters", "species"]

    # 获取input路径
    # 从降维聚类的task获取seurat对象，默认使用所有细胞的Seurat对象
    #
    ##读取rds路径
    rds_json = pd.read_json(
        f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json",
        orient='index', dtype={'id': 'str'})
    input_dir = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task', 'type']}/data_ob_v3.rds"
    gsea_param["--input"] = input_dir

    # outdir
    outdir = f'{workdir}/gsea'
    if not os.path.exists(f"{outdir}"):
        os.makedirs(f"{outdir}")
    gsea_param["--outdir"] = outdir

    # gmt
    global gmt
    if input.loc["parameters", "gmt"] == "":
        if species == "Human":
            gmt_dir = ["/data/database/GSEA_gmt/human/c2.cp.kegg.v7.2.symbols.gmt",
                    "/data/database/GSEA_gmt/human/c5.go.bp.v7.2.symbols.gmt"]
            gmt = ",".join(gmt_dir)
        if species == "Mouse":
            gmt_dir = ["/data/database/GSEA_gmt/mouse/c2.cp.kegg.v7.2.symbols.gmt",
                    "/data/database/GSEA_gmt/mouse/c5.go.bp.v7.2.symbols.gmt"]
            gmt = ",".join(gmt_dir)
    else:
        gmt_dir = f"{workdir}/input/" + input.loc["parameters", "gmt"]
        gmt = "".join(gmt_dir)

    gsea_param["--gmt"] = gmt

    # ident2use
    if input.loc["parameters", "ident2use"] != "":
        ident2use = input.loc["parameters", "ident2use"]
        gsea_param["--ident2use"] = ident2use

    # which_cells
    if input.loc["parameters", "which_cells"] != "":
        which_cells = input.loc["parameters", "which_cells"]
        gsea_param["--which_cells"] = which_cells

    # contrast group:xx:xx
    global contrast
    if input.loc["parameters", "contrast"] != "":
        contrast = input.loc["parameters", "contrast"]
        contrast = "".join(contrast)
        gsea_param["--contrast"] = contrast

    else:
        sys.exit()

    return gsea_param


def get_gsea_cmd_list(para):
    cmd_list = []

    for key in para:
        param = f'{key} {para[key]}'
        cmd_list.append(param)
    print(cmd_list)
    return cmd_list


def generate_output_tsv(workdir):
    factor, case, control = contrast.split(":")
    data = {}
    gmt_file = pd.read_csv(gmt, header=None, names=["termid"])
    prefix = (gmt.split("/")[-1]).split(".gmt")[0]
    downloadPath = [f"download/{case}-vs-{control}/{case}-vs-{control}-None/gsea.gsea.gene_set.report.csv"]
    input_name = ["gsea.gsea.gene_set.report.tsv"]
    # if (case != "all") and (control != "all"):
    #     for i in range(len(gmt_file)):
    #         num=i+1
    #         termid=gmt_file.loc[i,"termid"].split("\t")[0]
    #         downloadPath += [f"download/{case}-vs-{control}/{case}-vs-{control}-None/{prefix}/{num}_{termid}.gsea.png",
    #                          f"download/{case}-vs-{control}/{case}-vs-{control}-None/{prefix}/{num}_{termid}.gsea.pdf"]
    #         input_name += [f"{num}_{termid}.gsea.png",
    #                        f"{num}_{termid}.gsea.pdf"]

    data["task_type"] = ["gsea"]#["gsea"] * (len(gmt_file) +1) #* len(database)
    data["result_module"] =["data"] #(["diagram"] * len(gmt_file) + ["data"] * 1) #* len(database)
    data["input"] = input_name
    data["type"] ="tsv" #["png", "pdf"] * len(gmt_file) +["xls"] #* len(database)
    data["file"] = "gsea.gsea.gene_set.report.tsv"  #downloadPath
    data["title"] ="gsea.gsea.gene_set.report"  #[""] *len(gmt_file) + ["gsea.gsea.gene_set.report"]
    data["downloadName"] = input_name
    data["downloadPath"] = downloadPath

    df = pd.DataFrame(data)
    df.to_csv(f"{workdir}/gsea/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")
    return data


def make_links(workdir, projectid, taskid):
    factor, case, control = contrast.split(":")
    #for i in database:
    cmd_ln = f"ln -s {workdir}/gsea/GSEA/{case}-vs-{control}/{case}-vs-{control}-None/  {workdir}/output/download/"
    cmd_ln2 = f"ln -s {workdir}/gsea/GSEA/{case}-vs-{control}/{case}-vs-{control}-None/*/gsea.gsea.gene_set.report.csv {workdir}/output/gsea.gsea.gene_set.report.tsv"
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    with module_cmd(environment) as p:
        status2 = p(cmd_ln2, projectid, taskid)
    return status


def task_gsea(input, projectid, taskid, workdir):
    logger.info("task_gsea: 1.获得脚本的参数字典")
    gsea_param = set_gsea_params(input, projectid, taskid, workdir)

    logger.info("task_gsea: 2.获得脚本的命令行列表")
    gsea_cmd_line = " ".join(get_gsea_cmd_list(gsea_param))
    gsea_cmd = f"Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/GSEA_sc.R {gsea_cmd_line}"

    # 进入工作路径开始分析
    set_workdir_cmd = f"cd {workdir}/gsea"
    logger.info("task_gsea: 3.开始分析")
    start_gsea_analysis = f"{set_workdir_cmd} && module purge && module load {environment} && {gsea_cmd}"
    with module_cmd(environment) as p:
        status = p(start_gsea_analysis, projectid, taskid)
    logger.info("task_gsea: Down!")

    logger.info("task_gsea: 4.生成 output.json.tsv")
    generate_output_tsv(workdir)

    logger.info("task_gsea: 5.生成软链接")
    make_links(workdir, projectid, taskid)

