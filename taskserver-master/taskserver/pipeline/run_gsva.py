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

def set_gsva_params(input, projectid, taskid, workdir):
    '''
    从 input.json 文件获取 R 脚本的命令行参数
    '''
    global environment
    environment = "OESingleCell/2.0.0"

    if not os.path.exists(f'{workdir}/gsva'):
        os.makedirs(f'{workdir}/gsva')

    if not os.path.exists(f"{workdir}/gsva/output/download"):
        os.makedirs(f"{workdir}/output/download")


    #========================================== 1.GSVA_enrich.R 参数 ==========================================
    #判断上一步是否对数据进行了删选，若筛选，则读取cellmeta文件
    # cellmeta=input.loc['base', 'tasks'][0]['filename']
    # add_params = ""
    # if cellmeta != "":
    #     add_params = add_params + f" --metadata  {wkdir}/input/{cellmeta} "
    # gsva打分脚本参数字典
    gsva_enrich_param = {}
    #print(1)
    # 前端 获取物种，人和小鼠，下拉选择
    species = input.loc["parameters", "species"]
    # 前端 获取database，多选框 默认不传参，使用默认的kegg和bp（暂定）
    global database

    database = input.loc["parameters", "database"]



    ##读取rds路径
    rds_json = pd.read_json(f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient='index', dtype={'id': 'str'})
    input1=f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/{rds_json.loc['task','type']}/data_ob_v3.rds"

    # gmt and outdir
    gmt_outdir_dict = {} # gmt:outdir
    if input.loc["parameters", "gmt"] == "":
        if species == "Human" and (database == "kegg") :
            gmt1 = "/data/database/GSVA_gmt/hg38.kegg.gmt"
            outdir1 = f'{workdir}/gsva/kegg'
            gmt_outdir_dict[gmt1] = outdir1
            if not os.path.exists(f"{outdir1}"):
                os.makedirs(f"{outdir1}")
        if species == "Mouse" and database == "kegg":
            gmt2 = "/data/database/GSVA_gmt/mm10.kegg.gmt"
            outdir2 = f'{workdir}/gsva/kegg'
            gmt_outdir_dict[gmt2] = outdir2
        if species == "Human" and database == "bp":
            gmt3 = "/data/database/GSEA_gmt/human/c5.go.bp.v7.2.symbols.gmt"
            outdir3 = f'{workdir}/gsva/bp'
            if not os.path.exists(f"{outdir3}"):
                os.makedirs(f"{outdir3}")
            gmt_outdir_dict[gmt3] = outdir3
        if species == "Mouse" and database == "bp":
            gmt4 = "/data/database/GSEA_gmt/mouse/gene_go_bp.backgroud.gmt"
            outdir4 = f'{workdir}/gsva/bp'
            gmt_outdir_dict[gmt4] = outdir4
        # if species == "人" and ("MsigDB" in database):
        #     pass
        # if species == "小鼠" and ("MsigDB" in database):
        #     pass
        gsva_enrich_param["gmt_outdir"] = gmt_outdir_dict
    else:
        gmt = input.loc["parameters", "gmt"]
        outdir = f'{workdir}/gsva/{gmt}'
        if not os.path.exists(f"{outdir}"):
            os.makedirs(f"{outdir}")
        gmt_outdir_dict[gmt] = outdir
        gsva_enrich_param["gmt_outdir"] = gmt_outdir_dict
        print(gsva_enrich_param["gmt_outdir"])

    print(gsva_enrich_param,species,database)
    # abs_rank
    # if input.loc["parameters", "abs_rank"] != "":
    #     abs_rank = toupper(input.loc["parameters", "abs_rank"])
    # else:
    #     abs_rank = "F"
    # gsva_enrich_param["--abs_rank"] = abs_rank

    # which_group1
    if input.loc["parameters", "which_group1"] != "":
        which_group1 = input.loc["parameters", "which_group1"]
        gsva_enrich_param["--WHICH_GROUP"] = input.loc["parameters", "which_group1"]

    # which_cells1
    if input.loc["parameters", "which_cells1"] != "":
        which_cells1 = input.loc["parameters", "which_cells1"]
        gsva_enrich_param["--WHICH_CELLS"] = input.loc["parameters", "which_cells1"]

    gsva_enrich_param["-i"] = input1

    #===================================== 2.GSVA_pathway_diffxp.R 参数 =========================================
    # gsva画图脚本参数字典
    gsva_diffxp_param = {}

    # input_score_table = outdir/GSVA_enrichment_results.xls

    # 获取input2路径，如果不指定，默认为input1
    if input.loc["parameters", "input2"] != "":
        input2_dir = input.loc["parameters", "input2"]
    else:
        input2_dir = input1

    gsva_diffxp_param["--seurat"] = input2_dir

    # # 获取input2后缀
    # input2_suffix = input2_dir.rsplit(".", 1)[1]

    # # 如果传入input2后缀为h5seurat，则转换成rds
    # if input2_suffix == "rds":
    #     pass
    # else:
    #     pass

    # # metadata2
    # if input.loc["parameters", "metadata2"] != "":
    #     metadata2 = input.loc["parameters", "metadata2"]
    #     gsva_diffxp_param["--metadata"] = metadata2

    # contrast group:xx:xx
    global contrast
    if input.loc["parameters", "contrast"] != "":
        contrast = input.loc["parameters", "contrast"]
        contrast = " ".join(contrast)
    else:
        sys.exit()

    gsva_diffxp_param["--contrast"] = contrast

    # pval
    if input.loc["parameters", "pval"] != "":
        pval = float(input.loc["parameters", "pval"])
    else:
        pval = float(0.05)

    gsva_diffxp_param["--pval"] = pval

    # fdr
    if input.loc["parameters", "fdr"] != "":
        fdr = float(input.loc["parameters", "fdr"])
        gsva_diffxp_param["--fdr"] = fdr

    # topn
    if input.loc["parameters", "topn"] != "":
        topn = int(input.loc["parameters", "topn"])
    else:
        topn = int(10)

    gsva_diffxp_param["--topn"] = topn

    # matrix
    if input.loc["parameters", "matrix"] != "":
        matrix = input.loc["parameters", "matrix"]
        gsva_diffxp_param["--matrix"] = matrix

    # groupby
    if input.loc["parameters", "groupby"] != "":
        groupby = input.loc["parameters", "groupby"]
        gsva_diffxp_param["--groupby"] = groupby

    # termid_table
    if input.loc["parameters", "termid"] != "":
        termid_table = input.loc["parameters", "termid"]
        gsva_diffxp_param["--termid"] = termid_table

    # which_group2
    if input.loc["parameters", "which_group2"] != "":
        which_group2 = input.loc["parameters", "which_group2"]
        gsva_diffxp_param["--WHICH_GROUP"] = input.loc["parameters", "which_group2"]

    # which_cells2
    if input.loc["parameters", "which_cells2"] != "":
        which_cells2 = input.loc["parameters", "which_cells2"]
        gsva_diffxp_param["--WHICH_CELLS"] = input.loc["parameters", "which_cells2"]

    return gsva_enrich_param, gsva_diffxp_param


def get_gsva_cmd_dict(enrich, diffxp):
    enrich_cmd_dict = {} # "gmt1":"enri_cmd"
    diffxp_cmd_dict = {} # "gmt1":"diff_cmd"
    for gmt in enrich["gmt_outdir"].keys():
        gmt_dir = gmt
        out_dir = enrich["gmt_outdir"][gmt]
        input_score_table = enrich["gmt_outdir"][gmt]+"/GSVA_enrichment_results.xls"
        enrich["--gmt"] = gmt_dir
        enrich["--OUTDIR"] = out_dir
        enrich["-j"] = 20
        diffxp["--outdir"] = out_dir
        diffxp["--input"] = input_score_table
        enrich_cmd_list = []
        diffxp_cmd_list = []
        for key in diffxp:
            param1 = f'{key} {diffxp[key]}'
            diffxp_cmd_list.append(param1)
            diffxp_cmd_dict[gmt] = diffxp_cmd_list
        for key in enrich:
            if key != "gmt_outdir":
                param2 = f'{key} {enrich[key]}'
                enrich_cmd_list.append(param2)
                enrich_cmd_dict[gmt] = enrich_cmd_list
        diffxp_cmd_dict[gmt] = ' '.join(diffxp_cmd_list)
        enrich_cmd_dict[gmt] = ' '.join(enrich_cmd_list)
    return enrich_cmd_dict, diffxp_cmd_dict


def generate_output_tsv(workdir):
    factor, case, control = contrast.split(":")
    downloadPath = []
    input_name = []
    data = {}
    if (case != "all") and (control != "all"):
        #for i in database:
        downloadPath += [f"download/diffexp_genesets_GSVA_score4group_{case}-vs-{control}_barplot.png",
                          f"download/diffexp_genesets_GSVA_score4group_{case}-vs-{control}.xls"]
        input_name += [f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}_barplot.png",
                          f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}.xls"]

        data["task_type"] = ["gsva"] * 2
        data["result_module"] = (["diagram"] * 1 + ["data"] * 1)
        data["input"] = [f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}_barplot.png",
                          f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}.tsv"]
        data["type"] = ["image", "tsv"]
        data["file"] = [f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}_barplot.png",
                          f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}.tsv"]
        data["title"] = [f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}_barplot",
                          f"diffexp_genesets_GSVA_score4group_{case}-vs-{control}"]
        data["downloadPath"] = downloadPath
    else:
        #for i in database:
        downloadPath += ["download/top10_term_t_value.xls",
                          "download/top10_gsva_term.png",
                          "download/GSVA_top10_results.xls"]
        input_name += ["top10_term_t_value.xls",
                          "top10_gsva_term.png",
                          "GSVA_top10_results.xls"]

        data["task_type"] = ["gsva"] * 3
        data["result_module"] = (["data"] + ["diagram"]  + ["data"] )
        data["input"] = ["top10_term_t_value.tsv",
                          "top10_gsva_term.png",
                          "GSVA_top10_results.tsv"]
        data["type"] = ["tsv", "image", "tsv"]
        data["file"] = ["top10_term_t_value.tsv",
                          "top10_gsva_term.png",
                          "GSVA_top10_results.tsv"]
        data["title"] = ["top10_term_t_value",
                          "top10_gsva_term",
                          "GSVA_top10_results"
                          ]
        data["downloadName"] = input_name
        data["downloadPath"] = downloadPath

    df = pd.DataFrame(data)
    df.to_csv(f"{workdir}/gsva/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")
    return data


# def make_links(workdir, projectid, taskid):
#     #for i in database:
#     upload_lst = list(set(generate_output_tsv(workdir)["title"]))
#     for f in upload_lst:
#         cmd_ln = f"ln -s {workdir}/gsva/{database}/{f} {workdir}/output/download/{database}_{f}"
#         #cmd_ln2 = f"ln -s {workdir}/gsva/{database}/{f} {workdir}/{database}_{f}"
#         #cmd_ln3 = f"ln -s {workdir}/gsva/{database}/*.png {workdir}/output/download/"
#         cmd_ln4 = f"ln -s {workdir}/gsva/{database}/*.png {workdir}/output/"
#         with module_cmd(environment) as p:
#             status = p(cmd_ln, projectid, taskid)
#         #with module_cmd(environment) as p:
#         #    status2 = p(cmd_ln2, projectid, taskid)
#         #with module_cmd(environment) as p:
#           # status3 = p(cmd_ln3, projectid, taskid)
#         with module_cmd(environment) as p:
#             status4 = p(cmd_ln4, projectid, taskid)
#         return status


def task_gsva(input, projectid, taskid, workdir):
    logger.info("task_gsva: 1.获得两个r脚本的参数字典")
    gsva_enrich_param, gsva_diffxp_param = set_gsva_params(input, projectid, taskid, workdir)
    logger.info(gsva_enrich_param)
    logger.info("task_gsva: 获得两个r脚本的命令行列表")
    #print(gsva_enrich_param["gmt_outdir"])
    for gmt in gsva_enrich_param["gmt_outdir"].keys():
        enrich_cmd_dict, diffxp_cmd_dict = get_gsva_cmd_dict(gsva_enrich_param, gsva_diffxp_param)
        enrich_cmd_line = enrich_cmd_dict.get(gmt)
        diffxp_cmd_line = diffxp_cmd_dict.get(gmt)

        enrich_cmd = f"Rscript /home/luyao/10X_scRNAseq_v3/src/Enrichment/GSVA_enrich.R {enrich_cmd_line}"
        diffxp_cmd = f"Rscript /home/luyao/10X_scRNAseq_v3/src/Enrichment/GSVA_pathway_diffxp.R {diffxp_cmd_line}"

        # 进入工作路径开始分析
        set_workdir_cmd = f"cd {workdir}/gsva"
        logger.info("task_gsva: 2.开始分析")
        start_gsva_analysis = f"{set_workdir_cmd} && module purge && module load {environment} && {enrich_cmd} && {diffxp_cmd}"
        with module_cmd(environment) as p:
            status = p(start_gsva_analysis, projectid, taskid)
        logger.info("task_gsva: Down!")

    logger.info("task_gsva: 3.生成 output.json.tsv")
    generate_output_tsv(workdir)

    logger.info("task_gsva: 4.生成软链接...")
    #make_links(workdir, projectid, taskid)
    factor, case, control = contrast.split(":")
    if (case != "all") and (control != "all"):
        cmd3= f"ln -s {workdir}/gsva/{database}/*{{png,pdf,docx,xls}} {workdir}/output/download"+"\n" +\
          f"ln -s {workdir}/gsva/{database}/*.png  {workdir}/output/" +"\n" +\
          f"ln -s {workdir}/gsva/{database}/diffexp_genesets_GSVA_score4group_{case}-vs-{control}.xls  {workdir}/output/diffexp_genesets_GSVA_score4group_{case}-vs-{control}.tsv"
    else:
        cmd3= f"ln -s {workdir}/gsva/{database}/*{{png,pdf,docx,xls}} {workdir}/output/download"+"\n" +\
              f"ln -s {workdir}/gsva/{database}/*.png  {workdir}/output/" +"\n" +\
              f"ln -s {workdir}/gsva/{database}/top10_term_t_value.xls  {workdir}/output/top10_term_t_value.tsv" +"\n" +\
              f"ln -s {workdir}/gsva/{database}/GSVA_top10_results.xls  {workdir}/output/GSVA_top10_results.tsv"

    with module_cmd(environment) as p:
             status = p(cmd3, projectid, taskid)
    return status
