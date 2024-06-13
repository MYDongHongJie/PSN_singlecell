import os
import sys
import re
import json
import subprocess
import pandas as pd
import requests
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')


def task_diffexp(input, output_cfg, projectid, taskid, workdir):
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
    env_module = input.loc['project', 'environment']
    #filtered_h5seurat = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.h5seurat"
    ####进行输入对象rds或者h5seurat格式的识别和判断
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files

    rds_input = search_files(f'{base_wd}/{base_module_name}', "rds")
    h5seurat_input = search_files(f'{base_wd}/{base_module_name}', "h5seurat")
    seurat_ob = ""
    format_ob = ""
    if (len(rds_input) == 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    elif (len(rds_input) != 0 and len(h5seurat_input) == 0):
        seurat_ob = rds_input[0]
        format_ob = "rds"
    elif (len(rds_input) != 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    else:
        print(f"输入seurat对象为空，请检查！！！")
    print(seurat_ob)
    #seurat_ob_realpath=os.readlink(seurat_ob)
    ####以上完成输入对象的识别及赋值
    assay = input.loc['parameters', 'assay']
    contrast = input.loc['parameters', 'contrast']
    FC = input.loc['parameters', 'FC']
    DEG_Significance = input.loc['parameters', 'DEG_Significance']
    DEG_Significance_threshold = input.loc['parameters', 'DEG_Significance_threshold']
    DEGtest = input.loc['parameters', 'DEGtest']
    topn_diffgene = int(input.loc['parameters', 'topn_diffgene'])
    var2use = input.loc['parameters', 'var2use']
    levels4var = input.loc['parameters', 'levels4var']
    database = f"{input.loc['parameters', 'database']}/annotation/gene_annotation.xls"
    cellmeta = f"{input.loc['base', 'tasks'][0]['filename']}"

    add_params = ""
    # if var2use != "" and levels4var != "":
    #     add_params = f" --predicate {var2use} %in%  {levels4var}"

    if cellmeta != "":
        add_params = f" --cloudmeta  {workdir}/input/{cellmeta} "

    ## 2.获取命令行
    for i in range(len(contrast)):
        cmd1 = f"sctool " \
            f"-i {seurat_ob} " \
            f"-f {format_ob} " \
            f"-o {workdir}/{module} " \
            f"--assay {assay} " \
            f"--dataslot data,counts " \
            f"-j 10  " \
            f"{add_params} " \
            f"diffexp " \
            f"      -c {contrast[i]} " \
            f"      -k {FC}  " \
            f"      --{DEG_Significance} {DEG_Significance_threshold} " \
            f"      -e {DEGtest} "
        ## 3.执行分析
        with module_cmd(env_module) as p:
            status = p(cmd1, projectid, taskid)
    alldiff_files_list = [glob(e) for e in [f"{workdir}/{module}/*-all_diffexp_genes.tsv", f"{workdir}/{module}/*-diff-*FC*.tsv"]]
    print("test1")
    #diff_files_list = glob(f"{workdir}/{module}/*-diff-*FC*.tsv")
    print(alldiff_files_list)
    ## diff annotation
    alldiff_files_list2 = [i for k in alldiff_files_list for i in k]
    for alldiff_gene_file in alldiff_files_list2:
        cmd2 = f"sctool " \
            f"annotation " \
            f"      -g {alldiff_gene_file} " \
            f"      --anno {database}  "
        with module_cmd(env_module) as p:
            status = p(cmd2, projectid, taskid)

    ### diff heatmap
    diff_files_list2 = glob(f"{workdir}/{module}/*-diff-*FC*.tsv")
    group_all=[]
    for diff_gene_file in diff_files_list2:
        contrast_group = diff_gene_file.split('/')[-1].split('_')[0]
        case_name = diff_gene_file.split('/')[-1].split('-vs-')[0].split(contrast_group+'_')[-1]
        ctrl_name = diff_gene_file.split('/')[-1].split('-vs-')[-1].split('-diff-')[0]
        if case_name not in group_all:
            group_all.append(case_name)
        if ctrl_name not in group_all:
            group_all.append(ctrl_name)
    print(group_all)
    for diff_gene_file in diff_files_list2:
        contrast_group = diff_gene_file.split('/')[-1].split('_')[0]
        case_name = diff_gene_file.split('/')[-1].split('-vs-')[0].split(contrast_group+'_')[-1]
        ctrl_name = diff_gene_file.split('/')[-1].split('-vs-')[-1].split('-diff-')[0]
        groupcolor=f"{group_all.index(case_name)+1},{group_all.index(ctrl_name)+1}"
        cmd3 = f"scVis  " \
            f"-i {seurat_ob} " \
            f"-f {format_ob} " \
            f"-o {workdir}/{module} " \
            f"--assay {assay} " \
            f"--slot data,scale.data " \
            f"{add_params} " \
            f"--predicate '{contrast_group} %in% c(\"{case_name}\", \"{ctrl_name}\")' "\
            f"diff_heatmap " \
            f"      -d {diff_gene_file} " \
            f"      -n {topn_diffgene}  " \
            f"      -g {contrast_group}  " \
            f" --groupcolor {groupcolor}"\
            f"      --group_colors customecol2  " \
            f"      --sample_ratio 0.8   "
        with module_cmd(env_module) as p:
            status = p(cmd3, projectid, taskid)
    #删除Rplots.pdf文件
    cmd4=f"rm {workdir}/{module}/Rplots.pdf"
    with module_cmd(env_module) as p:
        status = p(cmd4, projectid, taskid)
    # =================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd = f"mv {workdir}/{module}/heatmap_cluster_*.txt {workdir}/output/ && ln -s {workdir}/{module}/* {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*.tsv"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    # logger.info("step3.3:生成output.json.xls")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    # df = output_df.loc[output_df['task_type'] == module]
    df = pd.DataFrame(columns=output_df.columns)
    df.loc[0] = ["diffexp", "summary", "diffexp_results_stat.tsv",
                 "tsv", "diffexp_results_stat.tsv", "数据统计", "", ""]
    diff_stat = pd.read_csv(f'{workdir}/{module}/diffexp_results_stat.tsv',
                            dtype=str, sep="\t", na_filter=False)
    info = {"num_diff": len(diff_stat), "diff_num": list(diff_stat.iloc[:, 4])}
    df.loc[1] = ["diffexp", "reportVars", info, "", "", "", "", ""]

    alldiff_files_list_anno = [glob(e) for e in [
        f"{workdir}/{module}/*-diff-*FC*_anno.tsv", f"{workdir}/{module}/*-all_diffexp_genes_anno.tsv"]]
    for i in range(len(alldiff_files_list_anno[0])):
        filename = alldiff_files_list_anno[0][i].split('/')[-1]
        df.loc[len(df)] = ["diffexp", "data", filename,
                           "gene", f"{filename}", re.sub('diff-pval-.*_anno.tsv', 'signification_DEG', filename), "", ""]
        # diff heatmap
        #heatmap_name = re.sub("-diff-*FC*.*_anno.tsv", "_heatmap.tsv",re.sub("^", "top"+str(topn_diffgene)+"_", filename))
        heatmap_name = re.sub("^", "top"+str(topn_diffgene)+"_",filename.split('-diff')[0]+"_heatmap.tsv")
        df.loc[len(df)] = ["diffexp", "diagram", heatmap_name,
                           "heatmap", heatmap_name, re.sub("_heatmap.tsv", "", re.sub('top\d+_', '', heatmap_name)), re.sub('.tsv', '', heatmap_name), ["download/"+heatmap_name.replace('tsv', 'pdf'), "download/" + heatmap_name.replace('tsv', 'png')]]
    for j in range(len(alldiff_files_list_anno[1])):
        filename = alldiff_files_list_anno[1][j].split('/')[-1]
        df.loc[len(df)] = ["diffexp", "data", filename,
                           "gene", f"{filename}", filename.replace("all_diffexp_genes_anno.tsv", "all_DEG"), "", ""]

    df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',
              index=False, header=True, encoding='utf-8')
    return status
