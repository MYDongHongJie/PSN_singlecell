import os, sys,re
import shutil
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
from glob import glob
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_diffexp(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作目录,pojectid层级
    :return: status
    """
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    #input = pd.read_json(f'{wkdir}/input/input.json', orient="index",dtype={"id":'str'})
    #input_meta=f"{wkdir}/input/cell.xls"
    cellmeta = input.loc['base', 'tasks'][0]['filename']
    add_params=""
    if cellmeta != "":
        add_params = add_params + f" --cloudmeta  {cellmeta} "
    seurat_ob = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/*/*.rds"
    #image=input.loc["parameters","image"]
    environment="OESingleCell/v_3.0.0_visium_produce"
    database=input.loc["parameters","database"]
    assay=input.loc["parameters","assay"]
    design=input.loc["parameters","design"]
    if input.loc["parameters", "foldchange"] != "":
        foldchange =  float(input.loc["parameters", "foldchange"])
    else:
        foldchange = 1.5
    DEGtest=input.loc["parameters","DEGtest"]
    DEG_Significance = input.loc['parameters', 'DEG_Significance']
    DEG_Significance_threshold = input.loc['parameters', 'DEG_Significance_threshold']
    contrast=input.loc["parameters","contrast"]

    ## 2.获取命令行
    for i in range(len(contrast)):
        cmd1=f"sctool "\
            f"  -i {seurat_ob} "\
            f"  -f rds "\
            f"  -o {workdir}/{module} "\
            f"  --image FALSE "\
            f"  --assay {assay} "\
            f"diffexp  "\
            f"  --design  ~{design} "\
            f"  --contrast  {contrast[i]} "\
            f"  --FC {foldchange} "\
            f"  --{DEG_Significance} {DEG_Significance_threshold} " \
            f"  -e {DEGtest} "\
            f"  --anno {database}annotation/gene_annotation.xls"
         

        logger.info(cmd1)
        logger.info("step1: 执行分析")
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
        with module_cmd(environment) as p:
            status=p(cmd1, projectid, taskid)
    
    for file_name in os.listdir(f"{workdir}/{module}"):
        if file_name.endswith('.xls'):
            new_file_name = os.path.splitext(file_name)[0] + '.tsv'
            shutil.move(os.path.join(f"{workdir}/{module}", file_name), os.path.join(f"{workdir}/{module}", new_file_name))
    
    logger.info("step2:生成output.json.xls")
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
        f"{workdir}/{module}/*-diff-pval-*_anno.tsv", f"{workdir}/{module}/*-all_diffexp_genes_anno.tsv"]]
    for i in range(len(alldiff_files_list_anno[0])):
        filename = alldiff_files_list_anno[0][i].split('/')[-1]
        df.loc[len(df)] = ["diffexp", "data", filename,
                           "gene", f"{filename}", re.sub('diff-pval-.*_anno.tsv', 'signification_DEG', filename), "", ""]
        # # diff heatmap
        # heatmap_name = re.sub("-diff-pval-.*_anno.tsv", "_heatmap.tsv",
                              # re.sub("^", "top"+str(topn_diffgene)+"_", filename))
        # df.loc[len(df)] = ["diffexp", "diagram", heatmap_name,
                           # "heatmap", heatmap_name, re.sub("_heatmap.tsv", "", re.sub('top\d+_', '', heatmap_name)), re.sub('.tsv', '', heatmap_name), ["download/"+heatmap_name.replace('tsv', 'pdf'), "download/" + heatmap_name.replace('tsv', 'png')]]
    for j in range(len(alldiff_files_list_anno[1])):
        filename = alldiff_files_list_anno[1][j].split('/')[-1]
        df.loc[len(df)] = ["diffexp", "data", filename,
                           "gene", f"{filename}", filename.replace("all_diffexp_genes_anno.tsv", "all_DEG"), "", ""]
    df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',
              index=False, header=True, encoding='utf-8')
    logger.info("step3:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd_ln= f"ln -s {workdir}/{module}/{{*_anno.tsv,*_stat.tsv}}  {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*.tsv"
    with module_cmd("OESingleCell/v_3.0.0_visium_produce") as p:
        status = p(cmd_ln, projectid, taskid)
    return status

