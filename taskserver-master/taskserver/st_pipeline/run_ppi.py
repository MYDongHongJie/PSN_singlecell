import os, sys
import json
import subprocess
import pandas as pd
import requests
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')
def task_ppi(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.tsv文件，用于生成output.json
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
    ## 1.获取diffexp路径及参数信息
    env_module="OESingleCell/v_3.0.0_visium_produce"
    input_dir = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}"
    diff_files_list = glob(f"{input_dir}/diffexp/*-diff-pvalue-*_anno.tsv")
    species_ppi = int(input.loc['parameters', 'species_ppi'])
    topn_ppi = int(input.loc['parameters', 'topn_ppi'])
    ## 2.获取命令行
    for diff_gene_file in diff_files_list:
        prefix = diff_gene_file.split("/")[-1].replace("_anno.tsv", "")
        cmd1 = f"/public/scRNA_works/works/guochy/ST_taskserver/scripts/ppi.py " \
              f"--inputfile {diff_gene_file} " \
              f"--prefix {prefix} " \
              f"--noft   {topn_ppi} " \
              f"--outputdir {workdir}/{module} " \
              f"--species {species_ppi} "
        ## 3.执行分析
        with module_cmd(env_module) as p:
            status = p(cmd1, projectid, taskid)
        cmd2 = f"Rscript /public/scRNA_works/works/guochy/ST_taskserver/scripts/ppi_circle.r " \
               f"-i {workdir}/{module}/{prefix}.ppi_network.tsv " \
               f"-d  {diff_gene_file} " \
               f"-o {workdir}/{module} " 
        with module_cmd("OESingleCell/v_3.0.0_visium_produce") as p:
            status = p(cmd2, projectid, taskid)
    # ==================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd3 = f"ln -s {workdir}/{module}/" \
          f"{{*-interaction.pdf," \
          f"*-interaction.png," \
          f"*-interaction.svg," \
          f"*_network.pdf," \
          f"*_network.png," \
          f"*-interaction.ne*.pdf," \
          f"*-interaction.ne*.png," \
          f"*-interaction.ne*.svg," \
          f"*-interaction.tsv," \
          f"*mapping_results.tsv}} " \
          f"{workdir}/output/download/"  
    with module_cmd(env_module) as p:
        status = p(cmd3, projectid, taskid)
    
    cmd4 = f"ln -s {workdir}/{module}/" \
    f"{{*-interaction.ne*.png,*_network.png}}" \
    f" {workdir}/output/"

    with module_cmd(env_module) as p:
        status = p(cmd4, projectid, taskid) 
    # ==================================================================================================================
    logger.info("step3.3:生成output.json.tsv")
    output_file = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    output_file = output_file.loc[output_file['task_type'] == module]
    output_json=  pd.DataFrame(columns=list(output_file.columns))
    for diff_gene_file in diff_files_list:
        prefix = diff_gene_file.split("/")[-1].replace("_anno.tsv", "")
        print(output_json)
        output_json = output_json.append(output_file.replace('{prefix}', prefix, regex=True),ignore_index=True)
        output_json.loc[len(output_json)] = ["PPI","diagram",f"{prefix}.top_25_ppi_network.png","image",f"download/{prefix}.top_25_ppi_network.png",f"{prefix}","",""]
    output_json.to_csv(f'{workdir}/{module}/output.json.tsv', sep='\t', index=False, header=True)

    return status

