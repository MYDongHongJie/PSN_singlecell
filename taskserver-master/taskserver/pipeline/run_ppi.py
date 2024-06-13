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
    ## 1.获取filterd.h5seurat路径及参数信息
    env_module=input.loc['project', 'environment']
    input_dir = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}"
    diff_files_list = glob(f"{input_dir}/diffexp/*-diff-*_anno.tsv")
    species_ppi = int(input.loc['parameters', 'species_ppi'])
    topn_ppi = int(input.loc['parameters', 'topn_ppi'])
    ## 2.获取命令行
    for diff_gene_file in diff_files_list:
        prefix = diff_gene_file.split("/")[-1].replace("_anno.tsv", "")
        cmd = f"oesc_ppi.py " \
              f"--inputfile {diff_gene_file} " \
              f"--prefix {prefix} " \
              f"--noft   {topn_ppi} " \
              f"--outputdir {workdir}/{module} " \
              f"--species {species_ppi} "
        ## 3.执行分析
        with module_cmd(env_module) as p:
            status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    network_files_list = glob(f"{workdir}/{module}/*-diff-*ppi_network.tsv")
    ## 2.获取命令行
    def pair_files_by_prefix(list1, list2):
        paired_files = {}
        for filename1 in list1:
            prefix1 = filename1[:filename1.rfind('_') ]
            for filename2 in list2:
                if filename2.startswith(prefix1):
                    paired_files[filename1] = filename2
                    break
        return paired_files

    network_files_list2=[x.split("/")[-1] for x in network_files_list ]
    diff_files_list2=[x.split("/")[-1] for x in diff_files_list ]
    pairs = pair_files_by_prefix( diff_files_list2,network_files_list2)
    print(network_files_list2,"\n",diff_files_list2,"\n",pairs)
    for diff_file,network_file in pairs.items():
        cmd2 = f"Rscript /opt/softwares/oe_packages/cloud_report_script/scripts_use/ppi_circle.r " \
              f"-i {workdir}/{module}/{network_file} " \
              f"-d {input_dir}/diffexp/{diff_file} " \
              f"-o {workdir}/{module} "
        ## 3.执行分析
        with module_cmd(env_module) as p:
            status = p(cmd2, projectid, taskid)
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd = f"ln -s {workdir}/{module}/" \
          f"{{*-interaction.pdf," \
          f"*-interaction.png," \
          f"*-interaction.svg," \
          f"*-interaction.ne*.pdf," \
          f"*-interaction.ne*.png," \
          f"*-interaction.ne*.svg," \
          f"*-interaction.tsv," \
          f"*mapping_results.tsv," \
          f"*top_{topn_ppi}_ppi_network.pdf," \
          f"*top_{topn_ppi}_ppi_network.png}} " \
          f"{workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    logger.info("step3.3:生成output.json.tsv")
    output_file = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    output_file = output_file.loc[output_file['task_type'] == module]
    output_json=  pd.DataFrame(columns=list(output_file.columns))
    for diff_gene_file in diff_files_list:
        prefix = diff_gene_file.split("/")[-1].replace("_anno.tsv", "")
        print(output_json)
        output_json = output_json.append(output_file.replace('{prefix}', prefix, regex=True),ignore_index=True)
    output_json.to_csv(f'{workdir}/{module}/output.json.tsv', sep='\t', index=False, header=True)

    return status
