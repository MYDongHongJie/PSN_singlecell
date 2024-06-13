import pandas as pd 
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger
import glob
import csv
logger = getLogger('oe.cloud.sc.qsub')

from taskserver.tools.module_cmd import module_cmd


def find_cellranger_local( project_name, sample):
    remote_host = "scrna@10.100.10.69"
    project_dir = f"/cds/scRNA-bambak/{project_name}"
    sample_dir = f"{project_dir}/{sample}"
    bam = f"{sample_dir}/*bam"
    cmd = f'''
    if ssh -t {remote_host} '[ -f {bam} ]'; then
        exit 0
    else
        exit 1
    fi
    '''
    result = os.system(cmd)
    if result == 0:
        logger.info(f"在本地找到bam：{sample_dir}")
        return sample_dir
    else:
        logger.warning(f"本地未找到bam,请重新运行cellranger：{sample_dir}")
        sys.exit(1)

def get_subdirectories(directory):
    subdirectories = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
    return subdirectories

def task_loom(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    wkdir =workdir
    module = input.loc['task', 'type']
    # cellranger_version = input.loc['parameters', 'cellranger_version']
    project_name = f"{input.loc['project','name']}"
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    if os.path.exists(f'{workdir}/../cellranger'):
        os.symlink(f'{workdir}/../cellranger',f'{workdir}/{module}/cellranger')
    else:
        if not os.path.exists(f'{workdir}/{module}/cellranger'):
            os.makedirs(f'{workdir}/{module}/cellranger')
        metadata = pd.DataFrame(input.loc["project","samples"])
        metadata = metadata.set_index('name',drop=True)
        #if any(metadata["group"] == ""): metadata["group"] = metadata["label"]
        metadata["group"] = metadata.apply(lambda row: row["label"] if row["group"].strip() == "" else row["group"], axis=1)
        metadata [["group"+str(i) if i!= 1 else "group"
            for i in range(1, 1+ metadata["group"].str.split(",",expand=True).shape[1])] ] = metadata["group"].str.split(",",expand=True)
        metadata["species"]= input.loc["project","species"]
        ## create soft link & write metadata.csv
        ## 加入判断bam是否存在，不存在则scp过来
        samples=[]
        for task in input.loc["base","tasks"]:
            if os.path.exists(f"{workdir}/../../../../20{task['taskId'][0:2]}/{task['taskId'][2:4]}/{task['taskId'][0:10]}/{task['taskId']}/cellranger"):
                if not os.path.exists(f'{workdir}/../cellranger'):
                    os.makedirs(f'{workdir}/../cellranger')
                sample = next(os.walk(f"{workdir}/../../../../20{task['taskId'][0:2]}/{task['taskId'][2:4]}/{task['taskId'][0:10]}/{task['taskId']}/cellranger"))[1][0]
                if sample in metadata.index:
                    # soft link
                    os.symlink(f"{workdir}/../../../../20{task['taskId'][0:2]}/{task['taskId'][2:4]}/{task['taskId'][0:10]}/{task['taskId']}/cellranger/{sample}",
                            f"{workdir}/{module}/cellranger/{metadata['label'][sample]}")
                    if not os.path.exists(f"{workdir}/../cellranger/{metadata['label'][sample]}"):
                        os.symlink(f"{workdir}/../../../../20{task['taskId'][0:2]}/{task['taskId'][2:4]}/{task['taskId'][0:10]}/{task['taskId']}/cellranger/{sample}",
                                f"{workdir}/../cellranger/{metadata['label'][sample]}")
                    samples.append(sample)
                else:
                    print("sample %s not found."%sample)
                    sys.exit(1)
                    #sys.exit(f"sample {sample} not found.")
            else:
                print("cellranger dir not found.")
                sys.exit(1)
                #sys.exit(f"sample {sample} not found.")
        metadata2 = metadata.loc[samples,:]
        metadata2.rename(columns={"label":"sampleid"},inplace=True)
        metadata2.to_csv(f"{workdir}/{module}/metadata.csv",sep=',',index=False,header=True)
    samples = get_subdirectories(f"{workdir}/{module}/cellranger")
    if not os.path.exists( f"{wkdir}/output/download"):
        os.makedirs( f"{wkdir}/output/download")
    if not os.path.exists( f"{wkdir}/loom"):
        os.makedirs( f"{wkdir}/loom")
    if not os.path.exists( f"{wkdir}/loom/loom"):
        os.makedirs( f"{wkdir}/loom/loom")
    species_dir = input.loc['parameters', 'species']
    if not os.path.exists(f"{species_dir}genes/genes.gtf"):
        logger.info(f"{species_dir}genes.gtf not found")
        sys.exit(1)
    velo_cmd = ''
    bam_cmd = ''
    ln_cmd = ''
    for sample in samples:
        bam_files = glob.glob(os.path.join(f'{workdir}/{module}/cellranger/{sample}/outs/', '*.bam'))
        loom_files = glob.glob(os.path.join(f'{workdir}/{module}/cellranger/{sample}/velocyto/', '*.loom'))
        if not bam_files:
            find_result = find_cellranger_local(project_name, sample)
            bam_cmd += f"scp scrna@10.100.10.69:{find_result}/* {workdir}/{module}/cellranger/{sample}/outs/ &&\n"
        else:
            logger.info(f"{sample} bam files found.")
        if not loom_files:
            velo_cmd += f"velocyto run10x  {workdir}/{module}/cellranger/{sample}/ {species_dir}genes/genes.gtf &&\n"
        else:
            logger.info(f"{sample} loom files found.")
        ln_cmd += f"ln -s {workdir}/{module}/cellranger/{sample}/velocyto/{sample}.loom {workdir}/{module}/loom/ &&\n"
    cmd=f"source /opt/softwares/oe_packages/script/velocity.sh &&"\
        f"{bam_cmd} {velo_cmd} {ln_cmd} echo 'loom已生成，可以开始运行velocity' > {wkdir}/loom/说明.txt "
    logger.info(cmd)
    logger.info("开始运行")
    environment="velocity/1.0.0"
    #with module_cmd(f"seurat/{input.loc['parameters', 'seurat_version']}") as p:####need to confirm the seurat module version
    with module_cmd(environment) as p:
        status=p(cmd, projectid, taskid)
    useless_file = f"{wkdir}/loom/说明.tsv"
    useless_file_data = [
        ["说明"],
        ["loom已生成，可以开始velocity分析"]
    ]
    with open(useless_file, 'w', newline='') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerows(useless_file_data)
    # 定义文件路径
    tsv_file_path = f"{wkdir}/loom/output.json.tsv"
    # 定义数据
    data = [
        ["task_type", "result_module", "input", "type", "file", "title", "downloadName", "downloadPath"],
        ["loom", "data", "说明.tsv", "tsv", "说明.tsv", "说明", "", ""]
    ]
    # 写入到 TSV 文件
    with open(tsv_file_path, 'w', newline='') as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerows(data)
    cmd_ln = f"ln -s {wkdir}/loom/说明.txt {wkdir}/output/download/ "
    logger.info("完成")
    with module_cmd(environment) as p:
        status = p(cmd_ln, projectid, taskid)
    return status ###run_cmd use the default function from xiufeng.yang
