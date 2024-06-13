import os, sys, yaml
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
import time
from oebio.utils.log import getLogger
import datetime

logger = getLogger('oe.cloud.sc.qsub')
from taskserver.tools.module_cmd import module_cmd

def find_cellranger( project_name, sample):
    project_dir = f"/public/cloud_scRNA/chenhaoruo/cellranger_dir/{project_name}"
    sample_dir = f"{project_dir}/result/cellranger/{sample}"
    config_dir = f"{project_dir}/cellranger_config.csv"
    if not os.path.exists(config_dir):
        logger.warning("未找到已完成的cellranger结果")
        return None
    elif not os.path.exists(sample_dir):
        logger.warning("未找到已完成的cellranger结果")
        return None
    else:
        config_file = pd.read_csv(config_dir)
        if len(config_file.loc[ config_file['Sample'].isin([sample]),])==0 :
            logger.warning("未找到已完成的cellranger结果")
            return None
        while (config_file.loc[ config_file['Sample'].isin([sample]), 'status'] == "running").any() :
            logger.info(f"抄作业模式启动，本地的样本{sample}正在运行，等待10分钟后检查，如长时间未完成，请检查路径以下的运行情况。异常请联系陈皓若。")
            logger.info(f"{sample_dir}")
            time.sleep(600)
            config_file = pd.read_csv(config_dir)
        if (config_file.loc[ config_file['Sample'].isin([sample]), 'status'] == "finish").all() :
            return sample_dir
        else:
            logger.warning("华为云未找到已完成的cellranger结果")
            return None

def find_cellranger_obs( project_name, sample):
    project_dir = f"obs://oe-scrna/works/scRNA/PROJECT/{project_name}"
    sample_dir = f"{project_dir}/cellranger/{sample}"
    html = f"{sample_dir}/outs/web_summary.html"
    cmd = f"obsutil stat {html}"
    result = os.system(cmd)
    if result == 0:
        logger.info(f"在OBS找到cellranger结果：{sample_dir}")
        return sample_dir
    else:
        logger.warning("OBS未找到已完成的cellranger结果，启动正常cellranger")
        return None

def find_cellranger_local( project_name, sample):
    remote_host = "scrna@10.100.10.69"
    project_dir = f"/cds/obs-scrnabak/auto_cellranger/{project_name}"
    sample_dir = f"{project_dir}/cellranger/{sample}"
    html = f"{sample_dir}/outs/web_summary.html"
    cmd = f'''
    if ssh -t {remote_host} '[ -f {html} ]'; then
        exit 0
    else
        exit 1
    fi
    '''
    result = os.system(cmd)
    if result == 0:
        logger.info(f"在本地找到cellranger结果：{sample_dir}")
        return sample_dir
    else:
        logger.warning("本地未找到已完成的cellranger结果，启动正常cellranger")
        return None

def update_cellranger_record_file(cellranger_record_file_name, samples, ddh ,status):
    #记录cellranger自动运行情况
    record_file = pd.read_csv(cellranger_record_file_name)
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    record_file.loc[(record_file['sample'].isin(samples)) & (record_file['project'].isin(ddh)), 'use_time'] = now
    record_file.loc[(record_file['sample'].isin(samples)) & (record_file['project'].isin(ddh)), 'status'] = status
    record_file.to_csv(cellranger_record_file_name, index=False)

def task_cellranger(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.tsv文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    env=input.loc['project', 'environment']
    cellranger_version = input.loc['parameters', 'cellranger_version']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    logger.info("step3.1:下载原始数据")
    sample = f"{input.loc['base', 'sample']['name']}"
    project_name = f"{input.loc['project','name']}"
    find_result = None
    record_file = "/public/cloud_scRNA/chenhaoruo/cellranger_record_file.csv"
    rl1 = f"--r1-length={input.loc['parameters', 'R1_length']}" if input.loc['parameters', 'R1_length'] else ""
    rl2 = f"--r2-length={input.loc['parameters', 'R2_length']} " if input.loc['parameters', 'R2_length'] else ""
    chemistry = f"--chemistry={input.loc['parameters', 'Chemistry']} " if input.loc['parameters', 'Chemistry'] else ""
    bam=f" " if input.loc['parameters', 'bam'] else "--no-bam  "
    force_cells = f"--force-cells={int(input.loc['parameters', 'force_cell'])} " if input.loc['parameters', 'force_cell'] !="" else ""
    if cellranger_version=="7.0.1":
        include_introns_param = "--include-introns=true" if input.loc['parameters', 'include_intron'] else "--include-introns=false"
    else:
        include_introns_param = "--include-introns" if input.loc['parameters', 'include_intron'] else ""
    if input.loc['parameters', 'copy']  :
        if rl1 == "" and rl2 == "" and chemistry =="" and force_cells == "" and cellranger_version == "7.0.1" and include_introns_param == "--include-introns=true" :
            find_result = find_cellranger( project_name, sample)
            if find_result is None:
                logger.info(f"在OBS中查找已完成的cellranger")
                find_result = find_cellranger_obs(project_name, sample)
                if find_result is None:
                    find_result = find_cellranger_local(project_name, sample)
                    if find_result is None:
                        update_cellranger_record_file( record_file, [sample], [project_name], "run_but_not_used")
                        logger.info(f"已在文件记录: {record_file}")
        else:
            update_cellranger_record_file( record_file, [sample], [project_name], "run_but_not_used")
            logger.warning("参数设置非常规，抄作业模式不启用")
            logger.info(f"已在文件记录: {record_file}")
    else:
        update_cellranger_record_file( record_file, [sample], [project_name], "run_but_not_used")
        logger.info("不使用抄作业模式，正常进行cellranger")
        ##1 download fastq files
    if find_result is None:
        cmd1 = f"oesc_data_prepare.py --input {workdir}/input/input.json  --outdir {workdir}"
        with module_cmd("obsutil/5.2.12") as p:
            status = p(cmd1, projectid, taskid)
    else:
        logger.info("找到已运行过的cellranger结果，跳过下载原始数据步骤，进入抄作业模式")
        logger.info(f"本样本原始:{find_result}")
        logger.warning(f"本模式将跳过所有在云平台中设置的cellranger参数，请检查网页版报告中设置的参数是否正确，如有错误请不要勾选抄作业模式，重新运行")
    ##2 run cellranger for input sample
    logger.info("step3.2:运行cellranger")

    if find_result is None:
        cmd2 = f"export PATH=/opt/softwares/cellranger-{cellranger_version}:$PATH && " \
               f"cd {workdir}/cellranger/ && if [ -s {sample} ]; then rm -rf {sample} ;fi  &&  cellranger count " \
            f"--id={sample} " \
            f"--transcriptome={input.loc['parameters', 'database']} " \
            f"--fastqs={workdir}/raw_data/{sample + '_fastqs'} " \
            f"{rl1} " \
            f"{rl2} " \
            f"{bam} " \
            f"--localcores=10 " \
            f"--localmem=120 " \
            f"{chemistry} " \
            f"--description={input.loc['project', 'name']}_{input.loc['project', 'species']}_{input.loc['project', 'sampletype']}  " \
            f" {include_introns_param} {force_cells}"
        with module_cmd(f"cellRanger/{cellranger_version}") as p:
            status=p(cmd2, projectid, taskid)
    else:
        if find_result.startswith("obs"):
            cmd2 = f"cd {workdir}/cellranger/ && if [ -s {sample} ]; then rm -rf {sample} ;fi  &&  "\
            f"obsutil cp {find_result} ./ -r -f -vlength -vmd5 "
            # f"bash /public/cloud_scRNA/chenhaoruo/scripts/write_todelete.sh '{find_result}'"
            with module_cmd(f"obsutil/5.2.12") as p:
                status=p(cmd2, projectid, taskid)
        elif find_result.startswith("/cds"):
            cmd2 = f"cd {workdir}/cellranger/ && if [ -s {sample} ]; then rm -rf {sample} ;fi  &&  "\
            f"scp -r scrna@10.100.10.69:{find_result} ./ "
            with module_cmd(f"OESingleCell/3.0.d") as p:
                status=p(cmd2, projectid, taskid)
        else:
            cmd2 = f"cd {workdir}/cellranger/ && if [ -s {sample} ]; then rm -rf {sample} ;fi  &&  "\
            f"cp -r {find_result} ./ && "\
            f"bash /public/cloud_scRNA/chenhaoruo/scripts/write_todelete.sh '{find_result}'"
            with module_cmd(f"OESingleCell/3.0.d") as p:
                status=p(cmd2, projectid, taskid)
        update_cellranger_record_file( record_file, [sample], [project_name], "used")
        logger.info(f"已在文件记录: {record_file}")
    if not os.path.exists(f"{workdir}/../cellranger/"):
        os.mkdir(f"{workdir}/../cellranger/")
    if not os.path.exists(f"{workdir}/../cellranger/{sample}"):
        os.symlink(f"{workdir}/cellranger/{sample}",f"{workdir}/../cellranger/{sample}")
    # upload bam
    if input.loc['parameters', 'bam']:
        logger.info("step3.21:上传bam文件")
        cmdx = f"bash /public/cloud_scRNA/chenhaoruo/scripts/auto-cellranger-txy/upload_bam.sh '{workdir}/cellranger/{sample}' '{project_name}'"
    elif os.path.isfile(f"{workdir}/cellranger/{sample}/outs/possorted_genome_bam.bam"):
        logger.info("选择不上传，删除bam文件")
        cmdx = f"rm {workdir}/cellranger/{sample}/outs/*.bam* "
    else :
        logger.info("选择不上传，无bam文件")
        cmdx = f"echo 'no need to upload bam file'"
    with module_cmd(f"OESingleCell/3.0.d") as p:
        status=p(cmdx, projectid, taskid)
    ##screenshot
    cmd3= f"oesc_screenshot.py -i {workdir}/{module}/{sample}/outs/web_summary.html -n {sample} -o {workdir}/{module}/{sample}/outs"
    #with module_cmd("chromedriver/88","oebio/1.2.15") as p:
    with module_cmd(env) as p:
        status=p(cmd3, projectid, taskid)
    # ==================================================================================================================
    # logger.info("删除中间文件fastq及SC_RNA_COUNTER_CS：")
    # cmd4=f"rm -rf {workdir}/{module}/{sample}/SC_RNA_COUNTER_CS && rm -rf {workdir}/raw_data/download_from_obs/"
    # with module_cmd(f"cellRanger/{cellranger_version}") as p:
    #     status=p(cmd4, projectid, taskid)
    # ==================================================================================================================
    logger.info("step3.3:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download/{sample}'):
        os.makedirs(f'{workdir}/output/download/{sample}')
    cmd3 = f"ln -s {workdir}/{module}/{sample}/outs/  {workdir}/output/download/{sample} && ln -s {workdir}/{module}/{sample}/outs/{sample}.png  {workdir}/output/download/"
    with module_cmd(f"cellRanger/{cellranger_version}") as p:
        status = p(cmd3, projectid, taskid)
    # ==================================================================================================================
    logger.info("step3.3:生成output.json.tsv")
    output_file = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)  # encoding='utf-8')
    output_file = output_file.loc[output_file['task_type'] == module]
    output_file = output_file.replace('{sample}', sample, regex=True)
    with open(f"{input.loc['parameters', 'database']}/config.yaml") as f:
        reference_config = yaml.safe_load(f)
    genome_linkage = f"{reference_config['genome_linkage']}"
    annotation_linkage = f"{reference_config['annotation_linkage']}"
    genome_version = f"{reference_config['genome_version']}"
    annotation_version = f"{reference_config['annotation_version']}"
    d_index = list(output_file.columns).index('input')
    output_file.iloc[2, d_index] =str({"Cell_Ranger_version": cellranger_version,
                                    "Reference_genome_link": genome_linkage,
                                    "Reference_genome_annotation_link": annotation_linkage,
                                    "Reference_genome_version": genome_version,
                                    "Reference_genome_annotation_version": annotation_version})
    output_file.to_csv(f'{workdir}/{module}/output.json.tsv', sep='\t', index=False, header=True)
    return status
