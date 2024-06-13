import os, sys, yaml
import json
import subprocess
import pandas as pd
import requests
from shutil import copyfile
import numpy as np
from oebio.utils.log import getLogger

logger = getLogger('oe.cloud.sc.qsub')
from taskserver.tools.module_cmd import module_cmd

def task_spaceranger(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.tsv文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    if not os.path.exists(f"{workdir}/output/download"):
        os.makedirs(f"{workdir}/output/download")
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    if not os.path.exists(f'{workdir}/config'):
        os.makedirs(f'{workdir}/config')
    env_module = "OESingleCell/v_3.0.0_visium_produce"
    #========================================================================
    logger.info("step1:解析命令行，执行分析")
    logger.info("step1:下载原始数据")
    ##0 activate environment
    cmd0 = f"cp -r /public/scRNA_works/works/guochy/ST_taskserver/latest_ST_snakemake_pipeline/st_pipeline/envs {workdir}"
    with module_cmd(env_module) as p:
        status = p(cmd0, projectid, taskid)
    ##1 download fastq files and images
    cmd1 = f"cd {workdir} && oest_data_prepare.py --input {workdir}/input/input.json --outdir {workdir}"
    with module_cmd("obsutil/5.2.12") as p:
        status = p(cmd1, projectid, taskid)
    ##2 run spaceranger for input samples
    logger.info("step1.2:运行spaceranger")
    sample_info = pd.read_csv(f'{workdir}/config/samples.csv', index_col=0)
    probe_dict = {"小鼠":"/public/dev_scRNA/software/spaceranger-2.0.0/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv","人":"/public/dev_scRNA/software/spaceranger-2.0.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"}
    sample = input.loc["base", "sample"]["name"]
    transcriptome = input.loc["parameters","transcriptome"]
    fastq = sample_info.loc[sample, "fastq"]
    image_file = sample_info.loc[sample, "image"]
    cyta_image = sample_info.loc[sample, "cytaimage"]
    slide = sample_info.loc[sample, "slide"]
    slide_area = sample_info.loc[sample, "slide_area"]
    species = sample_info.loc[sample, "species"]
    with open(f'{workdir}/config/config.yaml','r',encoding='utf-8') as y:
        config_info = yaml.load(y,yaml.Loader)
        library_type = config_info['params']['library_type'] 

    if library_type == "fresh":
        probe_pram = " "
    else:
        for k,v in probe_dict.items():
            if k in species:
                probe_set = probe_dict[k]
                probe_pram = f"--probe-set={probe_set}"
       

    if pd.isnull(sample_info.loc[sample, "image_json_file"]):
        json_param = "  "
    else:
        json_param = f"--loupe-alignment={sample_info.loc[sample, 'image_json_file']}"
    if input.loc["parameters","localcores"] != "":
        localcores = int(input.loc["parameters","localcores"])
    else:
        localcores = 10
    if input.loc["parameters","localmem"] != "":
        localmem = int(input.loc["parameters","localmem"])
    else:
        localmem = 120
    spaceranger_version = input.loc['parameters','spaceranger_version']
    if input.loc["parameters","generate_bam"] == "no":
        generate_bam = "--no-bam"
    else:
        generate_bam = " "
    if library_type == "cytassist":
        cytaimage_pram = f"--cytaimage={cyta_image}"
        if image_file != "":
            image_pram = f"--image={image_file}"
        else:
            image_pram = " "
    else:
        cytaimage_pram = " "
        image_pram = f"--image={image_file}"
    
    cmd2 = f"module load {workdir}/envs/spaceranger/{spaceranger_version} && cd {workdir}/spaceranger/ && if [ -s {sample} ]; then rm -rf {sample} ;fi  && spaceranger count "    \
           f"--id={sample} "    \
           f"--transcriptome={transcriptome} "   \
           f"--fastqs={fastq} "    \
           f"--sample={sample} "   \
           f"{image_pram} " \
           f"{cytaimage_pram} " \
           f"--slide={slide} "   \
           f"--slidefile={workdir}/raw_data/{slide}.gpr "   \
           f"--area={slide_area} "    \
           f"--r1-length=28 "    \
           f"--r2-length=91 "   \
           f"--localcores={localcores} "    \
           f"--localmem={localmem} "    \
           f"{probe_pram} " \
           f"--description={input.loc['project', 'name']}_{input.loc['project', 'species']}_{input.loc['project', 'sampletype']} "  \
           f"{generate_bam} {json_param}"

    with module_cmd(env_module) as p:
        status=p(cmd2, projectid, taskid)
    ##screenshot
    cmd3= f"oest_html_screenshot.py -i {workdir}/{module}/{sample}/outs/web_summary.html -n {sample} -o {workdir}/{module}/{sample}/outs"
    with module_cmd(env_module) as p:
        status=p(cmd3, projectid, taskid)
    cmd4= f"cp {image_file} {workdir}/{module}/{sample}/outs/spatial/"
    with module_cmd(env_module) as p:
        status=p(cmd4, projectid, taskid)
    ## transfer bam
    if input.loc["parameters","generate_bam"] == "yes":
        sample_dir = f'{workdir}/{module}/{sample}'
        project_id = input.loc["project", "name"]
        cmd5 = f"bash /public/scRNA_works/works/guochy/ST_taskserver/ST_snakemake_template/st_pipeline/scripts/upload_bam.sh {sample_dir} {project_id}"
        with module_cmd(env_module) as p:
            status = p(cmd5, projectid, taskid)
    
    logger.info("step1.3:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download/{sample}'):
        os.makedirs(f'{workdir}/output/download/{sample}')
    cmd6 = f"ln -s {workdir}/{module}/{sample}/outs/*  {workdir}/output/download/{sample} && ln -s {workdir}/{module}/{sample}/outs/{sample}.png  {workdir}/output/download/"
    with module_cmd(env_module) as p:
        status = p(cmd6, projectid, taskid)
    

#==================================================================================================================
    logger.info("step1.4:生成output.json.tsv")
    output_file = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)  # encoding='utf-8')
    output_file = output_file.loc[output_file['task_type'] == module]
    output_file = output_file.replace('{sample}', sample, regex=True)
    with open(f"{input.loc['parameters', 'transcriptome']}/config.yaml") as f:
        reference_config = yaml.safe_load(f)
    genome_linkage = f"{reference_config['genome_linkage']}"
    annotation_linkage = f"{reference_config['annotation_linkage']}"
    genome_version = f"{reference_config['genome_version']}"
    annotation_version = f"{reference_config['annotation_version']}"
    d_index = list(output_file.columns).index('input')
    output_file.iloc[2, d_index] =str({"spaceranger_version": spaceranger_version,
                                    "Reference_genome_link": genome_linkage,
                                    "Reference_genome_annotation_link": annotation_linkage,
                                    "Reference_genome_version": genome_version,
                                    "Reference_genome_annotation_version": annotation_version})
    output_file.to_csv(f'{workdir}/{module}/output.json.tsv', sep='\t', index=False, header=True)
    return status