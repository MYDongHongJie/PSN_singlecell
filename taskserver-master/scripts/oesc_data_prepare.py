#!/usr/bin/env python3
# encoding: utf-8
import os, sys
import time
import click
import subprocess
import pandas as pd
from glob import glob

def raw_data_download(fastqs, sample, outdir):
    ##===========================================Step1: download fastqs ===========================================
    for batch, fastq in fastqs.items():
        for file in fastq:
            split_path = file.split("/")
            year = split_path[-3]
            id = split_path[-2]
            print(file)
            subprocess.run(
                f"obsutil cp  {file}  {outdir}/raw_data/download_from_obs/{year}_{id}_{split_path[-1]}  -f -flat -vlength -vmd5",
                shell=True, check=True)
    ##===========================================Step2: formart fastqs name=========================================
    if not os.path.exists(f'{outdir}/raw_data/{sample}_fastqs'):
        os.mkdir(f'{outdir}/raw_data/{sample}_fastqs')
        ########################
        fq1_list = [item for sublist in
                    [glob(f'{outdir}/raw_data/download_from_obs' + ext) for ext in [f"/*R1.fastq.gz", "/*_1.fq.gz","/*R1_001.fastq.gz","/*.R1.fq.gz"]]
                    for item in sublist]
        fq1_list.sort(key=os.path.getmtime)
        j = 1
        for fastq1 in fq1_list:
            ####R1==================================================================================================
            des1 = f'{outdir}/raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R1_001.fastq.gz'
            fastq_show1 = fastq1.replace(f"{outdir}", "")
            des_show1 = f"raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R1_001.fastq.gz"
            print(des1)
            print(des_show1)
            if not os.path.exists(des1):
                os.symlink(os.path.abspath(fastq1), des1)
                print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show1}  {des_show1} ")
            else:
                print("\033[0;31m%s\033[0m" % f'Warning: {des_show1} has existed!!!')
            ####R2==================================================================================================
            fastq2 = fastq1.replace("R1.fastq.gz", "R2.fastq.gz").replace("_1.fq.gz", "_2.fq.gz").replace("R1_001.fastq.gz", "R2_001.fastq.gz").replace(".R1.fq.gz",".R2.fq.gz")
            des2 = f'{outdir}/raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R2_001.fastq.gz'
            fastq_show2 = fastq2.replace(f"{outdir}", "")
            des_show2 = f"raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R2_001.fastq.gz"
            print(des1)
            print(des_show1)
            if not os.path.exists(des2):
                os.symlink(os.path.abspath(fastq2), des2)
                print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show2}  {des_show2} ")
            else:
                print("\033[0;31m%s\033[0m" % f'Warning: {des_show2} has existed!!!')
            j = j + 1

@click.command()
@click.option("--input", "-i", type=click.Path("r"), help="input files with gene list")
@click.option("--outdir", "-o", type=click.Path(exists=False), help="output directory.")
def data_prepare(input="input.json文件绝对路径", outdir="结果输出目录"):
    input = pd.read_json(f'{input}', orient="index")
    fastqs = input.loc['base', "sample"]['dataFiles']
    #将数据桶路径由上海二切换到上海一
    fastqs= {key: [url.replace('obs://scrna-cloud/', 'obs://scrna-cloud-sh1/') for url in urls] for key, urls in fastqs.items()}
    sample = input.loc['base', 'sample']['name']
    raw_data_download(fastqs, sample, outdir)

if __name__ == "__main__":
    data_prepare()
