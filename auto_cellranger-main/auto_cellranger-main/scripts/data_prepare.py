#!/usr/bin/env python3
# encoding: utf-8
"""
Project : snake_scRNAseq_pipeline
Author  : Jiaoyang Dong
Contact : jiaoyang.dong@oebiotech.com
File   : data_prepare_.py
Time   : 2021-08-19 14:16:04
"""
import os,sys,yaml,subprocess,click
from glob import glob
from copy import deepcopy

exec(open('/data/software/modules/modules-v4.2.1/init/python.py').read())


def get_obs_sample(address):
    '''
    return all the sample name in a list, given a obs address 
    '''
    folders = os.popen("obsutil ls -d -s %s | grep obs  | grep -v '\(md5\|bam_files\)' | sed 's#%s##'"%(address,address)).readlines()
    return [i.strip().replace("/","") for i in folders if i.strip() != "" ]

# deprecated 
# def get_obs_subdir(address):
#     folders = os.popen("obsutil ls -d -s %s | grep obs | grep -v md5 | grep '%s.*/' "%(address,address)).readlines()
#     return [i.strip() for i in folders]
def get_obs_subdir(address,sample_list):
    '''
    return VALID obs subdirectory list (at sample level), given a obs address and a sample list. 
    '''
    remote = [i.strip() for i in os.popen("obsutil ls -d -s %s | grep obs | grep -v '\(md5\|bam_files\)' | grep '%s.*/' "%(address,address)).readlines() ] 
    local = [address + sample_list[i] + "/" for i in range(len(sample_list))]
    return list(set(remote) & set(local))

def download(sub_addresses,outdir = "raw_data"): # sub_addresses in a list
    for sub_address in sub_addresses:
        subprocess.run(
            f"/data/software/obsutil/5.2.12/obsutil cp {sub_address}  {outdir} -r -f -vlength -vmd5", shell=True, check=True)
        print("\033[0;32m%s\033[0m" % f"{sub_address} download successfully.")

def rename (outdir, batch):
    '''
    use after download a new # use every time
    '''
    all_names = [item 
                for sublist in [glob(f'{outdir}/*' + ext) for ext in [f"/*.gz", "/*.md5"]]
                for item in sublist]
    batch_names = [item for sublist in
                  [glob(f'{outdir}/*' + ext) for ext in [f"/batch_*.gz"]]
                  for item in sublist]
    filtered_names = [x for x in all_names if
                     all(y not in x for y in batch_names)]
    for name in filtered_names:
        file_path, file_name = os.path.split(name)
        os.rename(name, os.path.join(file_path, "batch_" + str(batch) + "_" + file_name))

##===========================================Step1: dowload fastqs name=========================================
def raw_data_download(addresses,sample_list,outdir):
    # skip existing folders.
    sample_list2download = deepcopy(sample_list)
    for sample in sample_list :
        print(sample)
        if os.path.exists(f'{outdir}/{sample}/'):
            print(f"\033[0;33mWaring: Folder {outdir}/{sample} exists. Skipping...\033[0m")
            sample_list2download.remove(sample) 
            #os.system(f'rm -r {outdir}/{sample}')
    # download
    if len(sample_list2download) >0:
        print("\033[0;32m%s%s\033[0m" % ("Samples to download: " , ', '.join(sample_list2download)))
        batch = 1
        for address in addresses:
            sub_addresses = get_obs_subdir(address,sample_list2download) # valid sub directories for each run
            if len(sub_addresses) != 0 :
                download(sub_addresses,outdir)
                rename (outdir, batch)
            else:
                sys.exit("\033[31;1m"+"ERROR: sample not found in obs address "+ address +"."+"\033[0m")
            batch += 1

##===========================================Step2: formart fastqs name=========================================
def raw_data_formart(sample_list,outdir):
    for sample in sample_list:
        if sample in os.listdir(f"{outdir}/"): 
            if not os.path.exists(f'{outdir}/{sample}/{sample}_fastqs'):
                os.mkdir(f'{outdir}/{sample}/{sample}_fastqs')
            ########################
            for i in {"1", "2"}:
                #fq1_list = glob(f'{outdir}/{sample}/*R{i}.fastq,fq.gz')
                fq_list = sorted([item for sublist in [glob( f'{outdir}/{sample}'+ ext) for ext in [f"/*R{i}.fastq.gz", f"/*_{i}.fq.gz", f"/*R{i}_001.fastq.gz"]] for item in sublist])
                j = 1
                for fastq in fq_list:
                    src =  fastq.replace(f"{outdir}/{sample}","..")
                    dst = f'{outdir}/{sample}/{sample}_fastqs/{sample}_S{j}_L001_R{i}_001.fastq.gz'
                    fastq_show=fastq.replace(f"{outdir}/","")
                    des_show=f"{sample}/{sample}_fastqs/{sample}_S{j}_L001_R{i}_001.fastq.gz"
                    if not os.path.exists(dst):
                        os.symlink(src, dst)
                        print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show}  {des_show} ")
                    else:
                        print("\033[0;33m%s\033[0m" % f'Waring: {des_show} exists. Skiping...')
                    j = j + 1
        else:
            print("\033[0;31m%s\033[0m" % f"The fastq files of {sample} does not exist, please check it manually!!")

                   
@click.command()
@click.option('-c', 'configfile', type=str,default = 'config/config.yaml', 
    help='config yaml files with raw_data_obs_address, samples information. df:config/config.yaml')
@click.option('-o', 'outdir', type=str, default = 'raw_data',
    help='Output directory, default: raw_data')

def run(configfile,outdir):
    """example: python data_prepare.py -c config/config.yaml -o raw_data """
    with open(configfile) as f:
        config = yaml.safe_load(f)
    addresses = config["cellranger_params"]["raw_data_obs_address"]
    samples_config = config["cellranger_params"]["samples"]
    ## get sample names
    if type(samples_config) == list:
        if all(isinstance(i, str) for i in samples_config):
            sample_list = samples_config
        elif all(isinstance(i, dict) for i in samples_config):
            sample_list = [v for i in samples_config for v,k in i.items()]
    elif type(samples_config) == str and samples_config == "all":
        print(addresses[0])
        if addresses is None or addresses[0] is None or addresses[0] == "":
            sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
        else:
            sample_list = []
            for address in addresses:
                # sample_list = get_obs_sample(address)
                # sample_list.extend(get_obs_sample(address))
                # sample_list = list(set(sample_list))
                [sample_list.append(i) for i in get_obs_sample(address) if not i in sample_list  ]
    else:
        sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples can only be string 'all' or a list of sample names."+"\033[0m")
    print("\033[0;32m%s%s\033[0m" % ("Samples to prepare: " , ', '.join(sample_list)))

    ## Step1: download raw data
    raw_data_download(addresses,sample_list,outdir)

    ## Step2: formart fastqs name
    raw_data_formart(sample_list,outdir)

if __name__ == '__main__':
    run()