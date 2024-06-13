#!/usr/bin/env python3
# encoding: utf-8
"""
Project : snake_multimodal_pipeline
Author  : Sun Kun
Contact : kun.sun@oebiotech.com
File   : data_prepare_.py
Time   : 2021-11-01
"""
import os, sys, yaml, subprocess, click, shutil
from glob import glob
from copy import deepcopy
import pandas as pd

##===================  from obs download rawdata ===============================
def download(addresses,outdir = "raw_data"): # download and adjust filename
    for address in addresses:
        subprocess.run(
            f"/data/software/obsutil/5.2.12/obsutil cp {address} {outdir}/download -r -f -vlength -vmd5", shell=True, check=True)
        print("\033[0;32m%s\033[0m" % f"{address} download successfully.")

##======================= rename rawdata ====================================
def raw_data_formart(sample, outdir):
    if not os.path.exists(f'{outdir}/{sample}_fastqs'):
        os.mkdir(f'{outdir}/{sample}_fastqs')
        os.mkdir(f'{outdir}/{sample}_fastqs/gex')
        os.mkdir(f'{outdir}/{sample}_fastqs/atac')
    if len(glob(f'{outdir}/download/*/{sample}/')):
        for i in {"1", "2"}:
            j=1
            fq_name_refile= glob(f'{outdir}/download/*/{sample}/*R{i}.fastq.gz')
            fq_name_refile.sort()
            for num in range(0,len(fq_name_refile)):
                fq_name = os.path.abspath(fq_name_refile[num])
                fq_newname = f'{outdir}/{sample}_fastqs/gex/{sample}_S1_L00{j}_R{i}_001.fastq.gz'
                os.symlink(fq_name,fq_newname)
                j = j + 1
    if len(glob(f'{outdir}/download/*/{sample}_ATAC/')):
        for i in {"1", "2", "3"}:
            if i == "1":
                j=1
                fq_R_oldname_refile = glob(f'{outdir}/download/*/{sample}_ATAC/*R1*.fastq.gz')
                fq_R_oldname_refile.sort()
                fq_I_oldname_refile = glob(f'{outdir}/download/*/{sample}_ATAC/*I1*.fastq.gz')
                fq_I_oldname_refile.sort()
                for num in range(0,len(fq_R_oldname_refile)):
                    fq_R_oldname = os.path.abspath(fq_R_oldname_refile[num])
                    fq_I_oldname = os.path.abspath(fq_I_oldname_refile[num])
                    fq_R_newname = f'{outdir}/{sample}_fastqs/atac/{sample}_S1_L00{j}_R1_001.fastq.gz'
                    fq_I_newname = f'{outdir}/{sample}_fastqs/atac/{sample}_S1_L00{j}_I1_001.fastq.gz'
                    os.symlink(fq_R_oldname, fq_R_newname)
                    os.symlink(fq_I_oldname, fq_I_newname)
                    j=j+1
            else:
                r=1
                fq_name_refile = glob(f'{outdir}/download/*/{sample}_ATAC/*R{i}*.fastq.gz')
                fq_name_refile.sort()
                for num in range(0,len(fq_name_refile)):
                    fq_name = os.path.abspath(fq_name_refile[num])
                    fq_newname = f'{outdir}/{sample}_fastqs/atac/{sample}_S1_L00{r}_R{i}_001.fastq.gz'
                    os.symlink(fq_name, fq_newname)
                    r=r+1

##====================== get and adjust data ===================================
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
        #download
        download(addresses,outdir)
        for sample in sample_list2download:
            raw_data_formart(sample,outdir)
            #creat libraries file
            fastqs = [os.path.abspath(f'{outdir}/{sample}_fastqs/gex'),os.path.abspath(f'{outdir}/{sample}_fastqs/atac')]
            library_type = ["Gene Expression","Chromatin Accessibility"]
            lib_sample = [sample,sample]
            libraries = pd.DataFrame({'fastqs':fastqs,'sample':lib_sample,'library_type':library_type})
            libraries.to_csv(f'config/library/{sample}_libraries.csv',index=False,sep=",",encoding="UTF-8")



#==================================================================================================================
@click.command()
@click.option('-i','metadatafile',type=str,default = 'config/metadata.csv',
    help='metadata file with sampleid, group, bach informaation, df:config/metadata.csv')
@click.option('-c', 'configfile', type=str,default = 'config/config.yaml', 
    help='config yaml files with raw_data_obs_address, samples information. df:config/config.yaml')
@click.option('-o', 'outdir', type=str, default = 'raw_data',
    help='Output directory, default: raw_data')

def run(metadatafile,configfile,outdir):
    """example: python data_prepare.py -i metadata.csv-c config/config.yaml -o raw_data """
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
            samples = pd.read_csv(metadatafile, dtype=str).set_index("sampleid", drop=False)
            samples["addition.anno"]=[os.path.abspath("result/cellranger/"+sample+"/outs/per_barcode_metrics.csv") for sample in samples.index]
            samples.to_csv("config/metadata.csv",encoding='utf-8', sep=",",index=False)
            samples.index.names = ["sample_id"]
            sample_list = list(samples.index)
    else:
        sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples can only be string 'all' or a list of sample names."+"\033[0m")
    print("\033[0;32m%s%s\033[0m" % ("Samples to prepare: " , ', '.join(sample_list)))
    raw_data_download(addresses,sample_list,outdir)
   ##=====================================================================
    libraries_aggr = pd.DataFrame(columns=['library_id','atac_fragments','per_barcode_metrics','gex_molecule_info'])

    for sample in list(samples.index):
        libraries_aggr.loc[sample,'library_id'] = sample
        libraries_aggr.loc[sample,'atac_fragments'] = os.path.abspath("result/cellranger/"+sample+"/outs/atac_fragments.tsv.gz")
        libraries_aggr.loc[sample,'per_barcode_metrics'] = os.path.abspath("result/cellranger/"+sample+"/outs/per_barcode_metrics.csv")
        libraries_aggr.loc[sample,'gex_molecule_info'] = os.path.abspath("result/cellranger/"+sample+"/outs/gex_molecule_info.h5")
        libraries_aggr.to_csv( f'config/library/libraries_aggr.csv', encoding='utf-8',   sep=",",index=False)

##===============================================================
if __name__ == '__main__':
    run()