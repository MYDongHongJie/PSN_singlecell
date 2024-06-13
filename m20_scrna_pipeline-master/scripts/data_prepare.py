import os,sys,yaml,subprocess,click
from glob import glob
from copy import deepcopy
import pandas as pd
##function =============================================================================================================
def get_obs_sample(address):
    '''
    return all the sample name in a list, given a obs address
    '''
    folders = os.popen("/data/software/obsutil/5.2.12/obsutil ls -d -s %s | grep obs  | grep -v '\(md5\|bam_files\)' | sed 's#%s##'"%(address,address)).readlines()
    return [i.strip().replace("/","") for i in folders if i.strip() != "" ]

def get_obs_subdir(address,sample_list):
    '''
    return VALID obs subdirectory list (at sample level), given a obs address and a sample list.
    '''
    remote = [i.strip() for i in os.popen("/data/software/obsutil/5.2.12/obsutil ls -d -s %s | grep obs | grep -v '\(md5\|bam_files\)' | grep '%s.*/' "%(address,address)).readlines() ]
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
## data download =======================================================================================================
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
## rename data =========================================================================================================
def raw_data_formart(sample_list,outdir,raw2clean):
    for sample in sample_list:
        if sample in os.listdir(f"{outdir}/"):
            if not os.path.exists(f'{outdir}/{sample}/{sample}_fastqs'):
                os.mkdir(f'{outdir}/{sample}/{sample}_fastqs')
            if not os.path.exists(f'config/raw2clean'):
                os.mkdir(f'config/raw2clean')
            ########################
            raw2clean_outfile = os.path.abspath(f'{outdir}/{sample}/{sample}_fastqs')
            raw2clean_file = ['pathOfAnchorbcm=/home/sunkun/m20/VITASEER',
                              f'outdir={raw2clean_outfile}',
                              f'sampleID={sample}']
            for i in {"1", "2"}:
                fq_list = sorted([item for sublist in [glob( f'{outdir}/{sample}'+ ext) for ext in [f"/*R{i}.f*q.gz", f"/*_{i}.fq.gz", f"/*R{i}_001.f*q.gz"]] for item in sublist])
                print(fq_list)
                if len(fq_list) == 1 :
                    src = fq_list[0].replace(f"{outdir}/{sample}", "..")
                    dst = f'{outdir}/{sample}/{sample}_fastqs/{sample}_R{i}.fq.gz'
                    if raw2clean:
                        fastq_file=os.path.abspath(fq_list[0])
                        raw2clean_file.append(f'fastq{i}={fastq_file}')
                    else:
                        os.symlink(src, dst)
                else:
                    dir_path=os.path.abspath(f'{outdir}/{sample}')
                    subprocess.run(['zcat']+fq_list,
                                   stdout=open(f'{outdir}/{sample}/{sample}_merger_R{i}.fastq', 'w'))
                    subprocess.call(f'pigz -f {outdir}/{sample}/{sample}_merger_R{i}.fastq',shell=True)
                    src_merger = f'{outdir}/{sample}/{sample}_merger_R{i}.fastq.gz'
                    src_file = src_merger.replace(f"{outdir}/{sample}", "..")
                    dst_file = f'{outdir}/{sample}/{sample}_fastqs/{sample}_R{i}.fq.gz'
                    if raw2clean:
                        fastq_file=os.path.abspath(src_merger)
                        raw2clean_file.append(f'fastq{i}={fastq_file}')
                    else:
                        os.symlink(src_file, dst_file)
            if raw2clean:
                raw2clean_export = pd.DataFrame(raw2clean_file)
                raw2clean_export.to_csv(f'config/raw2clean/{sample}_config.info',sep="\t",header=False,index=False)

        else:
            print("\033[0;31m%s\033[0m" % f"The fastq files of {sample} does not exist, please check it manually!!")

##run ==================================================================================================================
@click.command()
@click.option('-c', 'configfile', type=str,default = 'config/config.yaml',
    help='config yaml files with raw_data_obs_address, samples information. df:config/config.yaml')
@click.option('-o', 'outdir', type=str, default = 'raw_data',
    help='Output directory, default: raw_data')

def run(configfile,outdir):
    """example: python data_prepare.py -c config/config.yaml -o raw_data """
    with open(configfile) as f:
        config = yaml.safe_load(f)
    addresses = config["params"]["STARsolo"]['raw_data_obs_address']
    samples_config = config["params"]["STARsolo"]['samples']
    raw2clean = config['params']['STARsolo']['raw_data']
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
                [sample_list.append(i) for i in get_obs_sample(address) if not i in sample_list  ]
    else:
        sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples can only be string 'all' or a list of sample names."+"\033[0m")
    print("\033[0;32m%s%s\033[0m" % ("Samples to prepare: " , ', '.join(sample_list)))

    ## Step1: download raw data
    raw_data_download(addresses,sample_list,outdir)

    ## Step2: formart fastqs name
    raw_data_formart(sample_list,outdir,raw2clean)

if __name__ == '__main__':
    run()
