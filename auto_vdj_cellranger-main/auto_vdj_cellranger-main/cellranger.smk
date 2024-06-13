from snakemake.utils import validate
import pandas as pd #
import yaml,time,sys,os,glob
from scripts.yaml_prepare import update_local_file

exec(open('/data/software/modules/modules-v4.2.1/init/python.py').read())
###########################################################
configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")
project = config["report"]["Project_Num"]
# reference = config["report"]["Reference"]
Species = config["report"]["Species"]
test_samples = config["cellranger_params"]["test_samples"]
###########################################################
## ==================== 1. load config ====================
###########################################################
def ordered_yaml_load(stream, Loader=yaml.SafeLoader,
                      object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def _construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        _construct_mapping)
    return yaml.load(stream, OrderedLoader)

def ordered_yaml_dump(data, stream=None, Dumper=yaml.SafeDumper,
                      object_pairs_hook=OrderedDict, **kwds):
    class OrderedDumper(Dumper):
        pass

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())
    OrderedDumper.add_representer(object_pairs_hook, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwds)


## checking raw_data address
import re
config_info="config/config.yaml"
with open(config_info) as f:
    config = ordered_yaml_load(f)

project = config['report']['Task_Num']
if config["cellranger_params"]["raw_data_obs_address"] is None or config["cellranger_params"]["raw_data_obs_address"][0] is None or config["cellranger_params"]["raw_data_obs_address"][0] == "":
    print("\033[33;1m"+"Info: No obs address provided. ","\033[0m")
    if config['report']['Task_Num'] is None or config['report']['Task_Num'] == "":
        sys.exit("\033[31;1m"+"ERROR: You need to specify Task_Num in config.yaml if No obs address is provided."+"\033[0m")
    else:
        match = re.search(r'\d{4}', project)
        if match:
            year = match.group()
        else:
            sys.exit("\033[33;1m"+"ERROR: Task_Num format is not supported, check it." + "\033[0m")
        obs_address = "obs://scrna-cloud/samples/" + str(year) + "/" + str(config['report']['Task_Num'] + "/")
        config["cellranger_params"]["raw_data_obs_address"][0] = obs_address
        #temp_config = "config/config_temp.yaml"
        # with open(temp_config,'w',encoding='utf-8') as yw:
        #     yaml.dump(config,yw)
        with open(config_info, 'w', encoding='utf-8') as fp:
            ordered_yaml_dump(config, fp, default_flow_style=False,encoding='utf-8',allow_unicode=True)
        print("\033[33;1m"+"Info: Auto create obs address: "+ obs_address + "\033[0m")
    # if config["cellranger_params"]["samples"] == "all" :
    #     sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
else:
    for addr in config["cellranger_params"]["raw_data_obs_address"]:
        if addr[0:6] == "obs://":
            if addr[-1] !="/":
                sys.exit("\033[31;1m"+"ERROR: cellranger_params:raw_data_obs_address: "+addr+" need to end with '/' !!!"+"\033[0m")
        else: sys.exit("\033[31;1m"+"ERROR: cellranger_params:raw_data_obs_address: "+addr+" need to start with 'obs://' !!!"+"\033[0m")

## get sample names ##  
## from EITHER obs sub directory OR config.yaml input
import os,sys
import pandas as pd
work_dir = os.getcwd()
if type(config["cellranger_params"]["samples"]) == str :
    # type1. "all"
    if config["cellranger_params"]["samples"] == "all" :
        samples = pd.Series()
        for addr in config["cellranger_params"]["raw_data_obs_address"]:
            folders = os.popen("obsutil ls -d -s %s | grep obs  | grep -v '\(md5\|bam_files\)' | sed 's#%s##'"%(addr,addr)).readlines()
            for i in folders:
                if i.strip() != "":
                    samples[i.strip().replace("/","")] = i.strip().replace("/","")
    else :  sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples can only be strings 'all' or a list of sample."+"\033[0m")
else :
    # type2. list of sample:cell_number dictionary
    if all(isinstance(i, str) for i in config["cellranger_params"]["samples"]):
        samples = config["cellranger_params"]["samples"]
    # type3. list of sample name strings
    elif all(isinstance(i, dict) for i in config["cellranger_params"]["samples"]):
        # FC={}
        # for i in config["cellranger_params"]["samples"]:
        #     FC.update(i)
        # samples_list = [str(key)+"_FC"+str(value) for key,value in FC.items()]
        samples = [str(v)+"_"+str(k) for i in config["cellranger_params"]["samples"] for v,k in i.items()]
        FC = {str(v)+"_"+str(k):k for i in config["cellranger_params"]["samples"] for v,k in i.items()}
    else: sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples list can only be strings or <name:cell_number> dictionary."+"\033[0m")
    #samples = pd.Series(samples_list,index = samples_list) # a list of sample name


def get_FC(wildcards): 
    # return FC[widcards.sample]
    try: FC
    except NameError: return ""
    else: return FC[wildcards.sample]

def get_raw_data_dir(wildcards):
    try: FC
    except NameError: return "/raw_data/"+wildcards.sample+"/"+wildcards.sample+"_fastqs"
    else: return "/raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"_".join(wildcards.sample.split("_")[:-1])+"_fastqs"

def get_raw_data_fastqs(wildcards):
    try: FC
    except NameError: return "raw_data/"+wildcards.sample+"/"+wildcards.sample+"_fastqs/"+wildcards.sample+"_S1_L001_R?_001.fastq.gz"
    else: return "raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"_".join(wildcards.sample.split("_")[:-1])+"_fastqs/"+wildcards.sample+"_S1_L001_R?_001.fastq.gz"

##  mutil sample process  
if type(samples) == list:
    Samples = (pd.DataFrame(samples)[0].str.replace("_TCR|_BCR","")).drop_duplicates()
else:
    Samples = (samples.str.replace("_TCR$|_BCR$","")).drop_duplicates()
wildcard_constraints:
    sample="|".join(Samples)

################################################################
## ==================== 2. define all input ====================
################################################################
# def get_sample(wildcards):
#     return wildcards.sample
# # for cellranger
# def get_fastqs(wildcards):
#     """Get raw FASTQ files directory from sheet."""
#     fastqs = samples.loc[(wildcards.sample), ["fastq"]].dropna()
#     return fastqs
# def all_input_cellranger(wildcards):
#     wanted_input = []
#     # cellranger
#     if config['cellranger_params']['module']["cellranger"]:
#         for (sample) in Samples:
#             wanted_input.extend(
#                 expand(
#                     [
#                         "result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",
#                         "result/cellranger/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv",      
#                         "result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_cloupe.cloupe",
#                         "result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix.h5",
#                         #"result/cellranger_sum/{sample}_summary.html",
#                         "logs/{sample}.fastqc.check"
#                     ],
#                     sample=sample
#                 )
#             )
#             wanted_input.extend(
#                 [
#                     #"result/cellranger_sum/metrics_summary_stat.csv",
#                     "result/cellranger_sum/cellranger_summary_email.html",
#                     # "logs/upload_cellranger_sor"
#                 ])
    
#     if config['cellranger_params']['module']["cellranger_report"]:
#         wanted_input.extend(
#             [
#                 f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/质控说明.html",
#                 #"result/report/Cellranger_report_upload.check"
#             ])

#     if config['cellranger_params']['module']["upload_bam"]:
#       for (sample) in Samples:
#         wanted_input.extend(
#             expand(
#                 ["logs/{sample}_upload_bam.log"],
#                 sample=sample
#             )
#       )
#     return wanted_input
rule all:
    input: "vdj_email.check"
        #all_input_cellranger

#############################################################
## ==================== 3. define rules  ====================
#############################################################
rule data_prepare:
    """
    prepare fastqs for each sample
    """
    input:
        # metadata = "config/samples.csv",
    params:
        config = "config/config.yaml"
    output:
        temp_file=temp("raw_data.check")
    benchmark:
        "benchmarks/data_prepare.benchmark.txt"
    log:
        "logs/data_prepare.log"
    resources:
        qsub_mem=20,
        qsub_p=config['cpu']['data_prepare'],
        qsub_n=1
    shell:
        '''
python scripts/data_prepare.py -c {params.config}  -o raw_data/  |& tee {log} && touch raw_data.check
        '''

localrules: cellranger_config_csv_prep
rule cellranger_config_csv_prep:
    """
    getting cellranger config csv
    """ 
    input:
        temp_file="raw_data.check"
    output:
        config_csv_file="config/config_csv/{sample}_config.csv"
    params:
        scrna_reference = config["report"]["scrna_Reference"],
        include_introns = config["cellranger_params"]["include_introns"],
        generate_bam = config["cellranger_params"]["generate_bam"],
        chemistry = config["report"]["Others"],
        vdj_reference = config["report"]["vdj_Reference"],
        force_cells_num = get_FC,
        fastqs_load = get_raw_data_dir
    run:
        config_csv = pd.DataFrame(columns=['[gene-expression]'])
        config_csv.loc['scrna_reference'] = f"reference,{params.scrna_reference}"
        # include_introns
        if params.include_introns == True:
            config_csv.loc['include_introns_param'] = f"include-introns,TRUE"
        else:
            config_csv.loc['include_introns_param'] = f"include-introns,FALSE"

        # generate_bam
        if params.generate_bam == False:
            config_csv.loc['generate_bam_param'] = "no-bam,TRUE"
        else:
            config_csv.loc['generate_bam_param'] = "no-bam,FALSE"

        # chemistry
        if params.chemistry == "5' scRNA":
            config_csv.loc['chemistry'] = "chemistry,fiveprime"
        elif chemistry == "3' scRNA":
            config_csv.loc['chemistry'] = "chemistry,threeprime"
        else:
            sys.exit("\033[31;1m"+"ERROR: chemistry params (report:Others) can only be either 3' scRNA or 5' scRNA."+"\033[0m")
        #force cell
        try:FC
        except NameError: ""
        else:config_csv.loc['force-cells'] = f"force-cells,{params.force_cells_num}"
        
        config_csv.loc['vdj_flag'] = "[vdj]"
        config_csv.loc['vdj_reference'] = f"reference,{params.vdj_reference}"
        config_csv.loc['libraries_flag'] = f"[libraries]"
        config_csv.loc['libraries_header'] = f"fastq_id,fastqs,lanes,feature_types,subsample_rate"
        config_csv.loc['scrna_library'] = f"{wildcards.sample},{work_dir}{params.fastqs_load},any,gene expression,"
        try:FC
        except NameError: 
            config_csv.loc['scrna_library'] = f"{wildcards.sample},{work_dir}{params.fastqs_load},any,gene expression,"
            if f"{wildcards.sample}_TCR" in samples:
                sam=f"{wildcards.sample}_TCR"
                config_csv.loc['TCR_library'] = f"{sam},{work_dir}/raw_data/{sam}/{sam}_fastqs,any,VDJ-T,"
            if f"{wildcards.sample}_BCR" in samples:
                sam=f"{wildcards.sample}_BCR"
                config_csv.loc['BCR_library'] = f"{sam},{work_dir}/raw_data/{sam}/{sam}_fastqs,any,VDJ-B,"
            if f"{wildcards.sample}_TCR" not in samples and f"{wildcards.sample}_BCR" not in samples:
                sys.exit("\033[31;1m"+"ERROR: {wildcards.sample} 没有 vdj data, 因此不适合进行 cellrnager mutil 分析!!!"+"\033[0m")
            # sep="\t", 去除多余的双引号 
        else:
            raw_fq_name = "_".join(wildcards.sample.split("_")[:-1])
            config_csv.loc['scrna_library'] = f"{raw_fq_name},{work_dir}{params.fastqs_load},any,gene expression,"
            TCR = wildcards.sample.rsplit("_", 1)[0] +"_TCR_" +wildcards.sample.rsplit("_", 1)[1]
            BCR = wildcards.sample.rsplit("_", 1)[0] +"_BCR_" +wildcards.sample.rsplit("_", 1)[1]
            if TCR in samples:
                sam=wildcards.sample.rsplit("_", 1)[0] + "_TCR"
                config_csv.loc['TCR_library'] = f"{sam},{work_dir}/raw_data/{sam}/{sam}_fastqs,any,VDJ-T,"
            if BCR in samples:
                sam=wildcards.sample.rsplit("_", 1)[0] + "_BCR"
                config_csv.loc['BCR_library'] = f"{sam},{work_dir}/raw_data/{sam}/{sam}_fastqs,any,VDJ-B,"
            if TCR not in samples and BCR not in samples:
                sys.exit("\033[31;1m"+"ERROR: {wildcards.sample} 没有 vdj data, 因此不适合进行 cellrnager mutil 分析!!!"+"\033[0m")
        
        config_csv.to_csv(output.config_csv_file, sep="\t", encoding='utf-8', index=False)
# ========================================================================================================================
rule cellranger:
    """
    run cellranger for each samples 
    """
    input:
        config_csv_file="config/config_csv/{sample}_config.csv"
    output:
        reference_json="result/cellranger/{sample}/outs/vdj_reference/reference.json",
        regions_fa="result/cellranger/{sample}/outs/vdj_reference/fasta/regions.fa",
        web_summary="result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",
        metrics_summary="result/cellranger/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv",      
        cloupe="result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_cloupe.cloupe",
        filtered="result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix.h5",
        check="{sample}_cellranger.check",
        tmpfiles = temp(["result/cellranger/{sample}/_cmdline",
                         "result/cellranger/{sample}/_filelist",
                         "result/cellranger/{sample}/_finalstate",
                         "result/cellranger/{sample}/_invocation",
                         "result/cellranger/{sample}/_jobmode",
                         "result/cellranger/{sample}/_log",
                         "result/cellranger/{sample}/_mrosource",
                         "result/cellranger/{sample}/_perf",
                         "result/cellranger/{sample}/_perf._truncated_",
                         "result/cellranger/{sample}/_sitecheck",
                         "result/cellranger/{sample}/_tags",
                         "result/cellranger/{sample}/_timestamp",
                         "result/cellranger/{sample}/_uuid",
                         "result/cellranger/{sample}/_vdrkill",
                         "result/cellranger/{sample}/_versions",
                         directory("result/cellranger/{sample}/SC_MULTI_CS"),
                         "result/cellranger/{sample}/{sample}.mri.tgz"])
    params:
        fastq = get_raw_data_dir,
        force_cells = get_FC,
        cellranger_version = config['cellranger_params']['envmodules'].strip().replace("cellRanger/","")
    benchmark:
        "benchmarks/QC/{sample}.cellranger.benchmark.txt"
    log:
        "logs/{sample}_cellranger.log"
    resources:
        qsub_mem=64,
        qsub_p=config['cpu']['cellranger'],
        qsub_n=1
    envmodules:
        config["cellranger_params"]['envmodules']
    shell:
        '''
        export PATH=/opt/softwares/cellranger-{params.cellranger_version}:$PATH ;
        cd result/cellranger/ &&
        rm -rf {wildcards.sample} &&
        cellranger multi --id={wildcards.sample} \
        --csv=../../{input.config_csv_file} \
        --localmem={resources.qsub_mem} \
        --localcores={resources.qsub_p} \
        --description={project}_{Species} |& tee ../../{log} && 
          cd ../../ && touch {wildcards.sample}_cellranger.check && 
          if [ -d raw_data/{wildcards.sample} ]; then rm -r raw_data/{wildcards.sample} ;fi          
          if [ -d raw_data/{wildcards.sample}_BCR ]; then rm -r raw_data/{wildcards.sample}_BCR ;fi          
          if [ -d raw_data/{wildcards.sample}_TCR ]; then rm -r raw_data/{wildcards.sample}_TCR ;fi          
        '''
#====================================================================================================
if config["cellranger_params"]["module"]["upload_bam"]:
    
    #localrules: bam_prepare
    rule bam_prepare:
        """
        prepare bam files for each sample
        """
        input:
            check = "{sample}_bam_moved.check",
            #web_summary = "result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",
            #bam_bai = "result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/possorted_genome_bam.bam.bai"
        output:
            check = "{sample}_upload_bam.check"
        log:
            "logs/{sample}_upload_bam.log"
        params:
            Task_Num = config['report']['Task_Num']
        resources:
            qsub_mem=20,
            qsub_p=config['cpu']['bam_prepare'],
            qsub_n=1
        shell:
            '''
bash /public/cloud_scRNA/chenhaoruo/scripts/auto-email/upload_bam.sh "$PWD/result/cellranger_sum/{wildcards.sample}" {params.Task_Num} |& tee {log}
touch {output.check}
            '''
# ========================================================================================================================
#localrules: cellranger_sample_stat
my_array = " ".join(samples)
try:FC
except NameError : my_array = " ".join(samples)
else: my_array = " ".join([item.replace('_BCR', '').replace('_TCR', '') + '_BCR' if '_BCR' in item else item.replace('_BCR', '').replace('_TCR', '') + '_TCR' if '_TCR' in item else item for item in samples])
rule cellranger_sample_stat:
    input:
        web_summary="result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",
        metrics_summary="result/cellranger/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv"
    output:
        sample_summary="result/cellranger_sum/{sample}_summary.html",
        sample_metrics="result/cellranger_sum/{sample}_metrics_summary.csv",
        check = "{sample}_bam_moved.check"
    params:
        bam_dir="result/cellranger_sum/{sample}/outs/"
    benchmark:
        "benchmarks/QC_sum/{sample}.cellranger_sample_stat.benchmark.txt"
    log:
       "logs/{sample}.cellranger_sample_stat.log"
    envmodules:
         config["cellranger_params"]['envmodules']
    resources:
        qsub_mem=4,
        qsub_p=1,
        qsub_n=1
    shell:
        '''
            cp {input.web_summary}  {output.sample_summary}
            cp {input.metrics_summary}  {output.sample_metrics}
            python scripts/selenium_screenshot.py -i {input.web_summary} -n {wildcards.sample} -o result/cellranger/ |& tee {log}
            if echo {my_array} | grep -q "{wildcards.sample}_TCR" ; then
                python scripts/selenium_screenshot.py -i {input.web_summary} -n "{wildcards.sample}_TCR" -c TCR -o result/cellranger/ |& tee {log}
            fi            
            if echo {my_array} | grep -q "{wildcards.sample}_BCR" ; then
                python scripts/selenium_screenshot.py -i {input.web_summary} -n "{wildcards.sample}_BCR" -c BCR -o result/cellranger/ |& tee {log}
            fi  
            if [[ ! -d {params.bam_dir} ]]; then mkdir -p {params.bam_dir} ;fi
            if [[ -f "result/cellranger/{wildcards.sample}/outs/per_sample_outs/{wildcards.sample}/count/sample_alignments.bam" ]]; then
                    mv "result/cellranger/{wildcards.sample}/outs/per_sample_outs/{wildcards.sample}/count/sample_alignments.bam" {params.bam_dir} &&
                    mv "result/cellranger/{wildcards.sample}/outs/per_sample_outs/{wildcards.sample}/count/sample_alignments.bam.bai" {params.bam_dir} 
            fi && touch {output.check} |& tee -a {log}
         
        '''

# ========================================================================================================================
#localrules: cellranger_check
rule cellranger_check:
    input:
        samples_summary=expand("result/cellranger_sum/{sample}_summary.html",sample=Samples),
        metrics_summary=expand("result/cellranger_sum/{sample}_metrics_summary.csv",sample=Samples)
    output:
        check = "result/cellranger_sum/cellranger_summary_email.html",
        tempdir = temp(directory(project + "_cellranger_summary")),
        tempzip = temp(project + "_cellranger_summary.zip"),
    log:
        "logs/cellranger_check.log"
    params:
        outdir = "result/cellranger_sum/"
    resources:
        qsub_mem=10,
        qsub_p=1,
        qsub_n=1
    run:
        total_size = 0
        for file in input.samples_summary:
            total_size += os.path.getsize(file)
        if total_size > 45 * 1024 * 1024:
            print(f"附件的总大小为 {total_size/(1024*1024):.2f} MB，超过了45MB。")
            import shutil
            sum_dir = project + "_cellranger_summary"
            os.mkdir(sum_dir)
            for file in input.samples_summary:
                shutil.copy(file, sum_dir)
            shutil.make_archive(sum_dir,'zip',root_dir="./", base_dir= sum_dir)
            samples_html = sum_dir + ".zip"
            if os.path.getsize(samples_html) > 45 * 1024 * 1024:
                print(f"压缩包的总大小为 {os.path.getsize(samples_html)/(1024*1024):.2f} MB，超过了45MB。")
                samples_html = input.metrics_summary_stat
        else:
            print(f"附件的总大小为 {total_size/(1024*1024):.2f} MB，未超过50MB。")
            samples_html = ",".join(input.samples_summary)
            shell("touch {project}_cellranger_summary.zip && mkdir {project}_cellranger_summary")
        samples_metrics = ",".join(input.metrics_summary)
        shell("python scripts/cellranger_check.py --config config/config.yaml --stat {samples_metrics} --html {samples_html}  --outdir {params.outdir} |& tee {log}")

# ========================================================================================================================
folder_list = glob.glob(f"result/report/" + config['report']['Task_Num'] + f"_QC_Report_*")
if folder_list:
    report_dir = folder_list[0]
else:
    report_dir = f"result/report/" + config['report']['Task_Num'] + f"_QC_Report_" + time.strftime('%Y_%m_%d')

if config["cellranger_params"]["module"]["cellranger_report"]:
    rule cellranger_report:
        """
        cellranger report
        """
        input:
            web_summary=expand("result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",sample=Samples),
            cloupe=expand("result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_cloupe.cloupe",sample=Samples),
            check = expand("{sample}_bam_moved.check", sample=Samples)
        output:
            html = report_dir+"/质控说明.html"
        params:
            report = report_dir+"/",
            config = "config/config.yaml",
            write_todelete = "/public/cloud_scRNA/chenhaoruo/scripts/write_todelete.sh"
        log:
            "logs/Cellranger_report.log"
        benchmark:
            "benchmarks/Cellranger.report.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=1
        shell:
            '''
export CLOUD=True
export CLOUD_URL=http://192.168.20.201:9999
export CLOUDUSER=qingzhen.zi@oebiotech.com
export PASSWORD=quanzhilong0818
            rm -rf {params.report}
            python scripts/report/cellranger_mutil_vdj_report.py \
            -i result \
            -c {params.config} |& tee {log} && bash {params.write_todelete} $PWD;
            '''

# # ========================================================================================================================
# # ========================================================================================================================
# rule fastqc_check: 
#     input:
#         #stat =  "result/cellranger/{sample}/outs/metrics_summary.csv",
#         stat= "result/cellranger/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv"
#         #fastq = "raw_data/{sample}/" # get_fastqs
#     output:
#         "logs/{sample}.fastqc.check"
#     params:
#         output_dir = "./result/Fastqc/{sample}/",
#         fastqs = get_raw_data_fastqs 
#     resources:
#         qsub_mem=20,
#         qsub_p=config['cpu']['fastqc_check'],
#         qsub_n=1
#     run:
#         import pandas as pd
#         table = pd.read_csv(input.stat,index_col=0,header = 0)
#         Q30 = table[table['Metric Name']=='Q30 RNA read']['Metric Value'][0]
#         #Mapping = table.loc[0,"Reads Mapped Confidently to Transcriptome"]
#         Mapping = table[table['Metric Name']=='Confidently mapped to transcriptome' ]['Metric Value'][0]

#         if any([float(Q30[:-1]) < 80, float(Mapping[:-1]) < 30 ]):
#             rerun = True
#             shell("mkdir -p {params.output_dir}")
#             shell("module load FastQC/0.11.8 && fastqc  {params.fastqs} -t {resources.qsub_p}  -o {params.output_dir} && touch {output}")
            
#         else:
#             rerun = False
#             shell("echo 'Q30 in RNA read is {Q30}. [√] \nMapping rate to transcriptome is: {Mapping}. [√]' > {output}")

#         shell('rm -rf raw_data/{wildcards.sample}*')
# ========================================================================================================================
if not config["cellranger_params"]["module"]["cellranger_report"] and len(test_samples) == 0 :
    rule cellranger_upload:
        """
        cellranger upload
        """
        input:
            web_summary=expand("result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",sample=Samples),
            check = expand("{sample}_bam_moved.check", sample=Samples)
        output:
            check = "cellranger_processing.check"
        params:
            project = project,
            write_todelete = "/public/cloud_scRNA/chenhaoruo/scripts/write_todelete.sh"
        log:
            "logs/cellranger_upload.log"
        benchmark:
            "benchmarks/cellranger_upload.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=1
        shell:
            '''
bash /public/cloud_scRNA/chenhaoruo/scripts/auto-email/upload_cellranger.sh result/cellranger {params.project} |& tee {log} && bash {params.write_todelete} $PWD && touch {output.check}
            '''

if len(test_samples) > 0 :
    rule test_samples_processing:
        """
        processing test samples
        """
        input:
            web_summary=expand("result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",sample=test_samples),
            check = expand("{sample}_bam_moved.check", sample=test_samples)
        output:
            check = "cellranger_processing.check"
        log:
            "logs/test_samples_processing.log"
        benchmark:
            "benchmarks/test_samples_processing.benchmark.txt"
        run:
            for I in test_samples:
                update_local_file([I],"cellranger_config.csv","Sample","finish")
                update_local_file([I + f"_TCR"],"cellranger_config.csv","Sample","wait")
                update_local_file([I + f"_BCR"],"cellranger_config.csv","Sample","wait")
            shell("touch {output.check}")
        
# ========================================================================================================================
localrules: update_status
rule update_status: 
    input:
        cellranger_result=expand("{sample}_cellranger.check",sample=Samples),
        #metrics_summary_stat="result/cellranger_sum/metrics_summary_stat.csv",
        check = "result/cellranger_sum/cellranger_summary_email.html",
        bam_upload = expand("{sample}_upload_bam.check",sample=Samples) if config["cellranger_params"]["module"]["upload_bam"] else "config/config.yaml",
        html = report_dir+"/质控说明.html" if config["cellranger_params"]["module"]["cellranger_report"] else "cellranger_processing.check",
        sample_stat = expand("logs/{sample}.cellranger_sample_stat.log", sample=Samples),
        #upload_cellranger_sor = "logs/cellranger_sor.log",
        #update_sql = "scripts/update_sql.py"
    output:
        check = "vdj_email.check"
    params:
        all_samples = [config["cellranger_params"]["samples"]],
        project = project
    log:
        "logs/update_status.log"
    resources:
        qsub_mem=10,
        qsub_n=1
    run:
        test_samples_all = test_samples + [s + "_TCR" if not s.endswith("_TCR") else s for s in test_samples] + [s + "_BCR" if not s.endswith("_BCR") else s for s in test_samples]
        for I in samples:
            if I not in test_samples_all:
                update_local_file([I],"cellranger_config.csv","Sample","finish")
        update_local_file([params.project],"email_config.csv","rwdh","finish")
        shell("touch {output.check}")





