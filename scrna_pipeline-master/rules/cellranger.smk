###########################################################
## ==================== 1. load config ====================
###########################################################
## checking raw_data address

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


import re
if config["cellranger_params"]["raw_data_obs_address"] is None or config["cellranger_params"]["raw_data_obs_address"][0] is None or config["cellranger_params"]["raw_data_obs_address"][0] == "":
    print("\033[33;1m"+"Info: No obs address provided. ","\033[0m")
    if config['report']['Task_Num'] is None or config['report']['Task_Num'] == "":
        sys.exit("\033[31;1m"+"ERROR: You need to specify Task_Num in config.yaml if No obs address is provided."+"\033[0m")
    else:
        project = config['report']['Task_Num']
        match = re.search(r'\d{4}', project)
        if match:
            year = match.group()
        else:
            sys.exit("\033[33;1m"+"ERROR: Task_Num format is not supported, check it." + "\033[0m")
        obs_address = "obs://scrna-cloud/samples/" + str(year) + "/" + str(config['report']['Task_Num'] + "/")
        config["cellranger_params"]["raw_data_obs_address"][0] = obs_address
        print("\033[33;1m"+"Info: Auto create obs address: "+ obs_address + "\033[0m")
    # if config["cellranger_params"]["samples"] == "all" :
    #     sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
else:
    for addr in config["cellranger_params"]["raw_data_obs_address"]:
        if addr[0:6] == "obs://":
            if addr[-1] !="/":
                sys.exit("\033[31;1m"+"ERROR: cellranger_params:raw_data_obs_address: "+addr+" need to end with '/' !!!"+"\033[0m")
        else: sys.exit("\033[31;1m"+"ERROR: cellranger_params:raw_data_obs_address: "+addr+" need to start with 'obs://' !!!"+"\033[0m")
temp_config = "config/config_temp.yaml"
with open(temp_config,'w',encoding='utf-8') as yw:
    ordered_yaml_dump(config, yw, default_flow_style=False,encoding='utf-8',allow_unicode=True)

## get sample names ##  
## from EITHER obs sub directory OR config.yaml input
import os,sys
import pandas as pd
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
    else: return "--force-cells=" + str(FC[wildcards.sample])

def get_raw_data_dir(wildcards):
    try: FC
    except NameError: return "../../raw_data/"+wildcards.sample+"/"+wildcards.sample+"_fastqs"
    else: return "../../raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"_".join(wildcards.sample.split("_")[:-1])+"_fastqs"

def get_raw_data_fastqs(wildcards):
    try: FC
    except NameError: return "raw_data/"+wildcards.sample+"/"+wildcards.sample+"_fastqs/"+wildcards.sample+"_S1_L001_R?_001.fastq.gz"
    else: return "raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"_".join(wildcards.sample.split("_")[:-1])+"_fastqs/"+wildcards.sample+"_S1_L001_R?_001.fastq.gz"

wildcard_constraints:
    sample="|".join(samples)

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
def all_input_cellranger(wildcards):
    wanted_input = []
    # cellranger
    if config['cellranger_params']['module']["cellranger"]:
        for (sample) in samples:
            wanted_input.extend(
                expand(
                    [
                        "result/cellranger/{sample}/outs/web_summary.html",
                        "result/cellranger/{sample}/outs/metrics_summary.csv",
                        "result/cellranger/{sample}/outs/cloupe.cloupe",
                        "result/cellranger/{sample}/outs/raw_feature_bc_matrix.h5",
                        "result/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
                        "result/cellranger_sum/{sample}_summary.html",
                        "logs/{sample}.fastqc.check"
                    ],
                    sample=sample
                )
            )
            wanted_input.extend(
                [
                #    "result/cellranger_sum/metrics_summary_stat.csv",
                    "result/cellranger_sum/cellranger_check_log.txt",
                    "logs/upload_cellranger_sor"
                ])
    if config['cellranger_params']['module']["cellranger_report"]:
        wanted_input.extend(
            [
                "logs/Cellranger_report.log"
            ])
    if config['cellranger_params']['module']["upload_bam"]:
      for (sample) in samples:
        wanted_input.extend(
            expand(
                ["logs/{sample}_upload_bam.log"],
                sample=sample
            )
      )
    if config['cellranger_params']['module']["cellranger_upload"]:
      for (sample) in samples:
        wanted_input.extend(
            expand(
                ["logs/cellranger_upload.log"],
                sample=sample
            )
      )
    return wanted_input

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
        config = "config/config_temp.yaml"
        # obs_address = config["cellranger_params"]["raw_data_obs_address"][0]
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
python scripts/data_prepare.py -c {params.config}  -o raw_data/ |& tee {log} && touch raw_data.check
        '''
# ========================================================================================================================
# param1
include_introns = config["cellranger_params"]["include_introns"]
cellranger_version = config['cellranger_params']['envmodules'].strip().replace("cellRanger/","")
if include_introns==True :
    if cellranger_version == "7.0.1":
        include_introns_param = "--include-introns=true"
    else:
        include_introns_param = "--include-introns"
else:
    if cellranger_version == "7.0.1":
        include_introns_param = "--include-introns=false"
    else:
        include_introns_param = ""

# param2
chemistry = config["report"]["Others"]
if chemistry == "5' scRNA":
    chemistry_params = "fiveprime"
elif chemistry == "3' scRNA":
    chemistry_params = "threeprime"
else:
    sys.exit("\033[31;1m"+"ERROR: chemistry params (report:Others) can only be either 3' scRNA or 5' scRNA."+"\033[0m")
# param3
generate_bam = config["cellranger_params"]["generate_bam"]
if generate_bam == False:
    generate_bam_param = "--no-bam"
else:
    generate_bam_param = ""



rule cellranger:
    """
    run cellranger for each samples 
    """
    input:
        temp_file="raw_data.check"
    output:
        web_summary="result/cellranger/{sample}/outs/web_summary.html",
        metrics_summary="result/cellranger/{sample}/outs/metrics_summary.csv",
        cloupe="result/cellranger/{sample}/outs/cloupe.cloupe",
        raw="result/cellranger/{sample}/outs/raw_feature_bc_matrix.h5",
        filtered="result/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
        tmpfiles = temp(["result/cellranger/{sample}/_cmdline",
                         "result/cellranger/{sample}/_filelist",
                         "result/cellranger/{sample}/_finalstate",
                         "result/cellranger/{sample}/_invocation",
                         "result/cellranger/{sample}/_jobmode",
                         "result/cellranger/{sample}/_log",
                         "result/cellranger/{sample}/_mrosource",
                         "result/cellranger/{sample}/_perf",
                         "result/cellranger/{sample}/_sitecheck",
                         "result/cellranger/{sample}/_tags",
                         "result/cellranger/{sample}/_timestamp",
                         "result/cellranger/{sample}/_uuid",
                         "result/cellranger/{sample}/_vdrkill",
                         "result/cellranger/{sample}/_versions",
                         directory("result/cellranger/{sample}/SC_RNA_COUNTER_CS"),
                         "result/cellranger/{sample}/{sample}.mri.tgz"])
    params:
        fastq = get_raw_data_dir,
        force_cells = get_FC
    benchmark:
        "benchmarks/QC/{sample}.cellranger.benchmark.txt"
    resources:
        qsub_mem=64,
        qsub_p=config['cpu']['cellranger'],
        qsub_n=1
    envmodules:
        config["cellranger_params"]['envmodules']
    shell:
        '''
cd result/cellranger/ &&
rm -rf {wildcards.sample} &&
cellranger count \
  --id={wildcards.sample}   \
  --transcriptome={reference}  \
  --fastqs={params.fastq}    \
  --localcores={resources.qsub_p}  \
  --localmem={resources.qsub_mem}   \
  --chemistry={chemistry_params} \
  --description={project}_{Species} {include_introns_param} {generate_bam_param} {params.force_cells}
        '''
# ========================================================================================================================
localrules: cellranger_sample_stat
rule cellranger_sample_stat:
    input:
        web_summary="result/cellranger/{sample}/outs/web_summary.html",
        metrics_summary="result/cellranger/{sample}/outs/metrics_summary.csv"
    output:
        sample_summary="result/cellranger_sum/{sample}_summary.html",
        sample_metrics=temp("result/cellranger_sum/{sample}_metrics_summary.csv"),
        sample_png="result/cellranger/{sample}.png"
    benchmark:
        "benchmarks/QC_sum/{sample}.cellranger_sample_stat.benchmark.txt"
    log:
       "logs/{sample}.cellranger_sample_stat.log"
    envmodules:
        "OESingleCell/v_3.0.0_visium_produce"
    resources:
        qsub_mem=4,
        qsub_p=1,
        qsub_n=1
    shell:
        '''
oest_html_screenshot.py  -i {input.web_summary} -n {wildcards.sample} -o result/cellranger/ |& tee {log}
cp {input.web_summary}  {output.sample_summary}
echo -e "sampleid\\n{wildcards.sample}" | paste - {input.metrics_summary} | sed 's/\\t/,/g' - >  {output.sample_metrics}
        '''
# ========================================================================================================================
localrules: cellranger_summary
rule cellranger_summary:
    """
    cellranger stats
    """
    input:
        metrics_summary=expand("result/cellranger_sum/{sample}_metrics_summary.csv",sample=samples)
    output:
        metrics_summary_stat="result/cellranger_sum/metrics_summary_stat.csv",
    benchmark:
        "benchmarks/QC_sum/cellranger_summary.benchmark.txt"
    log:
        "logs/cellranger_summary.log"
    resources:
        qsub_mem=1,
        qsub_p=1,
        qsub_n=1
    shell:
        """
cat {input.metrics_summary} | awk '!a[$0]++' > {output.metrics_summary_stat} |& tee {log}
        """
# ========================================================================================================================
localrules: cellranger_check
rule cellranger_check:
    input:
        metrics_summary_stat = "result/cellranger_sum/metrics_summary_stat.csv",
        samples_summary=expand("result/cellranger_sum/{sample}_summary.html",sample=samples),
        sample_png=expand("result/cellranger/{sample}.png",sample=samples)
    output:
        check = "result/cellranger_sum/cellranger_check_log.txt"
    #log:
    #    "logs/cellranger/cellranger_check.log"
    resources:
        qsub_mem=10,
        qsub_p=1,
        qsub_n=1
    run:
        samples_html = ",".join(input.samples_summary)
        samples_png = ",".join(input.sample_png)
        shell("python scripts/cellranger_check.py --config config/config.yaml --stat {input.metrics_summary_stat} --html {samples_html} --png {samples_png} --output {output.check}")
# ========================================================================================================================
localrules: cellranger_sor
rule cellranger_sor:
    """
    upload cellranger SOR
    """
    input:
        metrics_summary_stat="result/cellranger_sum/metrics_summary_stat.csv",
        #metadata = "config/samples.csv",
    output:
        directory("logs/upload_cellranger_sor")
    params:
        config = "config/config.yaml"
    benchmark:
        "benchmarks/cellranger_sor.benchmark.txt"
    resources:
        qsub_mem=4,
        qsub_p=1
    envmodules:
         config["report_params"]['envmodules']
    shell:
        """
python  scripts/scrna_sor.py \
  -i {params.config} \
  -l logs/upload_cellranger_sor \
  -s "Project,QC,Cpu_Mem,Ref_database,Software"
        """

# ========================================================================================================================
rule fastqc_check: 
    input:
        stat =  "result/cellranger/{sample}/outs/metrics_summary.csv",
        #fastq = "raw_data/{sample}/" # get_fastqs
    output:
        "logs/{sample}.fastqc.check"
    params:
        output_dir = "./result/Fastqc/{sample}/",
        fastqs = get_raw_data_fastqs 
    resources:
        qsub_mem=20,
        qsub_p=config['cpu']['fastqc_check'],
        qsub_n=1
    run:
        import pandas as pd
        table = pd.read_csv( input.stat,header = 0)
        Q30 = table.loc[0,"Q30 Bases in RNA Read"]
        Mapping = table.loc[0,"Reads Mapped Confidently to Transcriptome"]
        if any([float(Q30[:-1]) < 80, float(Mapping[:-1]) < 30 ]):
            rerun = True
            shell("mkdir -p {params.output_dir}")
            shell("module load FastQC/0.11.8 && fastqc  {params.fastqs} -t {resources.qsub_p}  -o {params.output_dir} && touch {output}")
            
        else:
            rerun = False
            shell("echo 'Q30 in RNA read is {Q30}. [√] \nMapping rate to transcriptome is: {Mapping}. [√]' > {output}")

        shell("rm -r raw_data/{wildcards.sample}")


# ========================================================================================================================
if config["cellranger_params"]["module"]["upload_bam"]:
    
    localrules: bam_prepare
    rule bam_prepare:
        """
        prepare bam files for each sample
        """
        input:
            web_summary=expand("result/cellranger/{sample}/outs/web_summary.html",sample=samples)
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
cd result/cellranger/ &&
sh  ../../scripts/upload_bam.sh {wildcards.sample} {params.Task_Num}
            '''
# ========================================================================================================================
if config["cellranger_params"]["module"]["cellranger_report"]:
    rule cellranger_report:
        """
        cellranger report
        """
        input:
            web_summary=expand("result/cellranger/{sample}/outs/web_summary.html",sample=samples)
        output:
            html = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/质控说明.html"
        params:
            report = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/",
            config = "config/config.yaml"
        log:
           "logs/Cellranger_report.log"
        benchmark:
            "benchmarks/Cellranger.report.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=1
        shell:
            '''
rm -rf {params.report}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/Cellranger_report.py \
  -i result \
  -c {params.config} |& tee {log};
            '''

# ========================================================================================================================
if config["cellranger_params"]["module"]["cellranger_upload"]:
    localrules: cellranger_upload
    rule cellranger_upload:
        """
        cellranger upload
        """
        input:
            web_summary=expand("result/cellranger/{sample}/outs/web_summary.html",sample=samples),
        params:
            Task_Num = config['report']['Task_Num']
        log:
            "logs/cellranger_upload.log"
        benchmark:
            "benchmarks/cellranger_upload.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=1
        shell:
            '''
sh  scripts/upload_cellranger_m.sh result/cellranger {params.Task_Num}
            '''
