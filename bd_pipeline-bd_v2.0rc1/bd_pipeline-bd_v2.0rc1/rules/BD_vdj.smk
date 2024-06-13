###########################################################
## ==================== 1. load config ====================
###########################################################
## checking raw_data address
import re
if config["BD_params"]["raw_data_obs_address"] is None or config["BD_params"]["raw_data_obs_address"][0] is None or config["BD_params"]["raw_data_obs_address"][0] == "":
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
        config["BD_params"]["raw_data_obs_address"][0] = obs_address
        print("\033[33;1m"+"Info: Auto create obs address: "+ obs_address + "\033[0m")
    # if config["BD_params"]["samples"] == "all" :
    #     sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
else:
    for addr in config["BD_params"]["raw_data_obs_address"]:
        if addr[0:6] == "obs://":
            if addr[-1] !="/":
                sys.exit("\033[31;1m"+"ERROR: BD_params:raw_data_obs_address: "+addr+" need to end with '/' !!!"+"\033[0m")
        else: sys.exit("\033[31;1m"+"ERROR: BD_params:raw_data_obs_address: "+addr+" need to start with 'obs://' !!!"+"\033[0m")
temp_config = "config/config_temp.yaml"
with open(temp_config,'w',encoding='utf-8') as yw:
    yaml.dump(config,yw)


#if config["BD_params"]["raw_data_obs_address"] is None or config["BD_params"]["raw_data_obs_address"][0] is None or config["BD_params"]["raw_data_obs_address"][0] == "":
#    print("\033[33;1m"+"WARNING: No obs address provided.","\033[0m")
#    if config["BD_params"]["samples"] == "all" :
#        sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
#else:
#    for addr in config["BD_params"]["raw_data_obs_address"]:
#        if addr[0:6] == "obs://":
#            if addr[-1] !="/":
#                sys.exit("\033[31;1m"+"ERROR: BD_params:raw_data_obs_address: "+addr+" need to end with '/' !!!"+"\033[0m")
#        else: sys.exit("\033[31;1m"+"ERROR: BD_params:raw_data_obs_address: "+addr+" need to start with 'obs://' !!!"+"\033[0m")

addresses = config["BD_params"]["raw_data_obs_address"]
vdj = config["BD_params"]["VDJ"]
## get sample names ##  
## from EITHER obs sub directory OR config.yaml input
import os,sys
import pandas as pd
if type(config["BD_params"]["samples"]) == str :
    # type1. "all"
    if config["BD_params"]["samples"] == "all" :
        samples = pd.Series()
        for addr in config["BD_params"]["raw_data_obs_address"]:
            folders = os.popen("obsutil ls -d -s %s | grep obs  | grep -v '\(md5\|bam_files\)' | sed 's#%s##'"%(addr,addr)).readlines()
            for i in folders:
                if i.strip() != "":
                    samples[i.strip().replace("/","")] = i.strip().replace("/","")
        if(config["BD_params"]["VDJ"] == "mRNA"):
            samples = samples
        else:
            samples_tmp = list(set([item.replace('_TCR', '') for item in samples]))
            samples_tmp1 = list(set([item.replace('_dan', '') for item in samples_tmp]))
            samples_tmp2 = list(set([item.replace('_BCR', '') for item in samples_tmp1]))
            samples_merge = pd.Series()
            for j in samples_tmp2:
                if j.strip() != "":
                    samples_merge[j] = j
            samples = samples_merge
    else :  sys.exit("\033[31;1m"+"ERROR: BD_params:samples can only be strings 'all' or a list of sample."+"\033[0m")
else :
    # type2. list of sample:cell_number dictionary
    if all(isinstance(i, str) for i in config["BD_params"]["samples"]):
        samples = config["BD_params"]["samples"]
    # type3. list of sample name strings
    elif all(isinstance(i, dict) for i in config["BD_params"]["samples"]):
        # FC={}
        # for i in config["BD_params"]["samples"]:
        #     FC.update(i)
        # samples_list = [str(key)+"_FC"+str(value) for key,value in FC.items()]
        samples = [str(v)+"_"+str(k) for i in config["BD_params"]["samples"] for v,k in i.items()]
        FC = {str(v)+"_"+str(k):k for i in config["BD_params"]["samples"] for v,k in i.items()}
    else: sys.exit("\033[31;1m"+"ERROR: BD_params:samples list can only be strings or <name:cell_number> dictionary."+"\033[0m")
    #samples = pd.Series(samples_list,index = samples_list) # a list of sample name

def get_FC(wildcards): 
    # return FC[widcards.sample]
    try: FC
    except NameError: return ""
    else: return "--force_cells=" + str(FC[wildcards.sample])

#def get_raw_data_dir(wildcards):
#    try: FC
#    except NameError: return "../../raw_data/"+wildcards.sample+"/"+wildcards.sample+"_fastqs"
#    else: return "../../raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"_".join(wildcards.sample.split("_")[:-1])+"_fastqs"
#
#def get_raw_data_fastqs(wildcards):
#    try: FC
#    except NameError: return "raw_data/"+wildcards.sample+"/"+wildcards.sample+"_fastqs/"+wildcards.sample+"_S1_L001_R?_001.fastq.gz"
#    else: return "raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"_".join(wildcards.sample.split("_")[:-1])+"_fastqs/"+wildcards.sample+"_S1_L001_R?_001.fastq.gz"

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
def all_input_bd(wildcards):
    wanted_input = []
    # cellranger
    if config['BD_params']['module']["BD_analysis"]:
        for (sample) in samples:
            wanted_input.extend(
                expand(
                    [
                        "result/BD_Analysis/{sample}/{sample}_Pipeline_Report.html",
                        "result/BD_Analysis/{sample}/{sample}_Metrics_Summary.csv",
                        "result/BD_Analysis/{sample}/{sample}_RSEC_MolsPerCell_MEX.zip",
                        #"result/BD_Analysis/{sample}/{sample}_Cell_Label_Filter.png",
                        #"result/BD_Analysis/{sample}/{sample}_Cell_Label_Second_Derivative_Curve.png",
                        "result/BD_Analysis/{sample}/{sample}_Bioproduct_Stats.csv",
                        "logs/{sample}.upload_bam_log.txt"
                        #"logs/{sample}.fastqc.check"
                    ],
                    sample=sample
                )
            )
    if config['BD_params']['VDJ'] != "mRNA" and config['BD_params']['VDJ'] != "dan":
        wanted_input.extend(
            expand(
                [
                    "result/BD_Analysis/{sample}/{sample}_VDJ_metrics.csv"
                    #"logs/upload_cellranger_sor"
                ],
                sample=sample
            )
        )
    if config['BD_params']['VDJ'] == "dan":
        wanted_input.extend(
            expand(
                [
                    "result/BD_Analysis/{sample}/{sample}_DBEC_MolsPerCell_MEX.zip"
                    #"logs/upload_cellranger_sor"
                ],
                sample=sample
            )
        )
    if config['BD_params']['QC_report']:
        wanted_input.extend(
            expand(
                [
                    "logs/BD_qc.log",
                    #"logs/QC_report_upload.log"
                ])
        )
    wanted_input.extend(
        [
        #    "result/cellranger_sum/metrics_summary_stat.csv",
                "result/BD_Analysis/BD_check_log.txt"
        #    "logs/upload_cellranger_sor"
        ])
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
    output:
        temp_file="raw_data.check"
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

# ========================================================================================================================
## param1
#include_introns = config["BD_params"]["include_introns"]
#cellranger_version = config['BD_params']['envmodules'].strip().replace("cellRanger/","")
#if include_introns==True :
#    if cellranger_version == "7.0.1":
#        include_introns_param = "--include-introns=true"
#    else:
#        include_introns_param = "--include-introns"
#else:
#    if cellranger_version == "7.0.1":
#        include_introns_param = "--include-introns=false"
#    else:
#        include_introns_param = ""

# param2
bao = config["report"]["Others"]
if bao == "single":
    Basic_Algo_Only = "false"
else:
    sys.exit("\033[31;1m"+"ERROR: Basic_Algo_Only (report:Others) can only be false when multiple samples are run in batches."+"\033[0m")
## param3
generate_bam = config["BD_params"]["generate_bam"]
AbSeq_fa = config["report"]["AbSeq_Reference"]
if generate_bam == False:
    generate_bam_param = "false"
else:
    generate_bam_param = "true"

if  AbSeq_fa !="None":
    AbSeq_fa_param = AbSeq_fa
else :
    AbSeq_fa_param = "None"

if config['BD_params']['VDJ'] != "dan":
    rule BD_Analysis:
        """
        run BD_Analysis for each samples 
        """
        input:
            temp_file="raw_data.check"
        output:
            web_summary="result/BD_Analysis/{sample}/{sample}_Pipeline_Report.html",
            metrics_summary="result/BD_Analysis/{sample}/{sample}_Metrics_Summary.csv",
            MolsPerCell="result/BD_Analysis/{sample}/{sample}_RSEC_MolsPerCell_MEX.zip",
            Cell_h5ad="result/BD_Analysis/{sample}/{sample}.h5ad",
            Cell_rds="result/BD_Analysis/{sample}/{sample}_Seurat.rds",
            Bioproduct_Stats="result/BD_Analysis/{sample}/{sample}_Bioproduct_Stats.csv",
            VDJ_metrics_summary="result/BD_Analysis/{sample}/{sample}_VDJ_metrics.csv",
            log_file=directory("result/BD_Analysis/{sample}/Logs/")
            #temp_file=directory("result/BD_Analysis/{sample}temp")
        benchmark:
            "benchmarks/QC/{sample}.BD_analysis.benchmark.txt"
        resources:
            qsub_mem=64,
            qsub_p=16,
            qsub_n=1
        envmodules:
            config["BD_params"]['envmodules']
        shell:
            '''
    cd result/BD_Analysis/ &&
    rm -rf {wildcards.sample} &&
    python ../../scripts/BD_get_qsub.py \
      -i ../../raw_data/ \
      -s {wildcards.sample} \
      -d {reference}/genome.tar.gz \
      -t {Basic_Algo_Only} \
      --get_bam {generate_bam_param} \
      --AbSeq {AbSeq_fa_param} \
      -v {vdj} \
      -n {resources.qsub_p} 
            '''
if config['BD_params']['VDJ'] == "dan":
    rule BD_Analysis:
        """
        run BD_Analysis for each samples 
        """
        input:
            temp_file="raw_data.check"
        output:
            web_summary="result/BD_Analysis/{sample}/{sample}_Pipeline_Report.html",
            metrics_summary="result/BD_Analysis/{sample}/{sample}_Metrics_Summary.csv",
            MolsPerCell="result/BD_Analysis/{sample}/{sample}_RSEC_MolsPerCell_MEX.zip",
            Cell_h5ad="result/BD_Analysis/{sample}/{sample}.h5ad",
            Cell_rds="result/BD_Analysis/{sample}/{sample}_Seurat.rds",
            Bioproduct_Stats="result/BD_Analysis/{sample}/{sample}_Bioproduct_Stats.csv",
            DBEC="result/BD_Analysis/{sample}/{sample}_DBEC_MolsPerCell_MEX.zip",
            log_file=directory("result/BD_Analysis/{sample}/Logs/")
            #temp_file=directory("result/BD_Analysis/{sample}temp")
        benchmark:
            "benchmarks/QC/{sample}.BD_analysis.benchmark.txt"
        resources:
            qsub_mem=64,
            qsub_p=16,
            qsub_n=1
        envmodules:
            config["BD_params"]['envmodules']
        shell:
            '''
    cd result/BD_Analysis/ &&
    rm -rf {wildcards.sample} &&
    python ../../scripts/BD_get_qsub.py \
      -i ../../raw_data/ \
      -s {wildcards.sample} \
      -d {reference}/genome.tar.gz \
      -t {Basic_Algo_Only} \
      --get_bam {generate_bam_param} \
      --AbSeq {AbSeq_fa_param} \
      -v {vdj} \
      -n {resources.qsub_p} 
            '''
#======================================
#localrules: BD_sample_stat
if config['BD_params']['VDJ'] != "dan":
    rule BD_sample_stat:
        input:
            web_summary="result/BD_Analysis/{sample}/{sample}_Pipeline_Report.html",
            metrics_summary="result/BD_Analysis/{sample}/{sample}_Metrics_Summary.csv",
            vdj_summary="result/BD_Analysis/{sample}/{sample}_VDJ_metrics.csv"
        output:
            sample_summary=temp("result/BD_Analysis/{sample}_summary.gz"),
            temp_metrics=temp("result/BD_Analysis/{sample}_temp.csv"),
            sample_metrics=temp("result/BD_Analysis/{sample}_metrics_summary.csv")
        benchmark:
            "benchmarks/QC_sum/{sample}.BD_sample_stat.benchmark.txt"
        log:
           "logs/{sample}.BD_sample_stat.log"
        envmodules:
             config["report_params"]['envmodules']
        resources:
            qsub_mem=4,
            qsub_p=2,
            qsub_n=1
        shell:
            '''
    python scripts/BD_stat.py -i {input.metrics_summary} -o {output.temp_metrics}
    echo -e "sampleid\\n{wildcards.sample}" | paste - {output.temp_metrics} | sed 's/\\t/,/g' >  {output.sample_metrics}
    gzip -c  {input.web_summary} > {output.sample_summary} 
            '''
else :
    rule BD_sample_stat:
        input:
            web_summary="result/BD_Analysis/{sample}/{sample}_Pipeline_Report.html",
            metrics_summary="result/BD_Analysis/{sample}/{sample}_Metrics_Summary.csv"
        output:
            sample_summary=temp("result/BD_Analysis/{sample}_summary.gz"),
            temp_metrics=temp("result/BD_Analysis/{sample}_temp.csv"),
            sample_metrics=temp("result/BD_Analysis/{sample}_metrics_summary.csv")
        benchmark:
            "benchmarks/QC_sum/{sample}.BD_sample_stat.benchmark.txt"
        log:
           "logs/{sample}.BD_sample_stat.log"
        envmodules:
             config["report_params"]['envmodules']
        resources:
            qsub_mem=4,
            qsub_p=2,
            qsub_n=1
        shell:
            '''
    python scripts/BD_stat.py -i {input.metrics_summary} -o {output.temp_metrics} -t "dan"
    echo -e "sampleid\\n{wildcards.sample}\\n{wildcards.sample}_dan" | paste - {output.temp_metrics} | sed 's/\\t/,/g' >  {output.sample_metrics}
    gzip -c  {input.web_summary} > {output.sample_summary} 
            '''

# ========================================================================================================================
localrules: BD_summary
rule BD_summary:
    """
    BD stats
    """
    input:
        metrics_summary=expand("result/BD_Analysis/{sample}_metrics_summary.csv",sample=samples)
    output:
        metrics_summary_stat="result/BD_Analysis/metrics_summary_stat.csv",
    benchmark:
        "benchmarks/QC_sum/BD_summary.benchmark.txt"
    log:
        "logs/BD_summary.log"
    resources:
        qsub_mem=4,
        qsub_p=2,
        qsub_n=1
    shell:
        """
cat {input.metrics_summary} | awk '!a[$0]++' > {output.metrics_summary_stat} |& tee {log}
        """
#====================================
#localrules: BD_check
rule BD_check:
    input:
        metrics_summary_stat = "result/BD_Analysis/metrics_summary_stat.csv",
        samples_summary=expand("result/BD_Analysis/{sample}_summary.gz",sample=samples)
    output:
        check = "result/BD_Analysis/BD_check_log.txt"
    #log:
    #    "logs/cellranger/BD_check.log"
    resources:
        qsub_mem=10,
        qsub_p=2,
        qsub_n=1
    run:
        samples_html = ",".join(input.samples_summary)
        #samples_png = ",".join(input.sample_png)
        shell("python scripts/BD_check.py --config config/config.yaml --stat {input.metrics_summary_stat} --html {samples_html} --output {output.check}")

# ========================================================================================================================
if config["BD_params"]["generate_bam"]:
    rule BD_upload_bam:
        input:
            merge_summary="result/BD_Analysis/metrics_summary_stat.csv"
        output:
            check = "logs/{sample}.upload_bam_log.txt"
        params:
            task_num = config['report']['Task_Num'],
            result_dir="result/BD_Analysis"
        resources:
            qsub_mem=10,
            qsub_p=2,
            qsub_n=1
        shell:
            '''
    module purge
    bash scripts/auto_upload_bam.sh {params.result_dir}/{wildcards.sample} {params.task_num} |& tee {output.check}
            '''
else :
    rule BD_noupload_bam:
        input:
            merge_summary="result/BD_Analysis/metrics_summary_stat.csv"
        output:
            check = "logs/{sample}.upload_bam_log.txt"
        params:
            task_num = config['report']['Task_Num'],
            result_dir="result/BD_Analysis"
        resources:
            qsub_mem=10,
            qsub_p=2,
            qsub_n=1
        shell:
            '''
    module purge
    echo "Do not need to get bam" > {output.check}
            '''

#=====================================
if config["BD_params"]["QC_report"]:
    rule QC_report:
        """
        QC_report 
        """
        input:
            check = "result/BD_Analysis/BD_check_log.txt",
            bam_log=expand("logs/{sample}.upload_bam_log.txt",sample=samples)
        output:
            html = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/质控说明.html"
        params:
            report = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/",
            config = "config/config.yaml"
        log:
           "logs/BD_qc.log"
        benchmark:
            "benchmarks/BD_qc.report.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=2
        shell:
            '''
rm -rf {params.report}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report_BD/BD_QC_report.py \
  -i result -v {vdj} \
  -c {params.config} |& tee {log};
            '''