###########################################################
## ==================== 1. load config ====================
###########################################################
## checking raw_data address
if config["mobivision_params"]["raw_data_obs_address"] is None or config["mobivision_params"]["raw_data_obs_address"][0] is None or config["mobivision_params"]["raw_data_obs_address"][0] == "":
    print("\033[33;1m"+"WARNING: No obs address provided.","\033[0m")
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
        config["mobivision_params"]["raw_data_obs_address"][0] = obs_address
        print("\033[33;1m"+"Info: Auto create obs address: "+ obs_address + "\033[0m")
        update_config = "config/config_update.yaml"
        with open(update_config,'w',encoding='utf-8') as yw:
            yaml.dump(config,yw,allow_unicode=True)



#    if config["mobivision_params"]["samples"] == "all" :
#        sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
else:
    for addr in config["mobivision_params"]["raw_data_obs_address"]:
        if addr[0:6] == "obs://":
            if addr[-1] !="/":
                sys.exit("\033[31;1m"+"ERROR: mobivision_params:raw_data_obs_address: "+addr+" need to end with '/' !!!"+"\033[0m")
        else: sys.exit("\033[31;1m"+"ERROR: mobivision_params:raw_data_obs_address: "+addr+" need to start with 'obs://' !!!"+"\033[0m")

## get sample names ##  
## from EITHER obs sub directory OR config.yaml input
import os,sys
import pandas as pd
if type(config["mobivision_params"]["samples"]) == str :
    # type1. "all"
    if config["mobivision_params"]["samples"] == "all" :
        samples = pd.Series()
        for addr in config["mobivision_params"]["raw_data_obs_address"]:
            folders = os.popen("obsutil ls -d -s %s | grep obs  | grep -v '\(md5\|bam_files\)' | sed 's#%s##'"%(addr,addr)).readlines()
            for i in folders:
                if i.strip() != "":
                    samples[i.strip().replace("/","")] = i.strip().replace("/","")
    else :  sys.exit("\033[31;1m"+"ERROR: mobivision_params:samples can only be strings 'all' or a list of sample."+"\033[0m")
else :
    # type2. list of sample name strings
    if all(isinstance(i, str) for i in config["mobivision_params"]["samples"]):
        samples = config["mobivision_params"]["samples"]
    # type3. list of sample:cell_number dictionary    
    elif all(isinstance(i, dict) for i in config["mobivision_params"]["samples"]):
        # FC={}
        # for i in config["mobivision_params"]["samples"]:
        #     FC.update(i)
        # samples_list = [str(key)+"_FC"+str(value) for key,value in FC.items()]
        samples = [str(k)+"_"+str(v) for i in config["mobivision_params"]["samples"] for k,v in i.items()]
        FC = {str(k)+"_"+str(v):v for i in config["mobivision_params"]["samples"] for k,v in i.items()}
    else: sys.exit("\033[31;1m"+"ERROR: mobivision_params:samples list can only be strings or <name:cell_number> dictionary."+"\033[0m")
    #samples = pd.Series(samples_list,index = samples_list) # a list of sample name

def get_FC(wildcards): 
    # return FC[widcards.sample]
    try: FC
    except NameError: return ""
    else: return "--cellnumber " + str(FC[wildcards.sample])

def get_raw_data_dir(wildcards):
    try: FC
    except NameError: return "../../raw_data/"+wildcards.sample
    else: return "../../raw_data/"+"_".join(wildcards.sample.split("_")[:-1])

def get_raw_data_fastqs(wildcards):
    try: FC
    except NameError: return "raw_data/"+wildcards.sample+"/"+"*.gz"
    else: return "raw_data/"+"_".join(wildcards.sample.split("_")[:-1])+"/"+"*.gz"

wildcard_constraints:
    sample="|".join(samples)

################################################################
## ==================== 2. define all input ====================
################################################################
def all_input_mobivision(wildcards):
    wanted_input = []
    # mobivision
    if config['mobivision_params']['module']["mobivision"]:
        for (sample) in samples:
            wanted_input.extend(
                expand(
                    [
                        "result/mobivision/{sample}/outs/{sample}_Report.html",
                        "result/mobivision/{sample}/outs/{sample}_summary.csv",
                        "result/mobivision/{sample}/outs/{sample}_filtered.h5ad",
                        "result/mobivision_sum/{sample}_Report.html",
                        "logs/{sample}.fastqc.check"
                    ],
                    sample=sample
                )
            )
            wanted_input.extend(
                [
                #    "result/mobivision_sum/metrics_summary_stat.csv",
                    "result/mobivision_sum/mobivision_check_log.txt",
                    "logs/upload_mobivision_sor"
                ])
    if config['mobivision_params']['module']["upload_bam"]:
        wanted_input.extend(
            [
                "logs/upload_bam.log"
            ])
    if config['mobivision_params']['module']["mobivision_report"]:
        wanted_input.extend(
            [
                "logs/mobivision_report.log",
                "logs/report_upload.log"
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
        config = "config/config_update.yaml"
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
###################################################################
#param
if config["mobivision_params"]["merge_samples"]:
    merge_samples = "yes"
else:
    merge_samples = "no"
rule data_merge:
    """
    merge fastqs
    """
    input:
        temp_file="raw_data.check"
    output:
        tmpfiles=temp("merge_data.check")
    benchmark:
        "benchmarks/merge.benchmark.txt"
    log:
        "logs/data_merge.log" 
    resources:
        qsub_mem=120,
        qsub_p=config['cpu']['merge'],
        qsub_n=1
    shell:
        '''
if [ {merge_samples} == "yes" ];then bash scripts/merge.sh;fi && touch merge_data.check
        '''

# ========================================================================================================================
# param1
include_introns = config["mobivision_params"]["intron"]
if include_introns=="included":
    include_introns_param = "--intron included"
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
if config["mobivision_params"]["merge_samples"] == True:
    merge_param = "merged"
else:
    merge_param = ""
# param4
if config["mobivision_params"]["version"] == "v1.1":
    version_param = "/data/software/MobiVision/mobivision-v1.1/source.sh"
    method_param = "--ed"
else:
    version_param = "/public/scRNA_works/works/guochy/software/mobivision-v2.0/source.sh"
    method_param = " "

#

rule mobivision:
    """
    run mobivision for each samples 
    """
    input:
        temp_file="merge_data.check"
    output:
        web_summary="result/mobivision/{sample}/outs/{sample}_Report.html",
        metrics_summary="result/mobivision/{sample}/outs/{sample}_summary.csv",
        filtered="result/mobivision/{sample}/outs/{sample}_filtered.h5ad",
    params:
        fastq = get_raw_data_dir,
        force_cells = get_FC
    benchmark:
        "benchmarks/QC/{sample}.mobivision.benchmark.txt"
    resources:
        qsub_mem=500,
        qsub_p=config['cpu']['mobivision'],
        qsub_n=1        
    shell:
        '''
cd result/mobivision/ &&
rm -rf {wildcards.sample} &&
source {version_param}  &&
mobivision quantify -i {reference}  \
  -t {resources.qsub_p}  \
  -f {params.fastq}/{merge_param}    \
  -o {wildcards.sample}  \
  -s {wildcards.sample} \
  {include_introns_param} {params.force_cells}  {method_param} &&
mobivision rename -i {wildcards.sample}  \
  -d {wildcards.sample} &&
mv {wildcards.sample}/*outs {wildcards.sample}/outs &&
rm {wildcards.sample}/outs/*bam* &&
rm {wildcards.sample}/outs/run_analysis_cmds.txt

        '''

# ========================================================================================================================
localrules: mobivision_sample_stat
rule mobivision_sample_stat:
    input:
        web_summary="result/mobivision/{sample}/outs/{sample}_Report.html",
        metrics_summary="result/mobivision/{sample}/outs/{sample}_summary.csv"
    output:
        sample_summary="result/mobivision_sum/{sample}_Report.html",
        sample_metrics="result/mobivision_sum/{sample}_summary.csv",
        sample_png="result/mobivision/{sample}.png"
    benchmark:
        "benchmarks/QC_sum/{sample}.mobivision_sample_stat.benchmark.txt"
    log:
       "logs/{sample}.mobivision_sample_stat.log"
    envmodules:
         config["report_params"]['envmodules']
    resources:
        qsub_mem=4,
        qsub_p=1,
        qsub_n=1
    shell:
        '''
python scripts/selenium_screenshot.py -i {input.web_summary} -n {wildcards.sample} -o result/mobivision/ |& tee {log}
cp {input.web_summary}  {output.sample_summary}
sed -i 's/\\t/,/g'  {input.metrics_summary}
[ -n "$ (tail -c1  {input.metrics_summary})" ] && echo >> {input.metrics_summary}
cp {input.metrics_summary}  {output.sample_metrics}
        '''

# ========================================================================================================================
localrules: mobivision_summary
rule mobivision_summary:
    """
    mobivision stats
    """
    input:
        metrics_summary=expand("result/mobivision_sum/{sample}_summary.csv",sample=samples)
    output:
        metrics_summary_stat="result/mobivision_sum/summary_stat.csv",
    benchmark:
        "benchmarks/QC_sum/mobivision_summary.benchmark.txt"
    log:
        "logs/mobivision_summary.log"
    resources:
        qsub_mem=1,
        qsub_p=1,
        qsub_n=1
    shell:
        """
cat {input.metrics_summary} | awk '!a[$0]++' > {output.metrics_summary_stat} |& tee {log}
        """

# ========================================================================================================================
localrules: mobivision_check
rule mobivision_check:
    input:
        metrics_summary_stat = "result/mobivision_sum/summary_stat.csv",
        samples_summary=expand("result/mobivision_sum/{sample}_Report.html",sample=samples),
        sample_png=expand("result/mobivision/{sample}.png",sample=samples),
        sample_csv=expand("result/mobivision_sum/{sample}_summary.csv",sample=samples)

    output:
        check = "result/mobivision_sum/mobivision_check_log.txt"
    #log:
    #    "logs/mobivision/mobivision_check.log"
    resources:
        qsub_mem=10,
        qsub_p=1,
        qsub_n=1
    run:
        samples_html = ",".join(input.samples_summary)
        samples_png = ",".join(input.sample_png)
        samples_csv = ",".join(input.sample_csv)
        shell("python scripts/mobivision_check.py --config config/config.yaml --stat {input.metrics_summary_stat} --html {samples_html} --png {samples_png} --csv {samples_csv} --output {output.check}")

# ========================================================================================================================
localrules: mobivision_sor
rule mobivision_sor:
    """
    upload mobivision SOR
    """
    input:
        metrics_summary_stat="result/mobivision_sum/summary_stat.csv",
        #metadata = "config/samples.csv",
    output:
        directory("logs/upload_mobivision_sor")
    params:
        config = "config/config_update.yaml"
    benchmark:
        "benchmarks/mobivision_sor.benchmark.txt"
    resources:
        qsub_mem=4,
        qsub_p=1
    envmodules:
         config["report_params"]['envmodules']
    shell:
        """
python  scripts/scrna_sor.py \
  -i {params.config} \
  -l logs/upload_mobivision_sor \
  -s "Project,QC,Cpu_Mem,Ref_database,Software"
        """

# ========================================================================================================================
rule fastqc_check: 
    input:
        stat =  "result/mobivision/{sample}/outs/{sample}_summary.csv",
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
if config["mobivision_params"]["module"]["mobivision_report"]:
    rule mobivision_report:
        """
        mobivision_report
        """
        input:
            web_summary=expand("result/mobivision/{sample}/outs/{sample}_Report.html",sample=samples)
        output:
            html = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/质控说明.html"
        params:
            report = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/",
            config = "config/config.yaml"
        log:
           "logs/mobivision_report.log"
        benchmark:
            "benchmarks/mobivision.report.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=1
        shell:
            '''
rm -rf {params.report}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/mobivision_report.py \
  -i result \
  -c {params.config} |& tee {log};
            '''

# ========================================================================================================================
    rule mobivision_report_upload:
        """
        mobivision report upload
        """
        input:
            report = f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}/质控说明.html"
        output:
            protected(f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}.zip"),
            protected(f"result/report/{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}.zip.md5"),
            "result/report/mobivision_report_upload.check"
        params:
            report_num = config['report']['Project_Num'],
            report = f"{config['report']['Task_Num']}_QC_Report_{time.strftime('%Y_%m_%d')}"
        log:
           "logs/report_upload.log"
        benchmark:
            "benchmarks/report_upload.benchmark.txt"
        resources:
            qsub_mem=10,
            qsub_p=1
        shell:
            '''
cd result/report/ ; 
/public/scRNA_works/works/guokaiqi/software/zip -r {params.report}.zip {params.report} ;
obsutil cp {params.report}.zip   obs://oe-scrna/Analysis_Report/{params.report_num}/ -f -vlength &&
    ( echo upload succeed on : obs://oe-scrna/Analysis_Report/{params.report_num}/{params.report}.zip  |& tee mobivision_report_upload.check ) ||
    echo upload fails |& tee {log}
            '''

# ========================================================================================================================
if config["mobivision_params"]["module"]["upload_bam"]:
    # check format
    raw_data_dir = config["mobivision_params"]["raw_data_obs_address"][0]
    if raw_data_dir[-1] =="/":
        raw_data_dir = raw_data_dir[:-1]

    rule bam_prepare:
        """
        prepare bam files for each sample
        """
        input:
            web_summary=expand("result/mobivision/{sample}/outs/{sample}_Report.html",sample=samples)
        output:
            bam_md5=expand("bam_files/{sample}/possorted_genome_bam.bam.md5",sample=samples),
            bam_bai_md5=expand("bam_files/{sample}/possorted_genome_bam.bam.bai.md5",sample=samples)
        log:
           "logs/bam_prepare.log"
        benchmark:
            "benchmarks/bam_prepare.benchmark.txt"
        resources:
            qsub_mem=20,
            qsub_p=config['cpu']['bam_prepare'],
            qsub_n=1
        shell:
            '''
if [ ! -d "bam_files" ]; then mkdir bam_files ;fi
for i in `ls -d result/mobivision/*/ | grep -v  '\(md5\|.png\|aggr\)' | sed 's#result/mobivision/##'`; \
do  mv  result/mobivision/$i/outs/$i_Aligned.sort.bam* bam_files/$i ;done
for i in `ls  -d bam_files/*/`;do cd $i ;md5sum $i_Aligned.sort.bam > $i_Aligned.sort.bam.md5 ; \
md5sum $i_Aligned.sort.bam.bai > $i_Aligned.sort.bam.bai.md5 ;cd - ; done
            '''

# ========================================================================================================================
    rule bam_upload:
        """
        upload bam files to OBS
        """
        input:
            bam_md5=expand("bam_files/{sample}/{sample}_Aligned.sort.bam.md5",sample=samples),
            bam_bai_md5=expand("bam_files/{sample}/{sample}_Aligned.sort.bam.bai.md5",sample=samples),
            fastqc_check = expand("logs/{sample}.fastqc.check",sample=samples)
        output:
            "logs/upload_bam.log"
        params:
            raw_data_dir = raw_data_dir
        benchmark:
            "benchmarks/bam_upload.benchmark.txt"
        resources:
            qsub_mem=20,
            qsub_p=config['cpu']['bam_upload'],
            qsub_n=1
        envmodules:
            "obsutil/5.1.11"
        shell:
            '''
obsutil  cp  bam_files/  {params.raw_data_dir}  -r  -f -vlength -vmd5
if [ -d "bam_files" ]; then rm -r bam_files ;fi
if [ -d "raw_data" ]; then rm -r raw_data ;fi
touch logs/upload_bam.log
            '''

