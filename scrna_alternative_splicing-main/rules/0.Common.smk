
import pandas as pd
import numpy as np
import os

### load config file======================================================================================================
configfile: "config/config.yaml"
configfile: "config/cluster.yaml"

### tool and scripts =====================================================================================================
## tool 
mp = config["envmodules"]["mp"]
cellsnplite = config["envmodules"]["cellsnplite"]
samtools = config["envmodules"]["samtools"]
Rscript_splicing = config["envmodules"]["Rscript_splicing"]
star = config["envmodules"]["star"]
cellranger = config["envmodules"]["cellranger"]
enrichwrap = config["envmodules"]["enrichwrap"]
Rscript_splicing = config["envmodules"]["Rscript_splicing"]
sshpass = config["envmodules"]["sshpass"]
subsetBam = config["envmodules"]["subsetBam"]

## scripts
marvel_prepare = config["scripts"]["marvel_prepare"]
marvel_subsetBam = config["scripts"]["marvel_subsetBam"]
marvel_analysis = config["scripts"]["marvel_analysis"]
marvel_plot = config["scripts"]["marvel_plot"]
report_script = config["scripts"]["report_script"]
sashimi = config["scripts"]["sashimi"]

## bam path
# Raw_Task_Num = config["report"]["Raw_Task_Num"].split(',')
bam_path = config['bam_path']
password = config["report"]["password"]

## report
report_dir = str(config['report']['Project_Num']) + "_Report"
## library
if os.path.exists("result/cellranger"):
    if config["library"]["version_10X_chemistry"] == "V2":
        soloUMIlen = 10
    if config["library"]["version_10X_chemistry"] == "V3":
        soloUMIlen = 12
if os.path.exists("result/1.STARsolo"):
    soloUMIlen = 12


## reference
Species = config["report"]["Species"].capitalize()
if Species == "Human":
    if os.path.exists("result/cellranger"):
        genomeDir = config["database"]["Homo"]["10x"]["genomeDir_samewithm20"]
        gtf = config['database']['Homo']["m20"]['ref_genome'] + "/genome.gtf" 

    if os.path.exists("result/1.STARsolo"):
        genomeDir = config["database"]["Homo"]["10x"]["genomeDir_samewithm20"]
        gtf = config['database']['Homo']["m20"]['ref_genome'] + "/genome.gtf" 
        # 后面两个参数是为call snp 准备的
        ref_vcf = config["database"]["Homo"]["m20"]["ref_vcf"]
        ref_dir = config['database']['Homo']["m20"]['ref_genome']

if Species == "Mouse":
    if os.path.exists("result/cellranger"):
        genomeDir = config["database"]["Mus"]["10x"]["genomeDir_samewithm20"]
        gtf = config['database']['Mus']["m20"]['ref_genome'] + "/genome.gtf" 
    if os.path.exists("result/1.STARsolo"):
        genomeDir = config["database"]["Mus"]["m20"]["genomeDir_samewithm20"]
        gtf = config['database']['Mus']["m20"]['ref_genome'] + "/genome.gtf" 


## get sampleid=========================================================================================================
# 样本文件可以使用 “#” 进行注释
samplefile = config["metadatafile"]
difffile = config["diff_DEG_group"]
changenamefile = config["change_name"]

if os.path.exists(samplefile):
    samples = pd.read_csv(samplefile, dtype=str)
    samples.dropna(how='all', inplace=True) # 删除可能存在的空行
    samples = samples[ ~np.array([s.startswith("#") for s in samples.sampleid.to_list()])]
    if os.path.exists(changenamefile):
        change_name = pd.read_csv(changenamefile, dtype=str, sep="\t", header=None, names=["old_name", "sampleid"])
        change_name.dropna(how='all', inplace=True)
        change_name.old_name = "S" + change_name.old_name
        samples = pd.merge(samples, change_name, on='sampleid' )
    samples = samples.set_index("sampleid",drop=False)
    samples.index.names = ["sample_id"]
    wildcard_constraints:
        sample="|".join(samples.index)
    # print(samples)
else:
    print(f"您尚未准备好样本文件{samplefile}")
    os._exit(0)






##read diff_group file =================================================================================================
if config["module"]["splicing"]:
    ###########
    # 差异分析文件可以用 “#” 进行注释
    if os.path.exists(difffile):
        diff_groups = pd.read_csv(difffile, dtype=str, delimiter=",")
        diff_groups.dropna(how='all', inplace=True) # 删除可能存在的空行
        diff_groups = diff_groups[ ~np.array([s.startswith("#") for s in diff_groups.treatment.to_list()]) ]
        diff_groups.apply(lambda x: x.str.strip() if x.dtype == "object" else x) # 删除每个字符串两头的空格
        diff_groups.index = diff_groups["type"] + "_" + diff_groups["treatment_name"].replace(" ", "_", regex=True) + "-vs-" + diff_groups["control_name"].replace(" ", "_", regex=True)
        diff_groups['groupname'] = diff_groups.index
        wildcard_constraints:
            diff_group="|".join(diff_groups.index)
    else:
        print(f"您尚未准备好差异分组设置文件{difffile}")
        os._exit(0)

    betweencells = [x for x in diff_groups.index if not x.startswith('group') and not x.startswith('sampleid')]
    betweensamples = [x for x in diff_groups.index if x.startswith('sampleid')]
    ############
    # prepare
    ############

    dimension_type = config["params"]["splicing"]["prepare"]["dimension_type"]
    use_routine_genearray = config["params"]["splicing"]["prepare"]["use_routine_genearray"]
    if os.path.exists("result/cellranger"):
        # 10X 数据
        dimension_type_file = "result/Clustering/" + f"{dimension_type}" + "_Dimension_Reduction/" + f"{dimension_type}" + "_Dimension_Reduction_coordination.csv"
        dimension_type_cluster = "result/Clustering/clustering_results.csv"
        # if use_routine_genearray:
        #     Gene_datapath = "result/cellranger"
        # else:
        #     Gene_datapath = "result/Splicing/STARSolo/"
    if os.path.exists("result/1.STARsolo"):
        # M20 数据
        dimension_type_file = "result/3.Clustering/" + f"{dimension_type}" + "_Dimension_Reduction/" + f"{dimension_type}" + "_Dimension_Reduction_coordination.csv"
        dimension_type_cluster = "result/3.Clustering/" + f"{dimension_type}" + "_Dimension_Reduction/clustering_results.csv"
        # if use_routine_genearray:
        #     Gene_datapath = "result/1.STARsolo/"
        # else:
        #     Gene_datapath = "result/Splicing/STARSolo/"

    ############
    # analysis
    ############

    # 按照在细胞中表达百分比筛选基因，若 filter_hard 设置为 TRUE 直接定死，若为FALSE则使用marvel自带的方式进行筛选
    filter_hard = config["params"]["splicing"]["analysis"]["filter_hard"]
    # if filter_hard:
    #     filter_gene_in_percent_cell = "--filter_gene_in_percent_cell " + str(config["params"]["splicing"]["analysis"]["filter_gene_in_percent_cell"])
    #     filter_sj_in_percent_cell = "--filter_sj_in_percent_cell " + str(config["params"]["splicing"]["analysis"]["filter_sj_in_percent_cell"])
    #     downsample = False
    # else:
    #     filter_gene_in_percent_cell = ""
    #     filter_sj_in_percent_cell = ""
    #     downsample = False

    ############
    min_gene_expression = config["params"]["splicing"]["analysis"]["min_gene_expression"]
    min_gene_expression_name = "min.log2EP-" + str(min_gene_expression)
    # iterations = config["params"]["splicing"]["analysis"]["iterations"]

    ############
    # plot
    ############

    pvalue_sj = config["params"]["splicing"]["plot"]["pvalue_sj"]
    if "delta_sj" in config["params"]["splicing"]["plot"].keys() & "log2fc_sj" in config["params"]["splicing"].keys():
        print("Warning: 参数 delta_sj 和 log2fc_sj 二选一")
        os._exit(0)
    if "delta_sj" in config["params"]["splicing"]["plot"].keys():
        command_delta_or_log2fc = "--delta_sj " + str(config["params"]["splicing"]["plot"]["delta_sj"])
        delat_or_log2fc = "delta.sj-" + str(config["params"]["splicing"]["plot"]["delta_sj"])
    if "log2fc_sj" in config["params"]["splicing"]["plot"].keys():
        command_delta_or_log2fc = "--log2fc_sj " + str(config["params"]["splicing"]["plot"]["log2fc_sj"])
        delat_or_log2fc = "log2fc.sj-" + str(config["params"]["splicing"]["plot"]["log2fc_sj"])
    ############
    if "pvalue_gene" in config["params"]["splicing"]["plot"].keys() & "qvalue_gene" in config["params"]["splicing"]["plot"].keys():
        print("Warning: 参数 pvalue_gene 和 qvalue_gene 二选一")
        os._exit(0)
    if "pvalue_gene" in config["params"]["splicing"]["plot"].keys():
        pvalue_gene = config["params"]["splicing"]["plot"]["pvalue_gene"]
        # command_pval_or_padj = "--pvalue_gene " + str(pvalue_gene)
        pval_or_padj = "pvalue.gene-" + str(pvalue_gene)
    if "qvalue_gene" in config["params"]["splicing"]["plot"].keys():
        qvalue_gene = config["params"]["splicing"]["plot"]["qvalue_gene"]
        # command_pval_or_padj = "--qvalue_gene" + str(qvalue_gene)
        pval_or_padj = "qvalue.gene-" + str(qvalue_gene)
    log2fc_gene = config["params"]["splicing"]["plot"]["log2fc_gene"]
    volcan_long_name = "pvalue.sj-" + str(pvalue_sj) + "-" + delat_or_log2fc + "-" + min_gene_expression_name
    diff_long_name = "pvalue.sj-" + str(pvalue_sj) + "-" + delat_or_log2fc + "-" + pval_or_padj + "-" + "log2fc.gene-" + str(log2fc_gene) + "-" + min_gene_expression_name
    ############
    diff_exp_junction_only = config["params"]["splicing"]["plot"]["diff_exp_junction_only"]

    ############
    # top_DE_junctions = config["params"]["splicing"]["plot"]["top_DE_junctions"]

    report_or_not = config["params"]["splicing"]["report"]
    

## function =============================================================================================================

flatten = lambda x: [subitem for item in x for subitem in flatten(item)] if type(x) is list else [x]  




