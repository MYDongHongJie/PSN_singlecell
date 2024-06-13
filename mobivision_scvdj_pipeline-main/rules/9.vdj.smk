import glob
import os
VDJ_Types=[]
if  glob.glob(f"result/mobivision/*TCR/")  :
    VDJ_Types.append("TCR")
if  glob.glob(f"result/mobivision/*BCR/")  :
    VDJ_Types.append("BCR")
wildcard_constraints:
    vdj_type="|".join(VDJ_Types)

rule VDJ_plot:
    """
    VDJ result 
    """
    input:
        mobivision = "result/mobivision",
        metadata = config['report_params']['samples_file']
    output:
        clonotype="result/VDJ_aggr/{vdj_type}/Clonotypes/merged_vdj_contig_annotation.xls",
    params:
        outdir="result/VDJ_aggr/"
    log:
        "logs/{vdj_type}_vdj_plot.log"
    resources:
        qsub_mem=20,
        qsub_p=2
    shell:
        """
python scripts/VDJ_scripts/VDJ.py -v result/mobivision/ -o result/VDJ_aggr/ -q Sampleid -t {wildcards.vdj_type} -m {input.metadata} -p {resources.qsub_p} |& tee {log}
        """

rule VDJ_scRNA_Conjoint :
    input:
        filtered_rds = "result/Count_QC/filtered.h5seurat",
        clonotypes= "result/VDJ_aggr/{vdj_type}/Clonotypes/merged_vdj_contig_annotation.xls"
    output:
        "logs/{vdj_type}_vdj_scRNA_conjoint.log",
        "result/VDJ_aggr/{vdj_type}/vdj_scRNA_Conjoint_analysis/clonotype.metadata.xls"
    params:
        reduct = config['report_params']['reduct2'],
        topn = config['report_params']['topn']
    benchmark:
        "benchmarks/{vdj_type}_vdj_scRNA_conjoint.benchmark.txt"
    log:
        "logs/{vdj_type}_vdj_scRNA_conjoint.log"
    resources:
        qsub_mem=20,
        qsub_p=2
    envmodules:
        config['report_params']['envmodules'] 
    shell:
        """
    Rscript scripts/VDJ_scripts/vdj.R -i {input.filtered_rds} -c {input.clonotypes} -r {params.reduct} --topn {params.topn} --type {wildcards.vdj_type} -g clusters -y group,sampleid -a all -n 3 -o result/VDJ_aggr/{wildcards.vdj_type}/vdj_scRNA_Conjoint_analysis |& tee {log}
        """
