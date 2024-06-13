import glob
import os
VDJ_Types=[]
if  glob.glob(f"result/cellranger/*/outs/multi/vdj_t")  :
    VDJ_Types.append("TCR")
if  glob.glob(f"result/cellranger/*/outs/multi/vdj_b")  :
    VDJ_Types.append("BCR")
wildcard_constraints:
    vdj_type="|".join(VDJ_Types)

Project_Num=config['report']['Task_Num'].split("-")[0]

localrules: VDJ_report
rule VDJ_report:
    """
    report summary
    """
    input:
        config = "config/config.yaml",
        vdj_scRNA_log="logs/{vdj_type}_vdj_scRNA_conjoint.log" if config['report_params']['module']['VDJ_scRNA_Conjoint'] else "logs/{vdj_type}_vdj_plot.log"
    output:
        html = "result/VDJ_aggr/{vdj_type}/%s_{vdj_type}_Report_%s/分析说明.html" % (Project_Num,time.strftime('%Y_%m_%d')),
        #ll = "logs/{vdj_type}_report.log"
    params:
        report = "result/VDJ_aggr/{vdj_type}/%s_{vdj_type}_Report_%s/" % (Project_Num,time.strftime('%Y_%m_%d')),
        config = "config/config.yaml",
    benchmark:
        "benchmarks/{vdj_type}_VDJ_report.benchmark.txt"
    log:  "logs/{vdj_type}_report.log"
    resources:
        qsub_mem=10,
        qsub_p=1,
    envmodules:
        config['report_params']['envmodules']
        #"oebio/1.2.15"
    shell:
        """
rm -rf {params.report}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/VDJ_scripts/VDJ_report.py \
  -i  result/VDJ_aggr \
  -c {params.config} \
  -t {wildcards.vdj_type} |& tee {log};
        """
