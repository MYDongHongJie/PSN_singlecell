rule VDJ_report:
    input:
        vdj_check="logs/BD_vdj_plot.log",
        metadata = config['report_params']['samples_file']
    output:
        BCR_html = f"result/VDJ_aggr/BCR/{config['report']['Task_Num']}_BCR_Report_{time.strftime('%Y_%m_%d')}/分析说明.html",
        TCR_html = f"result/VDJ_aggr/TCR/{config['report']['Task_Num']}_TCR_Report_{time.strftime('%Y_%m_%d')}/分析说明.html"
    params:
        dir_TCR = f"result/VDJ_aggr/TCR/{config['report']['Task_Num']}_TCR_Report_{time.strftime('%Y_%m_%d')}/",
        dir_BCR = f"result/VDJ_aggr/BCR/{config['report']['Task_Num']}_BCR_Report_{time.strftime('%Y_%m_%d')}/",
        config = "config/config.yaml"
    log:
        bcr_log="logs/BCR_vdj_report.log",
        tcr_log="logs/TCR_vdj_report.log"
    resources:
        qsub_mem=20,
        qsub_p=3,
        qsub_n=1
    shell:
            '''
rm -rf {params.dir_TCR}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/vdj/BD_VDJ.report.py \
  -i result/VDJ_aggr -t TCR \
  -c {params.config} |& tee {log.tcr_log};
rm -rf {params.dir_BCR}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/vdj/BD_VDJ.report.py \
  -i result/VDJ_aggr -t BCR \
  -c {params.config} |& tee {log.bcr_log};
            '''
