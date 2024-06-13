rule BD_vdjplot:
    input:
        BD_check="result/BD_Analysis/BD_check_log.txt",
        metadata = config['report_params']['samples_file']
    output:
        "result/VDJ_aggr/TCR/Clonotypes/clonotypes_abundances.png",
        "result/VDJ_aggr/BCR/Clonotypes/clonotypes_abundances.png",
        temp_R1=temp("fancy_spectratype.r"),
        temp_R2=temp("join_venn.r"),
        temp_R3=temp("quantile_stats.r"),
        temp_R4=temp("rarefaction_curve.r"),
        temp_R5=temp("vexpr_plot.r"),
        temp_R6=temp("vj_pairing_plot.r")
    log:
        "logs/BD_vdj_plot.log"
    resources:
        qsub_mem=20,
        qsub_p=2,
        qsub_n=1
    shell:
        '''
    python scripts/vdj/BD_VDJ.analysis.py -v result/BD_Analysis/ -o result/VDJ_aggr/ -q Sampleid -t TCR,BCR -m {input.metadata} -p {resources.qsub_p} |& tee {log}
        '''
