##==================================================Quality Control =======================================================
rule qc:
    '''
    use sctool to quality control
    '''
    input:
        inputfile="result/seurat_ob.h5seurat"
    output:
        "result/2.Count_QC/seurat_ob.h5seurat",
        "result/2.Count_QC/statitics_before_after_QC.xls"
    benchmark:
        "benchmarks/qc/qc.benchmark.txt"
    log:
        "logs/qc/qc.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        vars2regress=config['params']['QC']["vars2regress"],
        filters=get_qc_filter_parameters,
        features2filter=config['params']['QC']["features2filter"],
        mincell4gene=config['params']['QC']["mincell4gene"],
        normmeth=config['params']['QC']["normmeth"],
        qsub_mem=config['qc']['mem'],
        qsub_p=config['qc']['cpu']
    shell:
        '''
        sctool --input result/seurat_ob.h5seurat \
              --informat h5seurat \
              --output result/2.Count_QC \
              --outformat h5seurat \
              --assay RNA \
              --dataslot counts \
              --prefix  seurat_ob \
              --ncores {params.qsub_p} \
              --update FALSE \
              qc --vars2regress "{params.vars2regress}" \
                --filters {params.filters[0]} \
                --lower {params.filters[1]} \
                --upper {params.filters[2]} \
                --cut.1 median \
                --cut.2 mad \
                --nfold 2 \
                --features2filter {params.features2filter} \
                --mincell4gene  {params.mincell4gene} \
                --ident NULL \
                --normmeth "{params.normmeth}" \
                --nvfeatures 2000 \
                --pointsize 0.1
        '''
