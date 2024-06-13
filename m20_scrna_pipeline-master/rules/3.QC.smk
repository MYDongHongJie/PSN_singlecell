##==================================================Quality Control =======================================================
rule qc:
    '''
    use sctool to quality control
    '''
    input:
        inputfile = "result/seurat_ob.h5seurat",
        donefile = "logs/Create.done"
    output:
        seurat = "result/2.Count_QC/seurat_ob.h5seurat",
        table = "result/2.Count_QC/statitics_before_after_QC.xls",
        donefile = "logs/QC.done"
    benchmark:
        "benchmarks/qc/qc.benchmark.txt"
    log:
        "logs/qc/qc.log"
    resources:
        qsub_mem = config['qc']['mem']
    threads:
        config['qc']['cpu']
    params:
        outdir = "result/2.Count_QC",
        filters = get_qc_filter_parameters,
        features2filter = get_features2filter,
        mincell4gene = mincell4gene,
        normmeth = normmeth,
        rmdoublets = rmdoublets,
        qc_method = qc_method
    shell:
        '''
        {mp}
        {m20sc_sctool}
        sctool --input {input.inputfile} \
              --informat h5seurat \
              --output {params.outdir} \
              --outformat h5seurat \
              --assay RNA \
              --dataslot counts \
              --prefix  seurat_ob \
              --ncores {threads} \
              --update FALSE \
              qc --vars2regress nCount_RNA \
                --filters {params.filters[0]} \
                --lower {params.filters[1]} \
                --upper {params.filters[2]} \
                --cut.1 median \
                --cut.2 mad \
                --nfold 2 \
                {params.features2filter} \
                --rmdoublets  {params.rmdoublets}\
                --method  {params.qc_method}\
                --mincell4gene  {params.mincell4gene} \
                --ident NULL \
                --normmeth {params.normmeth} \
                --nvfeatures 2000 \
                --pointsize 0.1 >&2 2>{log} && touch {output.donefile}
        '''

