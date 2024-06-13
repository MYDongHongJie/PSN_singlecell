rule qc:
    """
    xenium quality control
    """
    input:
        inputfile = "result/seurat_ob.rds"
    output:
        "result/2.Count_QC/filtered_seurat.rds",
        "result/2.Count_QC/statitics_for_QC.xls",
        "result/2.Count_QC/QC_featureplot_for_nCount_xenium.png",
        "result/2.Count_QC/QC_featureplot_for_nFeature_xenium.png"
    benchmark:
        "benchmarks/qc/qc.txt"
    log:
        "logs/qc/qc.log"
    params:
        vars2regress = config['params']['QC']['vars2regress'],
        filters = config['params']['QC']['filters'],
        lower = config['params']['QC']['lower'],
        upper = config['params']['QC']['upper'],
        normmeth = config['params']['QC']['normmeth'],
        rmdoublets =config['params']['QC']['rmdoublets'],
        crop=config['params']['QC']['crop'],
        location=config['params']['QC']['location']

    shell:
        '''
        sctool -i {input.inputfile} \
                -f rds \
                -o result/2.Count_QC \
                -d rds \
                --assay xenium \
                --prefix filtered_seurat \
                --update FALSE \
                qc -r {params.vars2regress} \
                    -c {params.filters} \
                    -l {params.lower} \
                    -L {params.upper} \
                    --loc {params.location} \
                    --rmdoublets {params.rmdoublets} \
                    --crop {params.crop} \
                    --normmeth {params.normmeth} \
                    --pointsize 0 >&2 2>{log}
        '''
