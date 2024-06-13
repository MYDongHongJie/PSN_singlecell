rule create:
    """
    create xenium object
    """
    input:
        summary_file = "result/1.Xenium/summary_all.csv"
    output:
        "result/seurat_ob.rds"
    benchmark:
        "benchmarks/create/create.txt"
    log:
        "logs/create/create.log"
    params:
        medatafile = config['database']['metadataFile']
    shell:
        '''
        sctool --input result/1.Xenium \
                --output result/ \
                --outformat rds \
                --assay xenium \
                --prefix seurat_ob \
                create --source xenium \
                       --metadata {params.medatafile} >&2 2>{log}
        '''
