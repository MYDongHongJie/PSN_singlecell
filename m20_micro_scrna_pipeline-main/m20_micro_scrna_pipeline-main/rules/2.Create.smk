##==================================================create seurat object ===============================================
rule create:
    '''
    run sctools to create multimodal seurat file
    '''
    input:
        metadatafile="config/samples.csv",
        cellranger_out=expand(
            [
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz"
            ],sample=samples.index)
    output:
        "result/seurat_ob.h5seurat"
    benchmark:
        "benchmarks/create/create.benchmark.txt"
    log:
        "logs/create/create.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        outformat="h5seurat",
        metrics=config['params']['Create']['metrics'],
        #reference=config['database']['ref_genome'],
        qsub_mem=config['create']['mem'],
        qsub_p=config['create']['cpu'],
        gset=get_gset_param
    shell:
        '''
        sctool --input  result/1.STARsolo/ \
               --output result/ \
               --outformat {params.outformat} \
               --assay RNA \
               --prefix  seurat_ob \
               --ncores {params.qsub_p} \
               create   --source  mtx \
                        --metadata {input.metadatafile} \
                        --metrics {params.metrics} \
                        {params.gset}  \
                        --gcolumn 2   >&2 2>{log} 
        '''
