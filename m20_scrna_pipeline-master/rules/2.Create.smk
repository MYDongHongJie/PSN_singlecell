##==================================================create seurat object ===============================================
rule create:
    '''
    run sctools to create multimodal seurat file
    '''
    input:
        samplefile = samplefile,
        cellranger_out = expand(
            [
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
            ],
                sample = samples.index
            ),
        starsolo_file= expand(
            [
                "logs/starsolo_{sample}.done"
            ],sample=samples.index) if filter_standard == "" else expand(
            [
                "logs/re_filter_{sample}.done"
            ],sample=samples.index)
    output:
        seurat = "result/seurat_ob.h5seurat",
        donefile = "logs/Create.done"
    benchmark:
        "benchmarks/create/create.benchmark.txt"
    log:
        "logs/create/create.log"
    resources:
        qsub_mem = config['create']['mem']
    params:
        indir = "result/1.STARsolo/",
        outdir = "result/",
        outformat = "h5seurat",
        metrics = get_metrics_param,
        reference = f"{ref_dir}",
        gset = get_gset_param
    threads:
        config['create']['cpu']
    shell:
        '''
        {mp}
        {m20sc_sctool}
        sctool --input {params.indir}  \
            --output  {params.outdir}\
            --outformat {params.outformat} \
            --assay RNA \
            --prefix  seurat_ob \
            --ncores {threads} \
            create   --source  mtx \
                        --metadata {input.samplefile} \
                        {params.metrics} \
                        {params.gset}  \
                        --gcolumn 2  >&2 2>{log} && touch {output.donefile}
        '''

        