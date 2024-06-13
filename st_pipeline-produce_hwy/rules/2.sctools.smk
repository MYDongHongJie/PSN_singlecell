rule sctools_create:
    """
    run sctools to create merge rds 
    """
    input:
        metadata = config['samples'],
        spaceranger_outfiles = expand(["result/spaceranger/{sample}/outs/web_summary.html",
                                    "result/spaceranger/{sample}/outs/metrics_summary.csv",
                                    "result/spaceranger/{sample}/outs/cloupe.cloupe",
                                    "result/spaceranger/{sample}/outs/spatial/tissue_lowres_image.png",
                                    "result/spaceranger/{sample}/outs/raw_feature_bc_matrix.h5",
                                    "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5"],sample=samples.index),
        #mt_gene = os.path.join(config['database']['reference'],"MT_genelist.gmt")
    output:
            f"result/sctools/spatial.rds"
    params:
            spaceranger_outdir = "result/spaceranger",
            outdir="result/sctools/",
            formart=config["params"]["object_formart"],
            qsub_mem=config['sctools_create']['mem'],
            qsub_p=config['sctools_create']['cpu'],
            mt = "" if not os.path.isfile(os.path.join(config['database']['reference'],'MT_genelist.gmt')) or not config["params"]["library_type"] == "fresh"  else f" --gset {os.path.join(config['database']['reference'],'MT_genelist.gmt')}"
    benchmark:
             "benchmarks/sctools_create.benchmark.txt"
    log:
       "logs/create/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        sctool  -i {params.spaceranger_outdir}  \
                -f {params.formart} \
                -o {params.outdir} \
                -d {params.formart} \
                -j {params.qsub_p} \
                --prefix spatial \
                --assay  Spatial \
                create --cell.meta FALSE  \
                       -s h5 \
                       -m  {input.metadata} \
                       {params.mt}  >&2 2>{log}                          
        """
