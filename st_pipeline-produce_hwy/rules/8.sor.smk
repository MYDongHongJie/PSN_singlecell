rule QC_sor_upload:
    """
    report QC spaceranger summary
    """
    input:
        all_input=all_input,
        metadata = "config/samples.csv",
        config = "config/config.yaml"
    output:
        "logs/sor/QC_upload_sor.check.done"
    benchmark:
        "benchmarks/QC_sor_upload.benchmark.txt"
    log:
        "logs/sor/run.snakemake.QC_sor_upload.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        python scripts/oest_sor_upload.py           \
                -i config/config.yaml \
                -m {input.metadata}   \
                -l logs/QC_upload_sor  >&2 2>{log} && \
        touch logs/sor/QC_upload_sor.check.done
        """
