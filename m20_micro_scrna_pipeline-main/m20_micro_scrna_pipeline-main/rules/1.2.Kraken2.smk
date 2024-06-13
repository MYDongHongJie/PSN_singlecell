# kraken2比对 ===========================================================================================================
rule kraken2:
    """
    run kraken2 for each sample 
    """
    input:
        CB_reads=get_fastqs_CB_reads,
        cDNA_reads=get_fastqs_cDNA_reads,
    output:
        "result/kraken_results/{sample}/{sample}_1.fq",
        "result/kraken_results/{sample}/{sample}_2.fq",
        "result/kraken_results/{sample}/{sample}.kraken.report.txt",
        "result/kraken_results/{sample}/{sample}.kraken.report.mpa.txt",
        "result/kraken_results/{sample}/{sample}.kraken.report.std.txt",
        "result/kraken_results/{sample}/{sample}.kraken.output.txt"
    params:
        sample=get_sample,
        outdir="result/kraken_results",
        ncbi_blast_path=config['tools']['ncbi_blast_path'],
        Kraken2Uniq_path=config['tools']['Kraken2Uniq_path'],
        kraken_database_path=config['tools']['kraken_database_path'],
        kreport2mpa_path=config['tools']['kreport2mpa_path'],
        qsub_p=int(config['kraken2']['cpu'])
    benchmark:
        "benchmarks/kraken2.{sample}.benchmark.txt"
    log:
        "logs/kraken2/run.kraken2.{sample}.output"
    envmodules:
        envmodules["dev"]
    shell:
        """ 
        Rscript scripts/run_kraken.r \
                --sample {params.sample} \
                --fq1 {input.CB_reads} \
                --fq2 {input.cDNA_reads} \
                --out_path {params.outdir}/{params.sample}/ \
                --ncbi_blast_path {params.ncbi_blast_path} \
                --Kraken2Uniq_path {params.Kraken2Uniq_path} \
                --kraken_database_path  {params.kraken_database_path} \
                --kreport2mpa_path {params.kreport2mpa_path} \
                --threads {params.qsub_p} \
                --paired T
        """
