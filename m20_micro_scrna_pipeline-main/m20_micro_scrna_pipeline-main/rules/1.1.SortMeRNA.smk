rule remove_sortmerna:
    input:
        CB_reads="result/clean_data/{sample}/{sample}.R1.gz",
        cDNA_reads="result/clean_data/{sample}/{sample}.R2.gz",
        rRNA_ref=config['params']['sortmerna']['ref_db']
    output:
        CB_reads="result/remove_sortmerna_pollution/after_remove/{sample}/{sample}.R1.fq.gz",
        cDNA_reads="result/remove_sortmerna_pollution/after_remove/{sample}/{sample}.R2.fq.gz"
    benchmark:
        "benchmarks/remove_sortmerna_pollution/{sample}.txt"
    log:
        "logs/remove_sortmerna_pollution/after_remove/{sample}.log"
    params:
        outdir="result/remove_sortmerna_pollution/after_remove/",
        prefix="{sample}"
    threads: 8
    shell:
        """
        sortmerna --ref   {input.rRNA_ref}  \
                  --reads   {input.cDNA_reads} \
                  --workdir result/sortmerna \
                  --num_alignments 1 \
                  --fastx -a 8 \
                  1>{log} 2>&1
        """
