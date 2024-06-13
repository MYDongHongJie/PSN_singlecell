
rule indexbam:
    input:
        bam = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        bai = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
    threads:
        config["indexbam"]["cpu"]
    params:
        threads = 2*config["indexbam"]["cpu"]
    log:
        "logs/indexbam.{sample}.txt"
    benchmark:
        "benchmarks/indexbam.{sample}.txt"
    shell:
        """
        {samtools}
        samtools index -b {input.bam} -@ {params.threads} >{log} 2>&1
        """
if ref_vcf != "":
    # 人类数据使用 1000G phase3 数据作为参考
    rule cellsnplite:
        input:
            vcf = f'{ref_vcf}',
            bam = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
            bai = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
            barcode = "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz"
        output:
            basevcf = "result/Cellsnp/{sample}/cellSNP.base.vcf.gz",
            cellvcf = "result/Cellsnp/{sample}/cellSNP.cells.vcf.gz",
            # tsv = "result/Cellsnp/{sample}/cellSNP.samples.tsv"
        threads:
            config["cellsnplite"]["cpu"]
        params:
            outpath = "result/Cellsnp/{sample}",
            minMAF = config["params"]["cellsnplite"]["minMAF"],
            minCOUNT = config["params"]["cellsnplite"]["minCOUNT"],
            minLEN = config["params"]["cellsnplite"]["minLEN"],
            minMAPQ = config["params"]["cellsnplite"]["minMAPQ"],
            threads = 2*config["cellsnplite"]["cpu"]
        log:
            "logs/cellsnplite.withref.{sample}.txt"
        benchmark:
            "benchmarks/cellsnplite.withref.{sample}.txt"
        shell:
            """
            {cellsnplite} \
                -s {input.bam} \
                -b {input.barcode} \
                -R {input.vcf} \
                -p {params.threads} \
                --minMAF {params.minMAF} \
                --minCOUNT {params.minCOUNT} \
                --gzip \
                --minLEN {params.minLEN} \
                --minMAPQ {params.minMAPQ} \
                --genotype \
                -O {params.outpath}  >{log} 2>&1
            """
else:
    # 人类以外的其他物种暂时没有参考，使用 cellsnp-lite model
    rule cellsnplite:
        input:
            bam = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
            bai = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
            barcode = "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz"
        output:
            basevcf = "result/Cellsnp/{sample}/cellSNP.base.vcf.gz",
            cellvcf = "result/Cellsnp/{sample}/cellSNP.cells.vcf.gz",
            # tsv = "result/Cellsnp/{sample}/cellSNP.samples.tsv"
        threads:
            config["cellsnplite"]["cpu"]
        params:
            outpath = "result/Cellsnp/{sample}",
            minMAF = config["params"]["cellsnplite"]["minMAF"],
            minCOUNT = config["params"]["cellsnplite"]["minCOUNT"],
            minLEN = config["params"]["cellsnplite"]["minLEN"],
            minMAPQ = config["params"]["cellsnplite"]["minMAPQ"],
            threads = 2*config["cellsnplite"]["cpu"]
        log:
            "logs/cellsnplite.withoutref.{sample}.txt"
        benchmark:
            "benchmarks/cellsnplite.withoutref.{sample}.txt"
        shell:
            """
            {cellsnplite} \
                -s {input.bam} \
                -b {input.barcode} \
                -p {params.threads} \
                --minMAF {params.minMAF} \
                --minCOUNT 100 \
                --gzip \
                --minLEN {params.minLEN} \
                --minMAPQ {params.minMAPQ} \
                --genotype \
                -O {params.outpath}  >{log} 2>&1
            """
