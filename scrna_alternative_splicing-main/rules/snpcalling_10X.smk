if ref_vcf != "":
    # 人类数据使用 1000G phase3 数据作为参考
    rule cellsnplite:
        input:
            vcf = f'{ref_vcf}',
            bam = "result/cellranger/{sample}/outs/possorted_genome_bam.bam",
            bai = "result/cellranger/{sample}/outs/possorted_genome_bam.bam.bai",
            barcode = "result/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
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
            bam = "result/cellranger/{sample}/outs/possorted_genome_bam.bam",
            bai = "result/cellranger/{sample}/outs/possorted_genome_bam.bam.bai",
            barcode = "result/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
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
