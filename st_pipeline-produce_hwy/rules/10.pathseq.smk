## 运行GATK.pathseq
rule GATK_pathseq:
    """
    运行GATK pathseq获取空间组织微生物含量 
    """
    input:
        barcodes="result/spaceranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    #bam="result/spaceranger/{sample}/outs/possorted_genome_bam.bam"
    output:
        "result/pathseq/GATK/{sample}.pathseq.complete.bam",
        "result/pathseq/GATK/{sample}.pathseq.complete.csv"
    params:
        outdir="result/pathseq/GATK/",
        sample=get_sample,
        pathseq_db=config['database']['pathseq_db'],
        qsub_mem=int(config['GATK_pathseq']['mem']),
        qsub_p=int(config['GATK_pathseq']['cpu'])
    benchmark:
        "benchmarks/GATK_pathseq.{sample}.benchmark.txt"
    log:
        "logs/GATK_pathseq/run.GATK_pathseq.{sample}.output"
    envmodules:
        envmodules["GATK"],
        envmodules["oesinglecell_path"]
    shell:
        """(
        mkdir -p {params.outdir}/tmp && 
        gatk --java-options "-Xmx{params.qsub_mem}g" PathSeqPipelineSpark \
             --tmp-dir={params.outdir}/tmp \
             --input  "result/spaceranger/{params.sample}/outs/possorted_genome_bam.bam" \
             --filter-bwa-image {params.pathseq_db}/pathseq_host.fa.img \
             --kmer-file {params.pathseq_db}/pathseq_host.bfi \
             --min-clipped-read-length 60 \
             --microbe-fasta {params.pathseq_db}/pathseq_microbe.fa \
             --microbe-bwa-image {params.pathseq_db}/pathseq_microbe.fa.img \
             --taxonomy-file {params.pathseq_db}/pathseq_taxonomy.db \
             --output {params.outdir}/{params.sample}.pathseq.complete.bam \
             --scores-output {params.outdir}/{params.sample}.pathseq.complete.csv \
             --is-host-aligned false \
             --filter-duplicates false \
             --min-score-identity .7    
        ) >{log} 2>&1      
	"""

## 可视化部分
rule pathseq_extract:
    """
    pathseq分析结果可视化
    """
    input:
        #"result/spaceranger/{sample}/outs/possorted_genome_bam.bam",
        "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "result/pathseq/GATK/{sample}.pathseq.complete.bam",
        "result/pathseq/GATK/{sample}.pathseq.complete.csv"
    output:
        "result/pathseq/pathseq_visualization/{sample}.visium.readname",
        "result/pathseq/pathseq_visualization/{sample}.visium.unmap_cbub.bam",
        "result/pathseq/pathseq_visualization/{sample}.visium.unmap_cbub.fasta",
        "result/pathseq/pathseq_visualization/{sample}.visium.list",
        "result/pathseq/pathseq_visualization/{sample}.visium.raw.readnamepath",
        "result/pathseq/pathseq_visualization/{sample}.visium.genus.cell",
        "result/pathseq/pathseq_visualization/{sample}.visium.genus.csv",
        "result/pathseq/pathseq_visualization/{sample}.visium.validate.csv"
    params:
        outdir="result/pathseq/pathseq_visualization",
        sample=get_sample,
        #pathseq_matrix=config['params']['pathseq_matrix'],
        qsub_mem=config['pathseq_extract']['mem'],
        qsub_p=config['pathseq_extract']['cpu']
    benchmark:
        "benchmarks/pathseq_extract.{sample}.benchmark.txt"
    log:
        "logs/GATK_pathseq/pathseq_extract.{sample}.output"
    envmodules:
        envmodules["oesinglecell_path"]
    shell:
        """
         python scripts/pathseq/UMI_annotator.py \
                  result/spaceranger/{params.sample}/outs/possorted_genome_bam.bam \
                  '' \
                  result/spaceranger/{params.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
                  result/pathseq/GATK/{params.sample}.pathseq.complete.bam \
                  result/pathseq/GATK/{params.sample}.pathseq.complete.csv \
                  {params.outdir}/{params.sample}.visium.readname \
                  {params.outdir}/{params.sample}.visium.unmap_cbub.bam \
                  {params.outdir}/{params.sample}.visium.unmap_cbub.fasta \
                  {params.outdir}/{params.sample}.visium.list \
                  {params.outdir}/{params.sample}.visium.raw.readnamepath \
                  {params.outdir}/{params.sample}.visium.genus.cell \
                  {params.outdir}/{params.sample}.visium.genus.csv \
                  {params.outdir}/{params.sample}.visium.validate.csv
        """

## 可视化部分
rule pathseq_visualization:
    """
    pathseq分析结果可视化
    """
    input:
        expand(["result/pathseq/pathseq_visualization/{sample}.visium.readname",
                "result/pathseq/pathseq_visualization/{sample}.visium.unmap_cbub.bam",
                "result/pathseq/pathseq_visualization/{sample}.visium.unmap_cbub.fasta",
                "result/pathseq/pathseq_visualization/{sample}.visium.list",
                "result/pathseq/pathseq_visualization/{sample}.visium.raw.readnamepath",
                "result/pathseq/pathseq_visualization/{sample}.visium.genus.cell",
                "result/pathseq/pathseq_visualization/{sample}.visium.genus.csv",
                "result/pathseq/pathseq_visualization/{sample}.visium.validate.csv",
                "result/cluster_seurat/singlecell_object.clustering_resolution0.4.rds"],
            sample=samples.index)
    output:
        "result/pathseq/pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution/summary_statistic.xls"
    params:
        outdir="result/pathseq/pathseq_visualization",
        qsub_mem=config['pathseq_visualization']['mem'],
        qsub_p=config['pathseq_visualization']['cpu']
    benchmark:
        "benchmarks/pathseq_visualization.benchmark.txt"
    log:
        "logs/GATK_pathseq/run.pathseq_visualization.output"
    envmodules:
        envmodules["oesinglecell_path"]
    shell:
        """ (
         scVis  -i  result/cluster_seurat/singlecell_object.clustering_resolution0.4.rds \
                    -f rds  \
                    -o result/pathseq/pathseq_visualization/ \
                   --assay SCT  \
                   --dataslot data \
                   pathseqvis \
                     --pathseq_dir result/pathseq/pathseq_visualization/ \
                     --groupby  sampleid,clusters  \
                     --palette  customecol2
        ) >{log} 2>&1
        """

## 网页报告
localrules: pathseq_report 
rule pathseq_report:
    """
    pathseq分析结果可视化
    """
    input:
        "result/pathseq/pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution/summary_statistic.xls",
        "config/config.yaml"
    output:
        f"result/report/{config['report']['Project_Num']}_PathSeq_Report/report.html"
    params:
        outdir="result/pathseq/pathseq_visualization",
        qsub_mem=config['pathseq_report']['mem'],
        qsub_p=config['pathseq_report']['cpu']
    benchmark:
        "benchmarks/pathseq_report.benchmark.txt"
    log:
        "logs/GATK_pathseq/run.pathseq_report.output"

    shell:
        """ (
         module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python  scripts/pathseq/report_pathseq.py  -i ./  -c  config/config.yaml
                 ) >{log} 2>&1
        """
