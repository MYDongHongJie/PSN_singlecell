# rule starsolo_genome_index:
#     input:
#         genome = join(config['database']['ref_genome'], 'fasta', 'genome.fa'),
#         gtf = join(config['database']['ref_genome'], 'genes', 'genes.gtf')
#     output:
#         join(config['database']['ref_genome'], 'starsolo', 'SA')
#     params:
#         index_dir =  join(config['database']['ref_genome'], 'starsolo'),
#         readlength = config.get('read_geometry', [28, 98])[-1]
#     log:
#         join(config['database']['ref_genome'], 'logs', 'STAR.index.log')
#     resources:
#         qsub_mem = 30,
#         qsub_p = config['cpu']['starsolo_genome_index'],
#         qsub_n = 1
#     shell:
#         '''
#         STAR
#         --runThreadN {threads}
#         --runMode genomeGenerate
#         --genomeDir {params.index_dir}
#         --genomeFastaFiles {input.genome}
#         --sjdbGTFfile {input.gtf}
#         --sjdbOverhang {params.readlength}
#         && mv Log.out {log}
#         '''

rule starsolo:
    input:
        CB_reads=get_fastqs_cb_reads,
        cDNA_reads=get_fastqs_cdna_reads,
        white_list=config['database']['white_list']
    output:
        barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
        matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
        raw_barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
        raw_features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
        raw_matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz",
        bam="result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        summary="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/Summary.csv"
    benchmark:
        "benchmarks/starsolo/{sample}.starsolo.benchmark.txt"
    log:
        "logs/starsolo/{sample}.starsolo.log"
    envmodules:
        config['envmodules']['star']
    params:
        sample=get_sample,
        STAR_index=config['database']['ref_genome'],
        qsub_mem=config['starsolo']['mem'],
        qsub_p=config['starsolo']['cpu']
    shell:
        '''
           STAR \
                --runThreadN     {params.qsub_p} \
                --soloType       CB_UMI_Simple \
                --soloCBwhitelist {input.white_list} \
                --genomeDir       {params.STAR_index}/STAR_index  \
                --soloCBstart      1  \
                --soloCBlen        20 \
                --soloUMIstart     21 \
                --soloUMIlen       8 \
                --readFilesCommand  'gzip -c -d'  \
                --readFilesIn     {input.cDNA_reads}   {input.CB_reads} \
                --outSAMtype BAM SortedByCoordinate \
                --outMultimapperOrder Random \
                --runRNGseed 1 \
                --outSAMattributes NH HI AS CB UB CR UR GX GN \
                --soloFeatures Gene GeneFull \
                --soloUMIdedup Exact \
                --outFileNamePrefix result/1.STARsolo/{params.sample}/{params.sample}_ \
                --soloStrand Reverse \
                --limitBAMsortRAM 41143265264  \
                --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 &&
            pigz -p  10 {params.sample}/{params.sample}_Solo.out/*/{{barcodes.tsv,features.tsv,matrix.mtx}}
            '''

if config['params']['STARsolo']['empty_drops'] is None:
    rule starsolo_filer:
        input:
            raw_barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
            raw_features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
            raw_matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz"
        output:
            barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/re-filtered/barcodes.tsv.gz",
            features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/re-filtered/features.tsv.gz",
            matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/re-filtered/matrix.mtx.gz"
        benchmark:
            "benchmarks/starsolo/{sample}.starsolo.re-filtered.benchmark.txt"
        log:
            "logs/starsolo/{sample}.starsolo.re-filtered.log"
        envmodules:
            config['envmodules']['star']
        params:
            sample=get_sample(),
            empty_drops=config["params"]["STARsolo"]["empty_drops"],
            qsub_mem=30,
            qsub_p=config['starsolo_filer']['cpu'],
        shell:
            '''(
            STAR --runMode soloCellFiltering  \
                 result/1.STARsolo/{params.sample}/{params.sample}_Solo.out/GeneFull/raw/ \
                 result/1.STARsolo/{params.sample}/{params.sample}_Solo.out/GeneFull/re-filtered/ \
                 --soloCellFilter EmptyDrops_CR  {params.empty_drops} && 
            pigz -p  10 result/1.STARsolo/{params.sample}/{params.sample}_Solo.out/GeneFull/*/*
           ) >{log} 2>&1'''

rule barcode_rank_gene_umi:
    input:
        barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
        matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
        raw_barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
        raw_features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
        raw_matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz"
    output:
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.png",
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.pdf",
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.png",
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.pdf",
    benchmark:
        "benchmarks/barcode_rank_gene_umi/{sample}.barcode_vs_umi.benchmark.txt"
    log:
        "logs/barcode_rank_gene_umi/{sample}.barcode_vs_umi.log"
    envmodules:
        config['envmodules']['oesinglecell_visium']
    params:
        sample=get_sample,
        qsub_mem=config['barcode_rank_gene_umi']['mem'],
        qsub_p=config['barcode_rank_gene_umi']['cpu'],
    params:
        sample=get_sample
    shell:
        '''(
        ./scripts/scVis \
              -o result/1.STARsolo/BarcodeRank_plot/ \
              barcode_rank_gene_umi \
              --starsolo result/1.STARsolo/{params.sample}/{params.sample}_Solo.out/GeneFull/ \
              --prefix  {params.sample}    
        ) >{log} 2>&1'''

rule reads_distribution:
    input:
        bam="result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        refFlat=f"{ref_dir}/refFlat.txt",
        rRNA_interval=f"{ref_dir}/rRNA.interval_list"
    output:
        rna_metrics="result/1.STARsolo/Read_Distribution/{sample}_rna_metrics.txt",
        png="result/1.STARsolo/Read_Distribution/reads_distribution_for_{sample}.pdf"
    benchmark:
        "benchmarks/read_distribution/{sample}.reads_distribution.benchmark.txt"
    log:
        "logs/read_distribution/{sample}.reads_distribution.log"
    envmodules:
        config['envmodules']['picard'],
        config['envmodules']['oesinglecell_visium']
    params:
        sample=get_sample,
        qsub_mem=config['reads_distribution']['mem'],
        qsub_p=config['reads_distribution']['cpu']
    shell:
        '''(
            java -Xmx{params.qsub_mem}G \
                 -jar ${{PicardCLASSPATH}}  \
                 CollectRnaSeqMetrics \
                    I={input.bam}\
                    O={output.rna_metrics}\
                    STRAND_SPECIFICITY=NONE \
                    REF_FLAT={input.refFlat} \
                    RIBOSOMAL_INTERVALS={input.rRNA_interval}  &&  \
            ./scripts/scVis \
                  -o  result/1.STARsolo/Read_Distribution/ \
                  read_distribution \
                    --rna_matrics {output.rna_metrics} \
                    --prefix {params.sample}
        ) >{log} 2>&1'''

localrules: starsolo_summary

rule starsolo_summary:
    '''
    汇总starsolos统计结果
    '''
    input:
        metadatafile="config/samples.csv",
        cellranger_out=expand(
            [
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/Summary.csv"
            ],sample=samples.index)
    output:
        "result/1.STARsolo/STARsolo.statistic.tsv"
    benchmark:
        "benchmarks/starsolo_summary/starsolo_summary.benchmark.txt"
    log:
        "logs/starsolo_summary/starsolo_summary.log"
    envmodules:
        config['envmodules']['oesinglecell_visium']
    shell:
        '''
         ./scripts/sctool -o result/1.STARsolo/ summary -s {input.metadatafile} -i result/1.STARsolo/  >&2 2>{log} 
        '''

localrules: QC_Report

rule QC_Report:
    '''
    report
    '''
    input:
        files=expand(["result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.png",
                      "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.pdf",
                      "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.png",
                      "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.pdf",
                      "result/1.STARsolo/Read_Distribution/{sample}_rna_metrics.txt",
                      "result/1.STARsolo/STARsolo.statistic.tsv",
                      "result/1.STARsolo/Read_Distribution/reads_distribution_for_{sample}.pdf",
                      "result/1.STARsolo/STARsolo.statistic.tsv"],sample=samples.index),
        config="config/config.yaml",
    output:
        html=f"result/report/{project}_QC_Report_{time.strftime('%Y_%m_%d')}/Report.html"
    log:
        "logs/report_qc.log"
    params:
        report=f"result/report/{project}_QC_Report_{time.strftime('%Y_%m_%d')}/"
    envmodules:
        config['envmodules']['oesinglecell_visium']
    shell:
        '''
        python scripts/report_QC_single_species.py -i result -c {input.config} |& tee {log}
        '''
