##数据准备
rule data_praper:
    """
    prepare fastqs for each sample
    """
    input:
        config = "config/config.yaml"
    output:
        praper_done = "logs/data_praper.done"
    benchmark:
        "benchmarks/data_prepare.benchmark.txt"
    log:
        "logs/data_praper/data_prepare.log"
    shell:
        '''
        python scripts/data_prepare.py -c {input.config}  -o raw_data/  > {log} 2>&1 && touch {output.praper_done}
        '''
## raw_data转clean_data
if config['params']['STARsolo']['raw_data']:
    rule raw2clean:
        input:
            praper_done ="logs/data_praper.done"
        output:
            donefile = "logs/raw2clean_{sample}.done"
        benchmark:
            "benchmarks/raw2clean_{sample}.benchmark.txt"
        log:
            "logs/data_praper/raw2clean_{sample}.log"
        params:
            config_file = "config/raw2clean/{sample}_config.info"
        shell:
            '''
            {mp}
            {m20sc_sctool}
            python /home/sunkun/m20/VITASEER/VITASEER.alpha.py --cfg {params.config_file} > {log} 2>&1 && touch {output.donefile}
            '''

##比对
rule starsolo:
    input:
        praper_done = expand(
            [
                "logs/raw2clean_{sample}.done"
            ],sample=samples.index) if config['params']['STARsolo']['raw_data'] else "logs/data_praper.done",
        gtf = ancient(f"{ref_dir}/genome.gtf"),
        rRNA = ancient(f"{ref_dir}/rRNA.interval_list")
    output:
        barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
        matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
        raw_barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
        raw_features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
        raw_matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz",
        bam="result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        summary="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/Summary.csv",
        donefile= "logs/starsolo_{sample}.done"
    benchmark:
        "benchmarks/starsolo/{sample}.starsolo.benchmark.txt"
    log:
        "logs/starsolo/{sample}.starsolo.log"
    params:
        STAR_index = ref_dir,
        CB_reads= "raw_data/{sample}/{sample}_fastqs/{sample}_R1.fq.gz",
        cDNA_reads="raw_data/{sample}/{sample}_fastqs/{sample}_R2.fq.gz"
    resources:
        qsub_mem = config['starsolo']['mem']
    threads:
        config['starsolo']['cpu']
    shell:
        '''
        {mp}
        {star}
           STAR \
                --runThreadN     {threads} \
                --soloType       CB_UMI_Simple \
                --genomeDir       {params.STAR_index}/STAR_index  \
                --soloCBstart      1  \
                --soloCBwhitelist None \
                --soloCBlen        20 \
                --soloUMIstart     21 \
                --soloUMIlen       8 \
                --readFilesCommand  'gzip -c -d'  \
                --readFilesIn     {params.cDNA_reads}   {params.CB_reads} \
                --outSAMtype BAM SortedByCoordinate \
                --outMultimapperOrder Random \
                --runRNGseed 1 \
                --outSAMattributes NH HI AS CB UB CR UR GX GN CY UY \
                --soloFeatures GeneFull \
                --soloUMIdedup Exact \
                --outFileNamePrefix result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_ \
                --soloStrand Reverse \
                --limitBAMsortRAM 41143265264  \
                --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 300 0.01 20000 0.01 10000 &&
            pigz -p {threads} -f result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/*/*/{{barcodes.tsv,features.tsv,matrix.mtx}} >{log} 2>&1 && touch {output.donefile}
            '''

if filter_standard != "":
    rule starsolo_filer:
        input:
            raw_barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
            raw_features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
            raw_matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz"
        output:
            donefile = "logs/re_filter_{sample}.done"
        benchmark:
            "benchmarks/starsolo/{sample}.starsolo.re-filtered.benchmark.txt"
        log:
            "logs/starsolo/{sample}.starsolo.re-filtered.log"
        params:
            empty_drops = filter_standard
        resources:
            qsub_mem = 30
        threads:
            config['starsolo_filer']['cpu']
        shell:
            '''(
            {mp}
            {star}
            gzip -d result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/raw/{{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz}} &&
            rm -rf result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/filtered
            STAR --runMode soloCellFiltering  \
                 result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/raw/ \
                 result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/filtered/ \
                 --soloCellFilter TopCells  {params.empty_drops} && 
            pigz -p  10 result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/{{raw,filtered}}/{{barcodes.tsv,features.tsv,matrix.mtx}}
           ) >{log} 2>&1 && touch {output.donefile}'''

### 拐点图及基因分布图 ###
rule barcode_rank_gene_umi:
    input:
        barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
        matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
        raw_barcode="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
        raw_features="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
        raw_matrix="result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz",
        starsolo_file= "logs/starsolo_{sample}.done" if filter_standard=="" else "logs/re_filter_{sample}.done"
    output:
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.png",
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.pdf",
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.png",
        "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.pdf",
        "result/1.STARsolo/Genetype_plot/{sample}_gene_type_totaldata.csv",
        "result/1.STARsolo/Genetype_plot/{sample}_gene_type_plotdata.csv",
        "result/1.STARsolo/Genetype_plot/{sample}_gene_type.png"
    benchmark:
        "benchmarks/barcode_rank_gene_umi/{sample}.barcode_vs_umi.benchmark.txt"
    log:
        "logs/barcode_rank_gene_umi/{sample}.barcode_vs_umi.log"
    resources:
        qsub_mem = config['barcode_rank_gene_umi']['mem']
    threads:
        config['barcode_rank_gene_umi']['cpu']
    params:
        STAR_index=ref_dir
    shell:
        '''(
        {mp}
        {oesinglecell_visium}
        ./scripts/scVis \
            -o result/1.STARsolo/BarcodeRank_plot/ \
            --ncores {threads} \
            barcode_rank_gene_umi \
            --starsolo result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/ \
            --prefix  {wildcards.sample}   \
            --filtered filtered && \
        ./scripts/scVis \
            -o result/1.STARsolo/Genetype_plot/ \
            --ncores {threads} \
            genetype \
            --starsolo result/1.STARsolo/{wildcards.sample}/{wildcards.sample}_Solo.out/GeneFull/ \
            --prefix  {wildcards.sample}   \
            --refgenome {params.STAR_index}
        ) >{log} 2>&1'''

### reads 分布 ###
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
    resources:
        qsub_men = config['reads_distribution']['mem']
    threads:
        config['reads_distribution']['cpu']
    shell:
        '''(
        {mp}
        {picard}
        {oesinglecell_visium}
            java -Xmx{resources.qsub_men}G \
                 -jar ${{PicardCLASSPATH}}  \
                 CollectRnaSeqMetrics \
                    I={input.bam} \
                    O={output.rna_metrics} \
                    STRAND_SPECIFICITY=NONE \
                    REF_FLAT={input.refFlat} \
                    RIBOSOMAL_INTERVALS={input.rRNA_interval}  &&  \
            ./scripts/scVis \
                  -o  result/1.STARsolo/Read_Distribution/ \
                  --ncores {threads} \
                  read_distribution \
                    --rna_matrics {output.rna_metrics} \
                    --prefix {wildcards.sample}
        ) >{log} 2>&1'''

## qualimap result ##
rule qualimap:
    '''
    qualimap result
    '''
    input:
        bam = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        'result/1.STARsolo/Qualimap_result/{sample}/rnaseq_qc_results.txt',
        'result/1.STARsolo/Qualimap_result/{sample}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt'
    benchmark:
        'benchmarks/qualimap/{sample}.qualimap.benchmark.txt'
    log:
        'logs/qualimap/{sample}.qualimap.log'
    resources:
        qsub_men = config['qualimap']['mem']
    params:
        STAR_index = ref_dir,
        outfile = "result/1.STARsolo/Qualimap_result/{sample}"
    shell:
        '''(
        {mp}
        {qualimap}
        qualimap rnaseq -outdir {params.outfile} -bam {input.bam} -gtf {params.STAR_index}/genome.gtf --java-mem-size={resources.qsub_men}G
        ) >{log} 2>&1
        '''

## qualimap_summary
rule qualimap_summary:
    input:
        qualimap_result=expand([
            'result/1.STARsolo/Qualimap_result/{sample}/rnaseq_qc_results.txt',
            'result/1.STARsolo/Qualimap_result/{sample}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt'
        ],sample=samples.index)
    output:
        'result/1.STARsolo/Qualimap_result/qualimap_plot.png'
    benchmark:
        'benchmarks/qualimap/qualimap_summary.benchmark.txt'
    log:
        'logs/qualimap/qualimap_summary.log'
    params:
        platte = config['params']['STARsolo']['platte'],
        inputfile= "result/1.STARsolo/Qualimap_result"
    shell:
        '''(
        {mp}
        {oesinglecell_visium}
        Rscript scripts/plot_qualimap.R -i {params.inputfile} -o {params.inputfile} -p {params.platte}
        ) >{log} 2>&1
        '''

### 总结数据质量 ###
rule starsolo_summary:
    '''
    汇总starsolos统计结果
    '''
    input:
        samplefile = samplefile,
        cellranger_out=expand(
            [
                "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/Summary.csv"
            ],sample=samples.index),
        starsolo_file= expand(
            [
                "logs/starsolo_{sample}.done"
            ],sample=samples.index) if filter_standard=="" else expand(
            [
                "logs/re_filter_{sample}.done"
            ],sample=samples.index)
    output:
        "result/1.STARsolo/STARsolo.statistic.tsv"
    benchmark:
        "benchmarks/starsolo_summary/starsolo_summary.benchmark.txt"
    params:
        re_filtered= "--filtered filtered" if filter_standard=="" else "--filtered re-filtered"
    log:
        "logs/starsolo_summary/starsolo_summary.log"
    shell:
        '''
        {mp}
        {oesinglecell_visium}
        ./scripts/sctool -o result/1.STARsolo/ summary -s {input.samplefile} -i result/1.STARsolo/  {params.re_filtered} >&2 2>{log} 
        '''

### 出具质控报告 ###
if config["module"]["QCreport"]:
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
                          ],
                          sample = samples.index
                          ),
            config="config/config.yaml",
        output:
            # html = f"result/report/{project}_QC_Report_{time.strftime('%Y_%m_%d')}/Report.html",
            html = f"result/report/{project}_QC_Report/Report.html",
            donefile = "logs/STARsolo.done"
        log:
            "logs/report_qc.log"
        shell:
            '''
            {mp}
            /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report_QC.py -i result -c {input.config}  >&2 2>{log} && touch {output.donefile}
            '''
