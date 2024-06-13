

localrules: white_list_prepare
rule white_list_prepare:
    input:
        barcode = "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz"
    output:
        white_list = ("result/Splicing/STARSolo/{sample}/white_list.txt")
    params:
        outpath = "result/Splicing/{sample}"
    shell:
        """
        gunzip -c {input.barcode} > {output.white_list}
        """


rule starsolo_extract_SJ:
    input:
        bam = "result/1.STARsolo/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        white_list = "result/Splicing/STARSolo/{sample}/white_list.txt"
    output:
        genebarcodes = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv"),
        genefeatures = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv"),
        genematrix = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx"),
        sjbarcodes = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv"),
        sjfeatures = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/features.tsv"),
        sjmatrix = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx"),
        sjout = "result/Splicing/STARSolo/{sample}/{sample}_SJ.out.tab",
        rawbarcodes = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv"),
        rawfeatures = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv"),
        rawmatrix = temp("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx"),
        sam = temp("result/Splicing/STARSolo/{sample}/{sample}_Aligned.out.sam"),
        logfinal = temp("result/Splicing/STARSolo/{sample}/{sample}_Log.final.out"),
        logout = temp("result/Splicing/STARSolo/{sample}/{sample}_Log.out"),
        logprogress = temp("result/Splicing/STARSolo/{sample}/{sample}_Log.progress.out")
    threads:
        config['starsolo_extract_SJ']['cpu']
    log:
        "logs/STARSolo/starsolo_extract_SJ.{sample}.log"
    benchmark:
        "benchmarks/STARSolo/starsolo_extract_SJ.{sample}.txt"
    params:
        STAR_index = ref_dir,
        outprefix = "result/Splicing/STARSolo/{sample}/{sample}_"
    resources:
        qsub_mem = config['starsolo_extract_SJ']['mem']
    shell:
        """
        {mp}
        {samtools}
        {star}
        STAR --runThreadN {threads} \
            --genomeDir {params.STAR_index}/STAR_index \
            --soloType CB_UMI_Simple \
            --readFilesIn {input.bam} \
            --soloCBstart      1 \
            --soloCBlen        20 \
            --soloUMIstart     21 \
            --soloUMIlen       8 \
            --readFilesCommand samtools view -F 0x100 \
            --readFilesType SAM SE  \
            --soloInputSAMattrBarcodeSeq CR UR  \
            --soloInputSAMattrBarcodeQual CY UY  \
            --soloCBwhitelist {input.white_list} \
            --outMultimapperOrder Random \
            --outSJfilterReads Unique \
            --soloFeatures GeneFull SJ \
            --soloStrand Reverse \
            --soloUMIdedup Exact \
            --runRNGseed 1 \
            --outFileNamePrefix {params.outprefix}  >{log} 2>&1
        """

localrules: deal_sj_features
rule deal_sj_features:
    input:
        genebarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv",
        genefeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv",
        genematrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx",
        sjbarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv",
        sjmatrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx",
        sjout = "result/Splicing/STARSolo/{sample}/{sample}_SJ.out.tab",
    output:
        genebarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        genefeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
        genematrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
        sjbarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv.gz",
        sjfeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/features.tsv.gz",
        sjmatrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx.gz",
    shell:
        """
        gzip -c {input.genebarcodes} > {output.genebarcodes}
        gzip -c {input.genefeatures} > {output.genefeatures}
        gzip -c {input.genematrix} > {output.genematrix}
        gzip -c {input.sjbarcodes} > {output.sjbarcodes}
        gzip -c {input.sjmatrix} > {output.sjmatrix}
        awk '{{print $1":"$2":"$3"\\ttest\\tSJ"}}' {input.sjout} | gzip > {output.sjfeatures}
        """


if use_routine_genearray:
    localrules: prepare_genearray
    rule prepare_genearray:
        input:
            genebarcodes = "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz", 
            genefeatures = "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz", 
            genematrix = "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz"
        output:
            genebarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/barcodes.tsv.gz",
            genefeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/features.tsv.gz", 
            genematrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/matrix.mtx.gz"
        shell:
            """
            cp {input.genefeatures} {output.genefeatures}
            cp {input.genematrix} {output.genematrix}
            cp {input.genebarcodes} {output.genebarcodes}
            """

    rule MARVEL_prepare:
        input:
            genebarcodes = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/barcodes.tsv.gz",  sample = samples.index),
            genefeatures = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/features.tsv.gz",  sample = samples.index),
            genematrix = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/matrix.mtx.gz",  sample = samples.index),
            sjbarcodes = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv.gz",  sample = samples.index),
            sjfeatures = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/features.tsv.gz",  sample = samples.index),
            sjmatrix = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx.gz",  sample = samples.index),
            coordinatefile = {dimension_type_file},
            coordubatecluster = {dimension_type_cluster},
            samplefile = {samplefile},
            h5seurat = "result/2.Count_QC/seurat_ob.h5seurat"
        output:
             MARVEL_data = "result/Splicing/MARVEL.RData"
        log:
            "logs/MARVEL_prepare.log"
        params:
            gtf = {gtf}
        threads:
            config['MARVEL_prepare']['cpu']
        resources:
            qsub_mem = config['MARVEL_prepare']['mem']
        benchmark:
            "benchmarks/MARVEL_prepare.txt"
        shell:
            """
            {mp}
            {Rscript_splicing}
            Rscript {marvel_prepare} \
                --samplefile {input.samplefile}  \
                --gtf {params.gtf} \
                --coordinate {input.coordinatefile} \
                --cluster {input.coordubatecluster} \
                --h5seurat {input.h5seurat} \
                --scale_factor 1e6 \
                --out_MARVEL_data {output.MARVEL_data} >{log} 2>&1
            """     

else:
    rule MARVEL_prepare:
        input:
            genebarcodes = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",  sample = samples.index),
            genefeatures = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz", sample = samples.index),
            genematrix = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",  sample = samples.index),
            sjbarcodes = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv.gz", sample = samples.index),
            sjfeatures = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/features.tsv.gz", sample = samples.index),
            sjmatrix = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx.gz", sample = samples.index),
            coordinatefile = {dimension_type_file},
            coordubatecluster = {dimension_type_cluster},
            samplefile = {samplefile},
            h5seurat = "result/2.Count_QC/seurat_ob.h5seurat"
        output:
            MARVEL_data = "result/Splicing/MARVEL.RData"
        log:
            "logs/MARVEL_prepare.log"
        params:
            gtf = {gtf},
            splincingpath = "result/Splicing/STARSolo"
        threads:
            config['MARVEL_prepare']['cpu']
        resources:
            qsub_mem = config['MARVEL_prepare']['mem']
        benchmark:
            "benchmarks/MARVEL_prepare.txt"
        shell:
            """
            {mp}
            {Rscript_splicing} 
            Rscript {marvel_prepare} \
                --samplefile {input.samplefile}  \
                --gtf {params.gtf} \
                --coordinate {input.coordinatefile} \
                --cluster {input.coordubatecluster} \
                --h5seurat {input.h5seurat} \
                --scale_factor 1e6 \
                --splincingpath {params.splincingpath} \
                --out_MARVEL_data {output.MARVEL_data} >{log} 2>&1
            """


