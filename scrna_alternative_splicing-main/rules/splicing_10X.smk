import os

fileexist = []
for s in samples.index:
    bamfile = "result/cellranger/" + s + "/outs/possorted_genome_bam.bam"
    fileexist.append(os.path.exists(bamfile))

if not all(fileexist):
    if 'old_name' in samples.columns:
        rule get_bam:
            output:
                bam = temp("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam"),
                bai = temp("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam.bai")
            params:
                bam_path = {bam_path},
                # Raw_Task_Num = {Raw_Task_Num},
                password = {password},
                Raw_Task_Num = lambda w: samples.loc["{}".format(w.sample)]["raw_task_number"],
                old_name = lambda w: samples.loc["{}".format(w.sample)]["old_name"],
                sampleid = lambda w: samples.loc["{}".format(w.sample)]["sampleid"]
            log:
                "logs/bam/{sample}.log"
            shell:
                """
                {sshpass}
                sshpass -p {params.password} scp {params.bam_path}/{params.Raw_Task_Num}/{params.old_name}/possorted_genome_bam.bam result/Splicing/STARSolo/{params.sampleid}/Bam/  >{log} 2>&1
                sshpass -p {params.password} scp {params.bam_path}/{params.Raw_Task_Num}/{params.old_name}/possorted_genome_bam.bam.bai result/Splicing/STARSolo/{params.sampleid}/Bam/ >>{log} 2>&1
                """
    else:
        rule get_bam:
            output:
                bam = temp("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam"),
                bai = temp("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam.bai")
            params:
                bam_path = {bam_path},
                # Raw_Task_Num = {Raw_Task_Num},
                Raw_Task_Num = lambda w: samples.loc["{}".format(w.sample)]["raw_task_number"],
                password = {password},
            log:
                "logs/bam/{sample}.log"
            shell:
                """
                {sshpass}
                sshpass -p {params.password} scp {params.bam_path}/{params.Raw_Task_Num}/{wildcards.sample}/possorted_genome_bam.bam result/Splicing/STARSolo/{wildcards.sample}/Bam/  >{log} 2>&1
                sshpass -p {params.password} scp {params.bam_path}/{params.Raw_Task_Num}/{wildcards.sample}/possorted_genome_bam.bam.bai result/Splicing/STARSolo/{wildcards.sample}/Bam/ >>{log} 2>&1
                """
else:
    localrules: get_bam_link
    rule get_bam_link:
        output:
            bam = ("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam"),
            bai = ("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam.bai")
        shell:
            """
            path=$(pwd)
            ln -f -s $path/result/cellranger/{wildcards.sample}/outs/possorted_genome_bam.bam  $path/result/Splicing/STARSolo/{wildcards.sample}/Bam/possorted_genome_bam.bam
            ln -f -s $path/result/cellranger/{wildcards.sample}/outs/possorted_genome_bam.bam.bai  $path/result/Splicing/STARSolo/{wildcards.sample}/Bam/possorted_genome_bam.bam.bai
            """

localrules: white_list_prepare
rule white_list_prepare:
    input:
        barcode = "result/cellranger/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz"
    output:
        white_list = temp("result/Splicing/STARSolo/{sample}/white_list.txt")
    shell:
        """
        gunzip -c {input.barcode} > {output.white_list}
        sed -i 's/-1//' {output.white_list}
        """

# 从cellranger比对的bam文件开始，提取gene和SJ
rule starsolo_extract_SJ:
    input:
        bam = "result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam",
        white_list = "result/Splicing/STARSolo/{sample}/white_list.txt"
    output:
        genebarcodes = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/barcodes.tsv"),
        genefeatures = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/features.tsv"),
        genematrix = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/matrix.mtx"),
        sjbarcodes = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv"),
        sjfeatures = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/features.tsv"),
        sjmatrix = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx"),
        sjout = ("result/Splicing/STARSolo/{sample}/{sample}_SJ.out.tab"),
        rawbarcodes = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/raw/barcodes.tsv"),
        rawfeatures = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/raw/features.tsv"),
        rawmatrix = ("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/raw/matrix.mtx"),
        sam = ("result/Splicing/STARSolo/{sample}/{sample}_Aligned.out.sam"),
        logfinal = ("result/Splicing/STARSolo/{sample}/{sample}_Log.final.out"),
        logout = ("result/Splicing/STARSolo/{sample}/{sample}_Log.out"),
        logprogress = ("result/Splicing/STARSolo/{sample}/{sample}_Log.progress.out")
    threads:
        config['starsolo_extract_SJ']['cpu']
    log:
        "logs/STARSolo/starsolo_extract_SJ.{sample}.log"
    benchmark:
        "benchmarks/STARSolo/starsolo_extract_SJ.{sample}.txt"
    params:
        genomeDir = genomeDir + "/STAR_index",
        outprefix = "result/Splicing/STARSolo/{sample}/{sample}_",
        soloUMIlen = {soloUMIlen}
    resources:
        qsub_mem = config['starsolo_extract_SJ']['mem']
    shell:
        """
        {mp}
        {samtools}
        {star}
        STAR --runThreadN {threads} \
            --genomeDir {params.genomeDir} \
            --soloType CB_UMI_Simple \
            --readFilesIn {input.bam} \
            --soloCBstart 1 \
            --soloCBlen 16 \
            --soloUMIstart 17 \
            --soloUMIlen {params.soloUMIlen} \
            --readFilesCommand samtools view -F 0x100 \
            --readFilesType SAM SE \
            --soloInputSAMattrBarcodeSeq CR UR \
            --soloInputSAMattrBarcodeQual CY UY \
            --soloCBwhitelist {input.white_list} \
            --outMultimapperOrder Random \
            --outSJfilterReads Unique \
            --soloFeatures Gene SJ \
            --soloUMIdedup Exact \
            --runRNGseed 1 \
            --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
            --outFileNamePrefix {params.outprefix} >{log} 2>&1
        """
# --soloStrand Reverse \
localrules: deal_sj_features
rule deal_sj_features:
    input:
        genebarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/barcodes.tsv",
        genefeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/features.tsv",
        genematrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/matrix.mtx",
        sjbarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv",
        sjmatrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx",
        sjout = "result/Splicing/STARSolo/{sample}/{sample}_SJ.out.tab",
    output:
        genebarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/barcodes.tsv.gz",
        genefeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/features.tsv.gz",
        genematrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/matrix.mtx.gz",
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
            genebarcodes = "result/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
            genefeatures = "result/cellranger/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz", 
            genematrix = "result/cellranger/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
        output:
            genebarcodes = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/barcodes.tsv.gz",
            genefeatures = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/features.tsv.gz", 
            genematrix = "result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene_upstream/filtered/matrix.mtx.gz"
        shell:
            """
            cp {input.genefeatures} {output.genefeatures}
            cp {input.genematrix} {output.genematrix}
            gunzip -c {input.genebarcodes} | sed 's/-1//' | gzip > {output.genebarcodes}
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
            h5seurat = "result/Count_QC/filtered.h5seurat",
            samplefile = {samplefile}
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
            genebarcodes = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/barcodes.tsv.gz",  sample = samples.index),
            genefeatures = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/features.tsv.gz", sample = samples.index),
            genematrix = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/Gene/filtered/matrix.mtx.gz",  sample = samples.index),
            sjbarcodes = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/barcodes.tsv.gz", sample = samples.index),
            sjfeatures = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/features.tsv.gz", sample = samples.index),
            sjmatrix = expand("result/Splicing/STARSolo/{sample}/{sample}_Solo.out/SJ/raw/matrix.mtx.gz", sample = samples.index),
            coordinatefile = {dimension_type_file},
            coordubatecluster = {dimension_type_cluster},
            h5seurat = "result/Count_QC/filtered.h5seurat",
            samplefile = {samplefile}
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
                # --SJ_datapath {params.SJ_datapath} \
                # --Gene_datapath {params.Gene_datapath} \