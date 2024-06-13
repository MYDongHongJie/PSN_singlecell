##============================================= step1: Differential expression analysis =======================================================
rule diffexp:
        '''
        use sctool to run differential expression analysis with RNA assay
        '''
        input:
            donefile = "logs/sub_cluster.done"
        output:
            difffile = "result/7.Diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.xls",
        params:
            inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
            outpath = "result/7.Diffexp",
            FC = foldchange,
            pvalue_or_fdr = "--fdr " if (padj) else "--pvalue ",
            sgin_thread =  padj if (padj) else pvalue,
            test = test,
            treatment = lambda w:diff_groups.loc["{}".format(w.diff_group)]["treatment"],
            control = lambda w:diff_groups.loc["{}".format(w.diff_group)]["control"],
            type = lambda w:diff_groups.loc["{}".format(w.diff_group)]["type"],
        log:
            "logs/diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.log"
        benchmark:
            "benchmarks/diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.benchmark.txt"
        resources:
            qsub_mem=config['diffexp']['mem']
        threads:
            config['diffexp']['cpu']
        shell:
            '''
            {mp}
            {m20sc_sctool}
            sctool -i {params.inputfile} \
            -f h5seurat \
            -o {params.outpath} \
            --assay SCT \
            -d h5seurat \
            --dataslot data,counts  \
            -j {threads} \
            diffexp --contrast {params.type}:{params.treatment}:{params.control} \
                --FC {params.FC} \
                {params.pvalue_or_fdr} {params.sgin_thread}  \
                --test {params.test} \
                -v {params.type} \
                -u {params.treatment},{params.control} >&2 2>{log}
            '''


##==================================================step2: plot diff_heatmap ==============================================
rule diffexp_heatmap:
    '''
    use sctool to plot diff_heatmap with RNA assay
    '''
    input:
        difffile = "result/7.Diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.xls"
    output:
        figurepdf = "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_heatmap.pdf",
        figurepng = "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_heatmap.png",
        table = "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_genes.xls"
    benchmark:
        "benchmarks/diffexp/diffexp_heatmap-top{topn}-{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}-top{topn}.benchmark.txt"
    log:
        "logs/diffexp/diffexp_heatmap-top{topn}-{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}-top{topn}.log"
    resources:
        qsub_mem = config['diffexp_heatmap']['mem']
    threads:
        config['diffexp_heatmap']['cpu']
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        outpath = "result/7.Diffexp/Heatmap",
        type = lambda w:diff_groups.loc["{}".format(w.diff_group)]["type"],
        treatment = lambda w:diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w:diff_groups.loc["{}".format(w.diff_group)]["control"],
        topn = topn,
        pval_or_padj = "padj" if (padj) else "pval",
        sgin_thread = padj if (padj) else pvalue,
        foldchange = foldchange,
    shell:
        '''
        {mp}
        {m20sc_sctool}
        sctool -i {params.inputfile} \
            -f h5seurat \
            -o  {params.outpath}\
            --assay SCT \
            --dataslot data,scale.data \
            -j {threads} \
            visualize --diffGene {input.difffile} \
                --topn {params.topn} \
                --vismethod diff_heatmap \
                --groupby {params.type} \
                -v {params.type} \
                -u {params.treatment},{params.control} >&2 2>{log} && \
        mv {params.outpath}/top{params.topn}_{wildcards.diff_group}_heatmap.pdf {output.figurepdf} && \
        mv {params.outpath}/top{params.topn}_{wildcards.diff_group}_heatmap.png {output.figurepng} && \
        mv {params.outpath}/top{params.topn}_{wildcards.diff_group}_genes.xls {output.table} >&2 2>>{log}
        '''

##================================================step3: GO and KEGG enrich=============================================
rule diff_gene_enrich:
    '''
    Add GO and KEGG enrich
    '''
    input:
        diff_DEG = "result/7.Diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.xls"
    output:
        done_check="logs/diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.done"
    log:
        "logs/diffexp/run.snakemake.enrich_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.log"
    resources:
        qsub_mem = config["diff_gene_enrich"]["mem"]
    params:
        reference_path = ref_dir,
        outdir = "result/7.Diffexp/enrichment",
        prefix="{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}"
    shell:
        '''(
        {mp}
        {enrich}
        enrichwrap -i {input.diff_DEG}  -g {params.reference_path}/annotation/  -o {params.outdir}/tmp/{params.prefix}   --fullname  --nomap  --minListHits 3 && \
        mkdir -p  {params.outdir}/KEGG_enrichment/{params.prefix}  {params.outdir}/GO_enrichment/{params.prefix}   && \
        mv {params.outdir}/tmp/{params.prefix}/1.GO_enrichment/{params.prefix}/*    {params.outdir}/tmp/{params.prefix}/1.GO_enrichment/enrichment_go.xls {params.outdir}/GO_enrichment/{params.prefix}/ && \
        mv {params.outdir}/tmp/{params.prefix}/2.KEGG_enrichment/{params.prefix}/*  {params.outdir}/tmp/{params.prefix}/2.KEGG_enrichment/enrichment_kegg.xls  {params.outdir}/KEGG_enrichment/{params.prefix}/ && \
        rm -rf  {params.outdir}/tmp/ && \
        touch {output.done_check}
        )>&2 2>{log}
        '''

##================================================step4: PPI ===========================================================
rule diff_gene_ppi:
    """
    add Annotation for diff expression genes
    """
    input:
        diff_DEG="result/7.Diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.xls"
    output:
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.pdf",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.png",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.svg",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.pdf",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.png",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.svg",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.tsv",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.mapping_results.tsv",
        "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.top_25_ppi_network.png"
    params:
        outdir="result/7.Diffexp/ppi",
        prefix="{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}" ,
        speciesid=config['params']['Diffexp']['speciesid']
    resources:
        mem=config['diff_gene_ppi']['mem'],
        cpu=config['diff_gene_ppi']['cpu']
    log:
        "logs/ppi/run.snakemake.{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.output"
    shell:
        """
        {mp}
        {oesinglecell_visium}
        python scripts/ppi/ppi.py  \
              -i {input.diff_DEG} \
              -s {params.speciesid} \
              -p {params.prefix} \
              -o {params.outdir} && \
        Rscript scripts/ppi/ppi_circle.r \
              -i {params.outdir}/{params.prefix}.ppi_network.tsv  \
              -d result/7.Diffexp/{params.prefix}.xls  \
              -o {params.outdir} >&2 2>{log}   
        """

rule diffexp_done:
    input:
        expand([
            "result/7.Diffexp/enrichment/1.GO_enrichment/enrichment_go.xls",
            "result/7.Diffexp/enrichment/2.KEGG_enrichment/enrichment_kegg.xls",
            "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_heatmap.pdf",
            "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_heatmap.png",
            "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_genes.xls",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.pdf",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.png",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.svg",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.pdf",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.png",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.svg",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.tsv",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.mapping_results.tsv",
            "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.top_25_ppi_network.png"
        ],
            topn = topn,
            diff_group = diff_groups.index.to_list(),
            pval_or_padj = "padj" if (padj) else "pval",
            sgin_thread = padj if (padj) else pvalue,
            foldchange = foldchange
        )
    output:
        donefile = "logs/Diffexp.done"
    shell:
        """
        touch {output.donefile}
        """  