##============================================= step1: Differential expression analysis =======================================================
rule diffexp:
    '''
    use sctool to run differential expression analysis with RNA assay
    '''
    input:
        clusterlog="logs/clustering/RNA_cluster.log",
        inputfile="result/2.Count_QC/seurat_ob.h5seurat"
    output:
        "result/7.Diffgene/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls"
    params:
        test=config['params']['Diffexp']["test"],
        foldchange=config['params']['Diffexp']["FC"],
        pvalue_or_fdr="--fdr " if (config['params']['Diffexp']['padj']) else "--pvalue ",
        pvalue=config['params']['Diffexp']['padj'] if (config['params']['Diffexp']['padj']) else
        config['params']['Diffexp']['pvalue'],
        treatment=lambda w: diff_groups.loc["{}".format(w.groupname)]["treatment"],
        control=lambda w: diff_groups.loc["{}".format(w.groupname)]["control"],
        type=lambda w: diff_groups.loc["{}".format(w.groupname)]["type"],
        qsub_mem=config['diffexp']['mem'],
        qsub_p=config['diffexp']['cpu']
    benchmark:
        "benchmarks/diffexp/diff_gene_for_{groupname}-{pval_or_padj}-{pvalue}-{foldchange}.benchmark.txt"
    log:
        "logs/diffexp/run.snakemake.{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.output"
    envmodules:
        config['envmodules']['oesinglecell']
    shell:
        '''(
        sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/7.Diffgene \
              --assay SCT \
              -d h5seurat \
              --dataslot data,counts  \
              -j {params.qsub_p} \
              diffexp --contrast {params.type}:{params.treatment}:{params.control} \
                     --FC {params.foldchange} \
                     {params.pvalue_or_fdr} {params.pvalue}  \
                     --test {params.test} \
                     -v {params.type} \
                     -u {params.treatment},{params.control}
        ) >{log} 2>&1'''

##==================================================step2: plot diff_heatmap ==============================================
rule diffexp_visualize:
    '''
    use sctool to plot diff_heatmap with RNA assay
    '''
    input:
        inputfile="result/2.Count_QC/seurat_ob.h5seurat",
        diff_DEG="result/7.Diffgene/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls"
    output:
        pdf="result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_heatmap.pdf",
        png="result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_heatmap.png",
        file="result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_genes.xls"
    benchmark:
        "benchmarks/diffexp/diffexp_RNA_heatmap_{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}.benchmark.txt"
    log:
        "logs/diffexp/diffexp_RNA_heatmap_{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        type=lambda w: diff_groups.loc["{}".format(w.groupname)]["type"],
        treatment=lambda w: diff_groups.loc["{}".format(w.groupname)]["treatment"],
        control=lambda w: diff_groups.loc["{}".format(w.groupname)]["control"],
        qsub_mem=config['diffexp_visualize']['mem'],
        qsub_p=config['diffexp_visualize']['cpu']
    shell:
        '''(
        sctool -i {input.inputfile} \
          -f h5seurat \
          -o result/7.Diffgene \
          --assay SCT \
          --dataslot data,scale.data \
          visualize --diffGene {input.diff_DEG} \
            --topn 25 \
            --vismethod diff_heatmap \
            --groupby {params.type} \
            -v {params.type} \
            -u {params.treatment},{params.control} && 
        mv result/7.Diffgene/top25_{params.type}_{params.treatment}-vs-{params.control}_heatmap.png {output.png} && 
        mv result/7.Diffgene/top25_{params.type}_{params.treatment}-vs-{params.control}_heatmap.pdf {output.pdf} && 
        mv result/7.Diffgene/top25_{params.type}_{params.treatment}-vs-{params.control}_genes.xls {output.file} 
        ) >{log} 2>&1'''

##==================================================step1: diffexp Enrich ==============================================
localrules: diffexp_enrich
rule diffexp_enrich:
    '''
    Diffexp gene enrichment
    '''
    input:
        log_annotation=expand(
            ["logs/diffexp/diffexp_RNA_heatmap_{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}.log",
             "result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_genes.xls"],
            groupname=diff_groups.index,
            pval_or_padj="padj" if (config['params']['Diffexp']['padj']) else "pval",
            pvalue=config['params']['Diffexp']['padj'] if (config['params']['Diffexp']['padj']) else
            config['params']['Diffexp']['pvalue'],
            foldchange=config['params']['Diffexp']['FC']),
        go_backfile="result/backfile/go.backgroud.xls",
        kegg_backfile="result/backfile/kegg.backgroud.xls",
        catagory_file="result/backfile/category.xls"
    output:
        "result/8.DiffgeneEnrich/GO_enrichment/enrichment_go.xls",
        "result/8.DiffgeneEnrich/KEGG_enrichment/enrichment_kegg.xls"
    benchmark:
        "benchmarks/diffexp/diff_gene_enrich.benchmark.txt"
    log:
        "logs/diffexp/diff_enrich.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        #reference=config['database']['ref_genome'],
        qsub_mem=config['diffexp_enrich']['mem'],
        qsub_p=config['diffexp_enrich']['cpu']
    shell:
        '''
        perl  scripts/anno/enrich/enrich_go_kegg_diff.pl -infile result/7.Diffgene/*-vs-*-diff-*.xls \
            -go_bg  {input.go_backfile} \
            -category {input.catagory_file}   \
            -kegg_bg {input.kegg_backfile} \
            -outdir result/8.DiffgeneEnrich  \
            -shelldir result/8.DiffgeneEnrich/enrichment_sh \
            -thread {params.qsub_p} \
            -queue big
        '''
