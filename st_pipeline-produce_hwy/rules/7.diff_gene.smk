rule  diff_gene:
    """
    diff expressed gene select
    """
    input:
         cluster_rds = "result/marker/seurat_object_markers.rds"
    output:
          "result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls",
          "result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}_anno.xls"
    params:
           outdir="result/diffexp/",
           foldchange=config['params']["foldchange"],
           pvalue_or_fdr="--fdr " if (config["params"]['qvalue']) else "--pvalue ",
           pvalue=config["params"]['qvalue'] if(config["params"]['qvalue']) else config["params"]['pvalue'],
           treatment=lambda w: diff_groups.loc["{}".format(w.groupname)]["treatment"],
           control=lambda w: diff_groups.loc["{}".format(w.groupname)]["control"],
           type=lambda w: diff_groups.loc["{}".format(w.groupname)]["type"],
           gene_annotation= os.path.join(config['database']['reference'],"annotation/gene_annotation.xls")
    # benchmark:
    #         "benchmarks/diff_gene_for_{groupname}-{pval_or_padj}-{pvalue}-{foldchange}.benchmark.txt"
    log:
        "logs/diffexp/run.snakemake.{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        sctool \
             -i {input.cluster_rds} \
             -f rds \
             -o {params.outdir} \
             --assay SCT \
             --image FALSE \
             diffexp \
                 --design ~{params.type} \
                 --contrast {params.type}:{params.treatment}:{params.control} \
                 --FC {params.foldchange} \
                 {params.pvalue_or_fdr} {params.pvalue} \
                 -e wilcox \
                 --anno {params.gene_annotation}   >&2 2> {log}   
         """

#========================================================================================================================
if config['params']['foldchange'] == 1.5:
    rule diff_gene_enrich:
        """
        add Annotation for diff expression genes
        """
        input:
            #diff_DEG = "result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls"
            diff_DEG = "result/diffexp/{groupname}-diff-pvalue-0.05-FC-1.5.xls"
        output:
            #mark="result/diffexp/enrichment/KEGG_enrichment/{groupname}/enrichment-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.log",
            go="result/diffexp/enrichment/GO_enrichment/{groupname}/enrichment-go-{groupname}-Total.xls",
            kegg="result/diffexp/enrichment/KEGG_enrichment/{groupname}/enrichment-kegg-{groupname}-Total.xls"
        params:
            outdir="result/diffexp/enrichemnt/{groupname}",
            prefix="{groupname}"

        # benchmark:
        #          "benchmarks/diff_gene_enrich_for_{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}.benchmark.txt"
        log:
            "logs/diffexp/run.snakemake.enrich.{groupname}.output"
        # envmodules:
        #     envmodules["oesinglecell"]
        shell:
            """
            module purge && module load OESingleCell/3.0.d 
            mkdir -p logs/diff_gene_enrich/{params.prefix}  &&  
            perl  scripts/enrichment/enrich_go_kegg.pl  \
                -infile {input.diff_DEG} \
                -go_bg {reference}/annotation/gene_go.backgroud.xls  \
                -kegg_bg {reference}/annotation/gene_kegg.backgroud.xls  \
                -category scripts/enrichment/category.xls  \
                -outdir  result/diffexp/enrichment \
                -shelldir logs/diff_gene_enrich/{params.prefix}    \
                -thread 1  \
                -queue big  >&2 2>{log}
            """
else:
    rule diff_gene_enrich:
        """
        add Annotation for diff expression genes
        """
        input:
            #diff_DEG = "result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls"
            diff_DEG = "result/diffexp/{groupname}-diff-pvalue-0.05-FC-1.2.xls"
        output:
            #mark="result/diffexp/enrichment/KEGG_enrichment/{groupname}/enrichment-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.log",
            go="result/diffexp/enrichment/GO_enrichment/{groupname}/enrichment-go-{groupname}-Total.xls",
            kegg="result/diffexp/enrichment/KEGG_enrichment/{groupname}/enrichment-kegg-{groupname}-Total.xls"
        params:
            outdir="result/diffexp/enrichemnt/{groupname}",
            prefix="{groupname}"

        # benchmark:
        #          "benchmarks/diff_gene_enrich_for_{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}.benchmark.txt"
        log:
            "logs/diffexp/run.snakemake.enrich.{groupname}.output"
        # envmodules:
        #     envmodules["oesinglecell"]
        shell:
            """
            module purge && module load OESingleCell/3.0.d 
            mkdir -p logs/diff_gene_enrich/{params.prefix}  &&  
            perl  scripts/enrichment/enrich_go_kegg.pl  \
                -infile {input.diff_DEG} \
                -go_bg {reference}/annotation/gene_go.backgroud.xls  \
                -kegg_bg {reference}/annotation/gene_kegg.backgroud.xls  \
                -category scripts/enrichment/category.xls  \
                -outdir  result/diffexp/enrichment \
                -shelldir logs/diff_gene_enrich/{params.prefix}    \
                -thread 1  \
                -queue big  >&2 2>{log}
            """

#========================================================================================================================
if config["module"]["diff_gene_ppi"]:
    rule diff_gene_ppi:
        """
        add Annotation for diff expression genes
        """
        input:
            diff_DEG="result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls"
        output:
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.new_colors.pdf",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.new_colors.png",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.new_colors.svg",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.pdf",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.png",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.svg",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.tsv",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.mapping_results.tsv",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.top_25_ppi_network.png",
             "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.top_25_ppi_network.pdf"

        
        params:
               outdir="result/diffexp/ppi",
               prefix="{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}",
               speciesid=config['params']["speciesid"],
               network_tsv="result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.ppi_network.tsv"
        # benchmark:
        #         "benchmarks/diff_gene_ppi_for_{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.benchmark.txt"
        log:
            "logs/ppi/run.snakemake.{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.output"
        envmodules:
            envmodules["oesinglecell"]
        shell:
            """ 
                python scripts/ppi/ppi.py \
                    -i {input.diff_DEG} \
                    -p {params.prefix} \
                    -o {params.outdir} \
                    -s {params.speciesid} >&2 2>{log} 

                Rscript scripts/ppi/ppi_circle.r \
                    -i {params.network_tsv} \
                    -d {input.diff_DEG} \
                    -o {params.outdir} 
            """




#========================================================================================================================


