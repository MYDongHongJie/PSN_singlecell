# 1. marker_findallmarkers =================================================================================================
rule marker_findallmarkers:
    """
    findallmarkers
    """
    input:
        donefile = "logs/Celltype.done"
    output:
        allmarker = "result/4.mRNA_marker/all_markers_for_each_cluster.xls",
        top10_marker = "result/4.mRNA_marker/top10_markers_for_each_cluster.xls"
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        outdir = "result/4.mRNA_marker",
        topn_marker = topn_marker,
        test_method = test_method,
    benchmark:
        "benchmarks/marker_findallmarkers.benchmark.txt"
    log:
        "logs/marker/marker_findallmarkers.log"
    resources:
        qsub_mem=config['marker_findallmarkers']['mem']
    threads:
        config['marker_findallmarkers']['cpu']
    shell:
        """
        {mp}
        {m20sc_sctool}
        sctool \
            -i {params.inputfile} \
            -f h5seurat \
            -o {params.outdir} \
            --assay SCT \
            --dataslot data,counts \
            -j {threads} \
            findallmarkers \
                --min_pct1 0.5 \
                --max_pct2 0.5 \
                --pct_fold 2 \
                --topn_marker {params.topn_marker} \
                --avg_log2FC 1  \
                --pvalue 0.05  \
                --strict F  \
                --cluster_name clusters \
                --test {params.test_method} >&2 2>{log}
        """

# 2. marker_visualize ====================================================================================================
rule marker_visualize:
    """
    Visualizing (vlnplot, featureplot) mRNA/lncRNA markers
    """
    input:
        allmarker = "result/4.mRNA_marker/all_markers_for_each_cluster.xls",
        top10_marker = "result/4.mRNA_marker/top10_markers_for_each_cluster.xls"
    output:
        all_featureplot = expand("result/4.mRNA_marker/markers_vis4cluster1/marker_gene_featureplot.{format}", format = figures),
        all_violinplot = expand("result/4.mRNA_marker/markers_vis4cluster1/marker_gene_violin_plot.{format}", format = figures),
        lncRNA_featureplot = expand("result/5.LncRNA_marker/markers_vis4cluster1/marker_gene_featureplot.{format}", format = figures),
        lncRNA_violinplot= expand("result/5.LncRNA_marker/markers_vis4cluster1/marker_gene_violin_plot.{format}",format=figures),
        donefile = "logs/Marker.done"
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        outdir = "result/4.mRNA_marker",
        lncRNA = "result/5.LncRNA_marker",
        reduct2 = reduct2,
        topn = RNA_topn,
        topby = topby,
        vismethod = vismethod,
        STAR_index= ref_dir
    benchmark:
        "benchmarks/marker_visualize.benchmark.txt"
    log:
        "logs/marker/marker_visualize.log"
    resources:
        qsub_mem=config['marker_visualize']['mem']
    threads:
        config['marker_visualize']['cpu']
    shell:
        """
        ({mp}
        {m20sc_sctool}
        sctool \
            -i {params.inputfile} \
            -f h5seurat \
            -o {params.outdir} \
            --assay SCT \
            --dataslot data,scale.data \
            -j {threads} \
            visualize \
                --markers {input.allmarker} \
                --selGene {params.STAR_index}/mRNA_genelist.gmt \
                --topn {params.topn} \
                --topby {params.topby} \
                --vismethod {params.vismethod} \
                --groupby clusters \
                --pointsize 0 \
                --reduct SCT_{params.reduct2} && \
        sctool \
            -i {params.inputfile} \
            -f h5seurat \
            -o {params.lncRNA} \
            --assay SCT \
            --dataslot data,scale.data \
            -j {threads} \
            visualize \
                --markers {input.allmarker} \
                --selGene {params.STAR_index}/lncRNA_genelist.gmt \
                --topn {params.topn} \
                --topby {params.topby} \
                --vismethod {params.vismethod} \
                --groupby clusters \
                --pointsize 0 \
                --reduct SCT_{params.reduct2}                       
        ) >&2 2>{log} && touch {output.donefile}    
        """

