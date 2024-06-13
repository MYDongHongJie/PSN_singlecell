##==================================================step1: find marker =======================================================
rule marker_findallmarkers:
    '''
    use sctool to find RNA marker
    '''
    input:
        cluster="logs/clustering/RNA_cluster.log",
        inputfile="result/2.Count_QC/seurat_ob.h5seurat",
        kegglog="logs/anno/kegg.log",
        swissprotlog="logs/anno/swissprot.log",
        golog="logs/anno/go.log",
        cardlog="logs/anno/card.log",
        cazlog="logs/anno/cazy.log"
    output:
        "result/5.Marker/all_markers_for_each_cluster.xls",
        "result/5.Marker/top10_markers_for_each_cluster.xls"
    benchmark:
        "benchmarks/marker/RNA_findmarker.benchmark.txt"
    log:
        "logs/marker/RNA_findmarker.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        topn_marker=config["params"]['Marker']["RNA_topn_marker"],
        test_method=config["params"]['Marker']["RNA_test_method"],
        reduct2=config['params']['Clustering']["RNA_reduct2"],
        resolution=config['params']['Clustering']["RNA_resolution"],
        qsub_mem=config['marker_findallmarkers']['mem'],
        qsub_p=config['marker_findallmarkers']['cpu']
    shell:
        '''(
        sctool -i {input.inputfile} \
               -f h5seurat \
               -o result/5.Marker \
               --assay SCT \
               --dataslot data,counts \
               -j {params.qsub_p} \
               findallmarkers --min_pct1 0.5 \
                              --max_pct2 0.5  \
                              --pct_fold 2 \
                              --topn_marker {params.topn_marker} \
                              --avg_log2FC 1 \
                              --pvalue 0.05 \
                              --strict F \
                              --cluster_name  SCT.{params.reduct2}.res.{params.resolution}\
                              --test {params.test_method}
        ) >{log} 2>&1'''

##==================================================step 2: marker visualization=============================================
rule marker_visualize:
    '''
    use sctool to RNA marker visualization
    '''
    input:
        log_marker="logs/marker/RNA_findmarker.log",
        inputfile="result/2.Count_QC/seurat_ob.h5seurat"
    output:
        "result/5.Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
        "result/5.Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf"
    benchmark:
        "benchmarks/marker/RNA_markervis.benchmark.txt"
    log:
        "logs/marker/RNA_markervis.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        topn=config["params"]['Marker']["RNA_topn"],
        topby=config["params"]['Marker']["RNA_topby"],
        vismethod=config["params"]['Marker']["RNA_vismethod"],
        reduct2=config['params']['Clustering']["RNA_reduct2"],
        resolution=config['params']['Clustering']["RNA_resolution"],
        qsub_mem=config['marker_visualize']['mem'],
        qsub_p=config['marker_visualize']['cpu']
    shell:
        '''(
        sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/5.Marker \
              --assay SCT \
              --dataslot data,scale.data \
              visualize --markers result/5.Marker/top10_markers_for_each_cluster.xls \
                        --topn {params.topn} \
                        --topby {params.topby} \
                        --vismethod {params.vismethod} \
                        --groupby SCT.{params.reduct2}.res.{params.resolution} \
                        --pointsize 0 \
                        --reduct SCT_{params.reduct2}
        ) >{log} 2>&1'''

##==================================================step1: Marker Enrich praper==========================================
localrules: marker_enrich_praper
rule marker_enrich_praper:
    '''
    Marker enrichment praper file
    '''
    input:
        markervis_log="logs/marker/RNA_markervis.log",
        go_anno="result/4.Annotation/GO/Unigene.GO.gene.anno.xls",
        kegg_anno="result/4.Annotation/KEGG/Unigene.KEGG.gene.anno.xls",
    output:
        "result/backfile/go.backgroud.xls",
        "result/backfile/kegg.backgroud.xls",
        "result/backfile/category.xls"
    benchmark:
        "benchmarks/marker/marker_enrich_praper.benchmark.txt"
    log:
        "logs/marker/marker_enrich_praper.log"
    envmodules:
        config['envmodules']['oesinglecell']
    shell:
        """
        perl scripts/anno/backgroud_prepare.pl -g {input.go_anno} \
                                            -k {input.kegg_anno} \
                                            -o result/backfile
        sed -i "s/_/-/g\" result/backfile/go.backgroud.xls
        sed -i "s/_/-/g\" result/backfile/kegg.backgroud.xls
        """

##==================================================step2: Marker Enrich ===============================================
localrules: marker_enrich
rule marker_enrich:
    '''
    Marker enrichment
    '''
    input:
        praper_file="logs/marker/marker_enrich_praper.log",
        go_backfile="result/backfile/go.backgroud.xls",
        kegg_backfile="result/backfile/kegg.backgroud.xls",
        cate_file="result/backfile/category.xls",
        marker_file="result/5.Marker/all_markers_for_each_cluster.xls"
    output:
        "result/6.MarkerEnrich/GO_enrichment/enrichment_go.xls",
        "result/6.MarkerEnrich/KEGG_enrichment/enrichment_kegg.xls"
    benchmark:
        "benchmarks/marker/marker_enrich.benchmark.txt"
    log:
        "logs/marker/marker_enrich.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        topn=config['params']['Marker']['RNA_topn'],
        qsub_p=config['marker_enrich']['cpu']
    shell:
        """
        Rscript scripts/anno/filter_markers.R -i {input.marker_file} -o result/6.MarkerEnrich -n {params.topn} -c gene_diff -s T && 
        perl scripts/anno/enrich/enrich_go_kegg.pl -infile result/6.MarkerEnrich/*markers_for_cluster*.xls \
                                    -go_bg {input.go_backfile} \
                                    -category {input.cate_file} \
                                    -kegg_bg {input.kegg_backfile} \
                                    -outdir result/6.MarkerEnrich \
                                    -shelldir result/6.MarkerEnrich/enrichment_sh \
                                    -thread {params.qsub_p} \
                                    -queue big
        """
