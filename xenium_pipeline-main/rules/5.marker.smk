##=========================================step1: marker基因鉴定 ========================================================
rule marker:
    """
    find marker
    """
    input:
        inputfile=expand(["result/3.Clustering/singlecell_object.clustering_resolution{resolution}.rds"],resolution=resolution)
    output:
        "result/5.Marker/all_markers_for_each_cluster_anno.xls",
        "result/5.Marker/top10_markers_for_each_cluster_anno.xls"
    benchmark:
        "benchmarks/marker/marker.txt"
    log:
        "logs/marker/marker.log"
    params:
        gene_anno=config['database']['gene_annotation'],
        test_method=config['params']['Marker']['test_method'],
        topn= config['params']['Marker']['topn']
    shell:
        '''
        sctool -i {input.inputfile} \
                -f rds \
                -o result/5.Marker \
                --assay SCT \
                findallmarkers -e {params.test_method} \
                        -n clusters \
                        --anno {params.gene_anno} \
                        -t 0.5 \
                        -T 0.5 \
                        -c 2 \
                        -N {params.topn} \
                        -k 1 \
                        -p 0.05 \
                        -s F    >&2 2>{log}
        '''

##=========================================step2: marker基因可视化 ======================================================
rule marker_vis:
    """
    Marker Visualization
    """
    input:
        inputfile=expand(["result/3.Clustering/singlecell_object.clustering_resolution{resolution}.rds"],resolution=resolution),
        top10_marker = "result/5.Marker/top10_markers_for_each_cluster_anno.xls"
    output:
        "result/5.Marker/topmarker_gene_heatmap.png"
    benchmark:
        "benchmarks/marker/marker_vis.txt"
    log:
        "logs/marker/marker_vis.log"
    params:
        vismethod = config['params']['Marker']['vismethod'],
        topby = config['params']['Marker']['topby'],
        topn_vis = config['params']['Marker']['topn_vis'],
        reduct=reduct2,
        location = config['params']['Marker']['location']
    shell:
        '''
        scVis -i {input.inputfile} \
                -o result/5.Marker \
                --assay SCT \
                markervis -l {input.top10_marker} \
                        -n {params.topn_vis} \
                        -c {params.topby} \
                        -m {params.vismethod} \
                        -s 0 \
                        --reduct SCT_{params.reduct} \
                        --loc {params.location} \
                        --crop TRUE  >&2 2>{log}
        '''
