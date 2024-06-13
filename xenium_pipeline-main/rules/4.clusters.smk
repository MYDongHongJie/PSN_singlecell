##=========================================step1: 降维聚类 ==============================================================
rule cluster:
    """
    Dimensionality reduction
    """
    input:
        inputfile = "result/2.Count_QC/filtered_seurat.rds"
    output:
        expand([
            "result/3.Clustering/singlecell_object.clustering_resolution{resolution}.rds",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/clusters_infor.csv",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png"
        ],resolution=resolution,reduct2=reduct2)
    benchmark:
        "benchmarks/cluster/cluster.txt"
    log:
        "logs/cluster/cluster.log"
    params:
        reduct1=reduct1,
        reduct2=reduct2,
        component=config['params']['Clustering']['component'],
        clusteringuse=config['params']['Clustering']['clusteringuse'],
        resolution=resolution,
        location = config['params']['Clustering']['location'],
        palette = config['params']['Clustering']['palette'],
        batchid = config['params']['Clustering']['batchid']
    shell:
        '''
        sctool -i {input.inputfile} \
               -f rds \
               --assay SCT \
               --update FALSE \
               -o result/3.Clustering \
               -d rds \
               --prefix singlecell_object.clustering_resolution{params.resolution} \
               bclust --reduct1 {params.reduct1} \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid {params.batchid} \
                      --component {params.component} \
                      --loc {params.location} >&2 2>{log}
        '''

##=========================================step2: 降维聚类可视化 =========================================================
rule cluster_vis:
    """
    Dimensionality reduction Visualization
    """
    input:
        inputfile=expand(["result/3.Clustering/singlecell_object.clustering_resolution{resolution}.rds"],resolution=resolution)
    output:
        expand([
            "result/3.Clustering/visualize_clusters_by_clusters/groupby-sampleid_resolution{resolution}_contrast_plot.png",
            "result/3.Clustering/visualize_clusters_by_clusters/groupby-sampleid_resolution-{resolution}_summary_plot.png",
            "result/3.Clustering/visualize_clusters_by_clusters/splitby-sampleid_resolution{resolution}_split_plot.png"
        ],resolution=resolution)
    benchmark:
        "benchmarks/cluster/cluster_vis.txt"
    log:
        "logs/cluster/cluster_vis.log"
    params:
        reduct2=reduct2,
        resolution=resolution,
        groupby=config['params']['Clustering']['groupby'],
        splitby=config['params']['Clustering']['splitby'],
        groups=config['params']['Clustering']['groups'],
        propby=config['params']['Clustering']['propby'],
        ptsize=config['params']['Clustering']['ptsize'],
        location= config['params']['Clustering']['location']
    shell:
        '''
        scVis -i {input.inputfile} \
                -f rds \
                -o result/3.Clustering \
                --reduct SCT_{params.reduct2} \
                --assay SCT \
                dimplot --resolution {params.resolution} \
                        --ptsize {params.ptsize} \
                        --groupby {params.groupby} \
                        --splitby {params.splitby} \
                        --groups {params.groups} \
                        --propby {params.propby} \
                        --crop TRUE \
                        --loc {params.location} >&2 2>{log}
        '''
