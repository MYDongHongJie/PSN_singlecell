##==========================================step 1: RNA Dimensionality reduction/ clustring =======================================================

rule cluster:
    '''
    use sctool to RNA cluster
    '''
    input:
        donefile = "logs/QC.done"
    output:
        files = expand([
            "result/3.Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.{formats}"
            ],
                reduct1 = reduct1,
                reduct2 = reduct2,
                resolution = resolution,
                formats = figures
        ),
        donefile = "logs/sub_cluster.done"
    log:
        "logs/clustering/cluster.log"
    benchmark:
        "benchmarks/clustering/cluster.benchmark.txt"
    resources:
        qsub_mem = config['cluster']['mem']
    threads:
        config['cluster']['cpu']
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        outdir = "result/3.Clustering",
        reduct1 = reduct1,
        reduct2 = reduct2,
        clusteringuse = clusteringuse,
        resolution = resolution,
        palette = palette,
        component = component,
        batchid = batchid
    shell:
        '''
        {mp}
        {m20sc_sctool}
        if  [ {params.reduct1} == "harmony" ]; then sctool -i {params.inputfile} \
               -f h5seurat \
               --assay SCT \
               --dataslot data,scale.data  \
               --update TRUE \
               -o {params.outdir}\
               -d h5seurat \
                --ncores {threads} \
               bclust --reduct1 pca \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid {params.batchid} \
                      --component {params.component} >&2 2>{log} && touch {output.donefile}; 
        fi
        
        sctool -i {params.inputfile} \
               -f h5seurat \
               --assay SCT \
               --dataslot data,scale.data  \
               --update TRUE \
               -o {params.outdir} \
               -d h5seurat \
               --ncores {threads} \
               bclust --reduct1 {params.reduct1} \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid {params.batchid} \
                      --component {params.component} >&2 2>{log} &&  touch {output.donefile}
        '''
##==========================================step 2: RNA cluster visualization ==================================================
rule RNA_vis:
    '''
    use sctool to RNA cluster visualization
    '''
    input:
        donefile = "logs/sub_cluster.done"
    output:
        "result/3.Clustering/visualize_cluster_by_clusters/clust_cond_freq_info.xls"
    benchmark:
        "benchmarks/clustering/RNA_vis.benchmark.txt"
    log:
        "logs/clustering/RNA_vis.log"
    resources:
        qsub_mem=config['RNA_vis']['mem']
    threads:
        config['RNA_vis']['cpu']
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        reduct2 = reduct2,
        groupby = groupby
    shell:
        '''
        {mp}
        {m20sc_sctool}
        sctool -i {params.inputfile} \
               -f h5seurat \
               --assay SCT \
               --dataslot data  \
               -o result/3.Clustering \
               --ncores {threads} \
               summarize --reduct SCT_{params.reduct2} \
                         --palette customecol2 \
                         -c clusters\
                         -b {params.groupby} \
                         --ptsize 0.5 \
                         --dosummary T >&2 2>{log}
        '''

rule clustering_coefficient:
    """
    Visualizing coefficient among clsuters after Dimension reduction & clustering
    """
    input:
        donefile = "logs/sub_cluster.done"
    output:
        expand([
            "result/3.Clustering/clusters_correlation/normalized_data_groupby_clusters.xls",
            "result/3.Clustering/clusters_correlation/coefficient_heatmap.{formats}"
        ],
            formats = figures
        )
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        outdir = "result/3.Clustering/clusters_correlation",
        reduct2 = reduct2
    benchmark:
        "benchmarks/clustering_coefficient.benchmark.txt"
    log:
        "logs/clustering/clustering_coefficient.log"
    resources:
        qsub_mem = config['clustering_coefficient']['mem']
    threads:
        config['clustering_coefficient']['cpu']
    shell:
        """
        {mp}
        {m20sc_sctool}
        scVis \
            -i {params.inputfile} \
            -f h5seurat \
            -o {params.outdir} \
            -t {threads} \
            --assay SCT \
            --dataslots data \
            --reduct SCT_{params.reduct2} \
            coefficient \
                -g clusters >&2 2>{log}
        """

rule cluster_done:
    input:
        expand([
            "result/3.Clustering/clusters_correlation/normalized_data_groupby_clusters.xls",
            "result/3.Clustering/clusters_correlation/coefficient_heatmap.{formats}",
            "result/3.Clustering/visualize_cluster_by_clusters/clust_cond_freq_info.xls"
        ],
            formats = figures
        )
    output:
        donefile = "logs/Clustering.done"
    shell:
        """
        touch {output.donefile}
        """