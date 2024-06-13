##==========================================step 1: RNA Dimensionality reduction/ clustring =======================================================
rule clustering_bclust:
    '''
    use sctool to RNA cluster
    '''
    input:
        inputfile="result/2.Count_QC/seurat_ob.h5seurat"
    output:
        expand([
            "result/3.Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
            "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf"
        ],
            reduct1=config['params']['Clustering']["RNA_reduct1"],
            reduct2=config['params']['Clustering']["RNA_reduct2"],
            resolution=config['params']['Clustering']["RNA_resolution"]
        )
    benchmark:
        "benchmarks/clustering/RNA_cluster.benchmark.txt"
    log:
        "logs/clustering/RNA_cluster.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        reduct1=config['params']['Clustering']["RNA_reduct1"],
        reduct2=config['params']['Clustering']["RNA_reduct2"],
        clusteringuse=config['params']['Clustering']["RNA_clusteringuse"],
        resolution=config['params']['Clustering']["RNA_resolution"],
        palette=config['params']['Clustering']["RNA_palette"],
        component=config['params']['Clustering']["RNA_component"],
        qsub_mem=config['clustering_bclust']['mem'],
        qsub_p=config['clustering_bclust']['cpu']
    shell:
        '''(
        if  [ {params.reduct1} == "harmony" ];then sctool -i {input.inputfile} \
               -f h5seurat \
               --assay SCT \
               --dataslot data,scale.data  \
               --update TRUE \
               -o result/3.Clustering\
               -d h5seurat \
               bclust --reduct1 pca \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid batchid \
                      --component {params.component} ; fi
        sctool -i {input.inputfile} \
               -f h5seurat \
               --assay SCT \
               --dataslot data,scale.data  \
               --update TRUE \
               -o result/3.Clustering\
               -d h5seurat \
               bclust --reduct1 {params.reduct1} \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid batchid \
                      --component {params.component}   
        ) >{log} 2>&1'''

##==========================================step 2: RNA cluster visualization ==================================================
rule clustering_summarize:
    '''
    use sctool to RNA cluster visualization
    '''
    input:
        inputfile="result/2.Count_QC/seurat_ob.h5seurat",
        arg="logs/clustering/RNA_cluster.log"
    output:
        expand(["result/3.Clustering/visualize_cluster_by_SCT.{reduct2}.res.{resolution}/clust_cond_freq_info.xls"],
            reduct2=config['params']['Clustering']["RNA_reduct2"],
            resolution=config['params']['Clustering']["RNA_resolution"])
    benchmark:
        "benchmarks/clustering/RNA_vis.benchmark.txt"
    log:
        "logs/clustering/RNA_vis.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        reduct2=config['params']['Clustering']["RNA_reduct2"],
        resolution=config['params']['Clustering']["RNA_resolution"],
        groupby=config['params']['Clustering']["RNA_groupby"],
        qsub_mem=config['clustering_summarize']['mem'],
        qsub_p=config['clustering_summarize']['cpu']
    shell:
        '''
        sctool -i {input.inputfile} \
               -f h5seurat \
               --assay SCT \
               --dataslot data,counts  \
               --update  TRUE \
               -o result/3.Clustering \
               -d h5seurat \
               summarize --reduct SCT_{params.reduct2} \
                         --palette blindless \
                         -c SCT.{params.reduct2}.res.{params.resolution}\
                         -b {params.groupby} \
                         --ptsize 0.5 \
                         --dosummary T   >&2  2>{log} 
        '''
