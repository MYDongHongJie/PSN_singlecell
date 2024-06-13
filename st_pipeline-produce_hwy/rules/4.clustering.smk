rule cluster_seurat:
    """
    Dimension reduction & clustering  
    """
    input:
         filtered_rds = "result/count_qc/filtered_seurat.rds"
    output:
        expand([
                "result/cluster_seurat/{reduct2_method}_Dimension_Reduction/{reduct2_method}_Dimension_Reduction_coordination.csv",
                "result/cluster_seurat/{reduct2_method}_Dimension_Reduction/{reduct2_method}_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/cluster_seurat/{reduct2_method}_Dimension_Reduction/{reduct2_method}_groupby_cluster_resolution{resolution}_plot.png",
               #"result/cluster_seurat/{reduct1_method}_Dimension_Reduction/{reduct1_method}_Dimension_Reduction_coordination.csv",
                "result/cluster_seurat/singlecell_object.clustering_resolution{resolution}.rds"
                ],
                        reduct2_method=config["params"]["reduct2_method"],
                        reduct1_method=config["params"]["reduct1_method"],
                        resolution=config["params"]["resolution"]
        )
    params:
            outdir="result/cluster_seurat",
           # reduct1_method=config["params"]["reduct1_method"],
            #lustering_method2=config["params"]["reduct2_method"],
           # batchid=config["params"]["batchid"] ,
            component= " "  if config["params"]["component"] is None  else f" --components  {config['params']['component']} ",
            resolution=config["params"]["resolution"]
            #spointsize=config['params']['spointsize']

    benchmark:
             "benchmarks/cluster_seurat.benchmark.txt"
    log:
        "logs/bclust/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        sctool \
          -i {input.filtered_rds} \
          -f rds \
          -o {params.outdir} \
          -d rds \
          --assay SCT \
          --image TRUE \
          --update FALSE \
          --prefix singlecell_object.clustering_resolution{params.resolution} \
          bclust \
              --reduct1 {config[params][reduct1_method]} \
              --reduct2 {config[params][reduct2_method]} \
              --batchid {config[params][batchid]} \
              {params.component} \
              --clusteringuse {config[params][clusteringuse]} \
              --resolution {params.resolution}  \
              --rerun  T >&2 2>{log}    
        """
#=======================================================================================================================
rule cluster_visualization:
    """
    Visualizing after Dimension reduction & clustering
    """
    input:
        cluster_rds =expand([
                            "result/cluster_seurat/singlecell_object.clustering_resolution{resolution}.rds"],
                            resolution=config["params"]["resolution"]
                            )

    output:
        expand([
            "result/cluster_seurat/visualize_cluster_by_clusters/clust_cond_freq_info.xls",
            "result/cluster_seurat/visualize_cluster_by_clusters/groupby-group_resolution{resolution}_contrast_plot.pdf",
            "result/cluster_seurat/visualize_cluster_by_clusters/groupby-group_resolution{resolution}_contrast_plot.png",
            "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution{resolution}_contrast_plot.pdf",
            "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution{resolution}_contrast_plot.png",
            "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution-{resolution}_summary_plot.pdf",
            "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution-{resolution}_summary_plot.png",
            "result/cluster_seurat/visualize_cluster_by_clusters/splitby-group_resolution{resolution}_split_plot.pdf",
            "result/cluster_seurat/visualize_cluster_by_clusters/splitby-group_resolution{resolution}_split_plot.png",
            "result/cluster_seurat/visualize_cluster_by_clusters/splitby-sampleid_resolution{resolution}_split_plot.pdf",
            "result/cluster_seurat/visualize_cluster_by_clusters/splitby-sampleid_resolution{resolution}_split_plot.png"
        ],
            resolution=config["params"]["resolution"]
        )
    params:
           outdir="result/cluster_seurat",
           reduct2_method=config["params"]["reduct2_method"],
           assay= "SCT" if config["params"]["normmeth"] == "sctransform" else "RNA",
           reduct = f"SCT_umap" if config["params"]["normmeth"] == "sctransform" else "RNA_umap",
           resolution=config["params"]["resolution"],
           crop= "TRUE" if config["params"]["library_type"] == "cytassist" else "FALSE"
    log:
        "logs/dimplot/run.snakemake.output"
    benchmark:
             "benchmarks/cluster_visualization.benchmark.txt"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        scVis \
              -i {input.cluster_rds} \
              -f rds \
              -o {params.outdir} \
              --reduct {params.reduct}  \
              --assay {params.assay}  \
              --image TRUE \
              dimplot \
                  --resolution {params.resolution} \
                  --ptsize 0.5 \
                  --spointsize  1.2 \
                  --groupby clusters \
                  --splitby sampleid,group \
                  --crop {params.crop} \
                  --common.legend FALSE  >&2 2>{log}    
        """
