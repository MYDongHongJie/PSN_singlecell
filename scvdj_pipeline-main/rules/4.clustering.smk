reduct1 = config["report_params"]["reduct1"]
if reduct1=="mnn":
    reduct1_params = " --batchid batchid --components 10 "
else:
    reduct1_params = ""

rule clustering_bclust:
    """
    Dimension reduction & clustering  
    """
    input:
         filtered_rds = "result/Count_QC/filtered.h5seurat"
    output:
        expand([
                "result/Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
                "result/Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
                "result/Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv"
                ],
                        reduct2=config["report_params"]["reduct2"],
                        reduct1=reduct1,
                        resolution=config["report_params"]["resolution"]
        ),
        temp_file="clustering_h5seurat_touched.check"
    params:
            outdir="result/Clustering",
            reduct1=reduct1,
            reduct1_params=reduct1_params,
            reduct2=config["report_params"]["reduct2"],
            pointsize_params = get_pointsize_params ,
            resolution=config["report_params"]["resolution"]
    benchmark:
             "benchmarks/clustering_bclust.benchmark.txt"
    #log:  "logs/clustering/clustering_bclust.log"
    resources:
             qsub_mem=30,
             qsub_p=config['cpu']['clustering_bclust']
    envmodules:
                config['report_params']['envmodules']
    shell:
        """
rm -r {params.outdir} &&
Rscript  {script_dir}/sctool \
  -i {input.filtered_rds} \
  -f h5seurat \
  -o {params.outdir} \
  -d h5seurat \
  --assay RNA \
  --dataslot counts,data,scale.data \
  --update T \
  bclust \
  --reduct1 {params.reduct1}  {params.reduct1_params} \
  --reduct2 {params.reduct2} \
  --clusteringuse snn \
  --resolution {params.resolution} \
  --rerun T \
  {params.pointsize_params} \
  --palette customecol2  && 
touch -r benchmarks/qc.benchmark.txt  result/Count_QC/filtered.h5seurat && 
touch clustering_h5seurat_touched.check && 
touch clustering_h5seurat_touched.`date +%Y%m%d-%H%M`
        """
#=======================================================================================================================
rule clustering_summarize:
    """
    Visualizing after Dimension reduction & clustering
    """
    input:
         temp_file="clustering_h5seurat_touched.check",
         cluster_rds = "result/Count_QC/filtered.h5seurat"
    output:
        "result/Clustering/visualize_cluster_by_clusters/clust_cond_freq_info.xls",
        "result/Clustering/visualize_cluster_by_clusters/groupby-clusters_summary_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/groupby-clusters_summary_plot.png",
        "result/Clustering/visualize_cluster_by_clusters/groupby-group_contrast_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/groupby-group_contrast_plot.png",
        "result/Clustering/visualize_cluster_by_clusters/groupby-group_summary_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/groupby-group_summary_plot.png",
        "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_contrast_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_contrast_plot.png",
        "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_summary_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_summary_plot.png",
        "result/Clustering/visualize_cluster_by_clusters/splitby-group_split_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/splitby-group_split_plot.png",
        "result/Clustering/visualize_cluster_by_clusters/splitby-sampleid_split_plot.pdf",
        "result/Clustering/visualize_cluster_by_clusters/splitby-sampleid_split_plot.png"
    params:
           outdir="result/Clustering",
           reduct2=config["report_params"]["reduct2"],
           pointsize_params = get_pointsize_params,
    benchmark:
             "benchmarks/clustering_summarize.benchmark.txt"
    # log:  "logs/clustering/clustering_summarize.log"
    resources:
             qsub_mem=50,
             qsub_p=config['cpu']['clustering_summarize']
    envmodules:
                config['report_params']['envmodules']
    shell:
        """
Rscript  {script_dir}/sctool \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  --assay RNA \
  --dataslot data \
  summarize \
  --reduct {params.reduct2} \
  --palette customecol2 \
  -c clusters \
  -b sampleid,group \
  {params.pointsize_params} \
  --dosummary T  
        """
#==============================================================================================
rule clustering_coefficient:
    """
    Visualizing coefficient among clsuters after Dimension reduction & clustering
    """
    input:
         temp_file="clustering_h5seurat_touched.check",
         cluster_rds = "result/Count_QC/filtered.h5seurat"
    output:
        "result/Clustering/clusters_correlation/normalized_data_groupby_clusters.xls",
        "result/Clustering/clusters_correlation/coefficient_heatmap.png",
        "result/Clustering/clusters_correlation/coefficient_heatmap.pdf",
    params:
        outdir="result/Clustering/clusters_correlation/",
        reduct2=config["report_params"]["reduct2"],
    benchmark:
        "benchmarks/clustering_coefficient.benchmark.txt"
    #log:  "logs/clustering/clustering_coefficient.log"
    resources:
        qsub_mem=30,
        qsub_p=config['cpu']['clustering_coefficient']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript  {script_dir}/scVis \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  -t 6 \
  --assay RNA \
  --slot data \
  --reduct {params.reduct2} \
  coefficient \
  -g clusters  
        """
