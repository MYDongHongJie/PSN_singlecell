# 1. marker_findallmarkers =================================================================================================
rule marker_findallmarkers:
    """
    findallmarkers
    """
    input:
         cluster_rds = "result/Count_QC/filtered.h5seurat",
         temp_file="clustering_h5seurat_touched.check"
    output:
        "result/Marker/all_markers_for_each_cluster.xls",
        "result/Marker/top10_markers_for_each_cluster.xls"
    params:
        outdir="result/Marker",
        marker_test=config['report_params']['marker_test']
    benchmark:
        "benchmarks/marker_findallmarkers.benchmark.txt"
#    log:  "logs/marker/marker_findallmarkers.log"
    resources:
        qsub_mem=50,
        qsub_p=config['cpu']['marker_findallmarkers']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
rm -r {params.outdir} &&
Rscript {script_dir}/sctool \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  --assay RNA \
  --dataslot data,counts \
  -j 10 \
  findallmarkers \
  -c 2 \
  -N 10 \
  -k 1 \
  -p 0.05 \
  -s F \
  -e {params.marker_test} \
  -n clusters 
        """
# 2. marker_heatmap ======================================================================================================
rule marker_heatmap:
    """
    Visualizing (heatmap) markers
    """
    input:
        cluster_rds = "result/Count_QC/filtered.h5seurat",
        temp_file="clustering_h5seurat_touched.check",
        top10_marker = "result/Marker/top10_markers_for_each_cluster.xls"
    output:
        "result/Marker/topmarker_gene_heatmap.png",
        "result/Marker/topmarker_gene_heatmap.pdf"
    params:
        outdir="result/Marker/",
    benchmark:
        "benchmarks/marker_heatmap.benchmark.txt"
    #log:  "logs/marker/marker_heatmap.log"
    resources:
        qsub_mem=50,
        qsub_p=config['cpu']['marker_heatmap']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript {script_dir}/scVis \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  -t 10 \
  --assay RNA \
  --slot data,scale.data \
  heatmap \
  -l {input.top10_marker} \
  -c gene_diff \
  -n 10 \
  -g clusters \
  --group_colors customecol2 \
  --sample_ratio 0.8 \
  --style seurat 
        """
# 3. marker_visualize ====================================================================================================
rule marker_visualize:
    """
    Visualizing (vlnplot, featureplot) markers
    """
    input:
         cluster_rds = "result/Count_QC/filtered.h5seurat",
         temp_file="clustering_h5seurat_touched.check",
         top10_marker = "result/Marker/top10_markers_for_each_cluster.xls"
    output:
        "result/Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
        "result/Marker/markers_vis4cluster1/marker_gene_featureplot.png",
        "result/Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
        "result/Marker/markers_vis4cluster1/marker_gene_violin_plot.png",
    params:
        outdir="result/Marker/",
        reduct2=config["report_params"]["reduct2"],
        pointsize_params = get_pointsize_params,
    benchmark:
        "benchmarks/marker_visualize.benchmark.txt"
    #log:  "logs/marker/marker_visualize.log"
    resources:
        qsub_mem=30,
        qsub_p=config['cpu']['marker_visualize']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript {script_dir}/sctool \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  -j 10 \
  --assay RNA \
  --dataslot data \
  visualize \
  -l {input.top10_marker} \
  -g clusters \
  --reduct {params.reduct2} \
  --topn  10  \
  --topby avg_log2FC \
  -m vlnplot,featureplot \
  --vcolors customecol2 \
  --ccolors spectral \
  {params.pointsize_params} \
  --dodge F
        """
# 4. marker_annotation ===================================================================================================
#localrules: marker_annotation
rule marker_annotation:
    """
    add Annotation for marker genes
    """
    input:
        all_DEG="result/Marker/all_markers_for_each_cluster.xls",
        top10_DEG="result/Marker/top10_markers_for_each_cluster.xls",
        gene_annotation = os.path.join(config['report']['scrna_Reference'],"annotation/gene_annotation.xls")
    output:
        "result/Marker/all_markers_for_each_cluster_anno.xls",
        "result/Marker/top10_markers_for_each_cluster_anno.xls"
    params:
        outdir="result/Marker",
    benchmark:
        "benchmarks/marker_annotation.benchmark.txt"
    #log:  "logs/marker/marker_annotation.log"
    resources:
        qsub_mem=10,
        qsub_p=1
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript {script_dir}/sctool annotation \
  -g {input.all_DEG} \
  --anno {input.gene_annotation}
Rscript {script_dir}/sctool annotation \
  -g {input.top10_DEG} \
  --anno {input.gene_annotation}
        """
