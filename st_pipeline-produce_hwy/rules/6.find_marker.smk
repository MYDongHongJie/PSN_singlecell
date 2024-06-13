#================================================find_marker============================================================
rule find_marker:
    """
    find_marker
    """
    input:
        cluster_rds=expand([
                        "result/cluster_seurat/singlecell_object.clustering_resolution{resolution}.rds"],
                        resolution=config["params"]["resolution"]
                        )
    output:
        "result/marker/all_markers_for_each_cluster_anno.xls",
        "result/marker/top10_markers_for_each_cluster_anno.xls",
        "result/marker/seurat_object_markers.rds"
    params:
           outdir="result/marker",
           gene_annotation = os.path.join(config['database']['reference'],"annotation/gene_annotation.xls"),
           resolution=config["params"]["resolution"]
    benchmark:
             "benchmarks/find_marker.benchmark.txt"
    log:
        "logs/findallmarkers/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        sctool \
          -i {input.cluster_rds} \
          -f rds \
          -o {params.outdir}\
          --image TRUE \
          --assay SCT \
          findallmarkers \
              -e bimod \
              -n clusters \
              --anno {params.gene_annotation} \
              -t 0.5 \
              -T 0.5 \
              -c 2 \
              -N 10 \
              -k 1 \
              -p 0.05 \
              -s F   >&2 2>{log}    

        """
#==========================================Visualizing markers==========================================================
rule visualizing_markers:
    """
    Visualizing markers
    """
    input:
        # cluster_rds=expand([
        #     "result/cluster_seurat/singlecell_object.clustering_resolution{resolution}.rds"],
        #     resolution=config["params"]["resolution"]
        # ),
        cluster_rds = "result/marker/seurat_object_markers.rds",
        top10_marker = "result/marker/top10_markers_for_each_cluster_anno.xls"
    output:
        "result/marker/topmarker_gene_heatmap.png",
        "result/marker/topmarker_gene_heatmap.pdf"
    params:
           outdir="result/marker/",
           assay = "SCT" if  config["params"]["normmeth"]=="sctransform" else "RNA",
           reduct = f"SCT_umap" if  config["params"]["normmeth"]=="sctransform" else "RNA_umap",
           crop= "TRUE" if config['params']['library_type'] == "cytassist" else "FALSE"
    benchmark:
             "benchmarks/visualizing_markers.benchmark.txt"
    log: "logs/markervis/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        scVis -i  {input.cluster_rds} \
               -o {params.outdir} \
              --assay {params.assay} \
              --image TRUE \
              markervis  \
                   -l {input.top10_marker}  \
                   -n  10 \
                   -c gene_diff  \
                    -m vlnplot,featureplot,heatmap \
                    -s 0 \
                    --reduct {params.reduct} \
                    --crop {params.crop} \
                    -p 1.2 \
                    --HE FALSE \
                    -a 1  >&2 2>{log}    
        """
