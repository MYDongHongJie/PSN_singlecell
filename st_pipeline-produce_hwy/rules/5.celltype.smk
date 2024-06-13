rule  celltype:
    """
    spotlight/RCTD  celltype
    """
    input:
        cluster_rds=expand([
                        "result/cluster_seurat/singlecell_object.clustering_resolution{resolution}.rds"],
                        resolution=config["params"]["resolution"]
                            )
    output:
         "result/celltype/spatial_seurat_celltype.rds"
    params:
           ref_rds = config['database']['ref_singlecell'],
           ref_marker = config['database']['ref_marker'],
           ref_celltype = config['database']['ref_celltype'],
           outdir="result/celltype/",
           celltype=config["params"]["celltype"],
           resolution=config["params"]["resolution"],
           qsub_mem=config['celltype']['mem'],
           qsub_p=config['celltype']['cpu']
    benchmark:
             "benchmarks/celltype.benchmark.txt"
    log:
        f"logs/{celltype}/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
            if [[ {params.celltype} == "spotlight" ]];then
                sctool \
                    -i {input.cluster_rds} \
                    -f rds \
                    -o {params.outdir} \
                    -d rds  \
                    -j {params.qsub_p}\
                    --assay Spatial \
                    --dataslot counts \
                    --prefix 'spatial_seurat_celltype' \
                    st_deconv \
                    --refexp {params.ref_rds} \
                    --refcelltype {params.ref_celltype} \
                    --refassay  RNA \
                    --refmarker {params.ref_marker}   >&2 2>{log}    
                    
            elif [ {params.celltype} == "RCTD" ];then        
                sctool \
                    -i  {input.cluster_rds}   \
                    -f rds  \
                    -o {params.outdir}  \
                    -d rds  \
                    -j {params.qsub_p} \
                    --assay Spatial  \
                    --dataslot counts  \
                    --prefix 'spatial_seurat_celltype' \
                     RCTD \
                     --doublet_mode FALSE \
                     --refexp {params.ref_rds}  \
                     --refcelltype {params.ref_celltype} \
                     --refassay  RNA  >&2 2>{log}    
            fi           
        """

rule  celltype_visual:
    """
    spotlight/RCTD  celltype visualization
    """
    input:
        celltype_rds ="result/celltype/spatial_seurat_celltype.rds"
    output:
        "result/celltype/combined_all_slice-celltype.xls",
        "result/celltype/2.celltype_interaction/combined_all_slice-interactions.xls",
        "result/celltype/2.celltype_interaction/all_celltype_interaction.pdf",
        "result/celltype/2.celltype_interaction/all_celltype_interaction.png"
    params:
        outdir="result/celltype/",
        miscname= "spotlight_results" if config["params"]["celltype"]=="spotlight" else "RCTD_results",
        crop= "TRUE" if config["params"]["library_type"] == "cytassist" else "FALSE"
    benchmark:
        "benchmarks/celltype_visual.benchmark.txt"
    log:
        "logs/spatialpieplot/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        module list
        which R
        which scVis
        which plot-spatialpieplot.R
        scVis    -i {input.celltype_rds} \
                 -f rds \
                 --output {params.outdir}  \
                 --assay Spatial \
                 spatialpieplot \
                     --group.by clusters \
                     --piescale 0.38 \
                     --misclist {params.miscname}  \
                     --crop {params.crop} \
                     --pt.alpha TRUE   >&2 2>{log}        
        """

