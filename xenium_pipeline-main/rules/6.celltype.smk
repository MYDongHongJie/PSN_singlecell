#========================================== step 1: celltype ===========================================================
rule celltype:
    """
    celltype
    """
    input:
        inputfile = expand(["result/3.Clustering/singlecell_object.clustering_resolution{resolution}.rds"],resolution=resolution)
    output:
        expand([
            "result/6.Reference_Celltype/{species}_celltyping_plot.png",
            "result/6.Reference_Celltype/{species}_celltyping_results.csv",
            "result/6.Reference_Celltype/{species}_celltyping_spatical_plot.png"
        ],species=config['report']['Species']),
        "result/6.Reference_Celltype/seurat_celltype.rds"
    benchmark:
        "benchmarks/celltype/celltype.txt"
    log:
        "logs/celltype/celltype.log"
    params:
        location=config['params']['Celltyping']['location'],
        refmarkerfile=config['database']['pannelcelltypeFile'],
        species=config['report']['Species'],
        clusterby=config['params']['Celltyping']['clusterby'],
        palette=config['params']['Celltyping']['palette'],
        reduct=reduct2
    shell:
        '''
        sctool -i {input.inputfile} \
                -d rds \
                -f rds \
                -o result/6.Reference_Celltype \
                --assay SCT \
                --update FALSE \
                --prefix seurat_celltype \
                celltyping --loc v \
                            --refmarkerfile {params.refmarkerfile} \
                            --refmethod sctype \
                            --species {params.species} \
                            --clusterby {params.clusterby} \
                            --palette {params.palette} \
                            --reduct SCT_{params.reduct}  >&2 2>{log}
        '''
#========================================== step 2: celltype vis =======================================================
rule celltype_vis:
    """
    celltype split by sampleid
    """
    input:
        inputfile = "result/6.Reference_Celltype/seurat_celltype.rds"
    output:
        expand([
            "result/6.Reference_Celltype/visualize_celltype_by_celltype/groupby-sampleid_resolution{resolution}_contrast_plot.png",
            "result/6.Reference_Celltype/visualize_celltype_by_celltype/groupby-sampleid_resolution-{resolution}_summary_plot.png",
            "result/6.Reference_Celltype/visualize_celltype_by_celltype/splitby-sampleid_resolution{resolution}_split_plot.png"
        ],resolution=resolution)
    benchmark:
        "benchmarks/celltype/celltype_vis.txt"
    log:
        "logs/celltype/celltype_vis.log"
    params:
        reduct2 = reduct2,
        resolution=resolution,
        ptsize = config['params']['Celltyping']['ptsize'],
        location = config['params']['Celltyping']['location']
    shell:
        '''
        scVis -i {input.inputfile} \
                -f rds \
                -o result/6.Reference_Celltype \
                --reduct SCT_{params.reduct2} \
                --assay SCT \
                dimplot --resolution {params.resolution} \
                        --ptsize {params.ptsize} \
                        --groupby celltype \
                        --splitby sampleid \
                        --groups sampleid \
                        --propby celltype \
                        --crop TRUE \
                        --loc {params.location} >&2 2>{log}
        '''
