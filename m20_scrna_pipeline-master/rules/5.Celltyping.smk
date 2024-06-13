##================================================== RNA celltyping=======================================================

rule celltyping:
    '''
    use sctool to identification of cell types
    '''
    input:
        donefile = "logs/sub_cluster.done"
    output:
        files = expand([
            "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_celltyping_heatmap.pdf",
            "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_celltyping_plot.pdf",
            "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_celltyping_results.xls",
            "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_simplified_celltyping_results.csv",
            "result/6.Reference_celltype/{Species}ref_{refbuiltin}_top.{annolevel}_celltyping_plot.pdf"
        ],
            Species = Species.lower(),
            refbuiltin = refbuiltin,
            annolevel = annolevel),
        donefile = "logs/Celltype.done"
    benchmark:
        "benchmarks/celltyping/RNA_celltyping.benchmark.txt"
    log:
        "logs/celltyping/RNA_celltyping.log"
    resources:
        qsub_mem = config['celltyping']['mem']
    threads:
        config['celltyping']['cpu']
    params:
        inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
        outpath = "result/6.Reference_celltype",
        refbuiltin = refbuiltin,
        annolevel = annolevel,
        demethod = demethod,
        reduct = reduct2,
        species = Species.lower()
    shell:
        '''
        {mp}
        {m20sc_sctool}
        sctool -i {params.inputfile} \
              -f h5seurat \
              -o {params.outpath} \
              -d h5seurat \
              --assay SCT \
              --dataslot data,counts \
              -j {threads} \
              celltyping    --refbuiltin {params.refbuiltin} \
                            --annolevel {params.annolevel} \
                            --demethod {params.demethod} \
                            --clusterby clusters \
                            -v 0.8 \
                            -n 25 \
                            --palette customecol2 \
                            --reduct SCT_{params.reduct} \
                            --species {params.species} \
                            --doplot TRUE >&2 2>{log} && touch {output.donefile} 
        '''

