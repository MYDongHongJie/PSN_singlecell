
def get_gset_param(wildcards):
    if MT_gene == True and HB_gene == True :
        return "--gset "+config['report']['scrna_Reference']+"/MT_genelist.gmt"+","+config['report']['scrna_Reference']+"/HB_genelist.gmt"
    elif MT_gene == True and HB_gene == False :
        return "--gset "+config['report']['scrna_Reference']+"/MT_genelist.gmt"
    elif MT_gene == False and HB_gene == True :
        return "--gset "+config['report']['scrna_Reference']+"/HB_genelist.gmt"
    else: 
        return ""

rule sctools_create:
    """
    run sctools to create merge rds 
    """
    input:
        metadata = config['report_params']['samples_file'],
        #MT_gene=f"{config['report']['scrna_Reference']}/MT_genelist.gmt",
        cellranger_outfiles = expand([
                                       "result/cellranger/{sample}/outs/per_sample_outs/{sample}/web_summary.html",
                                        "result/cellranger/{sample}/outs/per_sample_outs/{sample}/metrics_summary.csv",
                                        "result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_cloupe.cloupe",
                                        "result/cellranger/{sample}/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix.h5"], sample=samples.index)
    output:
        "result/create/seurat.h5seurat"
    params:
        cellranger_outdir = "result/cellranger",
        outdir = "result/create/",
        gset_param = get_gset_param
    benchmark:
        "benchmarks/sctools_create.benchmark.txt"
    #log:  "logs/create/sctools_create.log"
    resources:
        qsub_mem=30,                    ###  modify indentation
        qsub_p=config['cpu']['sctools_create']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript scripts/VDJ_scripts/VDJ_create.r \
-i {params.cellranger_outdir} \
-o {params.outdir} \
-d h5seurat \
--assay RNA  \
-s mtx \
-m {input.metadata} \
--gcolumn 2 \
{params.gset_param} 
        """
