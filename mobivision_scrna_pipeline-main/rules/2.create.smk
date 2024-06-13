
def get_gset_param(wildcards):
    if MT_gene == True and HB_gene == True :
        return "--gset "+config['report']['Reference']+"/MT_genelist.gmt"+","+config['report']['Reference']+"/HB_genelist.gmt"
    elif MT_gene == True and HB_gene == False :
        return "--gset "+config['report']['Reference']+"/MT_genelist.gmt"
    elif MT_gene == False and HB_gene == True :
        return "--gset "+config['report']['Reference']+"/HB_genelist.gmt"
    else: 
        return ""

rule sctools_create:
    """
    run sctools to create merge rds 
    """
    input:
        metadata = config['report_params']['samples_file'],
        #MT_gene=f"{config['report']['Reference']}/MT_genelist.gmt",
        mobivision_outfiles = expand(["result/mobivision/{sample}/outs/{sample}_Report.html",
        "result/mobivision/{sample}/outs/{sample}_summary.csv",
        "result/mobivision/{sample}/outs/{sample}_filtered.h5ad"], sample=samples.index)
    output:
        "result/create/seurat.h5seurat"
    params:
        mobivision_outdir = "result/mobivision",
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
Rscript {script_dir}/sctool \
-i {params.mobivision_outdir} \
-o {params.outdir} \
-d h5seurat \
--assay RNA  \
create \
-s mtx \
-m {input.metadata} \
--gcolumn 2 \
{params.gset_param} && Rscript scripts/R/set_ident.r -i {params.outdir} 
        """
