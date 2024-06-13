def get_gset_param(wildcards):
    if MT_gene == True and HB_gene == True :
        return "--gset "+config['report']['Reference']+"/MT_genelist.gmt"+","+config['report']['Reference']+"/HB_genelist.gmt"
    elif MT_gene == True and HB_gene == False :
        return "--gset "+config['report']['Reference']+"/MT_genelist.gmt"
    elif MT_gene == False and HB_gene == True :
        return "--gset "+config['report']['Reference']+"/HB_genelist.gmt"
    else: 
        return ""


localrules: file_link
rule file_link:
    """
    run cp csv file to create merge rds 
    """
    input:
        metrics_Mols="result/BD_Analysis/{sample}/{sample}_RSEC_MolsPerCell_MEX.zip"
    output:
        sample_metrics=("result/create/{sample}/barcodes.tsv.gz")
    params:
       sample_dir=("result/create/{sample}/")
    envmodules:
         config["report_params"]['envmodules']
    resources:
        qsub_mem=4,
        qsub_p=1,
        qsub_n=1
    shell:
        '''
#mkdir -p {params.sample_dir}
/usr/bin/unzip -d {params.sample_dir}/ {wkdir}/{input.metrics_Mols}
        '''

rule sctools_create:
    """
    run sctools to create merge rds 
    """
    input:
        metadata = config['report_params']['samples_file'],
        sample_metrics = expand(["result/create/{sample}/barcodes.tsv.gz"], sample=samples.index)
    output:
        "result/create/seurat.h5seurat"
    params:
        outdir = "result/create/",
        gset_param = get_gset_param
    benchmark:
        "benchmarks/sctools_create.benchmark.txt"
    #log:  "logs/create/sctools_create.log"
    resources:
        qsub_mem=30,                    ###  modify indentation
        #qsub_p=config['cpu']['sctools_create']
        qsub_p=4
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript scripts/sctool \
-i {params.outdir} \
-o {params.outdir} \
-d h5seurat \
--assay RNA \
create \
-s xsv \
-m {input.metadata} \
--gcolumn 2 \
{params.gset_param} \
--transpose true 
        """
