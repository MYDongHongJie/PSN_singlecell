def get_QC_params_1(wildcards):
    if MT_gene == True and HB_gene == True and percent == False:
        return f""" \
-c nFeature_RNA,nCount_RNA,log10GenesPerUMI,percent.mito,percent.HB \
-l 200,1000,{config['report_params']['log10GenesPerUMI']}," -Inf"," -Inf"  \
-L " Inf"," Inf"," Inf",{config['report_params']['percent_mito']},{config['report_params']['percent_HB']} \
-r nCount_RNA,percent.mito \
"""
    elif MT_gene == True and HB_gene == True and percent == True:
        return f""" \
-c nFeature_RNA,nCount_RNA,log10GenesPerUMI,percent.mito,percent.HB \
-l 200,1000,{config['report_params']['log10GenesPerUMI']}," -Inf"," -Inf"  \
-L " Inf"," Inf"," Inf",1,{config['report_params']['percent_HB']} \
-r nCount_RNA,percent.mito \
--gset {config['report']['Reference']}/MT_genelist.gmt \
--mito_four TRUE \
"""

    elif MT_gene == True and HB_gene == False and percent == False:
        return f""" \
-c nFeature_RNA,nCount_RNA,log10GenesPerUMI,percent.mito \
-l 200,1000,{config['report_params']['log10GenesPerUMI']}," -Inf"  \
-L " Inf"," Inf"," Inf",{config['report_params']['percent_mito']} \
-r nCount_RNA,percent.mito \
"""
    elif MT_gene == True and HB_gene == False and percent == True:
        return f""" \
-c nFeature_RNA,nCount_RNA,log10GenesPerUMI,percent.mito \
-l 200,1000,{config['report_params']['log10GenesPerUMI']}," -Inf"  \
-L " Inf"," Inf"," Inf",1 \
-r nCount_RNA,percent.mito \
--gset {config['report']['Reference']}/MT_genelist.gmt \
--mito_four TRUE \
"""

    elif MT_gene == False and HB_gene == True :
        return f""" \
-c nFeature_RNA,nCount_RNA,log10GenesPerUMI,percent.HB \
-l 200,1000,{config['report_params']['log10GenesPerUMI']}," -Inf"  \
-L "Inf","Inf","Inf",{config['report_params']['percent_HB']} \
-r nCount_RNA \
"""
    else: 
        return f""" \
-c nFeature_RNA,nCount_RNA,log10GenesPerUMI \
-l 200,1000,{config['report_params']['log10GenesPerUMI']}  \
-L "Inf","Inf","Inf" \
-r nCount_RNA \
"""
def get_QC_params_2(wildcards):
    if MT_gene == True and percent == False:
        return f""" \
-c nFeature_RNA,nCount_RNA,percent.mito \
-l NULL,NULL," -Inf"  \
-L NULL,NULL,{config['report_params']['percent_mito']} \
-r nCount_RNA,percent.mito \
--cut.1 {config['report_params']['cut1']} \
--cut.2 {config['report_params']['cut2']} \
--nfold {config['report_params']['ngene_numi_fold2sd']} \
{get_roboust_linear_model_params(wildcards)} \
"""

    elif MT_gene == True and percent == True:
        return f""" \
-c nFeature_RNA,nCount_RNA,percent.mito \
-l NULL,NULL," -Inf"  \
-L NULL,NULL,1 \
-r nCount_RNA,percent.mito \
--cut.1 {config['report_params']['cut1']} \
--cut.2 {config['report_params']['cut2']} \
--nfold {config['report_params']['ngene_numi_fold2sd']} \
--gset {config['report']['Reference']}/MT_genelist.gmt \
--mito_four TRUE \
{get_roboust_linear_model_params(wildcards)} \
"""

    else:
        return f""" \
-c nFeature_RNA,nCount_RNA \
-l NULL,NULL  \
-L NULL,NULL \
-r nCount_RNA \
--cut.1 {config['report_params']['cut1']} \
--cut.2 {config['report_params']['cut2']} \
--nfold {config['report_params']['ngene_numi_fold2sd']} \
{get_roboust_linear_model_params(wildcards)} \
"""
def get_QC_params(wildcards):
    try:
        ngene_numi_fold2sd = config['report_params']['ngene_numi_fold2sd'] 
    except KeyError:
        QC_params = get_QC_params_1(wildcards)
        print('#未启用ngene_numi_fold2sd')
    else:
        QC_params = get_QC_params_2(wildcards)
        print('#启用ngene_numi_fold2sd')
    finally:
        return QC_params


rule qc:
    """
    run  QC and qc vis 
    """
    input:
        raw_rds = "result/create/seurat.h5seurat",
    output:
        "result/Count_QC/QC_metrics_beforeQC.png",
        "result/Count_QC/QC_metrics_beforeQC.pdf",
        "result/Count_QC/QC_metrics_afterQC.png",
        "result/Count_QC/QC_metrics_afterQC.pdf",
        "result/Count_QC/cell_statitics_before_after_QC.xls",
        "result/Count_QC/filtered.h5seurat"
    params:
        outdir="result/Count_QC",
        rmdoublets = config['report_params']['rmdoublets'],
        rmdoublets_methods = config['report_params']['rmdoublets_methods'],
        QC_params = get_QC_params
        #method="nGene,nUMI,percent.mito,batch",
        #spointsize=config['params']['spointsize']
    benchmark:
        "benchmarks/qc.benchmark.txt"
    resources:
        qsub_mem=30,
        qsub_p=config['cpu']['qc']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript  {script_dir}/sctool \
-i {input.raw_rds} \
-f h5seurat \
-o {params.outdir} \
-d h5seurat \
--assay RNA \
--dataslot counts \
-j 6 \
--prefix filtered \
qc \
--QC TRUE \
--pointsize 0.1 \
--mincell4gene 0 \
{params.QC_params} \
-m LogNormalize \
--rmdoublets {params.rmdoublets} \
--method {params.rmdoublets_methods}  
"""
