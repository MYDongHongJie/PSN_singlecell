from snakemake.utils import validate
import pandas as pd
import os,time,sys

##load config file============================================================
configfile: "config/config.yaml"
project = config["report"]["Project_Num"]
reference = config["report"]["Reference"]
Species = config["report"]["Species"]

##judget MT-gene and HB-gene====================================================
## MT list exist
MT_gene = False
if os.path.exists(f"{config['report']['Reference']}/genes/MT_genelist.gmt") :
    with open(f"{config['report']['Reference']}/genes/MT_genelist.gmt") as gmt:
        MT_gene = any([1  if line.find("percent.mito") == 0 else 0 for line in gmt])
## HB list exist
HB_gene = False
if os.path.exists(f"{config['report']['Reference']}/genes/HB_genelist.gmt") :
    with open(f"{config['report']['Reference']}/genes/HB_genelist.gmt") as gmt:
        HB_gene = any([1  if line.find("percent.HB") == 0 else 0 for line in gmt])


##read diff_group file =========================================================
if config["module"]["diffexp"] and os.path.exists(config["diff_DEG_group"]):
    diff_groups = pd.read_csv(config["diff_DEG_group"], dtype=str,delimiter="\t")
    diff_groups.index=diff_groups["type"]+"_"+diff_groups["treatment_name"]+"-vs-"+diff_groups["control_name"]
    wildcard_constraints:
                    diff_group="|".join(diff_groups.index)
## checking raw_data address====================================================
if config["cellranger_params"]["raw_data_obs_address"] is None or config["cellranger_params"]["raw_data_obs_address"][0] is None or config["cellranger_params"]["raw_data_obs_address"][0] == "":
    print("\033[33;1m"+"WARNING: No obs address provided.","\033[0m")
    if config["cellranger_params"]["samples"] == "all" :
        sys.exit("\033[31;1m"+"ERROR: You need to specify sample name if No obs address is provided."+"\033[0m")
else:
    for addr in config["cellranger_params"]["raw_data_obs_address"]:
        if addr[0:6] == "obs://":
            if addr[-1] != "/":
                sys.exit("\033[31;1m"+"ERROR: cellranger_params:raw_data_obs_address: "+addr+" need to end with '/' !!!"+"\033[0m")
        else: sys.exit("\033[31;1m"+"ERROR: cellranger_params:raw_data_obs_address: "+addr+" need to start with 'obs://' !!!"+"\033[0m")

## get sampleid=============================================================
if type(config["cellranger_params"]["samples"]) == str :
    # type1. "all"
    if config["cellranger_params"]["samples"] == "all" :
        samples = pd.read_csv(config["metadatafile"],dtype=str).set_index("sampleid",drop=False)
        samples.index.names = ["sample_id"]
        wildcard_constraints:
            sample="|".join(samples.index)
    else :  sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples can only be strings 'all' or a list of sample."+"\033[0m")
else :
    # type2. list of sample:cell_number dictionary
    if all(isinstance(i, str) for i in config["cellranger_params"]["samples"]):
        samples = config["cellranger_params"]["samples"]
    # type3. list of sample name strings
    elif all(isinstance(i, dict) for i in config["cellranger_params"]["samples"]):
        samples = [str(v)+"_"+str(k) for i in config["cellranger_params"]["samples"] for v,k in i.items()]
        FC = {str(v)+"_"+str(k):k for i in config["cellranger_params"]["samples"] for v,k in i.items()}
    else:
        sys.exit("\033[31;1m"+"ERROR: cellranger_params:samples list can only be strings or <name:cell_number> dictionary."+"\033[0m")
        wildcard_constraints:
            sample="|".join(samples)
##defination function =========================================================
def get_sample(wildcards):
    return wildcards.sample

def get_libraries_path(wildcards):
    libraries_path=os.path.abspath("config/library/"+wildcards.sample+"_libraries.csv")
    return libraries_path

def get_FC(wildcards):
    try: FC
    except NameError: return ""
    else: return "--force-cells=" + str(FC[wildcards.sample])

def get_gset_param(wildcards):
    if MT_gene == True and HB_gene == True :
        return "--gset "+config['report']['Reference']+"/genes/MT_genelist.gmt"+","+config['report']['Reference']+"/genes/HB_genelist.gmt"
    elif MT_gene == True and HB_gene == False :
        return "--gset "+config['report']['Reference']+"/genes/MT_genelist.gmt"
    elif MT_gene == False and HB_gene == True :
        return "--gset "+config['report']['Reference']+"/genes/HB_genelist.gmt"
    else:
        return ""

def get_qc_filter_parameters(wildcards):
    filter=config['qc_params']['filters']
    lower=config['qc_params']['lower']
    upper=config['qc_params']['upper']
    if MT_gene == True and HB_gene == True:
        filter=filter+',percent.mito'+',percent.HB'
        lower=lower+',0'+',0'
        upper=upper+',0.15'+',0.05'
        return [filter,lower,upper]
    elif MT_gene == True and HB_gene == False:
        filter = filter + ',percent.mito'
        lower = lower + ',0'
        upper = upper + ',0.15'
        return [filter,lower,upper]
    elif MT_gene == False and HB_gene == True:
        filter = filter + ',percent.HB'
        lower = lower + ',0'
        upper = upper + ',0.05'
        return [filter,lower,upper]
    else:
        return [filter,lower,upper]

##all output file check ========================================================
def all_input(wildcards):
    wanted_input = []
    #=====================================================================
    if config["module"]["cellranger"]:
        for (sample) in samples.index:
            wanted_input.extend(
                expand([
                        ##data_prepare
                        "config/library/{sample}_libraries.csv" ,
                        os.path.abspath("config/library/libraries_aggr.csv") ,
                        # ##cellreanger_output
                        "result/cellranger/{sample}/outs/web_summary.html",
                        "result/cellranger/{sample}/outs/summary.csv",
                        "result/cellranger/{sample}/outs/atac_fragments.tsv.gz",
                        "result/cellranger/{sample}/outs/per_barcode_metrics.csv",
                        "result/cellranger/{sample}/outs/gex_molecule_info.h5",
                        ##report_prepare
                        "result/cellranger/{sample}/outs/{sample}_summary.csv",
                        "result/cellranger/{sample}/outs/{sample}_summary.html",
                        "result/cellranger/{sample}/{sample}.png"
                    ],
                    sample=sample
                )
            )
        ###
        if  config["cellranger_params"]["aggr"]:
             wanted_input.extend(["logs/cellranger/run_aggr.log"])
        ##report
        wanted_input.extend(
            [
                "logs/summary_report.log",
                "logs/upload_sor"
            ])
    if config["module"]["QCreport"]:
        wanted_input.extend([f"result/report/{project}_QCreport_{time.strftime('%Y_%m_%d')}/report.html"])
     ##====================================================================
    if config["module"]["create"]:
        wanted_input.extend(  [ "result/seurat_ob.h5seurat"])
    ##=====================================================================
    if config["module"]["qc"]:
        wanted_input.extend(
            [
                "result/2.Count_QC/seurat_ob.h5seurat",
                "result/2.Count_QC/statitics_before_after_QC.xls"
            ])
    ##=====================================================================
    if config["module"]["clustering"]:
        ##RNA
        wanted_input.extend(
            expand([
                "result/3.RNA_Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
                "result/3.RNA_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
                "result/3.RNA_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
                "result/3.RNA_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/3.RNA_Clustering/visualize_cluster_by_SCT.{reduct2}.res.{resolution}/clust_cond_freq_info.xls"
            ],
                reduct1=config["cluster_params"]["RNA_reduct1"],
                reduct2=config["cluster_params"]["RNA_reduct2"],
                resolution=config["cluster_params"]["RNA_resolution"]
            )
        )
        ##ATAC
        wanted_input.extend(
            expand([
                "result/6.ATAC_Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
                "result/6.ATAC_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
                "result/6.ATAC_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
                "result/6.ATAC_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/6.ATAC_Clustering/visualize_cluster_by_ATAC.{reduct2}.res.{resolution}/clust_cond_freq_info.xls"
            ],
                reduct1=config["cluster_params"]["ATAC_reduct1"],
                reduct2=config["cluster_params"]["ATAC_reduct2"],
                resolution=config["cluster_params"]["ATAC_resolution"]
            )
        )
        wanted_input.extend(["logs/clustering/ATAC_cluster.log","logs/clustering/ATAC_vis.log","logs/clustering/RNA_cluster.log","logs/clustering/RNA_vis.log"])
    ##========================================================================
    if config["module"]["celltyping"]:
        wanted_input.extend(
            expand([
                "result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_heatmap.pdf",
                "result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_plot.pdf",
                "result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_results.xls",
                "result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_simplified_celltyping_results.csv",
                "result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_top.{annolevel}_celltyping_plot.pdf"
            ],
                species=config["report"]["Species"],
                refbuiltin=config["celltyping_params"]["refbuiltin"],
                annolevel=config["celltyping_params"]["annolevel"])
        )
        wanted_input.extend(
            [
                "logs/celltyping/RNA_celltyping.log",
                "logs/celltyping/ATAC_celltyping.log"
            ])

    ##========================================================================
    if config["module"]["marker"]:
        wanted_input.extend(
            [
                "result/4.RNA_Marker/all_markers_for_each_cluster.xls",
                "result/4.RNA_Marker/top10_markers_for_each_cluster.xls",
                "result/4.RNA_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
                "result/4.RNA_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
                "logs/marker/RNA_findmarker.log",
                "logs/marker/RNA_markervis.log",
                "result/7.ATAC_Marker/all_markers_for_each_cluster.xls",
                "result/7.ATAC_Marker/top10_markers_for_each_cluster.xls",
                "result/7.ATAC_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
                "result/7.ATAC_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
                "logs/marker/ATAC_findmarker.log",
                "logs/marker/ATAC_markervis.log"
            ])
    ##========================================================================
    if config["module"]["wnn"]:
        wanted_input.extend(
            [
                "result/10.WNN_RNA_Marker/all_markers_for_each_cluster.xls",
		        "result/10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls",
                "result/10.WNN_RNA_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
                "result/10.WNN_RNA_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
                "result/11.WNN_ATAC_Marker/all_markers_for_each_cluster.xls",
                "result/11.WNN_ATAC_Marker/top10_markers_for_each_cluster.xls",
                "result/11.WNN_ATAC_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
                "result/11.WNN_ATAC_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
                "logs/wnn/wnn_RNA_marker.log",
                "logs/wnn/wnn_RNA_vis.log",
                "logs/wnn/wnn_ATAC_marker.log",
                "logs/wnn/wnn_ATAC_vis.log",
                "logs/wnn/wnn_clustering.log",
                "logs/wnn/wnn_vis.log",
                "logs/wnn/wnn_celltyping.log"
            ]
        )
        wanted_input.extend(
            expand([
                ##cluster
                "result/9.WNN_Clustering/visualize_cluster_by_wnn.res.{resolution}/clust_cond_freq_info.xls",
                "result/9.WNN_Clustering/wnn_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/9.WNN_Clustering/wnn_groupby_cluster_resolution{resolution}_plot.png",
                ##celltyping
                "result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_heatmap.pdf",
                "result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_plot.pdf",
                "result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_results.xls",
                "result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_simplified_celltyping_results.csv",
                "result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_top.{annolevel}_celltyping_plot.pdf"
            ],
                species=config["report"]["Species"],
                refbuiltin=config["wnn_params"]["refbuiltin"],
                annolevel=config["wnn_params"]["annolevel"],
                resolution=config["wnn_params"]["resolution"])
        )
        ##======================================================================
    if config["module"]["diffexp"]:
        wanted_input.extend(
            expand([
                ##run diffexp
                "result/13.Diffexp_RNA/{diff_group}-all_diffexp_genes.xls",
		        "result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5.xls",
                "result/14.Diffexp_ATAC/{diff_group}-all_diffexp_genes.xls",
                "result/14.Diffexp_ATAC/{diff_group}-diff-pval-0.05-FC-1.5.xls",
                ##run diff_heatmap
                "result/13.Diffexp_RNA/top25_{diff_group}_heatmap.pdf",
                "result/13.Diffexp_RNA/top25_{diff_group}_heatmap.png",
                "result/13.Diffexp_RNA/top25_{diff_group}_genes.xls",
                "logs/diffexp/diffexp_RNA_heatmap_{diff_group}.log",
                "result/14.Diffexp_ATAC/top25_{diff_group}_heatmap.pdf",
                "result/14.Diffexp_ATAC/top25_{diff_group}_heatmap.png",
                "result/14.Diffexp_ATAC/top25_{diff_group}_genes.xls",
                "logs/diffexp/diffexp_ATAC_heatmap_{diff_group}.log",
                ##run annotation
                "result/13.Diffexp_RNA/{diff_group}-all_diffexp_genes_anno.xls",
                "result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5_anno.xls",
                "logs/diffexp/diffgene_annotation_{diff_group}.log",
                ##run enrich
                "result/15.enrichment/GO_enrichment/enrichment_go.xls",
                "result/15.enrichment/KEGG_enrichment/enrichment_kegg.xls"
            ],diff_group=diff_groups.index)
        )
    ##=======================================================================
    if config["module"]["report"]:
        wanted_input.extend([f"result/report/{project}_Report_{time.strftime('%Y_%m_%d')}/report.html"])

    return wanted_input
