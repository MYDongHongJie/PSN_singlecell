import os
import pandas as pd
import numpy as np
import sys
import time
from snakemake.utils import validate

#=======================================================================================================================
##### load config and sample sheets #####
configfile: "config/config.yaml"
configfile: "config/cluster.yaml"

samplefile = config['database']['metadataFile']
samples_config = config['params']['Xenium']['samples']
xenium_method = config['params']['Xenium']['method']
reduct1 = config['params']['Clustering']['reduct1']
reduct2 = config['params']['Clustering']['reduct2']
resolution = config['params']['Clustering']['resolution']
project = config['report']['Project_Num']
localPath = config['params']['Xenium']['raw_data_local_address']
##======================================================================================================================
if os.path.exists(samplefile) and samples_config == "all":
    samples = pd.read_csv(samplefile, dtype=str).set_index("sampleid",drop=False)
    samples = samples[ ~np.array([s.startswith("#") for s in samples.sampleid.to_list()])]
    samples.index.names = ["sample_id"]
    wildcard_constraints:
        sample="|".join(samples.index)
elif type(samples_config) == list:
    samples = pd.DataFrame({"sampleid":samples_config})
    samples.index = samples_config
    samples.index.names = ["sample_id"]
    wildcard_constraints:
        sample="|".join(samples.index)

##======================================================================================================================
def get_sample(wildcards):
    return wildcards.sample


##======================================================================================================================
## 以下模块为针对整个snakemake rule动态地设定最终期待文件目标================================================================
def all_input(wildcards):
    wanted_input = []
    if config["module"]["Xenium"]:
        for (sample) in samples.index:
            wanted_input.extend(
                expand([
                    "result/1.Xenium/{sample}/{sample}.png",
                    "raw_data/{sample}/gene_panel.json",
                    "result/1.Xenium/{sample}/outs/cells.csv.gz",
                    "result/1.Xenium/{sample}/outs/transcripts.csv.gz",
                    "result/1.Xenium/{sample}/outs/cell_boundaries.csv.gz",
                    "result/1.Xenium/{sample}/outs/metrics_summary.csv",
                    "result/1.Xenium/{sample}/outs/analysis_summary.html",
                ],sample=sample)
            )
        wanted_input.extend([
            "result/1.Xenium/summary_all.csv"
        ])

    if config["module"]["QCreport"]:
        wanted_input.extend([f"result/report/{project}_QCReport/Report.html"])

    if config["module"]["Create"]:
        wanted_input.extend(
            ["result/seurat_ob.rds"]
        )

    if config["module"]["QC"]:
        wanted_input.extend(
            [
                "result/2.Count_QC/filtered_seurat.rds",
                "result/2.Count_QC/statitics_for_QC.xls",
                "result/2.Count_QC/QC_featureplot_for_nCount_xenium.png",
                "result/2.Count_QC/QC_featureplot_for_nFeature_xenium.png"
            ]
        )

    if config["module"]["Clustering"]:
        wanted_input.extend(
            expand([
                "result/3.Clustering/singlecell_object.clustering_resolution{resolution}.rds",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/clusters_infor.csv",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
                "result/3.Clustering/visualize_clusters_by_clusters/groupby-sampleid_resolution{resolution}_contrast_plot.png",
                "result/3.Clustering/visualize_clusters_by_clusters/groupby-sampleid_resolution-{resolution}_summary_plot.png",
                "result/3.Clustering/visualize_clusters_by_clusters/splitby-sampleid_resolution{resolution}_split_plot.png"
            ],
                resolution=resolution,
                reduct2=reduct2)
        )

    if config["module"]["Marker"]:
        wanted_input.extend(
            [
                "result/5.Marker/all_markers_for_each_cluster_anno.xls",
                "result/5.Marker/top10_markers_for_each_cluster_anno.xls",
                "result/5.Marker/topmarker_gene_heatmap.png"
            ]
        )

    if config["module"]["Celltype"]:
        wanted_input.extend(
            expand([
                "result/6.Reference_Celltype/{species}_celltyping_plot.png",
                "result/6.Reference_Celltype/{species}_celltyping_results.csv",
                "result/6.Reference_Celltype/{species}_celltyping_spatical_plot.png"
            ],species=config['report']['Species'])
        )
        wanted_input.extend(expand([
            "result/6.Reference_Celltype/visualize_celltype_by_celltype/groupby-sampleid_resolution{resolution}_contrast_plot.png",
            "result/6.Reference_Celltype/visualize_celltype_by_celltype/groupby-sampleid_resolution-{resolution}_summary_plot.png",
            "result/6.Reference_Celltype/visualize_celltype_by_celltype/splitby-sampleid_resolution{resolution}_split_plot.png"
        ],resolution=resolution))
        wanted_input.extend(["result/6.Reference_Celltype/seurat_celltype.rds"])
    return wanted_input
##======================================================================================================================
def report_input(wildcards):
    wanted_input = []
    if config["module"]["Report"]:
        wanted_input.extend(
            [f"result/report/{project}_Report/Report.html"]
            )
    return wanted_input
