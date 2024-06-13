import os
import pandas as pd
import sys
import time
from snakemake.utils import validate

##load config file======================================================================================================
configfile: "config/config.yaml"
configfile: "config/cluster.yaml"
project = config["report"]["Project_Num"]
Species = config["report"]["Species"]
ref_dir = config['database']['ref_genome']

## get sampleid=========================================================================================================
if os.path.exists(config["metadatafile"]):
    samples = pd.read_csv(config["metadatafile"],dtype=str).set_index("sampleid",drop=False)
    samples.index.names = ["sample_id"]
    wildcard_constraints:
        sample="|".join(samples.index)


##defination function ==================================================================================================
def get_sample(wildcards):
    return wildcards.sample


def get_fastqs_cb_reads(wildcards):
    sampleid = samples.loc[wildcards.sample, ["sampleid"]].dropna().iloc[0]
    fastqs_cb_reads = os.path.abspath(f'result/0.clean_data/{sampleid}/{sampleid}_1.fq.gz')
    return fastqs_cb_reads


def get_fastqs_cdna_reads(wildcards):
    sampleid = samples.loc[wildcards.sample, ["sampleid"]].dropna().iloc[0]
    fastqs_cdna_reads = os.path.abspath(f'result/0.clean_data/{sampleid}/{sampleid}_2.fq.gz')
    return fastqs_cdna_reads


#=======================================================================================================================
##judget MT-gene and HB-gene=
## rRNA list exist
rRNA_gene = False
if os.path.exists(f"{ref_dir}/rRNA_genelist.txt"):
    with open(f"{ref_dir}/rRNA_genelist.txt") as gmt:
        rRNA_gene = any([1 if line.find("percent.rRNA") == 0 else 0 for line in gmt])
## tRNA list exist
tRNA_gene = False
if os.path.exists(f"{ref_dir}/tRNA_genelist.txt"):
    with open(f"{ref_dir}/tRNA_genelist.txt") as gmt:
        tRNA_gene = any([1 if line.find("percent.tRNA") == 0 else 0 for line in gmt])


#=======================================================================================================================
def get_gset_param(wildcards):
    if rRNA_gene == True and tRNA_gene == True:
        file_list = f"--gset {ref_dir}/rRNA_genelist.txt,{ref_dir}/tRNA_genelist.txt"
    elif rRNA_gene == True and tRNA_gene == False:
        file_list = f"--gset  {ref_dir}/rRNA_genelist.txt"
    elif rRNA_gene == False and tRNA_gene == True:
        file_list = f"--gset {ref_dir}/tRNA_genelist.txt"
    else:
        file_list = ""
    return file_list


#=======================================================================================================================
def get_qc_filter_parameters(wildcards):
    filter = config["params"]['QC']['filters']
    lower = config["params"]['QC']['lower']
    upper = config["params"]['QC']['upper']
    #============================================
    if rRNA_gene == True and tRNA_gene == True:
        filter = filter + ',percent.rRNA' + ',percent.tRNA'
        lower = lower + ',0' + ',0'
        upper = upper + ',Inf' + ',Inf'
        return [filter, lower, upper]
    elif rRNA_gene == True and tRNA_gene == False:
        filter = filter + ',percent.rRNA'
        lower = lower + ',0'
        upper = upper + ',Inf'
        return [filter, lower, upper]
    elif rRNA_gene == False and tRNA_gene == True:
        filter = filter + ',percent.tRNA'
        lower = lower + ',0'
        upper = upper + ',Inf'
        return [filter, lower, upper]
    else:
        return [filter, lower, upper]


##read diff_group file =================================================================================================
if config["module"]["Diffgene"] and os.path.exists(config["diff_DEG_group"]):
    diff_groups = pd.read_csv(config["diff_DEG_group"],dtype=str,delimiter="\t")
    diff_groups.index = diff_groups["type"] + "_" + diff_groups["treatment_name"] + "-vs-" + diff_groups["control_name"]
    diff_groups['groupname'] = diff_groups.index
    diff_groups['foldchange'] = config['params']['Diffexp']['FC']
    diff_groups['pvalue'] = config['params']['Diffexp']['padj'] if (config['params']['Diffexp']['padj']) else \
        config['params']['Diffexp']['pvalue']
    diff_groups['pval_or_padj'] = "padj" if (config['params']['Diffexp']['padj']) else "pval"
    wildcard_constraints:
        diff_group="|".join(diff_groups.index)
    groupname = lambda w: diff_groups.loc["{}".format(w.diff_group)]["groupname"]
    pvalue = lambda w: diff_groups.loc["{}".format(w.diff_group)]["pvalue"]
    foldchange = lambda w: diff_groups.loc["{}".format(w.diff_group)]["foldchange"]
    pval_or_padj = lambda w: diff_groups.loc["{}".format(w.diff_group)]["pval_or_padj"]


##all output file check ================================================================================================
def all_input(wildcards):
    wanted_input = []
    #===================================================================================================================
    if config["module"]["STARsolo"]:
        for (sample) in samples.index:
            wanted_input.extend(
                expand([
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/Summary.csv",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.png",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.pdf",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.png",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.pdf",
                    "result/1.STARsolo/Read_Distribution/{sample}_rna_metrics.txt",
                    "result/1.STARsolo/Read_Distribution/reads_distribution_for_{sample}.pdf",
                    "result/1.STARsolo/STARsolo.statistic.tsv",
                    f"result/report/{project}_QC_Report_{time.strftime('%Y_%m_%d')}/Report.html"
                ],
                    sample=sample
                )
            )
    #===================================================================================================================
    if config["module"]["Create"]:
        wanted_input.extend(["result/seurat_ob.h5seurat"])
    #===================================================================================================================
    if config["module"]["QC"]:
        wanted_input.extend(
            [
                "result/2.Count_QC/seurat_ob.h5seurat",
                "result/2.Count_QC/statitics_before_after_QC.xls"
            ])
    #===================================================================================================================
    if config["module"]["Clustering"]:
        ##RNA
        wanted_input.extend(
            expand([
                "result/3.Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/3.Clustering/visualize_cluster_by_SCT.{reduct2}.res.{resolution}/clust_cond_freq_info.xls"
            ],
                reduct1=config['params']['Clustering']["RNA_reduct1"],
                reduct2=config['params']['Clustering']["RNA_reduct2"],
                resolution=config['params']['Clustering']["RNA_resolution"]
            )
        )
        wanted_input.extend(["logs/clustering/RNA_cluster.log", "logs/clustering/RNA_vis.log"])
    #===================================================================================================================
    if config["module"]["Anotation"]:
        wanted_input.extend([
            #kegg anno
            "result/4.Annotation/KEGG/KEGG.blast.xml",
            "result/4.Annotation/KEGG/KEGG.blastm8.xls",
            "result/4.Annotation/KEGG/KEGG.blast.best.xls",
            "result/4.Annotation/KEGG/Unigene.KEGG.gene.anno.xls",
            "logs/anno/kegg.log",
            #swissprot
            "result/4.Annotation/SWISSPROT/swissprot.blast.xml",
            "result/4.Annotation/SWISSPROT/swissprot.blast.best.xls",
            "result/4.Annotation/SWISSPROT/swissprot.blast.anno.xls",
            "logs/anno/swissprot.log",
            #go anno
            "result/4.Annotation/GO/Unigene.GO.classification.stat.pdf",
            "result/4.Annotation/GO/Unigene.GO.gene.anno.xls",
            "logs/anno/go.log",
            #CARD anno
            "result/4.Annotation/CARD/UnigeneCARD.anno.xls",
            "logs/anno/card.log",
            #CAZy anno
            "result/4.Annotation/CAZy/CAZy.blast.xml",
            "result/4.Annotation/CAZy/CAZy.blast.best.xls",
            "result/4.Annotation/CAZy/CAZy_class_plot.pdf",
            "logs/anno/cazy.log"
        ])
    ##==================================================================================================================
    if config["module"]["Marker"]:
        wanted_input.extend(
            [
                "result/5.Marker/all_markers_for_each_cluster.xls",
                "result/5.Marker/top10_markers_for_each_cluster.xls",
                "result/5.Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
                "result/5.Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
                "logs/marker/RNA_findmarker.log",
                "logs/marker/RNA_markervis.log",
                "result/6.MarkerEnrich/GO_enrichment/enrichment_go.xls",
                "result/6.MarkerEnrich/KEGG_enrichment/enrichment_kegg.xls",
                "logs/marker/marker_enrich_praper.log",
                "logs/marker/marker_enrich.log",
                "result/backfile/go.backgroud.xls",
                "result/backfile/kegg.backgroud.xls",
                "result/backfile/category.xls"
            ])

    ##==================================================================================================================
    if config["module"]["Diffgene"]:
        for groupname in diff_groups.index:
            wanted_input.extend(
                expand(
                    ["result/7.Diffgene/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls",
                     "result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_heatmap.pdf",
                     "result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_heatmap.png",
                     "result/7.Diffgene/{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}_top25_genes.xls",
                     "logs/diffexp/diffexp_RNA_heatmap_{groupname}-{pval_or_padj}-{pvalue}-FC-{foldchange}.log",
                     "result/8.DiffgeneEnrich/GO_enrichment/enrichment_go.xls",
                     "result/8.DiffgeneEnrich/KEGG_enrichment/enrichment_kegg.xls"
                     ],
                    groupname=groupname,
                    pval_or_padj="padj" if (config['params']['Diffexp']['padj']) else "pval",
                    pvalue=config['params']['Diffexp']['padj'] if (config['params']['Diffexp']['padj']) else
                    config['params']['Diffexp']['pvalue'],
                    foldchange=config['params']['Diffexp']['FC']
                )
            )

    return wanted_input


##======================================================================================================================
##======================================================================================================================
def report_input(wildcards):
    wanted_input = []
    if config["module"]["Report"]:
        wanted_input.extend(
            [f"result/report/{project}_Report_{time.strftime('%Y_%m_%d')}/Report.html"])
        return wanted_input
