import os
import pandas as pd
import numpy as np
import sys
import time
from snakemake.utils import validate
from oebio.utls import cd
import copy
from glob import glob

### load config file======================================================================================================
configfile: "config/config.yaml"
configfile: "config/cluster.yaml"

samplefile = config["metadatafile"]
difffile = config["diff_DEG_group"]

project = config["report"]["Project_Num"]
Species = config["report"]["Species"].capitalize()
figures = ["png", "pdf"]

# genome
ref_dir = config['database']['ref_genome']
white_list = config['database']['white_list']
filter_genelist = config['database']['filter_genelist']
# human_gtf_url = config['database_url']['Homo']['Ref_genome_annotation']['Linkage']
# human_genome_url = config['database_url']['Homo']['Ref_genome']['Linkage']
# mouse_gtf_url = config['database_url']['Mus']['Ref_genome_annotation']['Linkage']
# mouse_genome_url = config['database_url']['Mus']['Ref_genome']['Linkage']
# STARsolo
filter_standard = config['params']['STARsolo']['empty_drops']
samples_config = config["params"]["STARsolo"]['samples']
# Clustering
reduct1 = config['params']['Clustering']["RNA_reduct1"],
reduct2 = config['params']['Clustering']["RNA_reduct2"],
resolution = config['params']['Clustering']["RNA_resolution"],
# Celltype
refbuiltin = config["params"]["Celltyping"]["refbuiltin"]
annolevel = config["params"]["Celltyping"]["annolevel"]
# Diffexp
foldchange = config['params']['Diffexp']['FC']
topn = config['params']['Diffexp']['topn']
test = config["params"]["Diffexp"]["test"]
padj = config['params']['Diffexp']['padj']
pvalue = config['params']['Diffexp']['pvalue']
# qc
mincell4gene = config['params']['QC']['mincell4gene']
normmeth = config['params']['QC']['normmeth']
qc_filters = config['params']['QC']['filters']
qc_lower = config['params']['QC']['lower']
qc_upper = config['params']['QC']['upper']
features2filter = config['params']['QC']['features2filter']
rmdoublets = config['params']['QC']['rmdoublets']
qc_method = config['params']['QC']['method']

# cluster
clusteringuse = config["params"]["Clustering"]["RNA_clusteringuse"]
palette = config["params"]["Clustering"]["RNA_palette"]
component = config["params"]["Clustering"]["RNA_component"]
batchid = config["params"]["Clustering"]["RNA_batchid"]
groupby = config["params"]["Clustering"]["RNA_groupby"]
# marker
test_method = config['params']['Marker']['RNA_test_method']
topn_marker = config['params']['Marker']['RNA_topn_marker']
RNA_topn = config['params']['Marker']['RNA_topn']
topby = config['params']['Marker']['RNA_topby']
vismethod = config['params']['Marker']['RNA_vismethod']
# celltyping
refbuiltin = config["params"]["Celltyping"]["refbuiltin"]
annolevel = config["params"]["Celltyping"]["annolevel"]
demethod = config["params"]["Celltyping"]["demethod"]


### load env ===========================================================================================================
mp = "module purge"
star = "module load " + config['envmodules']['star']
picard = "module load " + config['envmodules']['picard']
obsutil = "module load " + config['envmodules']['obsutil']
oesinglecell_visium = "module load " + config['envmodules']['oesinglecell_visium']
oesinglecell = "module load " + config['envmodules']['oesinglecell']
cellranger = "module load " + config['envmodules']['CellRanger']
gtfToGenePred = "module load " + config['envmodules']['gtfToGenePred']
samtools = "module load " + config['envmodules']['samtools']
m20sc_sctool = "module load " + os.getcwd() + "/" + config['envmodules']['m20sc_sctool']
enrich = "module load " + os.getcwd() + "/" + config['envmodules']['enrich']
qualimap = "module load " + config['envmodules']['qualimap']


### =============================================================

## 处理白名单路径：如果在白名单上指定了路径，则使用对应路径下的文件，如果没有指定路径，则认为白名单位于参考基因组的文件夹下
## 若白名单不存在，直接退出
if len(white_list.split("/")) > 1:
    if not os.path.exists(f'{white_list}'):
        print(f'{white_list} not exists')
        os._exit(0)
else:
    white_list = f"{ref_dir}/" + white_list
    if not os.path.exists(f'{white_list}'):
        print(f'{white_list} not exists in {ref_dir}')
        os._exit(0)

## 处理 filter_genelist

if len(filter_genelist.split("/")) > 1:
    if not os.path.exists(f'{filter_genelist}'):
        print(f'{filter_genelist} not exists')
        os._exit(0)
else:
    filter_genelist = f"{ref_dir}/" + filter_genelist
    if not os.path.exists(f'{filter_genelist}'):
        print(f'{filter_genelist} not exists in {ref_dir}')
        os._exit(0)

def get_features2filter(wildcards):
    if features2filter == True and os.path.exists(f'{filter_genelist}'):
        return "--features2filter " + filter_genelist
    else:
        return ""


## get sampleid=========================================================================================================
# 样本文件可以使用 “#” 进行注释
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

##defination function ==================================================================================================
def get_sample(wildcards):
    return wildcards.sample


# def get_fastqs_cb_reads(wildcards):
#     sampleid = samples.loc[wildcards.sample, ["sampleid"]].dropna().iloc[0]
#     fastqs_cb_reads = os.path.abspath(f"raw_data/{sampleid}/{sampleid}_fastqs/{sampleid}_R1.fq.gz")
#     return fastqs_cb_reads
#
# def get_fastqs_cdna_reads(wildcards):
#     sampleid = samples.loc[wildcards.sample, ["sampleid"]].dropna().iloc[0]
#     fastqs_cdna_reads = os.path.abspath(f"raw_data/{sampleid}/{sampleid}_fastqs/{sampleid}_R2.fq.gz")
#     return fastqs_cdna_reads
#
# def get_raw2clean_file(wildcards):
#     sampleid = samples.loc[wildcards.sample, ["sampleid"]].dropna().iloc[0]
#     raw2clean_inputfile = os.path.abspath(f'config/raw2clean/{sampleid}_config.info')
#     return raw2clean_inputfile
#=======================================================================================================================
##judget MT-gene and HB-gene====================================================
## MT list exists
MT_gene = False
if os.path.exists(f"{ref_dir}/MT_genelist.gmt") :
    with open(f"{ref_dir}/MT_genelist.gmt") as gmt:
        MT_gene = any([1  if line.find("percent.mito") == 0 else 0 for line in gmt])
## HB list exists
HB_gene = False
if os.path.exists(f"{ref_dir}/HB_genelist.gmt") :
    with open(f"{ref_dir}/HB_genelist.gmt") as gmt:
        HB_gene = any([1  if line.find("percent.HB") == 0 else 0 for line in gmt])

## lncRNA list exists
lncRNA_gene = False
if os.path.exists(f"{ref_dir}/lncRNA_genelist.gmt") :
    with open(f"{ref_dir}/lncRNA_genelist.gmt") as gmt:
        lncRNA_gene = any([1  if line.find("percent.lncRNA") == 0 else 0 for line in gmt])

## lncRNA list exists
rRNA_gene = False
if os.path.exists(f"{ref_dir}/rRNA_genelist.gmt") :
    with open(f"{ref_dir}/rRNA_genelist.gmt") as gmt:
        rRNA_gene = any([1  if line.find("percent.rRNA") == 0 else 0 for line in gmt])

#=======================================================================================================================

logic_list = [MT_gene, HB_gene, lncRNA_gene, rRNA_gene]
file_list = ["MT_genelist.gmt", "HB_genelist.gmt", "lncRNA_genelist.gmt", "rRNA_genelist.gmt"]
percent_lsit = ["percent.mito", "percent.HB", "percent.lncRNA", "percent.rRNA"]

def get_gset_param(wildcards):
    if any(logic_list):
        return "--gset " + f"{ref_dir}/" + f",{ref_dir}/".join(np.asarray(file_list)[np.where(np.array(logic_list) == True)[0]].tolist())
    else:
        return ""

def get_metrics_param(wildcards):
    if any(logic_list):
        return "--metrics " + ",".join(np.asarray(percent_lsit)[np.where(np.array(logic_list) == True)[0]].tolist())
    else:
        return ""

#=======================================================================================================================
def get_qc_filter_parameters(wildcards):
    filter = ''.join(qc_filters.split()).split(",")
    lower = ''.join(qc_lower.split()).split(",")
    upper = ''.join(qc_upper.split()).split(",")
    #============================================
    qc_logic_list = [True] * 3 + logic_list
    filter = ",".join(np.asarray(filter)[np.where(np.array(qc_logic_list) == True)[0]].tolist())
    lower = ",".join(np.asarray(lower)[np.where(np.array(qc_logic_list) == True)[0]].tolist())
    upper = ",".join(np.asarray(upper)[np.where(np.array(qc_logic_list) == True)[0]].tolist())
    return filter, lower, upper

##read diff_group file =================================================================================================
if config["module"]["Diffexp"] and os.path.exists(difffile):
    diff_groups = pd.read_csv(difffile,dtype=str,delimiter="\t")
    diff_groups.index = diff_groups["type"] + "_" + diff_groups["treatment_name"] + "-vs-" + diff_groups["control_name"]
    diff_groups['groupname'] = diff_groups.index
    wildcard_constraints:
        diff_group="|".join(diff_groups.index)


### =============================================================
# report modules
if config["module"]["Report"]:
    report_modules = copy.deepcopy(config['module'])
    # del report_modules['genome']
    del report_modules['Report']
    del report_modules['QCreport']
    inputfiles_report=["logs/" + module + ".done" for module in [module for module, b in report_modules.items() if b == True]]


# if Species.capitalize() == 'human'.capitalize():
#     gtf_url = human_gtf_url
#     gtffile_gz = gtf_url.split("/")[-1]
#     gtffile = gtffile_gz.split(".gz")[0]
#     genome_url = human_genome_url
#     genomefile_gz = genome_url.split("/")[-1]
#     genomefile = genomefile_gz.split(".gz")[0]
#
#     if not os.path.exists(f'{ref_dir}/genome.fa') or not os.path.exists(f"{ref_dir}/rRNA.interval_list"):
#         if not os.path.exists("genome"): os.mkdir("genome")
#         with cd(f'{ref_dir}'):
#             os.system(f'wget {genome_url}')
#             os.system(f'gzip -d {genomefile_gz}')
#             os.system(f'ln -s {genomefile} genome.fa')
#         with cd(f'{ref_dir}'):
#             os.system(f'wget {gtf_url}')
#             os.system(f'gzip -d {gtffile_gz}')
#
# if Species.capitalize() == 'mouse'.capitalize():
#     gtf_url = mouse_gtf_url
#     gtffile_gz = gtf_url.split("/")[-1]
#     gtffile = gtffile_gz.split(".gz")[0]
#     genome_url = mouse_genome_url
#     genomefile_gz = genome_url.split("/")[-1]
#     genomefile = genomefile_gz.split(".gz")[0]
#
#     if not os.path.exists(f'{ref_dir}/genome.fa') or not os.path.exists(f"{ref_dir}/rRNA.interval_list"):
#         if not os.path.exists("genome"): os.mkdir("genome")
#         with cd(f'{ref_dir}'):
#             os.system(f'wget {genome_url}')
#             os.system(f'gzip -d {genomefile_gz}')
#             os.system(f'')
#         with cd(f'{ref_dir}'):
#             os.system(f'wget {gtf_url}')
#             os.system(f'gzip -d {gtffile_gz}')
#
# if config["module"]["genome"]:
#     rule preparing_genome:
#         input:
#             gtf = ancient(f"{ref_dir}/{gtffile}"),
#             fa = ancient(f"{ref_dir}/{genomefile}")
#         output:
#             gtf = f"{ref_dir}/genome.gtf",
#             refFlat = f"{ref_dir}/refFlat.txt",
#             rRNA = f"{ref_dir}/rRNA.interval_list"
#         threads:
#             config['genome']['cpu']
#         params:
#             path = f'{ref_dir}'
#         shell:
#             """
#             {mp}
#             {cellranger}
#             cellranger mkgtf {input.gtf} {output.gtf} \
#                     --attribute=gene_type:protein_coding \
#                     --attribute=gene_type:lncRNA \
#                     --attribute=gene_type:antisense \
#                     --attribute=gene_type:IG_LV_gene \
#                     --attribute=gene_type:IG_V_gene \
#                     --attribute=gene_type:IG_V_pseudogene \
#                     --attribute=gene_type:IG_D_gene \
#                     --attribute=gene_type:IG_J_gene \
#                     --attribute=gene_type:IG_J_pseudogene \
#                     --attribute=gene_type:IG_C_gene \
#                     --attribute=gene_type:IG_C_pseudogene \
#                     --attribute=gene_type:TR_V_gene \
#                     --attribute=gene_type:TR_V_pseudogene \
#                     --attribute=gene_type:TR_D_gene \
#                     --attribute=gene_type:TR_J_gene \
#                     --attribute=gene_type:TR_J_pseudogene \
#                     --attribute=gene_type:TR_C_gene \
#                     --attribute=gene_type:Mt_rRNA  \
#                     --attribute=gene_type:rRNA  \
#                     --attribute=gene_type:rRNA_pseudogene
#             {mp}
#             {star}
#             STAR --runMode genomeGenerate \
#                 --genomeDir {params.path}/STAR_index \
#                 --runThreadN {threads} \
#                 --genomeFastaFiles {params.path}/genome.fa \
#                 --sjdbGTFfile {output.gtf} \
#                 --sjdbOverhang 50
#             {mp}
#             {gtfToGenePred}
#             gtfToGenePred -genePredExt \
#                 -ignoreGroupsWithoutExons {params.path}/genome.gtf  \
#                 {params.path}/refFlat.tmp.txt
#             {mp}
#             paste <(cut -f 12 {params.path}/refFlat.tmp.txt) <(cut -f 1-10 {params.path}/refFlat.tmp.txt) > {output.refFlat}
#             grep \'gene_type \"rRNA\"\' {params.path}/genome.gtf | cut -f1,4,5,7,9 | perl -lane \'/gene_id \"([^\"]+)\"/ or die \"no transcript_id on $.\";print join \"\t\", (@F[0,1,2,3], $1)\' | sort -k1V -k2n -k3n |uniq  >> {params.path}/tmp.rRNA.interval_list
#             {samtools}
#             cat <(samtools dict {params.path}/genome.fa | cut -f 1-3) {params.path}/tmp.rRNA.interval_list > {output.rRNA}
#             rm {params.path}/tmp.rRNA.interval_list
#             """

##all output file check ================================================================================================
def all_input(wildcards):
    wanted_input = []
    #===================================================================================================================
    # preparing genome and gtf
    # if config["module"]["genome"]:
    #     if not os.path.exists(f'{ref_dir}/genome.gtf') or not os.path.exists(f"{ref_dir}/rRNA.interval_list"):
    #         wanted_input = [f"{ref_dir}/genome.gtf", f"{ref_dir}/rRNA.interval_list"]
    #===================================================================================================================
    if config["module"]["STARsolo"]:
        for (sample) in samples.index:
            wanted_input.extend(
                expand([
                    # starsolo
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/barcodes.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/features.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/filtered/matrix.mtx.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/barcodes.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/features.tsv.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/raw/matrix.mtx.gz",
                    "result/1.STARsolo/{sample}/{sample}_Solo.out/GeneFull/Summary.csv",
                    # barcode_rank_gene_umi
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.png",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_Gene_Counts_for_{sample}.pdf",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.png",
                    "result/1.STARsolo/BarcodeRank_plot/Barcode_vs_UMI_Counts_for_{sample}.pdf",
                    # genetype_plot
                    "result/1.STARsolo/Genetype_plot/{sample}_gene_type_totaldata.csv",
                    "result/1.STARsolo/Genetype_plot/{sample}_gene_type_plotdata.csv",
                    "result/1.STARsolo/Genetype_plot/{sample}_gene_type.png",
                    # reads_distribution
                    "result/1.STARsolo/Read_Distribution/{sample}_rna_metrics.txt",
                    # qualimap
                    "result/1.STARsolo/Qualimap_result/{sample}/rnaseq_qc_results.txt",
                    "result/1.STARsolo/Qualimap_result/{sample}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt",
                    "logs/qualimap/{sample}.qualimap.log",
                    # starsolo_summary
                    "result/1.STARsolo/Read_Distribution/reads_distribution_for_{sample}.pdf",
                    "result/1.STARsolo/STARsolo.statistic.tsv",
                    "logs/starsolo_{sample}.done"

                ],
                    sample=sample
                )
            )
            wanted_input.extend(["logs/data_praper.done","result/1.STARsolo/Qualimap_result/qualimap_plot.png"])

    if filter_standard != "":
        for (sample) in samples.index:
            wanted_input.extend(
                expand([
                    "logs/re_filter_{sample}.done"
                ],
                    sample = sample
                )
            )
    if config["module"]["QCreport"]:
        wanted_input.extend(
            expand(
                [f"result/report/{project}_QC_Report/Report.html","logs/STARsolo.done"]
            )
        )

    if config['params']['STARsolo']['raw_data']:
        for (sample) in samples.index:
            wanted_input.extend(
                expand([
                    "logs/raw2clean_{sample}.done"
                ],
                    sample = sample
                )
            )
    #===================================================================================================================
    if config["module"]["Create"]:
        wanted_input.extend([
            "result/seurat_ob.h5seurat", 
            "logs/Create.done"
            ])
    #===================================================================================================================
    if config["module"]["QC"]:
        wanted_input.extend(
            [
                "result/2.Count_QC/seurat_ob.h5seurat",
                "result/2.Count_QC/statitics_before_after_QC.xls",
                "logs/QC.done"
            ])
    #===================================================================================================================
    if config["module"]["Clustering"]:
        ##RNA
        wanted_input.extend(
            expand([
                
                "result/3.Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
                "result/3.Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.{formats}",
                "logs/sub_cluster.done",
                "result/3.Clustering/visualize_cluster_by_clusters/clust_cond_freq_info.xls",
                "result/3.Clustering/clusters_correlation/normalized_data_groupby_clusters.xls",
                "result/3.Clustering/clusters_correlation/coefficient_heatmap.{formats}",
                "logs/Clustering.done"
            ],
                reduct1 = reduct1,
                reduct2 = reduct2,
                resolution = resolution,
                formats = figures
            )
        )
        
    ##==================================================================================================================
    if config["module"]["Marker"]:
        wanted_input.extend(
            expand([
                "result/4.mRNA_marker/all_markers_for_each_cluster.xls",
                "result/4.mRNA_marker/top10_markers_for_each_cluster.xls",
                "result/4.mRNA_marker/markers_vis4cluster1/marker_gene_featureplot.{formats}",
                "result/4.mRNA_marker/markers_vis4cluster1/marker_gene_violin_plot.{formats}",
                "result/5.LncRNA_marker/markers_vis4cluster1/marker_gene_featureplot.{formats}",
                "result/5.LncRNA_marker/markers_vis4cluster1/marker_gene_violin_plot.{formats}",
                "logs/Marker.done"
            ],
                formats = figures
            ))

    ##==================================================================================================================

    if config["module"]["Celltype"]:
        wanted_input.extend(
            expand([
                "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_celltyping_heatmap.pdf",
                "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_celltyping_plot.pdf",
                "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_celltyping_results.xls",
                "result/6.Reference_celltype/{Species}ref_{refbuiltin}_{annolevel}_simplified_celltyping_results.csv",
                "result/6.Reference_celltype/{Species}ref_{refbuiltin}_top.{annolevel}_celltyping_plot.pdf",
                "logs/Celltype.done"
            ],
                Species = Species.lower(),
                refbuiltin = refbuiltin,
                annolevel = annolevel
            )
        )

    ##==================================================================================================================
    if config["module"]["Diffexp"]:
        wanted_input.extend(
            expand([
                ##run diffexp
                # "result/6.Diffexp/{diff_group}-all_diffexp_genes.xls",
                "result/7.Diffexp/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.xls",
                ##run diff_heatmap
                "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_heatmap.pdf",
                "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_heatmap.png",
                "result/7.Diffexp/Heatmap/top{topn}_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}_genes.xls",
                # ##run enrich
                "logs/diffexp/run.snakemake.enrich_{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.log",
                ## run diff_ppi
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.pdf",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.png",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.new_colors.svg",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.pdf",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.png",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.svg",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.string_protein-protein-interaction.tsv",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.mapping_results.tsv",
                "result/7.Diffexp/ppi/{diff_group}-diff-{pval_or_padj}-{sgin_thread}-FC-{foldchange}.top_25_ppi_network.png"
            ], 
            diff_group = diff_groups.index.to_list(),
            sgin_thread = padj if (padj) else pvalue,
            pval_or_padj = "padj" if (padj) else "pval",
            foldchange = foldchange,
            topn = topn
            )
        )

    ##==================================================================================================================

    if config["module"]["Report"]:
        wanted_input.extend(
            [f"result/report/{project}_Report/Report.html"]            
            )

    # print(wanted_input)
    return wanted_input
