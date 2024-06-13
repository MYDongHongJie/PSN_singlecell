from snakemake.utils import validate
import pandas as pd
import os,time


#exec(open('/data/software/modules/modules-v4.2.1/init/python.py').read())
#os.environ['MODULEPATH'] =os.path.abspath("./rules/envs/")
#print("Using MODULEPATH "+os.environ['MODULEPATH'])

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# container: "docker://continuumio/miniconda3"
#=======================================================================================================================
##### load config and sample sheets #####

configfile: "config/config.yaml"
configfile: "config/cluster.yaml"
envmodules = config['envmodules']
#validate(config, schema="../config/schemas/config.schema.yaml")
samples = pd.read_csv(config["samples"], dtype=str).set_index("sampleid", drop=False)
samples.index.names = ["sample_id"]
#validate(samples, schema="../config/schemas/samples.schema.yaml")

if config["module"]["diff_gene"] and os.path.exists(config["diff_DEG_group"]):
    diff_groups = pd.read_csv(config["diff_DEG_group"], dtype=str,delimiter="\t")
    diff_groups.index= diff_groups["type"]+"_"+diff_groups["treatment_name"]+"-vs-"+diff_groups["control_name"]
    diff_groups['groupname'] = diff_groups.index
    diff_groups['foldchange']=  config['params']['foldchange']
    diff_groups['pvalue'] =  config["params"]['qvalue'] if (config["params"]['qvalue']) else config["params"]['pvalue']
    diff_groups['pval_or_padj']= "qvalue" if(config["params"]['qvalue']) else "pvalue"
    wildcard_constraints:
                    diff_group="|".join(diff_groups.index)

groupname = lambda w: diff_groups.loc["{}".format(w.diff_group)]["groupname"]
pvalue = lambda w: diff_groups.loc["{}".format(w.diff_group)]["pvalue"]
foldchange = lambda w: diff_groups.loc["{}".format(w.diff_group)]["foldchange"]
pval_or_padj= lambda w: diff_groups.loc["{}".format(w.diff_group)]["pval_or_padj"]


reference = config["database"]["reference"]
project = config["report"]["Task_Num"]
probe_set=config["database"]["probe_set"]
project_num = config['report']['Project_Num']
project_end_time = config['report']['Project_End_Time']
#report_file = "result/report/"+ project_num+"_Report_"+project_end_time+"/report.html"

celltype=config["params"]["celltype"]
##### wildcard constraints #####
wildcard_constraints:
                    sample="|".join(samples.index)
#=======================================================================================================================
####### helpers ###########
def get_sample(wildcards):
    return wildcards.sample

def get_sample_type(wildcards):
    """Get  sample_type  from sheet."""
    sample_type = samples.loc[wildcards.sample, ["species"]].dropna()
    return sample_type.species

#### for HE json_file
def get_sample_image_json_file(wildcards):
    """Get  sample_jsaddress  from sheet."""
    if "image_json_file"  in samples.columns   :
        sample_json_file = samples.loc[wildcards.sample, ["image_json_file"]]
        if not pd.isna(sample_json_file.values):
            if os.path.exists( str(sample_json_file.image_json_file) ):
                he_json_parameter = f'--loupe-alignment={sample_json_file.image_json_file}'
            else:
               print(f"The {sample_json_file} is not exist ,please check !!!! ")
        else:
            he_json_parameter = f" "
    else:
        he_json_parameter = f" "
    return he_json_parameter

##for fastqc
def get_individual_fastq(wildcards):
    if wildcards.reads == "0" or wildcards.reads == "1":
        return samples.loc[wildcards.sample, "fq1"]
    elif wildcards.reads == "2":
        return samples.loc[wildcards.sample, "fq2"]
################################################ for spaceranger########################################################
def get_fastqs(wildcards):
    """Get raw FASTQ files directory from sheet."""
    fastqs= os.path.abspath((samples.loc[wildcards.sample, ["fastq"]].dropna()).iloc[0])
    return fastqs

def get_slide(wildcards):
    """Get  slide id  from sheet."""
    slide = samples.loc[wildcards.sample, ["slide"]].dropna()
    return slide.slide

def get_slidearea(wildcards):
    """Get  slide area from sheet."""
    slidearea = samples.loc[wildcards.sample, ["slide_area"]].dropna()
    return slidearea.slide_area

def get_slide_file(wildcards):
    """Get  slide file from sheet's slideid."""
    slide = samples.loc[wildcards.sample, ["slide"]].dropna()
    slide_file =os.path.abspath(f'raw_data/{slide.slide}.gpr')
    return slide_file

def get_slide_image(wildcards):
    """Get  slide image file  from sheet."""
    image = os.path.abspath((samples.loc[wildcards.sample, ["image"]].dropna()).iloc[0])
    if config["params"]["library_type"] == "cytassist":
        sample_image_file = samples.loc[wildcards.sample, ["image"]]
        if not pd.isna(sample_image_file.values):
            slide_image = f"--image={image}"
        else:
            slide_image = ""
    else:
        slide_image=f"--image={image}"
    return slide_image

def get_cytaimage_image(wildcards):
    """Get  cytaimage  image file  from sheet."""
    cytaimage = os.path.abspath((samples.loc[wildcards.sample, ["cytaimage"]].dropna()).iloc[0])
    return cytaimage

def get_origin_image(wildcards):
    """Get origin image  image file  from sheet."""
    origin_image = os.path.abspath((samples.loc[wildcards.sample, ["image"]].dropna()).iloc[0])
    return origin_image

################################################ all output file check##################################################
def all_input(wildcards):
    wanted_input = []
    #===================================================================================================================
    # spaceranger
    if config["module"]["spaceranger"]:
        for (sample) in samples.index:
            wanted_input.extend(
                expand(
                    [
                        "result/spaceranger/{sample}/outs/web_summary.html",
                        "result/spaceranger/{sample}/outs/metrics_summary.csv",
                        "result/spaceranger/{sample}/outs/cloupe.cloupe",
                        "result/spaceranger/{sample}/outs/raw_feature_bc_matrix.h5",
                        "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5",
                        "result/spaceranger/{sample}/outs/{sample}.png",
                        "result/spaceranger/{sample}/outs/spatial/tissue_lowres_image.png",
                        "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
                    ],
                    sample=sample
                )
            )

    if config["module"]["spaceranger_email"]:
        wanted_input.extend(
                    expand(
                        [
                            f"result/report/{project}-QC_spaceranger_report/summary.csv",
                            f"result/report/{project}-QC_spaceranger_report/summary.csv.email.html",
                            "logs/spaceranger_report.success.log",
                           # "logs/QC_upload_sor/QC_upload_sor.success.log"
                        ],project=project))

    if config["params"]["generate_bam"] and config['params']['bam_upload_rm']:
        for (sample) in samples.index:
            wanted_input.extend(
                expand(
                    [
                        "logs/{sample}_upload_bam_check"
                    ],
                    sample=sample
                )
            )

    if config["module"]["spaceranger_upload"]:
        wanted_input.extend([
            "logs/spaceranger_upload.check"
                ])

    if config["module"]["spaceranger_report"]:
        wanted_input.extend([
            "logs/spaceranger_QC_report.txt"
            ])
        return wanted_input

    #===================================================================================================================
    #spaceranger_aggr
    if config["module"]["spaceranger_aggr"]:
        wanted_input.extend(
                [
                        "result/spaceranger/aggr/outs/web_summary.html",
                        "result/spaceranger/aggr/outs/summary.json",
                        "result/spaceranger/aggr/outs/cloupe.cloupe",
                       # "result/spaceranger/aggr/outs/raw_feature_bc_matrix.h5",
                        "result/spaceranger/aggr/outs/filtered_feature_bc_matrix.h5"
                        ])
    #===================================================================================================================
    #seurat对象创建
    if config["module"]["sctools"]:
        wanted_input.extend(
            [ "result/sctools/spatial.rds"]
            #expand(["result/sctools/spatial.{object_formart}"],object_formart=config["params"]["object_formart"])
            )
    #===================================================================================================================
    ##质控
    if config["module"]["count_qc"]:
       wanted_input.extend(
            expand(
                [  #os.path.join(config['database']['reference'],"MT_genelist.txt"),
                    "result/count_qc/statitics_for_QC.xls",
                    "result/count_qc/filtered_seurat.rds",
                    "result/count_qc/QC_featureplot_for_nCount_Spatial.pdf",
                    "result/count_qc/QC_featureplot_for_nCount_Spatial.png",
                    "result/count_qc/QC_featureplot_for_nFeature_Spatial.pdf",
                    "result/count_qc/QC_featureplot_for_nFeature_Spatial.png"
                ]
                # sample=sample
            )
       )
    #===================================================================================================================
    #降维聚类
    if config["module"]["clustering"]:
        for sample in samples.index:
            wanted_input.extend(
                expand(
                    [
                        "result/cluster_seurat/{reduct2_method}_Dimension_Reduction/{reduct2_method}_Dimension_Reduction_coordination.csv",
                        "result/cluster_seurat/{reduct2_method}_Dimension_Reduction/{reduct2_method}_groupby_cluster_resolution{resolution}_plot.pdf",
                        "result/cluster_seurat/{reduct2_method}_Dimension_Reduction/{reduct2_method}_groupby_cluster_resolution{resolution}_plot.png",
                        ##"result/cluster_seurat/{reduct1_method}_Dimension_Reduction/{reduct1_method}_Dimension_Reduction_coordination.csv",
                        "result/cluster_seurat/singlecell_object.clustering_resolution{resolution}.rds"
                    ],
                    sample=sample,
                    reduct2_method=config["params"]["reduct2_method"],
                    reduct1_method=config["params"]["reduct1_method"],
                    resolution=config["params"]["resolution"]
                )
            )
            wanted_input.extend(
                expand(
                    [
                        "result/cluster_seurat/visualize_cluster_by_clusters/clust_cond_freq_info.xls",
                        "result/cluster_seurat/visualize_cluster_by_clusters/groupby-group_resolution{resolution}_contrast_plot.pdf",
                        "result/cluster_seurat/visualize_cluster_by_clusters/groupby-group_resolution{resolution}_contrast_plot.png",
                        "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution{resolution}_contrast_plot.pdf",
                        "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution{resolution}_contrast_plot.png",
                        "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution-{resolution}_summary_plot.pdf",
                        "result/cluster_seurat/visualize_cluster_by_clusters/groupby-sampleid_resolution-{resolution}_summary_plot.png",
                        "result/cluster_seurat/visualize_cluster_by_clusters/splitby-group_resolution{resolution}_split_plot.pdf",
                        "result/cluster_seurat/visualize_cluster_by_clusters/splitby-group_resolution{resolution}_split_plot.png",
                        "result/cluster_seurat/visualize_cluster_by_clusters/splitby-sampleid_resolution{resolution}_split_plot.pdf",
                        "result/cluster_seurat/visualize_cluster_by_clusters/splitby-sampleid_resolution{resolution}_split_plot.png"
                    ],
                    sample=sample,
                    resolution=config["params"]["resolution"]
                )
            )
    #===================================================================================================================
    ##marker基因生成
    if config["module"]["find_marker"]:
        wanted_input.extend(
                    expand(
                        [   # "result/marker/all_markers_for_each_cluster.xls",
                            # "result/marker/top10_markers_for_each_cluster.xls",
                            "result/marker/topmarker_gene_heatmap.png",
                            "result/marker/topmarker_gene_heatmap.pdf",
                            "result/marker/all_markers_for_each_cluster_anno.xls",
                            "result/marker/top10_markers_for_each_cluster_anno.xls",
                            "result/marker/seurat_object_markers.rds"
                        ],
                        sample=sample
                    )
                )
    #===================================================================================================================
    ##细胞类型
    if config["module"]["celltype"]:
        wanted_input.extend(
                    expand(
                        [
                            "result/celltype/spatial_seurat_celltype.rds",
                            "result/celltype/combined_all_slice-celltype.xls",
                            "result/celltype/2.celltype_interaction/combined_all_slice-interactions.xls",
                            "result/celltype/2.celltype_interaction/all_celltype_interaction.pdf",
                            "result/celltype/2.celltype_interaction/all_celltype_interaction.png"
                        ],
                        sample=sample
                    )
                )
    #===================================================================================================================
    ## 差异基因
    if config["module"]["diff_gene"]:
        for groupname in diff_groups.index:
            wanted_input.extend(
                expand(
                    [
                        #"result/diffexp/{groupname}-all_diffexp_genes.xls",
                        #"result/diffexp/{groupname}-all_diffexp_genes_anno.xls",
                        #"result/diffexp/{groupname}-all_diffexp_genes-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls",
                        "result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.xls",
                        #"result/diffexp/{groupname}-all_diffexp_genes-{pval_or_padj}-{pvalue}-FC-{foldchange}_anno.xls",
                        "result/diffexp/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}_anno.xls",
                        "result/diffexp/enrichment/GO_enrichment/{groupname}/enrichment-go-{groupname}-Total.xls",
                        "result/diffexp/enrichment/KEGG_enrichment/{groupname}/enrichment-kegg-{groupname}-Total.xls"
                    ],
                    groupname = groupname ,
                    pval_or_padj="qvalue" if (config["params"]['qvalue']) else "pvalue",
                    pvalue=config["params"]['qvalue'] if (config["params"]['qvalue']) else config["params"]['pvalue'],
                    foldchange=config['params']["foldchange"]
                )
            )
        #===================================================================================================================
        ## 差异基因PPI
        if config["module"]["diff_gene_ppi"]:
            for groupname in diff_groups.index:
                wanted_input.extend(
                    expand(
                        [
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.new_colors.pdf",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.new_colors.png",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.new_colors.svg",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.pdf",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.png",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.svg",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.string_protein-protein-interaction.tsv",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.mapping_results.tsv",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.top_25_ppi_network.png",
                            "result/diffexp/ppi/{groupname}-diff-{pval_or_padj}-{pvalue}-FC-{foldchange}.top_25_ppi_network.pdf"
                        ],
                        groupname=groupname,
                        pval_or_padj="qvalue" if (config["params"]['qvalue']) else "pvalue",
                        pvalue=config["params"]['qvalue'] if (config["params"]['qvalue']) else config["params"]['pvalue'],
                        foldchange=config['params']["foldchange"]
                    )
                )
    #===================================================================================================================
    ## pathseq分析
    if config["module"]["pathseq"]:
        for sample in samples.index:
            wanted_input.extend(
                expand(
                    [
                        "result/pathseq/GATK/{sample}.pathseq.complete.bam",
                        "result/pathseq/GATK/{sample}.pathseq.complete.csv",
                        "result/pathseq/pathseq_visualization/{sample}.visium.readname",
                        "result/pathseq/pathseq_visualization/{sample}.visium.unmap_cbub.bam",
                        "result/pathseq/pathseq_visualization/{sample}.visium.unmap_cbub.fasta",
                        "result/pathseq/pathseq_visualization/{sample}.visium.list",
                        "result/pathseq/pathseq_visualization/{sample}.visium.raw.readnamepath",
                        "result/pathseq/pathseq_visualization/{sample}.visium.genus.cell",
                        "result/pathseq/pathseq_visualization/{sample}.visium.genus.csv",
                        "result/pathseq/pathseq_visualization/{sample}.visium.validate.csv",
                        "result/pathseq/pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution/summary_statistic.xls",
                        f"result/report/{config['report']['Project_Num']}_PathSeq_Report/report.html"
                    ],
                    sample=sample
                )
            )
    return wanted_input


##======================================================================================================================
##======================================================================================================================
def report_input(wildcards):
    wanted_input = []
    if config["module"]["report"]:
        wanted_input.extend(
            [
                "logs/sor/QC_upload_sor.check.done",
                f"result/report/{config['report']['Project_Num']}_Report_{config['report']['Project_End_Time']}/分析说明.html"
            ])
        return wanted_input
