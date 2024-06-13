###########################################################
## ==================== 1. load config ====================
###########################################################
script_dir = "/public/scRNA_works/pipeline/oesinglecell3/exec"
#script_dir = "scripts/R"
# script_dir = "/home/dongjiaoyang/oesinglecell3/exec" # 脚本合并测试

# Reference_celltype
species = config["report"]["Species"]
if not isinstance(species, str) or "-" not in species:
    #if species in ["人","猪","猴","鸡"]:
    #   celltypingdb = "/data/database/celltype_refdata/logNorm_rds/hpca.rds"
    #elif species in ["小鼠","大鼠"]:
    #    celltypingdb = "/data/database/celltype_refdata/logNorm_rds/immgen.rds"
    #elif species in ["拟南芥"]:
    #   celltypingdb = "/data/database/celltype_refdata/Ath_wangjiawei_celltype_genename.rds"
    print("Error:Please fill in the sample type in the correct format.")
    exit(0)
else:
    if species.split("-")[0] in ["人"]:
        celltype_dict = {
            "血液": "/data/database/celltype_refdata/logNorm_rds/blueprint_encode.rds",
            "schcl": "/data/database/celltype_refdata/logNorm_rds/schcl.rds",
            "乳腺癌": "/gpfs/oe-database/celltype_refdata/report/GSE114727.rds",
            "肝内胆管癌": "/gpfs/oe-database/celltype_refdata/report/CHOL_GSE125449.rds",
            "结肠癌": "/gpfs/oe-database/celltype_refdata/report/GSE146771.rds",
            "胶质瘤": "/gpfs/oe-database/celltype_refdata/report/GSE131928.rds",
            "肝癌": "/gpfs/oe-database/celltype_refdata/report/LIHC_GSE125449.rds",
            "鼻咽癌": "/gpfs/oe-database/celltype_refdata/report/GSE150430.rds",
            "子宫内膜": "/gpfs/oe-database/celltype_refdata/report/GSE111976.rds"
        }

        if any(x in species for x in celltype_dict.keys()):
            celltypingdb = celltype_dict[species.split("-")[1]]
        else:
            celltypingdb = "/data/database/celltype_refdata/logNorm_rds/hpca.rds"

    elif species.split("-")[0] in ["小鼠"]:
        celltype_dict = {
            "scmca": "/data/database/celltype_refdata/logNorm_rds/scmca.rds",
            "肝脏": "/data/database/celltype_refdata/logNorm_rds/mouse.rnaseq.rds",
            "心脏": "/data/database/celltype_refdata/logNorm_rds/mouse.rnaseq.rds",
            "小脑": "/data/database/celltype_refdata/logNorm_rds/mouse.rnaseq.rds",
            "嗅球": "/data/database/celltype_refdata/logNorm_rds/mouse.rnaseq.rds",
            "干细胞": "/gpfs/oe-database/celltype_refdata/report/GSE116530.rds",
            "脑": "/gpfs/oe-database/celltype_refdata/report/allen_mouse_cortex_smart_v3.rds",
            "肾脏": "/gpfs/oe-database/celltype_refdata/report/GSE107585.rds"
        }

        if any(x in species for x in celltype_dict.keys()):
            celltypingdb = celltype_dict[species.split("-")[1]]
        else:
            celltypingdb = "/data/database/celltype_refdata/logNorm_rds/immgen.rds"

    elif species.split("-")[0] in ["猪","猴","鸡"]:
        celltypingdb = "/data/database/celltype_refdata/logNorm_rds/hpca.rds"

    elif species.split("-")[0] in ["大鼠"]:
        celltypingdb = "/data/database/celltype_refdata/logNorm_rds/immgen.rds"

    elif species.split("-")[0] in ["拟南芥"]:
        celltypingdb = "/gpfs/oe-database/celltype_refdata/report/Ath_PMID_31004836.rds"
    else:
        print("Special species do not require cell type identification.")
        exit(0)


if any(x in species for x in ["人", "猪", "猴", "鸡"]):
    ref_species = "Human"
else:
    ref_species = "Mouse"


## ==================== 1.1 load sample.csv ====================
samples = pd.read_csv(config["report_params"]["samples_file"], dtype=str).set_index("sampleid", drop=False)
samples.index.names = ["sample_id"]
redo_flag = False
for s in samples.index:
    if samples.loc[s, "species"] != config["report"]["Species"]:
        samples.loc[s, "species"] = config["report"]["Species"]
        redo_flag = True
    # if samples.loc[s, "fastq"] != f'raw_data/{s}/{s}_fastqs':
    #     samples.loc[s, "fastq"] = f'raw_data/{s}/{s}_fastqs'
    #     redo_flag = True
    #h5="%s/result/cellranger/%s/outs/molecule_info.h5" %( os.getcwd(),s)
if redo_flag:
    samples.to_csv(config["report_params"]["samples_file"],encoding='utf-8',index=False)
#validate(samples, schema="../schemas/samples.schema.yaml")
wildcard_constraints:
                    sample="|".join(samples.index)
## ==================== 2.2 load diffexp group ====================
if config["report_params"]["module"]["diffexp"] and os.path.exists(config["report_params"]["diffexp_file"]):
    diff_groups = pd.read_table(config["report_params"]["diffexp_file"], dtype=str)
    diff_groups.index=diff_groups["type"]+"_"+diff_groups["treatment"]+"-vs-"+diff_groups["control"]
    wildcard_constraints:
                    diff_group="|".join(diff_groups.index)

## MT list exist
MT_gene = False
if os.path.exists(f"{config['report']['Reference']}/MT_genelist.gmt") :
    with open(f"{config['report']['Reference']}/MT_genelist.gmt") as gmt:
        MT_gene = any([1  if line.find("percent.mito") == 0 else 0 for line in gmt])
## HB list exist
HB_gene = False
if os.path.exists(f"{config['report']['Reference']}/HB_genelist.gmt") :
    with open(f"{config['report']['Reference']}/HB_genelist.gmt") as gmt:
        HB_gene = any([1  if line.find("percent.HB") == 0 else 0 for line in gmt])
############################################################################
## ==================== 2. define functions & all input ====================
############################################################################
def get_sample(wildcards):
    return wildcards.sample

def get_sample_type(wildcards):
    """Get  sample_type  from sheet."""
    sample_type = samples.loc[(wildcards.sample), ["species"]].dropna()
    return sample_type.species

def get_pointsize_params(wildcards):
    if config['report_params']['pointsize'] is None:
        return ""
    else:
        return f"--pointsize {config['report_params']['pointsize']}"

def get_roboust_linear_model_params(wildcards):
    if config['report_params']['roboust_linear_model'] is None:
        return ""
    else:
        return f"--rlm.xy {config['report_params']['roboust_linear_model']}"

def all_input_report(wildcards):
    wanted_input = []
    #if config["report_params"]["module"]["cellranger_aggr"]:
    #    wanted_input.extend(
    #            [
    #                    "result/cellranger/aggr/outs/web_summary.html",
    #                    "result/cellranger/aggr/outs/count/cloupe.cloupe"
    #                    ])
    if config["report_params"]["module"]["create"]:
        wanted_input.extend(
            [
                "result/create/seurat.h5seurat"
            ])
    if config["report_params"]["module"]["qc"]:
        wanted_input.extend(
            [
                "result/Count_QC/QC_metrics_beforeQC.png",
                "result/Count_QC/QC_metrics_beforeQC.pdf",
                "result/Count_QC/QC_metrics_afterQC.png",
                "result/Count_QC/QC_metrics_afterQC.pdf",
                "result/Count_QC/cell_statitics_before_after_QC.xls",
                "result/Count_QC/filtered.h5seurat"
            ])
    if config["report_params"]["module"]["clustering"]:
        wanted_input.extend(
            expand([
                "result/Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
                "result/Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf",
                "result/Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
                "result/Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv"
                ],
                    reduct2=config["report_params"]["reduct2"],
                    reduct1=config["report_params"]["reduct1"],
                    resolution=config["report_params"]["resolution"]
            )
        )
        wanted_input.extend([
            "result/Clustering/clusters_correlation/normalized_data_groupby_clusters.xls",
            "result/Clustering/clusters_correlation/coefficient_heatmap.png",
            "result/Clustering/clusters_correlation/coefficient_heatmap.pdf"])
        if config["report"]["Sample_Num"] > 1:
            wanted_input.extend([
                    "result/Clustering/visualize_cluster_by_clusters/clust_cond_freq_info.xls",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-clusters_summary_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-clusters_summary_plot.png",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-group_contrast_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-group_contrast_plot.png",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-group_summary_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-group_summary_plot.png",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_contrast_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_contrast_plot.png",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_summary_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/groupby-sampleid_summary_plot.png",
                    "result/Clustering/visualize_cluster_by_clusters/splitby-group_split_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/splitby-group_split_plot.png",
                    "result/Clustering/visualize_cluster_by_clusters/splitby-sampleid_split_plot.pdf",
                    "result/Clustering/visualize_cluster_by_clusters/splitby-sampleid_split_plot.png"]
            )
    if config["report_params"]["module"]["marker"]:
        wanted_input.extend(
            [
                "result/Marker/all_markers_for_each_cluster.xls",
                "result/Marker/top10_markers_for_each_cluster.xls",
                "result/Marker/all_markers_for_each_cluster_anno.xls",
                "result/Marker/top10_markers_for_each_cluster_anno.xls",
                "result/Marker/topmarker_gene_heatmap.png",
                "result/Marker/topmarker_gene_heatmap.pdf",
                "result/Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
                "result/Marker/markers_vis4cluster1/marker_gene_featureplot.png",
                "result/Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf",
                "result/Marker/markers_vis4cluster1/marker_gene_violin_plot.png"
            ])
    if config["report_params"]["module"]["celltyping"]:
        wanted_input.extend(
            expand([
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_statistics.xls",
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_results.xls",
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_simplified_celltyping_results.csv",
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_heatmap.png",
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_heatmap.pdf",
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_plot.png",
                "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_plot.pdf",
                "result/Reference_celltype/{ref_species}ref_{db_name}_top.main_celltyping_plot.png",
                "result/Reference_celltype/{ref_species}ref_{db_name}_top.main_celltyping_plot.pdf"
                ],
                  #celltypingdb=config["report_params"]["celltypingdb"],
                  db_name=celltypingdb.split('/')[-1].split('.')[0],
                  ref_species=ref_species
            )
        )
    if config["report_params"]["module"]["diffexp"]:
        for diff_group in diff_groups.index:
            wanted_input.extend(
                expand([
                    "result/Diffexp/{diff_group}-all_diffexp_genes.xls",
                    "result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5.xls",
                    "result/Diffexp/{diff_group}-all_diffexp_genes_anno.xls",
                    "result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5_anno.xls",
                    "result/Diffexp/enrichment/KEGG_enrichment/enrichment_kegg.xls",
                    "result/Diffexp/enrichment/GO_enrichment/enrichment_go.xls",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.new_colors.pdf",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.new_colors.png",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.new_colors.svg",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.pdf",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.png",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.svg",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.tsv",
					"result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction-top25.tsv",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.ppi_network.tsv",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.top_25_ppi_network.pdf",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.top_25_ppi_network.png",
                    "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.mapping_results.tsv",
					"result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.mapping_results-top25.tsv",
                    "result/Diffexp/top20_{diff_group}_genes.xls",
                    "result/Diffexp/top20_{diff_group}_heatmap.png",
                    "result/Diffexp/top20_{diff_group}_heatmap.pdf"],
                diff_group=diff_group
                )
            )
        wanted_input.extend(["logs/diffexp_fc12.check"])
    if config["report_params"]["module"]["report"]:
        wanted_input.extend([f"result/report/{config['report']['Task_Num'].split('-')[0]}_Report_{time.strftime('%Y_%m_%d')}/分析说明.html",
        "logs/upload_report_sor"])
        # wanted_input.extend(["result/rds/data_ob_v3.rds"])
    if config["report_params"]["module"]["report_upload"]:
        wanted_input.extend([
            f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}.zip",
            f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}.zip.md5",
            "result/report/report_upload.check",
            #"logs/upload_report_sor"
        ]) 
    return wanted_input

###################################################################
## ==================== 3. include other rules ====================
###################################################################
#if config["report_params"]["module"]["cellranger_aggr"]:
#  include: "1.cellranger_aggr.smk"
if config["report_params"]["module"]["create"]:
    include: "2.create.smk"
if config["report_params"]["module"]["qc"]:
    include: "3.qc.smk"
if config["report_params"]["module"]["clustering"]:
    include: "4.clustering.smk"
if config["report_params"]["module"]["marker"]:
    include: "5.marker.smk"
if config["report_params"]["module"]["celltyping"]:
    include: "6.celltyping.smk"
if config["report_params"]["module"]["diffexp"]:
    include: "7.diffexp.smk"
if config["report_params"]["module"]["report"] or config["report_params"]["module"]["report_upload"]:
    include: "8.report.smk"
