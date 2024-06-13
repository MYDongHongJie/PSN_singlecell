# -*-coding:utf-8-*-
import pandas as pd
import numpy as np
import os
import sys
import requests
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger

logger = getLogger('oe.cloud.sc.qsub')

from taskserver.tools.module_cmd import module_cmd
#===================change all the input into standard str form===========================================================================
def params_setting(d):
    params = []
    for key in d:
        param = f'{key} {d[key]}'
        params.append(param)
    return params

#===================input SCENIC input and run the cmds===================================================================================
def task_SCENIC(input,projectid="projectid", taskid="taskid", workdir="workdir"):
    wkdir  = workdir
    if not os.path.exists(f"{wkdir}/output/download"):
        os.makedirs(f"{wkdir}/output/download")
    if not os.path.exists(f"{wkdir}/SCENIC"):
        os.makedirs(f"{wkdir}/SCENIC")
    environment = "OESingleCell/3.0.d"
    # step 1 scenic.R required params
    scenic_param = {}
    # step 2 RunRAS-RSS.R required params
    RunRAS_RSS_param = {}
    # step 3 RunCSI.R required params
    RunCSI_param = {}
    # All three scripts share the same input, so we add input param to all three directory
    #if input.loc["parameters", "input"] != "":
    #    input_rds = input.loc["parameters", "input"]
    #else:
    input_rds = f"{wkdir}/../{input.loc['base', 'tasks'][0]['taskId']}/bclust/data_ob_v3.rds"
    scenic_param["--input"] = input_rds
    #==============================setting parameters for scenic.R script================================================================
    # define the input format for scenic.R
    if input.loc["parameters", "informat"] != "":
        informat = input.loc["parameters", "informat"]
    else:
        informat = "seurat"
    scenic_param["--informat"] = informat
    # define cisTargetdb for scenic.R
    scenic_param["--cisTargetdb"] = input.loc["parameters", "cisTargetdb"]
    # define species for scenic.R
    scenic_param["--species"] = input.loc["parameters", "species"]
    # define minCell4gene for scenic.R default is 0.01
    if input.loc["parameters", "minCell4gene"] != "":
        minCell4gene = float(input.loc["parameters", "minCell4gene"])
    else:
        minCell4gene = 0.01
    scenic_param["--minCell4gene"] = minCell4gene
    # define tfdb for scenic.R default this tgdb is not set
    if input.loc["parameters", "tfs"] != "":
        tf_db = input.loc["parameters", "tfs"]
        scenic_param["--tfs"] = tf_db
    # define coexMethod for scenic.R, w001,w005,top50,top50perTarget,top10perTarget,top5perTarget can be cosen
    scenic_param["--coexMethod"] = input.loc["parameters", "coexMethod"]
    # define CPUs for scenic.R
    if input.loc["parameters", "scenic_ncores"] != "":
        scenic_ncores = int(input.loc["parameters", "scenic_ncores"])
    else:
        scenic_ncores = 10
    scenic_param["--ncores"] = scenic_ncores
    # define downsampling number of cells for scenic.R
    if input.loc["parameters", "downsample"] != "":
        downsample = input.loc["parameters", "downsample"]
    else:
        downsample = "30000"
    scenic_param["--downsample"] = downsample
    # whether to use hvg in scenic.R
    if input.loc["parameters", "hvg"] != "":
        hvg = input.loc["parameters", "hvg"]
    else:
        hvg = "F"
    scenic_param["--hvg"] = hvg
    # whether to extend the searching area for regulon in scenic.R,and this param is also used in RunRAS_RSS.R and RunCSI.R
    if input.loc["parameters", "extended"] != "":
        extended = input.loc["parameters", "extended"]
    else:
        extended = "F"
    scenic_param["--extended"] = extended
    #==============================setting parameters for RunRAS_RSS.R script================================================================
    RunRAS_RSS_param["--input"] = scenic_param["--input"]
    # set the path for the auc matrix input
    #if input.loc["parameters", "auc"] != "":
    #    auc = input.loc["parameters", "auc"]
    #else:
    auc = f"{wkdir}/SCENIC/int/3.4_regulonAUC.Rds"
    # set auc matrix input format, default is rds
    RunRAS_RSS_param["--auc"] = auc
    if input.loc["parameters", "aucformat"] != "":
        aucformat = input.loc["parameters", "aucformat"]
    else:
        aucformat = "rds"
    RunRAS_RSS_param["--aucformat"] = aucformat
    # whether to use top genes
    if input.loc["parameters", "topGenes"] != "":
        topGenes = int(input.loc["parameters", "topGenes"])
        RunRAS_RSS_param["--topGenes"] = topGenes
    # set groupby
    RunRAS_RSS_param["--groupby"] = input.loc["parameters", "groupby"]
    # set binary method
    if input.loc["parameters", "binmethod"] != "":
        binmethod = input.loc["parameters", "binmethod"]
    else:
        binmethod = "aucell"
    RunRAS_RSS_param["--binmethod"] = binmethod
    # define CPUs for RunRAS_RSS.R
    if input.loc["parameters", "RAS_RSS_ncores"] != "":
        RAS_RSS_ncores = int(input.loc["parameters", "RAS_RSS_ncores"])
    else:
        RAS_RSS_ncores = 10
    RunRAS_RSS_param["--ncores"] = RAS_RSS_ncores
    # set cutoff for AUC activity
    RunRAS_RSS_param["--threshold"] = int(input.loc["parameters", "threshold"])
    RunRAS_RSS_param["--extended"] = scenic_param["--extended"]
    if input.loc["parameters", "predicate"] != "":
        RunRAS_RSS_param["--predicate"] = input.loc["parameters", "predicate"]
#==============================setting parameters for RunCSI.R script================================================================
    RunCSI_param["--input"] = scenic_param["--input"]
    RunCSI_param["--auc"] = RunRAS_RSS_param["--auc"]
    RunCSI_param["--aucformat"] = RunRAS_RSS_param["--aucformat"]
    RunCSI_param["--groupby"] = RunRAS_RSS_param["--groupby"]
    if input.loc["parameters", "nclust"] != "":
        nclust = int(input.loc["parameters", "nclust"])
    else:
        nclust = 5
    RunCSI_param["--nclust"] = nclust
    RunCSI_param["--extended"] = scenic_param["--extended"]
    if RunRAS_RSS_param.get("--predicate") != None:
        RunCSI_param["--predicate"] = RunRAS_RSS_param["--predicate"]
#====================setting output path==============================================================================================
    scenic_param["--output"] = f"{wkdir}/SCENIC/"
    RunRAS_RSS_param["--output"] = f"{wkdir}/SCENIC/byclusters"
    RunCSI_param["--output"] = f"{wkdir}/SCENIC/byclusters"
#====================setting command line=============================================================================================
    scenic_param_line = ' '.join(params_setting(scenic_param))
    scenic_cmd = f"Rscript /public/scRNA_works/works/guochy/SCENIC_CloadPlatform_pipeline_test/scenic.R {scenic_param_line}"
    RunRAS_RSS_param_line = ' '.join(params_setting(RunRAS_RSS_param))
    RunRAS_RSS_cmd = f"Rscript /public/scRNA_works/works/guochy/SCENIC_CloadPlatform_pipeline_test/RunRAS-RSS.R {RunRAS_RSS_param_line}"
    RunCSI_param_line = ' '.join(params_setting(RunCSI_param))
    RunCSI_cmd = f"Rscript /public/scRNA_works/works/guochy/SCENIC_CloadPlatform_pipeline_test/RunCSI.R {RunCSI_param_line}"
#=============Enter wkdir and start analysis=========================================================================================
    set_wkdir_cmd = f"cd {wkdir}/SCENIC/"
    start_scenic_analysis = f"module purge && source /home/lipeng/miniconda3/bin/activate Scenic && {set_wkdir_cmd} && {scenic_cmd}"
    logger.info(start_scenic_analysis)
    with module_cmd(environment) as p:
        status = p(start_scenic_analysis, projectid, taskid)
    logger.info("SCENIC analysis finished!")
    start_ras_rss_analysis = f"module purge && source /home/lipeng/miniconda3/bin/activate Scenic && {set_wkdir_cmd} && {RunRAS_RSS_cmd}"
    logger.info(start_ras_rss_analysis)
    with module_cmd(environment) as p:
        status = p(start_ras_rss_analysis, projectid, taskid)
    logger.info("RAS_RSS analysis finished!")
    start_csi_analysis = f"module purge && source /home/lipeng/miniconda3/bin/activate Scenic && {set_wkdir_cmd} && {RunCSI_cmd}"
    logger.info(start_csi_analysis)
    with module_cmd(environment) as p:
        status = p(start_csi_analysis, projectid, taskid)
    logger.info("CSI analysis finished!")
# ============generate output.json.tsv so that the diagrams can be uploaded to obs====================================================
    d = {"task_type":["SCENIC"]*12,"result_module":["diagram"]* 6 + ["data"] * 6,
    "input":["1.1.regulon_activity_heatmap_groupby_cells_data.tsv", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "2.1.regulon_RSS_annotation.tsv", "3.2.regulon_csi_correlation_heatmap_data.tsv", "3.3.csi_module_activity_heatmap_data.tsv", "0.1.TF_target_enrichment_annotation.tsv", "0.2.regulon_annotation.tsv", "0.3.MotifEnrichment_preview.html", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "3.1.csi_module_annotation.tsv"],
    "type":["heatmap","heatmap","rankingplot", "heatmap", "heatmap", "heatmap"] + ["tsv", "tsv", "html"] + ["tsv"] * 3,
    "file":["1.1.regulon_activity_heatmap_groupby_cells_data.tsv", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "2.1.regulon_RSS_annotation.tsv", "3.2.regulon_csi_correlation_heatmap_data.tsv", "3.3.csi_module_activity_heatmap_data.tsv", "0.1.TF_target_enrichment_annotation.tsv", "0.2.regulon_annotation.tsv", "0.3.MotifEnrichment_preview.html", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "3.1.csi_module_annotation.tsv"],
    "title":["cell_regulon_activity_heatmap", "cluster_regulon_activity_heatmap", "RSS_ranking", "RSS_heatmap", "regulons_csi_correlation", "csi_module_activity", "TF_target_enrichment_annotation", "regulon_annotation", "MotifEnrichment_preview", "centered_regulon_activity_groupby_design", "regulon_RSS_annotation", "csi_module_annotation"],
    "downloadName":["1.1.regulon_activity_heatmap_groupby_cells.pdf", "1.3.regulon_activity_heatmap.pdf", "2.2.RSS_ranking_plot.pdf", "2.3.RSS_heatmap.pdf", "3.2.regulons_csi_correlation_heatmap.pdf", "3.3.csi_module_activity_heatmap.pdf", "0.1.TF_target_enrichment_annotation.tsv", "0.2.regulon_annotation.tsv", "0.3.MotifEnrichment_preview.html", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "3.1.csi_module_annotation.tsv"],
    "downloadPath":["download/1.1.regulon_activity_heatmap_groupby_cells.pdf", "download/1.3.regulon_activity_heatmap.pdf", "download/2.2.RSS_ranking_plot.pdf", "download/2.3.RSS_heatmap.pdf", "download/3.2.regulons_csi_correlation_heatmap.pdf", "download/3.3.csi_module_activity_heatmap.pdf", "download/0.1.TF_target_enrichment_annotation.tsv", "download/0.2.regulon_annotation.tsv", "download/0.3.MotifEnrichment_preview.html", "download/1.2.centered_regulon_activity_groupby_design.tsv", "download/2.1.regulon_RSS_annotation.tsv", "download/3.1.csi_module_annotation.tsv"]}
    df = pd.DataFrame(d)
    df.to_csv(f"{wkdir}/SCENIC/output.json.tsv", index=False, sep="\t", header=True, encoding="utf-8")
#============link necessary results to output/download directory====================================================================
    res_lst = ["1.1.regulon_activity_heatmap_groupby_cells.pdf", "1.1.regulon_activity_heatmap_groupby_cells_data.tsv", "1.3.regulon_activity_heatmap.pdf", "2.2.RSS_ranking_plot.pdf", "2.3.RSS_heatmap.pdf", "3.2.regulons_csi_correlation_heatmap.pdf", "3.2.regulon_csi_correlation_heatmap_data.tsv", "3.3.csi_module_activity_heatmap.pdf", "3.3.csi_module_activity_heatmap_data.tsv", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "3.1.csi_module_annotation.tsv"]
    for f in res_lst:
        cmd_cp = f"mv {wkdir}/SCENIC/byclusters/{f} {wkdir}/SCENIC/"
        with module_cmd(environment) as p:
            status = p(cmd_cp, projectid, taskid)

    upload_lst = ["1.1.regulon_activity_heatmap_groupby_cells.pdf", "1.3.regulon_activity_heatmap.pdf", "2.2.RSS_ranking_plot.pdf", "2.3.RSS_heatmap.pdf", "3.2.regulons_csi_correlation_heatmap.pdf", "3.3.csi_module_activity_heatmap.pdf", "1.2.centered_regulon_activity_groupby_design.tsv", "2.1.regulon_RSS_annotation.tsv", "3.1.csi_module_annotation.tsv", "0.1.TF_target_enrichment_annotation.tsv", "0.2.regulon_annotation.tsv", "0.3.MotifEnrichment_preview.html"]
    for f in upload_lst:
        cmd_ln = f"ln -s {wkdir}/SCENIC/{f} {wkdir}/output/download/"
        with module_cmd(environment) as p:
            status = p(cmd_ln, projectid, taskid)
    return status
