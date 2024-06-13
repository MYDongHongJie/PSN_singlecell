import os, sys, re
import json
import subprocess
import pandas as pd
import requests
from glob import glob
from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

logger = getLogger('oe.cloud.sc.qsub')

def task_enrichment(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    ## 1.获取filterd.h5seurat路径及参数信息
    infile = f"{input.loc['parameters', 'enrich_gene_file']}"
    if infile == "":
        infile = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/diffexp/*-diff-pvalue-*_anno.tsv"
    go_bg = f"{input.loc['parameters', 'database']}/annotation/gene_go.backgroud.xls"
    kegg_bg = f"{input.loc['parameters', 'database']}/annotation/gene_kegg.backgroud.xls"

    ## 2.获取命令行
    cmd = f"module purge && perl  /public/scRNA_works/works/ziqingzhen/test/report_ht2023/scrna_pipeline/scripts/enrichment/enrich_go_kegg.pl " \
          f"-infile {infile} " \
          f"-go_bg {go_bg} " \
          f"-category /public/scRNA_works/works/ziqingzhen/test/report_ht2023/scrna_pipeline/scripts/enrichment/category.xls " \
          f"-kegg_bg {kegg_bg} " \
          f"-outdir {workdir}/{module} " \
          f"-shelldir {workdir}/{module}/enrichment_sh  " \
          f"-thread 4 " \
          f"-queue big "
    ## 3.执行分析
    with module_cmd(f"OESingleCell/v_3.0.0_visium_produce") as p:
        status = p(cmd, projectid, taskid)

    ## change xls to tsv
    os.rename(f"{workdir}/{module}/GO_enrichment/enrichment_go.xls",
              f"{workdir}/{module}/GO_enrichment/enrichment_go.tsv")
    os.rename(f"{workdir}/{module}/KEGG_enrichment/enrichment_kegg.xls",
              f"{workdir}/{module}/KEGG_enrichment/enrichment_kegg.tsv")
    enrich_files = glob(
        f'{workdir}/{module}/*_enrichment/*/*.xls')
    for file in enrich_files:
        os.rename(file, file.replace('.xls', '.tsv'))
        
    # =================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd = f"rm -r {workdir}/{module}/enrichment_sh && cp -r {workdir}/{module}/ {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*/*/*  {workdir}/output/download/*/*/*/*"  
    with module_cmd(f"OESingleCell/v_3.0.0_visium_produce") as p:
        status = p(cmd, projectid, taskid)

    #cmd = f"mv {workdir}/{module}/*/*/*.Total.tsv {workdir}/{module}/*/*/*.Down.tsv {workdir}/{module}/*/*/*.Up.tsv {workdir}/output/  "
    #with module_cmd(f"OESingleCell/v_3.0.0_visium_produce") as p:
    #    status = p(cmd, projectid, taskid)    
    
    # ==================================================================================================================
    # logger.info("step3.3:生成output.json.xls")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    # df = output_df.loc[output_df['task_type'] == module]
    df= pd.DataFrame(columns=output_df.columns)
    i = 0
    go_groups = os.listdir(f"{workdir}/{module}/GO_enrichment/")
    kegg_groups = os.listdir(f"{workdir}/{module}/KEGG_enrichment/")
    for go_group_file in go_groups:
        if not go_group_file.endswith('.tsv'):
            go_summary_filename = "enrichment-go-"+go_group_file+"-Total.tsv"
            go_summary_prefix = "enrichment-go-"+go_group_file 
            df.loc[len(df)] = ["enrichment", "summary", go_summary_filename,
                            "tsv", go_summary_prefix+"-Total.tsv", "go-"+go_group_file, "", ""]
            df.loc[len(df)] = ["enrichment", "diagram", f"{go_summary_prefix}-Total.circos.png",
                            "image", f"{go_summary_prefix}-Total.circos.png", "go-"+go_group_file, "go-"+go_group_file, f"download/enrichment/GO_enrichment/{go_group_file}/{go_summary_prefix}-Total.circos.pdf,download/enrichment/GO_enrichment/{go_group_file}/{go_summary_prefix}-Total.circos.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.top.Up.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.top.Up.tsv',
                                "bar_group", f'{go_group_file}-GO.top.Up.tsv', f'{go_group_file}-GO.top.Up', f"{go_group_file}-GO.top.Up", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Up.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Up.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.chord.top.Up.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.chord.top.Up.png',
                                "image", f'{go_group_file}-GO.chord.top.Up.png', f'{go_group_file}-GO.chord.top.Up', f"{go_group_file}-GO.chord.top.Up", f"download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Up.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Up.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.top.Total.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.top.Total.tsv',
                                "bar_group", f'{go_group_file}-GO.top.Total.tsv', f'{go_group_file}-GO.top.Total', f"{go_group_file}-GO.top.Total", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.chord.top.Total.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.chord.top.Total.png',
                                "image", f'{go_group_file}-GO.chord.top.Total.png', f'{go_group_file}-GO.chord.top.Total', f"{go_group_file}-GO.chord.top.Total", f"download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Total.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Total.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.top.Down.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.top.Down.tsv',
                                "bar_group", f'{go_group_file}-GO.top.Down.tsv', f'{go_group_file}-GO.top.Down', f"{go_group_file}-GO.top.Down", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Down.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Down.png"]
            
            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.chord.top.Down.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.chord.top.Down.png',
                                "image", f'{go_group_file}-GO.chord.top.Down.png', f'{go_group_file}-GO.chord.top.Down', f"{go_group_file}-GO.chord.top.Down", f"download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Down.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Down.png"]
            while i < len(df):
                if not os.path.exists(f"{workdir}/output/{df.loc[i, 'file']}"):
                    cmd1 = f"cp -rf {workdir}/{module}/GO_enrichment/{go_group_file}/{df.loc[i, 'input']} " \
                    f" {workdir}/output/{df.loc[i, 'file']}"
                    with module_cmd("obsutil/5.2.12") as p:
                        status=p(cmd1, projectid, taskid)
                i += 1


    for kegg_group_file in kegg_groups:
        if not kegg_group_file.endswith('.tsv'):
            kegg_summary_filename = "enrichment-kegg-"+kegg_group_file+"-Total.tsv"
            kegg_summary_prefix = "enrichment-kegg-"+kegg_group_file 
            df.loc[len(df)] = ["enrichment", "summary", kegg_summary_filename,
                            "tsv", "enrichment-kegg-"+kegg_group_file+"-Total.tsv", "kegg-"+kegg_group_file, "", ""]
            df.loc[len(df)] = ["enrichment", "diagram", f"{kegg_summary_prefix}-Total.circos.png",
                            "image", f"{kegg_summary_prefix}-Total.circos.png", "kegg-"+kegg_group_file, "kegg-"+kegg_group_file, f"download/enrichment/KEGG_enrichment/{kegg_group_file}/{kegg_summary_prefix}-Total.circos.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/{kegg_summary_prefix}-Total.circos.png"]

            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.top.Up.tsv',
                                "bar_group", f'{kegg_group_file}-KEGG.top.Up.tsv', f'{kegg_group_file}-KEGG.top.Up', f"{kegg_group_file}-KEGG.top.Up", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.png"]

            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Up.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.chord.top.Up.png',
                                "image", f'{kegg_group_file}-KEGG.chord.top.Up.png', f'{kegg_group_file}-KEGG.chord.top.Up', f"{kegg_group_file}-KEGG.chord.top.Up", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Up.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Up.png"]


            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.top.Total.tsv',
                                "bar_group", f'{kegg_group_file}-KEGG.top.Total.tsv', f'{kegg_group_file}-KEGG.top.Total', f"{kegg_group_file}-KEGG.top.Total", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.png"]

            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Total.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.chord.top.Total.png',
                                "image", f'{kegg_group_file}-KEGG.chord.top.Total.png', f'{kegg_group_file}-KEGG.chord.top.Total', f"{kegg_group_file}-KEGG.chord.top.Total", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Total.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Total.png"]


            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.top.Down.tsv',
                                "bar_group", f'{kegg_group_file}-KEGG.top.Down.tsv', f'{kegg_group_file}-KEGG.top.Down', f"{kegg_group_file}-KEGG.top.Down", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.png"]
            
            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Down.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.chord.top.Down.png',
                                "image", f'{kegg_group_file}-KEGG.chord.top.Down.png', f'{kegg_group_file}-KEGG.chord.top.Down', f"{kegg_group_file}-KEGG.chord.top.Down", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Down.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Down.png"]
            while i < len(df):
                if not os.path.exists(f"{workdir}/output/{df.loc[i, 'file']}"):
                    cmd1 = f"cp -rf {workdir}/{module}/KEGG_enrichment/{kegg_group_file}/{df.loc[i, 'input']} " \
                    f" {workdir}/output/{df.loc[i, 'file']}"
                    with module_cmd("obsutil/5.2.12") as p:
                        status=p(cmd1, projectid, taskid)
                i += 1


    df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',
              index=False, header=True, encoding='utf-8')


    
    return status

