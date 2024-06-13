import os, sys, re
import json
import subprocess
import pandas as pd
import requests
from glob import glob
import shutil
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
    if input.loc["parameters", "enrich_gene_file"] != "":
        infile = f"{workdir}/input/{input.loc['parameters', 'enrich_gene_file']}"
    else :
        infile = f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/diffexp/*-diff-*_anno.tsv"
    go_bg = f"{input.loc['parameters', 'database']}/annotation/gene_go.backgroud.xls"
    kegg_bg = f"{input.loc['parameters', 'database']}/annotation/gene_kegg.backgroud.xls"

    ## 2.获取命令行
    cmd = 'cat /proc/1/cgroup | grep -qi docker     && echo "Docker" || echo "Not Docker"'
    output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    if os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
        enrichment_script_dir = "/public/scRNA_works/works/ziqingzhen/test/sc_cloud/enrichment/"
        cmd = f"perl  {enrichment_script_dir}/enrich_go_kegg.pl " \
              f"-infile {infile} " \
              f"-go_bg {go_bg} " \
              f"-category {enrichment_script_dir}/category.xls " \
              f"-kegg_bg {kegg_bg} " \
              f"-outdir {workdir}/{module} " \
              f"-shelldir {workdir}/{module}/enrichment_sh  " \
              f"-thread 4 " \
              f"-queue big "
    else:
        cmd=f"enrichwrap -i {infile}  -g {input.loc['parameters', 'database']}/annotation/  -o {workdir}/tmp   --nomap  --minListHits 3   && " \
            f"mv  {workdir}/tmp/1.GO_enrichment    {workdir}/{module}/GO_enrichment && " \
            f"mv  {workdir}/tmp/2.KEGG_enrichment  {workdir}/{module}/KEGG_enrichment &&" \
            f"rm -rf {workdir}/tmp  "

    ## 3.执行分析
    with module_cmd(f"OESingleCell/3.0.d") as p:
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
        
    ##background_files
    backgroundfiles = {
                       f"{input.loc['parameters', 'database']}/annotation/gene_annotation.xls": f"{workdir}/{module}/background_files/gene_annotation.xls",
                       f"{input.loc['parameters', 'database']}/annotation/gene_kegg.backgroud.xls": f"{workdir}/{module}/background_files/gene_kegg.background.xls",
                       f"{input.loc['parameters', 'database']}/annotation/gene_go.backgroud.xls": f"{workdir}/{module}/background_files/gene_go.background.xls"
                      }
    os.makedirs(f"{workdir}/{module}/background_files")
    for source_file, target_file in backgroundfiles.items():
        shutil.copyfile(source_file, target_file)
    # =================================================================================================================
    logger.info("step3.2:构建output/download/分析结果")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    cmd = f" if [ -d {workdir}/{module}/enrichment_sh ] ;then rm -r {workdir}/{module}/enrichment_sh ;fi &&" \
          f" cp -r {workdir}/{module}/ {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*/*/*  {workdir}/output/download/*/*/*/* "
        # f" mv {workdir}/{module}/*/*/*_Total.tsv {workdir}/{module}/*/*/*_Down.tsv {workdir}/{module}/*/*/*_Up.tsv {workdir}/output/ &&" \
    with module_cmd(f"OESingleCell/3.0.d") as p:
        status = p(cmd, projectid, taskid)
    # ==================================================================================================================
    # logger.info("step3.3:生成output.json.xls")
    # output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)
    # # df = output_df.loc[output_df['task_type'] == module]
    # df= pd.DataFrame(columns=output_df.columns)
    # go_groups = os.listdir(f"{workdir}/{module}/GO_enrichment/")
    # kegg_groups = os.listdir(f"{workdir}/{module}/KEGG_enrichment/")
    # for go_group_file in go_groups:
        # if not go_group_file.endswith('.tsv'):
            # go_summary_filename = "GO_enrichment/"+go_group_file+"/enrichment-go-"+go_group_file+"-Total.tsv"
            # df.loc[len(df)] = ["enrichment", "summary", go_summary_filename,
                            # "tsv", "enrichment-go-"+go_group_file+"-Total.tsv", "go-"+go_group_file, "", ""]
            # # go_group = go_group_file
            # # if os.path.exists(f"{workdir}/output/GO_{go_group}_Up.tsv"):
            # #     df.loc[len(df)] = ["enrichment", "diagram", 'GO_'+go_group+'_Up.tsv',
            # #                     "bar_group", 'GO_'+go_group+'_Up.tsv', 'GO_'+go_group+'_Up', "GO.top.Up", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Up.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Up.png"]
            # # if os.path.exists(f"{workdir}/output/GO_{go_group}_Total.tsv"):
            # #     df.loc[len(df)] = ["enrichment", "diagram", 'GO_'+go_group+'_Total.tsv',
            # #                     "bar_group", 'GO_'+go_group+'_Total.tsv', 'GO_'+go_group+'_Total', "GO.top.Total", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.png"]
            # # if os.path.exists(f"{workdir}/output/GO_{go_group}_Down.tsv"):
            # #     df.loc[len(df)] = ["enrichment", "diagram", 'GO_'+go_group+'_Down.tsv',
            # #                     "bar_group", 'GO_'+go_group+'_Down.tsv', 'GO_'+go_group+'_Down', "GO.top.Down", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Down.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Down.png"]
    # for kegg_group_file in kegg_groups:
        # if not kegg_group_file.endswith('.tsv'):
            # kegg_summary_filename = "KEGG_enrichment/"+kegg_group_file + \
                # "/enrichment-kegg-"+kegg_group_file+"-Total.tsv"
            # df.loc[len(df)] = ["enrichment", "summary", kegg_summary_filename,
                            # "tsv", "enrichment-kegg-"+kegg_group_file+"-Total.tsv", "kegg-"+kegg_group_file, "", ""]
    # #         kegg_group = kegg_group_file
    # #         if os.path.exists(f"{workdir}/output/KEGG_{kegg_group}_Up.tsv"):
    # #             df.loc[len(df)] = ["enrichment", "diagram", 'KEGG_'+kegg_group+'_Up.tsv',
    # #                             "bubble", 'KEGG_'+kegg_group+'_Up.tsv', 'KEGG_'+kegg_group+'_Up', "KEGG.top.Up", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.png"]
    # #         if os.path.exists(f"{workdir}/output/KEGG_{kegg_group}_Total.tsv"):
    # #             df.loc[len(df)] = ["enrichment", "diagram", 'KEGG_'+kegg_group+'_Total.tsv',
    # #                             "bubble", 'KEGG_'+kegg_group+'_Total.tsv', 'KEGG_'+kegg_group+'_Total', "KEGG.top.Total", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.png"]
    # #         if os.path.exists(f"{workdir}/output/KEGG_{kegg_group}_Down.tsv"):
    # #             df.loc[len(df)] = ["enrichment", "diagram", 'KEGG_'+kegg_group+'_Down.tsv',
    # #                             "bubble", 'KEGG_'+kegg_group+'_Down.tsv', 'KEGG_'+kegg_group+'_Down', "KEGG.top.Down", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.png"]
    # df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',
              # index=False, header=True, encoding='utf-8')
    df = pd.DataFrame(columns=['task_type', 'result_module', 'input', 'type', 'file', 'title','downloadName', 'downloadPath'])
    i = 0
    go_groups = os.listdir(f"{workdir}/{module}/GO_enrichment/")
    kegg_groups = os.listdir(f"{workdir}/{module}/KEGG_enrichment/")
    for go_group_file in go_groups:
        if not go_group_file.endswith('.tsv'):
            go_summary_filename = "enrichment-go-"+go_group_file+"-Total.tsv"
            go_summary_prefix = "enrichment-go-"+go_group_file 
            df.loc[len(df)] = ["enrichment", "data", go_summary_filename,
                            "tsv", go_summary_prefix+"-Total.tsv", "go-"+go_group_file, "", ""]
            df.loc[len(df)] = ["enrichment", "diagram", f"{go_summary_prefix}-Total.circos.png",
                            "image", f"{go_summary_prefix}-Total.circos.png", "go-"+go_group_file, "go-"+go_group_file, f"download/enrichment/GO_enrichment/{go_group_file}/{go_summary_prefix}-Total.circos.pdf,download/enrichment/GO_enrichment/{go_group_file}/{go_summary_prefix}-Total.circos.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.top.Up.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.top.Up.png',
                                "image", f'{go_group_file}-GO.top.Up.png', f'{go_group_file}-GO.top.Up', f"{go_group_file}-GO.top.Up", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Up.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Up.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.chord.top.Up.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.chord.top.Up.png',
                                "image", f'{go_group_file}-GO.chord.top.Up.png', f'{go_group_file}-GO.chord.top.Up', f"{go_group_file}-GO.chord.top.Up", f"download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Up.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Up.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.top.Total.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.top.Total.png',
                                "image", f'{go_group_file}-GO.top.Total.png', f'{go_group_file}-GO.top.Total', f"{go_group_file}-GO.top.Total", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.chord.top.Total.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.chord.top.Total.png',
                                "image", f'{go_group_file}-GO.chord.top.Total.png', f'{go_group_file}-GO.chord.top.Total', f"{go_group_file}-GO.chord.top.Total", f"download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Total.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.chord.top.Total.png"]

            if os.path.exists(f"{workdir}/{module}/GO_enrichment/{go_group_file}/GO.top.Down.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'GO.top.Down.png',
                                "image", f'{go_group_file}-GO.top.Down.png', f'{go_group_file}-GO.top.Down', f"{go_group_file}-GO.top.Down", f"download/enrichment/GO_enrichment/{go_group_file}/GO.top.Down.pdf,download/enrichment/GO_enrichment/{go_group_file}/GO.top.Down.png"]
            
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
            df.loc[len(df)] = ["enrichment", "data", kegg_summary_filename,
                            "tsv", "enrichment-kegg-"+kegg_group_file+"-Total.tsv", "kegg-"+kegg_group_file, "", ""]
            df.loc[len(df)] = ["enrichment", "diagram", f"{kegg_summary_prefix}-Total.circos.png",
                            "image", f"{kegg_summary_prefix}-Total.circos.png", "kegg-"+kegg_group_file, "kegg-"+kegg_group_file, f"download/enrichment/KEGG_enrichment/{kegg_group_file}/{kegg_summary_prefix}-Total.circos.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/{kegg_summary_prefix}-Total.circos.png"]

            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.top.Up.png',
                                "image", f'{kegg_group_file}-KEGG.top.Up.png', f'{kegg_group_file}-KEGG.top.Up', f"{kegg_group_file}-KEGG.top.Up", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Up.png"]

            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Up.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.chord.top.Up.png',
                                "image", f'{kegg_group_file}-KEGG.chord.top.Up.png', f'{kegg_group_file}-KEGG.chord.top.Up', f"{kegg_group_file}-KEGG.chord.top.Up", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Up.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Up.png"]


            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.top.Total.png',
                                "image", f'{kegg_group_file}-KEGG.top.Total.png', f'{kegg_group_file}-KEGG.top.Total', f"{kegg_group_file}-KEGG.top.Total", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.png"]

            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Total.png"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.chord.top.Total.png',
                                "image", f'{kegg_group_file}-KEGG.chord.top.Total.png', f'{kegg_group_file}-KEGG.chord.top.Total', f"{kegg_group_file}-KEGG.chord.top.Total", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Total.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.chord.top.Total.png"]


            if os.path.exists(f"{workdir}/{module}/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.tsv"):
                df.loc[len(df)] = ["enrichment", "diagram", 'KEGG.top.Down.png',
                                "image", f'{kegg_group_file}-KEGG.top.Down.png', f'{kegg_group_file}-KEGG.top.Down', f"{kegg_group_file}-KEGG.top.Down", f"download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.pdf,download/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Down.png"]
            
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

    df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',index=False, header=True, encoding='utf-8')
    return status
