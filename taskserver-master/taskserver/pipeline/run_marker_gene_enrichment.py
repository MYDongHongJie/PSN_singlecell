import pandas as pd 
import numpy as np
import os
import sys
import requests
from glob import glob
import json
import subprocess
from shutil import copyfile
from oebio.utils.log import getLogger
logger = getLogger('oe.cloud.sc.qsub')


from taskserver.tools.module_cmd import module_cmd

def task_marker_gene_enrichment(input,projectid="项目ID", taskid="任务ID", workdir="分析目录"):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.xls文件，用于生成output.json
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    wkdir =workdir
    os.makedirs( f"{wkdir}/output/download")
    os.makedirs( f"{wkdir}/marker_gene_enrichment")

    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')

    infile = f"{input.loc['parameters', 'gene_file']}"
    subset = f"{input.loc['parameters', 'subset']}"

    if infile == "":
        infile = f"{wkdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/findallmarkers/all_markers_for_each_cluster_anno.tsv"

    go_bg = f"{input.loc['parameters', 'database']}/annotation/gene_go.backgroud.xls"
    kegg_bg = f"{input.loc['parameters', 'database']}/annotation/gene_kegg.backgroud.xls"

    ### 2. 整理文件 且 获取命令行 
    print(infile)
    file = pd.read_table(infile)
    if ("cluster" in file.columns) and len(file.columns) >1: 
        cluster_number = file.cluster.drop_duplicates().tolist()
        if infile.split("/")[-1] == 'all_markers_for_each_cluster_anno.tsv':
            prefix = f"all_markers_for"
        else:
            prefix = ""
        
        if subset == "":
            for index in cluster_number: 
                df = file[file.cluster == index ]
                df.to_csv(f"{wkdir}/{module}/{prefix}_cluster{index}.xls", sep="\t", encoding='utf-8', index=False) 
        else:
            list = subset.split(",")
            if type(cluster_number[1]) == int:
                numbers = [ int(float(x) ) for x in list ]
                set(numbers).issubset(set(cluster_number))
                for index in numbers: 
                    df = file[file.cluster == index ]
                    df.to_csv(f"{wkdir}/{module}/{prefix}_cluster{index}.xls", sep="\t", encoding='utf-8', index=False) 
                print(f"{wkdir}/{module}/{prefix}_cluster{index}.xls")
            else:
                set(list).issubset(set(cluster_number))
                for index in list: 
                    df = file[file.cluster == index ]
                    df.to_csv(f"{wkdir}/{module}/{prefix}_cluster{index}.xls", sep="\t", encoding='utf-8', index=False) 
            
        gene_files=f"{wkdir}/{module}/*_cluster*.xls"       

    elif ("gene_module" in file.columns) and len(file.columns) >1:
        gene_module_number = file.gene_module.drop_duplicates()
        if infile.split("/")[-1] == 'pseudotime_heatmap_gene_module_anno.xls':
            prefix = f"pseudotime_heatmap_"
        else:
            prefix = ""
        for index in cluster_number: 
            df = file[file.cluster == index ]
            df.to_csv(f"{wkdir}/{module}/{prefix}gene_module{index}.xls", sep="\t", encoding='utf-8', index=False) 

        
        gene_files=f"{workdir}/{module}/*_cluster*.xls" 
        
    else:
        gene_files=infile

    print(gene_files)
    cmd = f"module purge && module load OESingleCell/3.0.d && perl /public/scRNA_works/works/liuhongyan/Test/oecloude/marker_enrichment/enrich_go_kegg.pl " \
          f"-infile {gene_files} " \
          f"-go_bg {go_bg} " \
          f"-category /home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls  "  \
          f"-kegg_bg {kegg_bg} " \
          f"-outdir {workdir}/{module}/enrichment " \
          f"-shelldir {workdir}/{module}/enrichment_sh " \
          f"-thread 1 "\
          f"-queue big "

    logger.info(cmd)
    logger.info("step1: 开始 marker gene 富集分析")
    ### env ####
    ret = subprocess.run(cmd,
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        encoding="utf-8")
    logger.info(ret.stdout)  
    logger.info(ret.stderr)  


    ## change xls to tsv
    os.rename(f"{workdir}/{module}/enrichment/GO_enrichment/enrichment_go.xls",
              f"{workdir}/{module}/enrichment/GO_enrichment/enrichment_go.tsv")
    os.rename(f"{workdir}/{module}/enrichment/KEGG_enrichment/enrichment_kegg.xls",
              f"{workdir}/{module}/enrichment/KEGG_enrichment/enrichment_kegg.tsv")
    enrich_files = glob(f'{workdir}/{module}/enrichment/*/*/*.xls')
    for file in enrich_files:
        os.rename(file, file.replace('.xls', '.tsv'))

# ==================================================================================================================
    logger.info("step2: 构建output/download/分析结果")
    ## change xls to tsv
    cmd = f"cp {workdir}/{module}/*/*/*/*-Total.tsv {workdir}/output/ && cp -r {workdir}/{module}/ {workdir}/output/download/ && rename tsv xls {workdir}/output/download/*/*/*  {workdir}/output/download/*/*/*/* && rm -r {workdir}/output/download/*/enrichment_sh"

    with module_cmd(f"OESingleCell/3.0.d") as p:
        status = p(cmd, projectid, taskid)

    logger.info("step3:生成output.json.xls")
    df = pd.DataFrame(columns=['task_type', 'result_module', 'input', 'type', 'file', 'title','downloadName', 'downloadPath'])
    go_groups = os.listdir(f"{workdir}/{module}/enrichment/GO_enrichment/")
    kegg_groups = os.listdir(f"{workdir}/{module}/enrichment/KEGG_enrichment/")
    for go_group_file in go_groups:
        if not go_group_file.endswith('.tsv'):
            go_summary_filename = "GO_enrichment/"+go_group_file+"/enrichment-go-"+go_group_file+"-Total.tsv"
            cmd = f"cp {workdir}/{module}/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.png {workdir}/output/{go_group_file}.GO.top.Total.png"
            os.system(cmd)
            df.loc[len(df)] = ["marker_gene_enrichment", "data", go_summary_filename,
                            "tsv", "enrichment-go-"+go_group_file+"-Total.tsv", "go-"+go_group_file, "", ""]      
            df.loc[len(df)] = ["marker_gene_enrichment", "diagram", go_group_file+'.GO.'+'top.Total.png',
                                "image", go_group_file+'.GO.'+'top.Total.png', go_group_file+'.GO.'+'top.Total', go_group_file+".GO.top.Total", f"download/{module}/enrichment/GO_enrichment/{go_group_file}/GO.top.Total.png"]   
    for kegg_group_file in kegg_groups:
        if not kegg_group_file.endswith('.tsv'):
            kegg_summary_filename = "KEGG_enrichment/"+kegg_group_file + \
                "/enrichment-kegg-"+kegg_group_file+"-Total.tsv"
            df.loc[len(df)] = ["enrichment", "summary", kegg_summary_filename,
                            "tsv", "enrichment-kegg-"+kegg_group_file+"-Total.tsv", "kegg-"+kegg_group_file, "", ""]
            cmd = f"cp {workdir}/{module}/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.png {workdir}/output/{kegg_group_file}.KEGG.top.Total.png"
            os.system(cmd)
            df.loc[len(df)] = ["marker_gene_enrichment", "diagram", kegg_group_file+'.KEGG.'+'top.Total.png',
                                "bar_group", kegg_group_file+'.KEGG.'+'top.Total.png', kegg_group_file+'.KEGG.'+'top.Total', kegg_group_file+".KEGG.top.Total", f"download/{module}/enrichment/KEGG_enrichment/{kegg_group_file}/KEGG.top.Total.png"] 
    df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t',index=False, header=True, encoding='utf-8')
