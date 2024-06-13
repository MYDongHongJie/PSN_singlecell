import os,sys,shutil
import pandas as pd
# import json,subprocess,requests
# from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd
from glob import glob
logger = getLogger('oe.cloud.sc.qsub')
# for testing
# rscript_dir = "/home/dongjiaoyang/git/oesinglecell3-cloud"

def task_QC(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    proj_workdir = f"{workdir}/../"
    module = input.loc['task', 'type']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}')
    # alternative for testing:
    # if os.path.exists(f'{workdir}/{module}'):
    #     shutil.rmtree(f'{workdir}/{module}')
    #     shutil.rmtree(f'{workdir}/output',ignore_errors=True)
    # os.makedirs(f'{workdir}/{module}')
    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    ## =================================================
    ## 1. qc
    ## =================================================
    logger.info("step3.1解析运行命令，运行sctool qc。")
    ## parsing qc_params_df:
    qc_params_df=pd.DataFrame(columns=("c","l","L"))
    vars2vis=[]
    if input.loc['parameters']['QC_default'] == "default1":
        logger.info("选用第一套默认QC模板进行过滤……")
        vars2vis=["nFeature_RNA","nCount_RNA","percent.mito","percent.HB","log10GenesPerUMI"]
        qc_params_df.loc[len(qc_params_df.index)] = "nFeature_RNA","200",'Inf'
        qc_params_df.loc[len(qc_params_df.index)] = "nCount_RNA","1000",'Inf'
        qc_params_df.loc[len(qc_params_df.index)] = "percent.mito","\' -Inf\'","None"
        qc_params_df.loc[len(qc_params_df.index)] = "percent.HB", "\' -Inf\'", "0.05"
        qc_params_df.loc[len(qc_params_df.index)] = "log10GenesPerUMI", "0.7", "Inf"
        nUMI_nGene_mad_parmas = f"  "
    elif input.loc['parameters']['QC_default'] == "default2":
        logger.info("选用第二套默认QC模板进行过滤……")
        vars2vis=["nFeature_RNA","nCount_RNA","percent.mito"]
        qc_params_df.loc[len(qc_params_df.index)] = "nFeature_RNA","NULL",'NULL'
        qc_params_df.loc[len(qc_params_df.index)] = "nCount_RNA","NULL",'NULL'
        qc_params_df.loc[len(qc_params_df.index)] = "percent.mito","0","None"
        nUMI_nGene_mad_parmas = f"  --nfold 2 "
    else:
        logger.info("根据自选参数进行过滤……")
        if input.loc['parameters']['nGene'] == True:
            vars2vis.append("nFeature_RNA")
            if input.loc['parameters']["nGene_max"]=="" and input.loc['parameters']["nGene_min"]=="":
                qc_params_df.loc[len(qc_params_df.index)] = "nFeature_RNA","NULL","NULL"
            else:
                qc_params_df.loc[len(qc_params_df.index)] = ["nFeature_RNA",
                '0' if input.loc['parameters']['nGene_min'] == '' else str(input.loc['parameters']['nGene_min']),
                'Inf' if input.loc['parameters']['nGene_max'] == '' else str(input.loc['parameters']['nGene_max'])]
        if input.loc['parameters']['nUMI'] == True:
            vars2vis.append("nCount_RNA")
            if input.loc['parameters']["nUMI_max"]=="" and input.loc['parameters']["nUMI_min"]=="":
                qc_params_df.loc[len(qc_params_df.index)] = "nCount_RNA","NULL","NULL"
            else:
                qc_params_df.loc[len(qc_params_df.index)] = ["nCount_RNA",
                '0' if input.loc['parameters']['nUMI_min'] == '' else str(input.loc['parameters']['nUMI_min']),
                'Inf' if input.loc['parameters']['nUMI_max'] == '' else str(input.loc['parameters']['nUMI_max'])]
        if input.loc['parameters']['percent_mito'] == True :
            vars2vis.append("percent.mito")
            qc_params_df.loc[len(qc_params_df.index)] = "percent.mito","\' -Inf\'",str(input.loc['parameters']['percent_mito_threshold'])
        if input.loc['parameters']['percent_HB'] == True :
            vars2vis.append("percent.HB")
            qc_params_df.loc[len(
                qc_params_df.index)] = "percent.HB", "\' -Inf\'", str(input.loc['parameters']['percent_HB_threshold'])
        if input.loc['parameters']['log10GenesPerUMI'] == True :
            vars2vis.append("log10GenesPerUMI")
            qc_params_df.loc[len(
                qc_params_df.index)] = "log10GenesPerUMI", str(input.loc['parameters']['log10GenesPerUMI_threshold']), "Inf"
        nUMI_nGene_mad_parmas = f"  --nfold {str(input.loc['parameters', 'nUMI_nGene_mad'])}  " if input.loc['parameters',
                'nUMI_nGene_mad'] !='' else f"  "

    qc_params = ""
    for i in qc_params_df.columns:
        qc_params += f"-{i}  {','.join(qc_params_df[i])}   "
    rlm_param = "-u nCount_RNA:nFeature_RNA " if input.loc['parameters', 'use_rlm'] else ""
    genes2filter_param = f"--genes2filter {workdir}/input/{input.loc['parameters', 'genes2filter_file']}" if input.loc['parameters', 'genes2filter_file'] else ""
    cloud_meta_param = f"--cloudmeta {workdir}/input/cell.tsv" if os.path.exists(f"{workdir}/input/cell.tsv") else ""
    ## pasting cmd
    # cmd_qc =  f"Rscript {rscript_dir}/exec/sctool  "  \
    ####进行输入对象rds或者h5seurat格式的识别和判断
    base_taskid = input.loc["base", "tasks"][0]["taskId"]
    base_wd = f'/public/cloud_scRNA/20{base_taskid[0:2]}/{base_taskid[2:4]}/{base_taskid[0:10]}/{base_taskid}'
    base_json = pd.read_json(f'{base_wd}/input/input.json', orient="index", dtype={"id": 'str'})
    base_module_name = base_json.loc['task', 'type']
    def search_files(directory, extension):
        pattern = f"{directory}/*.{extension}"
        files = glob(pattern)
        return files

    rds_input = search_files(f'{base_wd}/{base_module_name}', "rds")
    h5seurat_input = search_files(f'{base_wd}/{base_module_name}', "h5seurat")
    seurat_ob = ""
    format_ob = ""
    if (len(rds_input) == 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    elif (len(rds_input) != 0 and len(h5seurat_input) == 0):
        seurat_ob = rds_input[0]
        format_ob = "rds"
    elif (len(rds_input) != 0 and len(h5seurat_input) != 0):
        seurat_ob = h5seurat_input[0]
        format_ob = "h5seurat"
    else:
        print(f"输入seurat对象为空，请检查！！！")
    #seurat_ob_realpath=os.readlink(seurat_ob)
    min_cell_for_gene = '0' if input.loc['parameters']['min_cell_for_gene'] == '' else str(input.loc['parameters']['min_cell_for_gene'])
    ####以上完成输入对象的识别及赋值
    cmd_qc =  f"sctool  "  \
    f"-i {seurat_ob} "  \
    f"-o {workdir}/{module}  " \
    f"-f {format_ob}  -d {format_ob}  --assay RNA  --dataslot counts  -j 6  --prefix filtered qc " \
    f"--QC TRUE  " \
    f"{qc_params}  " \
    f"-r {','.join(input.loc['parameters', 'regress'])}  " \
    f"--cut.1 median  " \
    f"--cut.2 mad  " \
    f"{nUMI_nGene_mad_parmas}  " \
    f"--rmdoublets {bool(input.loc['parameters', 'rmdoublets'])}  " \
    f"--method {input.loc['parameters', 'rmdoublets_method']}  " \
    f"{rlm_param}  " \
    f"--mincell4gene {min_cell_for_gene}  " \
    f"-m {input.loc['parameters', 'normmeth']}  " \
    f"--pointsize 0.1  " \
    f"{genes2filter_param}  " \
    f"--qc_groupby {input.loc['parameters', 'qc_groupby']}  " \
    #f"--qc_groupby sampleid  " \
    f"{cloud_meta_param}"
    ## run cmd
    with module_cmd(f"{input.loc['project', 'environment']}") as p:
        status=p(cmd_qc, projectid, taskid)
    ## =================================================
    ## 2. amending output.json
    ## =================================================
    logger.info("step3.2 生成output.json.tsv。")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False) #, encoding='gb2312')
    # df = output_df.loc[output_df['task_type'] == module]
    df= pd.DataFrame(columns=output_df.columns)
    df.loc[len(df.index)] = ["QC","data","metadata.tsv","cell","metadata.tsv","metadata","",""]
    df.loc[len(df.index)] = ["QC","summary","cell_count_before_after_QC.tsv","tsv","cell_count_before_after_QC.tsv","细胞统计","",""]
    df.loc[len(df.index)] = ["QC","summary","cell_statitics_before_after_QC.tsv","tsv","cell_statitics_before_after_QC.tsv","数据统计","",""]
    for i in vars2vis:
        df.loc[len(df.index)] = ["QC","diagram",f"QC_{i}_afterQC.tsv","violin",
            f"QC_{i}_afterQC.tsv",i,f"QC_{i}_afterQC",f"download/QC_{i}_afterQC.png,download/QC_{i}_afterQC.pdf"]
    # reportVars
    cell_count_df = pd.read_csv(f"{workdir}/{module}/cell_count_before_after_QC.tsv", dtype=str, sep="\t", na_filter=False)
    cell_count = {
    "sample_num": cell_count_df.shape[0],
    "min_beforeQC_cell": min(cell_count_df["Total_cells_beforeQC"]),
    "max_beforeQC_cell": max(cell_count_df["Total_cells_beforeQC"]),
    "min_afterQC_cell": min(cell_count_df["Total_cells_afterQC"]),
    "max_afterQC_cell": max(cell_count_df["Total_cells_afterQC"])
    }
    if 'nFeature_RNA' in vars2vis:
        cell_count["min_afterQC_gene"]=min(cell_count_df["nFeature_RNA_min"])
        cell_count["max_afterQC_gene"]=max(cell_count_df["nFeature_RNA_max"])

    if 'nCount_RNA' in vars2vis:
        cell_count["min_afterQC_umi"]=min(cell_count_df["nCount_RNA_min"])
        cell_count["max_afterQC_umi"]=max(cell_count_df["nCount_RNA_max"])

    if 'percent_mito' in vars2vis:
        cell_count["min_afterQC_mito"]=min(cell_count_df["percent.mito_min"])
        cell_count["max_afterQC_mito"]=max(cell_count_df["percent.mito_max"])
    if 'percent_HB' in vars2vis:
        cell_count["min_afterQC_HB"]=min(cell_count_df["percent.HB_min"])
        cell_count["max_afterQC_HB"]=max(cell_count_df["percent.HB_max"])
    if 'log10GenesPerUMI' in vars2vis:
        cell_count["min_afterQC_log10GenesPerUMI"]=min(cell_count_df["log10GenesPerUMI_min"])
        cell_count["max_afterQC_log10GenesPerUMI"]=max(cell_count_df["log10GenesPerUMI_max"])

    df.loc[len(df.index)] = ["QC","reportVars","metadata.tsv",str(cell_count),"","","",""]
    # write output.json.tsv
    df.to_csv(f"{workdir}/{module}/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')

    ## =================================================
    ## 4. link image to download
    ## =================================================
    logger.info("step3.3 构建output/download/分析结果。")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    os.symlink(f"{workdir}/{module}/cell_count_before_after_QC.tsv",f"{workdir}/output/download/cell_count_before_after_QC.xls")
    os.symlink(f"{workdir}/{module}/cell_statitics_before_after_QC.tsv",f"{workdir}/output/download/cell_statitics_before_after_QC.xls")
    cmd_ln = f"ln -s {workdir}/{module}/{{*.pdf,*png}}  {workdir}/output/download/"
    with module_cmd(f"{input.loc['project', 'environment']}") as p:
        status = p(cmd_ln, projectid, taskid)

    return status
