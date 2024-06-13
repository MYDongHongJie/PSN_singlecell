import os,sys,shutil
import pandas as pd
# import json,subprocess,requests
# from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd

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
    # env_module = input.loc['project', 'environment']
    env_module = input.loc['parameters', 'environment']
    # env_module = "/home/xfyang/modulefiles/OESingleCell/dev"
    ## =================================================
    ## 1. aggr 
    ## =================================================
    logger.info("step3.1 解析样本，运行spaceranger aggr。")
    spaceranger_version = input.loc['parameters', 'spaceranger_version']

    if input.loc["parameters","aggr"] == True and len(input.loc["base","tasks"]) > 1:
    # # create library.csv
        libraries = pd.DataFrame(columns=['library_id', 'molecule_h5', 'cloupe_file', 'spatial_folder'])
        for task in input.loc["base","tasks"]:
            spaceranger_dir = f"{workdir}/../../{task['projectId']}/{task['taskId']}/spaceranger"
            sample = next(os.walk(f"{workdir}/../../{task['projectId']}/{task['taskId']}/spaceranger"))[1][0]
             # estimate tissue_positions_list.csv
    #         if os.path.exists(f"{spaceranger_dir}/{sample}/outs/spatial/tissue_positions.csv") and os.path.exists(f"{spaceranger_dir}/{sample}/outs/spatial/tissue_positions_list.csv"):
    #             os.remove(f"{spaceranger_dir}/{sample}/outs/spatial/tissue_positions_list.csv")
            libraries.loc[sample, 'library_id'] = sample
            libraries.loc[sample, 'molecule_h5'] = f"{spaceranger_dir}/{sample}/outs/molecule_info.h5"
            libraries.loc[sample, 'cloupe_file'] = f"{spaceranger_dir}/{sample}/outs/cloupe.cloupe"
            libraries.loc[sample, 'spatial_folder'] = f"{spaceranger_dir}/{sample}/outs/spatial"
        libraries.to_csv(f"{workdir}/{module}/libraries.csv", encoding='utf-8', index=False)
        # run aggr
        with module_cmd(env_module) as p:
            status=p(f"module load /public/scRNA_works/works/guochy/ST_taskserver/latest_ST_snakemake_pipeline/st_pipeline/envs/spaceranger/{spaceranger_version} && cd  {workdir}/{module}/ &&  rm -rf aggr &&  spaceranger aggr --id=aggr --csv={workdir}/{module}/libraries.csv --normalize=none ", projectid, taskid)
    else:
        logger.info("跳过 aggr（aggr选择为False或只存在一个样本）")

    ## =================================================
    ## 2. create 
    ## =================================================
    logger.info("step3.2 解析样本，生成metadata.csv文件，运行sctool create。")
    # conscruct metadata
    metadata = pd.DataFrame(input.loc["project","samples"])
    metadata=metadata.set_index('name',drop=True)
    if any(metadata["batchid"] == ""): metadata["batchid"] = 1
    if any(metadata["group"] == ""): metadata["group"] = metadata["label"]
    metadata [["group"+str(i) if i!= 1 else "group" 
        for i in range(1, 1+ metadata["group"].str.split(",",expand=True).shape[1])] ] = metadata["group"].str.split(",",expand=True)
    metadata["species"]= input.loc["project","species"]


    ## create soft link & write metadata.csv
    samples=[]
    for task in input.loc["base","tasks"]:
        spaceranger_dir = f"{workdir}/../../{task['projectId']}/{task['taskId']}/spaceranger"
        sample = next(os.walk(f"{workdir}/../../{task['projectId']}/{task['taskId']}/spaceranger"))[1][0]
        if sample in metadata.index:
            # soft link
            if not os.path.exists(f'{workdir}/{module}/spaceranger'):
                os.makedirs(f'{workdir}/{module}/spaceranger')
            os.symlink(f"{spaceranger_dir}/{sample}",
                       f"{workdir}/{module}/spaceranger/{metadata['label'][sample]}")
            # if not os.path.exists(f"{spaceranger_dir}/{sample}/outs/spatial/tissue_positions_list.csv"):
            #      tissue_positions = pd.read_csv(f"{spaceranger_dir}/{sample}/outs/spatial/tissue_positions.csv", dtype=str, sep=",", na_filter=False)
            #      tissue_positions.to_csv(f"{spaceranger_dir}/{sample}/outs/spatial/tissue_positions_list.csv",sep=',',index=False,header=False,encoding='utf-8')
            samples.append(sample)
        else: 
            print("sample %s not found."%sample)
            sys.exit(1)
            #sys.exit(f"sample {sample} not found.")
    metadata = metadata.loc[samples,:]
    metadata.rename(columns={"label":"sampleid"},inplace=True)
    metadata.to_csv(f"{workdir}/{module}/metadata.csv",sep=',',index=False,header=True)
    metadata.to_csv(f"{workdir}/{module}/metadata.tsv",sep='\t',index=False,header=True)
    library = input.loc['parameters', 'library_size']
    print(f"文库类型为：{library}")

    if library != "fresh":
        gset_param = " "
    
    if input.loc['parameters', 'vars2vis'] != []:
        metrics = input.loc['parameters', 'vars2vis']
        print(f"质控可视化指标： {', '.join(metrics)}")
    else:
        metrics = ""
    if 'percent.mito' in metrics:
        if len(input.loc['parameters', 'mito_gene_file'])!=0:
            gset_param = f" --metrics percent.mito --gset {workdir}/input/{input.loc['parameters', 'mito_gene_file']}"
        else:
            gset_param = f" --metrics percent.mito --gset {input.loc['parameters', 'database']}/MT_genelist.gmt"
    else:
        gset_param = " --metrics NULL "
    
    ## run cmd
    cmd_create = f"sctool  " \
        f"-i {workdir}/{module}/spaceranger --prefix spatial " \
        f"-f rds " \
        f"-o {workdir}/{module}/create  -d rds --assay Spatial  create  " \
        f"-s h5 -m {workdir}/{module}/metadata.csv --gcolumn 2  " \
        f"--cell.meta FALSE " \
        f"{gset_param}  "
    with module_cmd(env_module) as p:
        status=p(cmd_create, projectid, taskid)

    # ## =================================================
    # ## 3. qc
    # ## =================================================
    logger.info("step3.3:解析运行命令，运行sctool qc。")
    ## parsing qc_params_df:
    # env_module = "/home/xfyang/modulefiles/OESingleCell/dev"
    qc_params_df=pd.DataFrame(columns=("c","l","L"))
    vars2vis=[]
    if input.loc['parameters']['nGene'] == True:
        vars2vis.append("nFeature_Spatial")
        # if input.loc['parameters']["nGene_max"]=="" and input.loc['parameters']["nGene_min"]=="":
        qc_params_df.loc[len(qc_params_df.index)] = "nFeature_Spatial","NULL","NULL"
        # else:
            # qc_params_df.loc[len(qc_params_df.index)] = ["nFeature_Spatial",
            # '0' if input.loc['parameters']['nGene_min'] == '' else str(input.loc['parameters']['nGene_min']),
            # 'Inf' if input.loc['parameters']['nGene_max'] == '' else str(input.loc['parameters']['nGene_max'])]
    if input.loc['parameters']['nUMI'] == True:
        vars2vis.append("nCount_Spatial")
        # if input.loc['parameters']["nUMI_max"]=="" and input.loc['parameters']["nUMI_min"]=="":
        qc_params_df.loc[len(qc_params_df.index)] = "nCount_Spatial","NULL","NULL"
        # else:
            # qc_params_df.loc[len(qc_params_df.index)] = ["nCount_Spatial",
            # '0' if input.loc['parameters']['nUMI_min'] == '' else str(input.loc['parameters']['nUMI_min']),
            # 'Inf' if input.loc['parameters']['nUMI_max'] == '' else str(input.loc['parameters']['nUMI_max'])]
    if input.loc['parameters']['percent_mito'] == True :
        vars2vis.append("percent.mito")
        qc_params_df.loc[len(qc_params_df.index)] = "percent.mito","NULL",str(input.loc['parameters']['percent_mito_threshold'])
    
    if input.loc['parameters']['library_size'] == 'fresh' :
        vars2regress = "percent.mito"
    else:
        vars2regress = "NULL"
    
    if input.loc['parameters']['library_size'] == 'cytassist' :
        crop = "--crop TRUE" 
    else:
        crop = " "    


    qc_params = ""
    for i in qc_params_df.columns:
        qc_params += f"-{i}  {','.join(qc_params_df[i])}   "
    rlm_param = "-u nCount_Spatial,nFeature_Spatial " if input.loc['parameters', 'use_rlm'] else ""
    genes2filter_param = f"--features2filter {workdir}/input/{input.loc['parameters', 'genes2filter_file']}" if input.loc['parameters', 'genes2filter_file'] else ""
    sct_split= "--sct_split sampleid" if {input.loc['parameters', 'sct_split']} else " "
    cloud_meta_param = f"--cloudmeta {workdir}/input/cell.tsv" if os.path.exists(f"{workdir}/input/cell.tsv") else ""
    ## pasting cmd
    # cmd_qc =  f"Rscript {rscript_dir}/exec/sctool  "  \
    cmd_qc =  f"/public/scRNA_works/works/guochy/ST_taskserver/exec/sctool  "  \
    f"-i {workdir}/{module}/create/spatial.rds  "  \
    f"-o {workdir}/{module}  " \
    f"--informat rds --outformat rds  --assay Spatial  --prefix filtered_seurat " \
    f"qc {qc_params}  " \
    f"-r {vars2regress}  " \
    f"--rmdoublets F " \
    f"{rlm_param}  " \
    f"--mincell4gene {str(input.loc['parameters', 'min_cell_for_gene'])}  " \
    f"--normmeth {input.loc['parameters', 'normmeth']}  " \
    f"--pointsize 0.1  " \
    f"{genes2filter_param}  " \
    f"{sct_split}  " \
    f"--nvfeatures {str(input.loc['parameters', 'nvfeatures'])}  " \
    f"{crop} " \
    f"{cloud_meta_param}" 
    ## run cmd
    with module_cmd(env_module) as p:
        status=p(cmd_qc, projectid, taskid)

    ## =================================================
    ## 2. amending output.json
    ## =================================================
    logger.info("step3.3 生成output.json.tsv。")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False) #, encoding='gb2312')
    # df = output_df.loc[output_df['task_type'] == module]
    df= pd.DataFrame(columns=output_df.columns)
    df.loc[len(df.index)] = ["Count_QC","data","metadata.tsv","cell","metadata.tsv","metadata","",""]
    # df.loc[len(df.index)] = ["Count_QC","summary","cell_count_before_after_QC.tsv","tsv","cell_count_before_after_QC.tsv","细胞统计","",""]
    df.loc[len(df.index)] = ["Count_QC","summary","statitics_for_QC.tsv","tsv","statitics_for_QC.tsv","数据统计","",""]
    for i in vars2vis:
        df.loc[len(df.index)] = ["Count_QC","diagram",f"QC_featureplot_for_{i}.tsv","violin",
            f"QC_featureplot_for_{i}.tsv",i,f"QC_featureplot_for_{i}",f"download/QC_featureplot_for_{i}.png,download/QC_featureplot_for_{i}.pdf"]
    # reportVars
    cell_count_df = pd.read_csv(f"{workdir}/{module}/statitics_for_QC.xls", dtype=str, sep="\t", na_filter=False)
    cell_count_df.to_csv(f"{workdir}/{module}/statitics_for_QC.tsv",sep='\t',index=False,header=True,encoding='utf-8')
    # write output.json.tsv
    df.to_csv(f"{workdir}/{module}/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')

    ## =================================================
    ## 4. link image to download
    ## =================================================
    logger.info("step3.3 构建output/download/分析结果。")
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    # os.symlink(f"{workdir}/{module}/cell_count_before_after_QC.tsv",f"{workdir}/output/download/cell_count_before_after_QC.xls")
    os.symlink(f"{workdir}/{module}/statitics_for_QC.tsv",f"{workdir}/output/download/statitics_for_QC.xls")
    if input.loc["parameters","aggr"] == True and len(input.loc["base","tasks"]) > 1:
        cmd_ln = f"""ln -s {workdir}/{module}/{{*.pdf,*png}}  {workdir}/output/download/
                 ln -s {workdir}/{module}/aggr/outs  {workdir}/output/download/aggr"""
    else:
        cmd_ln = f"ln -s {workdir}/{module}/{{*.pdf,*png}}  {workdir}/output/download/"

    with module_cmd(env_module) as p:
        status = p(cmd_ln, projectid, taskid)
    
    return status
