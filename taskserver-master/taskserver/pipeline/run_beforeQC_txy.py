import os,sys,shutil
import pandas as pd
# import json,subprocess,requests
# from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd
from selenium import webdriver
from glob import glob
logger = getLogger('oe.cloud.sc.qsub')
# for testing
# rscript_dir = "/home/dongjiaoyang/git/oesinglecell3-cloud"

def open_html(file):
    options = webdriver.ChromeOptions()
    options.add_argument("headless")
    options.add_argument('no-sandbox')
    options.add_argument('dns-prefetch-disable')
    options.add_argument('disable-gpu')
    options.add_argument('disable-dev-shm-usage')
    options.add_argument('disable-features=VizDisplayCompositor')
    options.add_argument('disable-features=NetworkService')
    options.add_argument('window-size=1920x945')
    options.add_experimental_option("excludeSwitches",["ignore-certificate-errors"])
    browser = webdriver.Chrome('/home/luyao/chromedriver', chrome_options=options)
    browser.get(f"file://{file}")
    # chemistry = browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(3) > td.tableMetric-value-chopped")
    try:
        reference_path =  browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(5) > td.tableMetric-value")
        return reference_path.text
    except:
        sys.exit("An error occur while finding reference path element.")
    # include-introns
    # if browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(7) > td.tableMetric-value-chopped").text == "cellranger-5.0.0":
    #     include_introns = browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(4) > td.tableMetric-value").text
    # else:
    #     include_introns = "unknown in v3"


def task_beforeQC_txy(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    # cellranger_version = input.loc['parameters', 'cellranger_version']
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}/cellranger')
    # alternative for testing:
    # if os.path.exists(f'{workdir}/{module}'):
    #     shutil.rmtree(f'{workdir}/{module}')
    #     shutil.rmtree(f'{workdir}/output',ignore_errors=True)
    # os.makedirs(f'{workdir}/{module}/cellranger')
    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    ## =================================================
    ## 1. create
    ## =================================================
    ## conscruct metadata
    logger.info("step3.1 解析样本，生成metadta.csv文件，运行sctool create。")
    metadata = pd.DataFrame(input.loc["project","samples"])
    metadata=metadata.set_index('name',drop=True)
    #if any(metadata["group"] == ""): metadata["group"] = metadata["label"]
    metadata["group"] = metadata.apply(lambda row: row["label"] if row["group"].strip() == "" else row["group"], axis=1)
    metadata = metadata.rename(columns={'group': 'group_raw'})
    #metadata [["group"+str(i) if i!= 1 else "group" for i in range(1, 1+ metadata["group"].str.split(",",expand=True).shape[1])] ] = metadata["group"].str.split(",",expand=True)
    metadata [["group"+str(i)  if i!= 0 else "group"  for i in range(0, metadata["group_raw"].str.split(",",expand=True).shape[1])] ] = metadata["group_raw"].str.split(",",expand=True)
    metadata["species"]= input.loc["project","species"]

    ## create soft link & write metadata.csv
    samples=[]
    for task in input.loc["base","tasks"]:
        sample = next(os.walk(f"{workdir}/../../../../20{task['projectId'][0:2]}/{task['projectId'][2:4]}/{task['projectId']}/{task['taskId']}/cellranger"))[1][0]
        if sample in metadata.index:
            # soft link
            os.symlink(f"{workdir}/../../../../20{task['projectId'][0:2]}/{task['projectId'][2:4]}/{task['projectId']}/{task['taskId']}/cellranger/{sample}",
                       f"{workdir}/{module}/cellranger/{metadata['label'][sample]}")
            samples.append(sample)
        else:
            print("sample %s not found."%sample)
            sys.exit(1)
            #sys.exit(f"sample {sample} not found.")
    metadata = metadata.loc[samples,:]
    metadata.rename(columns={"label":"sampleid"},inplace=True)
    metadata.to_csv(f"{workdir}/{module}/metadata.csv",sep=',',index=False,header=True)

    metrics = input.loc['parameters', 'vars2vis']
    print(f"质控可视化指标： {', '.join(metrics)}")

    if {'percent.mito','percent.HB'}.issubset(metrics):
        if len(input.loc['parameters', 'mito_gene_file'])!=0:
            gset_param = f"--gset {workdir}/input/{input.loc['parameters', 'mito_gene_file']},{input.loc['parameters', 'database']}/HB_genelist.gmt"
        else:
            gset_param = f"--gset {input.loc['parameters', 'database']}/MT_genelist.gmt,{input.loc['parameters', 'database']}/HB_genelist.gmt"
    elif 'percent.mito' in metrics:
        if len(input.loc['parameters', 'mito_gene_file'])!=0:
            gset_param = f"--gset {workdir}/input/{input.loc['parameters', 'mito_gene_file']}"
        else:
            gset_param = f"--gset {input.loc['parameters', 'database']}/MT_genelist.gmt"
    elif 'percent.HB' in metrics:
            gset_param = f"--gset {input.loc['parameters', 'database']}/HB_genelist.gmt"
    else:
        gset_param = ""

    ## run cmd
    # cmd_create =  f"Rscript {rscript_dir}/exec/sctool  " \
    cmd_create = f"mkdir $PWD/tmp \n" \
        f"export TMPDIR=$PWD/tmp \n" \
        f"sctool  " \
        f"-i {workdir}/{module}/cellranger  " \
        f"-o {workdir}/{module}  -d h5seurat --assay RNA create  " \
        f"-s mtx -m {workdir}/{module}/metadata.csv --gcolumn 2  " \
        f"--vars2vis {','.join(metrics)}  " \
        f"--qc_groupby {input.loc['parameters', 'qc_groupby']}  " \
        f"{gset_param}  "
    with module_cmd(f"{input.loc['project', 'environment']}") as p:
        status=p(cmd_create, projectid, taskid)

    ## =================================================
    ## 2. amending output.json
    ## =================================================
    logger.info("step3.3 生成output.json.tsv。")
    output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False) #, encoding='gb2312')
    # df = output_df.loc[output_df['task_type'] == module]
    df= pd.DataFrame(columns=output_df.columns)
    df.loc[len(df.index)] = ["beforeQC_txy","data","metadata.tsv","cell","metadata.tsv","metadata","",""]
    for i in input.loc['parameters', 'vars2vis']:
        df.loc[len(df.index)] = ["beforeQC_txy","diagram",f"QC_{i}_beforeQC.tsv","violin",
            f"QC_{i}_beforeQC.tsv",i,f"QC_{i}_beforeQC",f"download/{module}/QC_{i}_beforeQC.png,download/{module}/QC_{i}_beforeQC.pdf"]
    df.to_csv(f"{workdir}/{module}/output.json.tsv",sep='\t',index=False,header=True,encoding='utf-8')

    ## =================================================
    ## 4. link image to download
    ## =================================================
    logger.info("step3.4 构建output/download/分析结果。")
    if not os.path.exists(f'{workdir}/output/download/beforeQC_txy/'):
        os.makedirs(f'{workdir}/output/download/beforeQC_txy/')
    cmd_ln = f"ln -s {workdir}/{module}/{{*.pdf,*png}}  {workdir}/output/download/beforeQC_txy/"
    with module_cmd(f"{input.loc['project', 'environment']}") as p:
        status = p(cmd_ln, projectid, taskid)


    return status
