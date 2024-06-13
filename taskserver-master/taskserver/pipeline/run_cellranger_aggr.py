import os, sys, shutil
import pandas as pd

# import json,subprocess,requests
# from shutil import copyfile
from oebio.utils.log import getLogger
from taskserver.tools.module_cmd import module_cmd
from selenium import webdriver
from taskserver.tools.parse_json import parse_input
from taskserver.tools.api_upload import api_upload

logger = getLogger('oe.cloud.sc.qsub')
def task_cellranger_aggr(input, output_cfg, projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """
    module = input.loc['task', 'type']
    proj_workdir = f"{workdir}/../"
    if not os.path.exists(f'{workdir}/{module}'):
        os.makedirs(f'{workdir}/{module}/cellranger')

    # ==================================================================================================================
    logger.info("step3:解析命令行，执行分析")
    ## =================================================
    ## 1. aggr
    ## =================================================
    logger.info("step3.1:cellranger aggr")
    cellranger_version = input.loc['parameters', 'cellranger_version']
    if input.loc["parameters", "aggr"]:
        # create library.csv
        libraries_index = 'sample_id'  if cellranger_version>="7.0.1" else  "library_id"
        libraries = pd.DataFrame(columns=[libraries_index, 'molecule_h5'])
        ## 基于前一步beforeQC_txy中cellranger来构建
        input_before = pd.read_json(f"{workdir}/../../{input.loc['base', 'tasks'][0]['projectId']}/{input.loc['base', 'tasks'][0]['taskId']}/input/input.json", orient="index", dtype={"id": 'str'})
        if input_before.loc["task", "type"]=="beforeQC_txy":
            if len(input_before.loc["base", "tasks"]) > 1:
                for task in input_before.loc["base", "tasks"]:
                    sample = next(os.walk(f"{workdir}/../../../../20{task['projectId'][0:2]}/{task['projectId'][2:4]}/{task['projectId']}/{task['taskId']}/cellranger"))[1][0]
                    libraries.loc[sample, libraries_index] = sample
                    libraries.loc[sample, 'molecule_h5'] = f"{workdir}/../../../../20{task['projectId'][0:2]}/{task['projectId'][2:4]}/{task['projectId']}/{task['taskId']}/cellranger/{sample}/outs/molecule_info.h5"
                    libraries.to_csv(f"{workdir}/{module}/libraries.csv", encoding='utf-8', index=False)
                # run aggr
                with module_cmd(f"cellRanger/{input.loc['parameters', 'cellranger_version']}") as p:
                    # with module_cmd(f"cellRanger/5.0.0") as p:
                    status = p( f"export PATH=/opt/softwares/cellranger-{cellranger_version}:$PATH && "
                                f"cd  {workdir}/{module}/ &&  rm -rf aggr &&  cellranger aggr --id=aggr --csv={workdir}/{module}/libraries.csv --normalize=none",
                        projectid, taskid)
                # status=p(f"cd  {workdir}/{module}/ &&  touch QC_nCount_RNA_beforeQC.tsv QC_nFeature_RNA_beforeQC.tsv QC_percent.mito_beforeQC.tsv metadata.tsv ", projectid, taskid)
                # =================================================
                # 3. amending output.json
                # =================================================
                logger.info("step3.3 生成output.json.tsv。")
                output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)  # , encoding='gb2312')
                df = output_df.loc[output_df['task_type'] == module]
                # df = pd.DataFrame(columns=output_df.columns)
                # df.loc[len(df.index)] = ["beforeQC", "data", "metadata.tsv", "cell", "metadata.tsv", "metadata", "", ""]
                # for i in input.loc['parameters', 'vars2vis']:
                #     df.loc[len(df.index)] = ["beforeQC", "diagram", f"QC_{i}_beforeQC.tsv", "violin",
                #                              f"QC_{i}_beforeQC.tsv", i, f"QC_{i}_beforeQC",
                #                              f"download/{module}/QC_{i}_beforeQC.png,download/{module}/QC_{i}_beforeQC.pdf"]
                df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t', index=False, header=True, encoding='utf-8')

                ## =================================================
                ## 4. link image to download
                ## =================================================
                logger.info("step3.4 构建output/download/分析结果。")
                if not os.path.exists(f'{workdir}/output/download/aggr/'):
                    os.makedirs(f'{workdir}/output/download/aggr/')
                # cmd_ln = f"ln -s {workdir}/{module}/{{*.pdf,*png}}  {workdir}/output/download/beforeQC/"
                # with module_cmd(f"{input.loc['project', 'environment']}") as p:
                #     status = p(cmd_ln, projectid, taskid)
                cmd_ln_aggr = f"ln -s {workdir}/{module}/aggr/outs  {workdir}/output/download/aggr"
                with module_cmd(f"{input.loc['project', 'environment']}") as p:
                    status = p(cmd_ln_aggr, projectid, taskid)
            else:
                logger.info("跳过 aggr（只存在一个样本）")
                api_upload(taskid, "success", 5)
                os._exit()
        elif input_before.loc["task", "type"] == "cellranger":
            if len(input.loc["base", "tasks"]) > 1:
                for task in input.loc["base", "tasks"]:
                    sample = next(os.walk(f"{workdir}/../../../../20{task['projectId'][0:2]}/{task['projectId'][2:4]}/{task['projectId']}/{task['taskId']}/cellranger"))[1][0]
                    libraries.loc[sample, libraries_index] = sample
                    libraries.loc[sample, 'molecule_h5'] = f"{workdir}/../../../../20{task['projectId'][0:2]}/{task['projectId'][2:4]}/{task['projectId']}/{task['taskId']}/cellranger/{sample}/outs/molecule_info.h5"
                    libraries.to_csv(f"{workdir}/{module}/libraries.csv", encoding='utf-8', index=False)
                # run aggr
                with module_cmd(f"cellRanger/{input.loc['parameters', 'cellranger_version']}") as p:
                    # with module_cmd(f"cellRanger/5.0.0") as p:
                    status = p( f"export PATH=/opt/softwares/cellranger-{cellranger_version}:$PATH && "
                                f"cd  {workdir}/{module}/ &&  rm -rf aggr &&  cellranger aggr --id=aggr --csv={workdir}/{module}/libraries.csv",
                        projectid, taskid)
                # status=p(f"cd  {workdir}/{module}/ &&  touch QC_nCount_RNA_beforeQC.tsv QC_nFeature_RNA_beforeQC.tsv QC_percent.mito_beforeQC.tsv metadata.tsv ", projectid, taskid)
                # =================================================
                # 3. amending output.json
                # =================================================
                logger.info("step3.3 生成output.json.tsv。")
                output_df = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)  # , encoding='gb2312')
                df = output_df.loc[output_df['task_type'] == module]
                # df = pd.DataFrame(columns=output_df.columns)
                # df.loc[len(df.index)] = ["beforeQC", "data", "metadata.tsv", "cell", "metadata.tsv", "metadata", "", ""]
                # for i in input.loc['parameters', 'vars2vis']:
                #     df.loc[len(df.index)] = ["beforeQC", "diagram", f"QC_{i}_beforeQC.tsv", "violin",
                #                              f"QC_{i}_beforeQC.tsv", i, f"QC_{i}_beforeQC",
                #                              f"download/{module}/QC_{i}_beforeQC.png,download/{module}/QC_{i}_beforeQC.pdf"]
                df.to_csv(f"{workdir}/{module}/output.json.tsv", sep='\t', index=False, header=True, encoding='utf-8')

                ## =================================================
                ## 4. link image to download
                ## =================================================
                logger.info("step3.4 构建output/download/分析结果。")
                if not os.path.exists(f'{workdir}/output/download/aggr/'):
                    os.makedirs(f'{workdir}/output/download/aggr/')
                # cmd_ln = f"ln -s {workdir}/{module}/{{*.pdf,*png}}  {workdir}/output/download/beforeQC/"
                # with module_cmd(f"{input.loc['project', 'environment']}") as p:
                #     status = p(cmd_ln, projectid, taskid)
                cmd_ln_aggr = f"ln -s {workdir}/{module}/aggr/outs  {workdir}/output/download/aggr"
                with module_cmd(f"{input.loc['project', 'environment']}") as p:
                    status = p(cmd_ln_aggr, projectid, taskid)
            else:
                logger.info("跳过 aggr（只存在一个样本）")
                api_upload(taskid, "success", 5)
                os._exit()
    else:
        logger.info("跳过 aggr（aggr选择为False）")
        api_upload(taskid, "success", 5)
        os._exit()
    return status
