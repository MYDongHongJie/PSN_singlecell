import os, sys
import json
import pandas as pd
from taskserver.tools.module_cmd import module_cmd
from taskserver.tools.api_upload import api_upload
from oebio.utils.log import getLogger
logger = getLogger('oe.cloud.sc.qsub')
obs_dir = os.environ.get('SCC_OBS_ROOT')  #"obs://scrna-cloud/projects"
def parse_input(projectid, taskid, workdir):
    """
    :param projectid:项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: input.json解析dataframe
    """
    # ===================================================================================================================
    logger.info(f"step1:传回项目开始执行信息-项目号:{projectid}-任务号:{taskid}")
    api_upload(taskid, "start",5)
    # ===================================================================================================================
    logger.info("step2:从obs下拉input.json并解析")

    if not os.path.exists(f'{workdir}/input'):
        os.makedirs(f'{workdir}/input')
    cmd = f"obsutil cp {obs_dir}/{str(projectid)[0:4]}/{projectid}/tasks/{taskid}/input/ {workdir}/input -r -f -flat -vlength -vmd5 -o /public/cloud_scRNA/txy_cp_log "
    with module_cmd("obsutil/5.2.12") as p:
        status=p(cmd, projectid, taskid)
    input = pd.read_json(f'{workdir}/input/input.json', orient="index", dtype={"id": 'str'})
    return input

def save_output(input,output_cfg,projectid, taskid, workdir):
    """
    :param input: input.json解析dataframe
    :param output_cfg: output.json.tsv对应路径
    :param projectid: 项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: status
    """

    logger.info("step4:生成outputfile")  # =============================================================================
    module = input.loc['task', 'type']
    output_file = pd.read_csv(output_cfg, dtype=str, sep="\t", na_filter=False)# encoding='utf-8')
    output_file = output_file.loc[output_file['task_type'] == module]
    output_file_cp = output_file.loc[output_file['result_module'] != "reportVars"]
    #output_file_cp = output_file_cp.loc[output_file['type'] != "image"]
    output_file_cp = output_file_cp.loc[output_file['type'] != "link"]
    if not os.path.exists(f'{workdir}/output/download'):
        os.makedirs(f'{workdir}/output/download')
    if len(output_file_cp) > 0:
        for i in output_file_cp.index :
            if not os.path.exists(f"{workdir}/output/{output_file_cp.loc[i, 'file']}"):
                cmd1 = f"cp -rf {workdir}/{module}/{output_file_cp.loc[i, 'input']} " \
                    f" {workdir}/output/{output_file_cp.loc[i, 'file']}"
                with module_cmd("obsutil/5.2.12") as p:
                    status=p(cmd1, projectid, taskid)
            else:
                status=0
    else:
        status=0

    logger.info("step5:生成output.json")  # =============================================================================
    project = json.loads(input.loc["project",].dropna().to_json(orient='index'))
    task = json.loads(input.loc["task",].dropna().to_json(orient='index'))
    task['success'] = 1 if status == 1 else 0
    ##添加project,task信息
    output_json = {}
    output_json["project"] = project
    output_json["task"] = task
    # ==================================================================================================================

    output_diagram = output_file.loc[output_file['result_module'] == "diagram"]
    if len(output_diagram) > 0:
        logger.info("添加diagram信息")
        diagram = []
        for index in output_diagram.index:
            diagram.append({"type": f"{output_diagram.loc[index, 'type']}",
                            "file": f"{output_diagram.loc[index, 'file']}",
                            "title": f"{output_diagram.loc[index, 'title']}",
                            "downloadName": f"{output_diagram.loc[index, 'downloadName']}",
                            "downloadPath": f"{output_diagram.loc[index, 'downloadPath']}".split(",")})
        output_json["diagram"] = diagram

    # ==================================================================================================================
    output_summary = output_file.loc[output_file['result_module'] == "summary"]
    if len(output_summary) > 0:
        logger.info("添加summary信息")
        summary = []
        for index in output_summary.index:
            summary.append({"type": f"{output_summary.loc[index, 'type']}",
                            "file": f"{output_summary.loc[index, 'file']}",
                            "title": f"{output_summary.loc[index, 'title']}"})
        output_json["summary"] = summary
    # ==================================================================================================================
    output_data = output_file.loc[output_file['result_module'] == "data"]
    if len(output_data) > 0:
        logger.info("添加data信息")
        data = []
        for index in output_data.index:
            data.append({"type": f"{output_data.loc[index, 'type']}",
                         "file": f"{output_data.loc[index, 'file']}",
                         "title": f"{output_data.loc[index, 'title']}"})
        output_json["data"] = data
    # ==================================================================================================================
    output_reportVars = output_file.loc[output_file['result_module'] == "reportVars"]
    if len(output_reportVars) > 0:
        logger.info("添加reportVars信息")
        reportVars = []
        for index in output_reportVars.index:
            reportVars.append(output_reportVars.loc[index, 'input'])
        output_json["reportVars"] = reportVars
    # ==================================================================================================================
    # json.dumps(output_json,indent=4)
    with open(f"{workdir}/output/output.json", "w") as outfile:
        json.dump(output_json, outfile, indent=4, ensure_ascii=False)
    # ===================================================================================================================
    logger.info("step6: 上传数据 and 返回运行成功信息")
    cmd2 = f"obsutil cp {workdir}/output  {obs_dir}/{str(projectid)[0:4]}/{projectid}/tasks/{taskid}/output  -link -r -f -flat -vlength -vmd5 -o /public/cloud_scRNA/txy_cp_log "
    with module_cmd("obsutil/5.2.12") as p:
        status = p(cmd2, projectid, taskid)
    if status == 0 :
        api_upload(taskid, "success",5)
