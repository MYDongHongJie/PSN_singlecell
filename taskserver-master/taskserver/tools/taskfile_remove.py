from taskserver.tools.parse_json import parse_input
from taskserver.tools.module_cmd import module_cmd
from oebio.utils.log import getLogger
logger = getLogger('oe.cloud.sc.qsub')

def taskfile_remove(module, projectid, taskid, workdir):
    """
    :param module:项目类型
    :param projectid:项目ID
    :param taskid: 任务ID
    :param workdir: 工作根目录
    :return: input.json解析dataframe
    """
    # ===================================================================================================================
    logger.info("移除非必要中间文件")
    if module=="cellranger":
        ## step1
        input = parse_input(projectid, taskid, workdir)
        #
        sample = f"{input.loc['base', 'sample']['name']}"
        cmd =f"rm -rf {workdir}/{module}/{sample}/SC_RNA_COUNTER_CS && " \
             f"rm -rf {workdir}/raw_data/download_from_obs/  "
        with module_cmd(f"OESingleCell/v_3.0.0_cloud") as p:
            status=p(cmd, projectid, taskid)
    elif module=="spaceranger":
        ## step1
        input = parse_input(projectid, taskid, workdir)
        #
        sample = f"{input.loc['base', 'sample']['name']}"
        cmd =f"rm -rf {workdir}/{module}/{sample}/SPATIAL_RNA_COUNTER_CS && " \
             f"rm -rf {workdir}/raw_data/download_from_obs/ " 
        with module_cmd(f"OESingleCell/v_3.0.0_visium") as p:
            status=p(cmd, projectid, taskid)
    elif module=="cellranger_email":
        ## step1
        input = parse_input(projectid, taskid, workdir)
        Projectid = input.loc['project', 'name']
        if input.loc['parameters', 'generate_QC_report']:
            cmd =f"rm -rf {workdir}/{module}/{Projectid}_QC_Report_* " 
            with module_cmd(f"OESingleCell/v_3.0.0_cloud") as p:
                status=p(cmd, projectid, taskid)
    else:
        status=0
    return status 
