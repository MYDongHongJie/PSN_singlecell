import logging
from fastapi import FastAPI
from subprocess import run

from .items import TaskAddItem, TaskRemoveItem, ProjectItem, Sfs2OBSItem

logger = logging.getLogger(__name__)
app = FastAPI(title='scRNA qsub Server', description="欧易单细胞云平台qsub服务接口",
              version='0.2.0', )


@app.on_event('startup')
async def load_schedule():
    from apscheduler.jobstores.memory import MemoryJobStore
    from apscheduler.schedulers.asyncio import AsyncIOScheduler
    from taskserver.utils import scheduler_jobs
    from taskserver.settings import JOB_INTERVAL

    schedule = AsyncIOScheduler(jobstores={'default': MemoryJobStore()})
    schedule.add_job(scheduler_jobs.check_sge_job, trigger='interval', minutes=JOB_INTERVAL, id='check_sge_job')
    schedule.start()


@app.post('/tasks/add')
async def add_task(task: TaskAddItem):
    return task.run()


@app.post('/tasks/kill')
async def remove_task(task: TaskRemoveItem):
    logger.debug(f"CMD: kill {task.qsubId}")
    if '-' in task.qsubId:
        output = run(['kubectl', 'delete', 'jobs', '-n', 'snakemake-scrna', task.qsubId])
    else:
        output = run(['qdel', task.qsubId], text=True, capture_output=True)
    logger.debug(f'output: {output}')
    return {'success': 1 if output.returncode == 0 else 1, 'error': output.stderr}


@app.post('/project/stat',)
async def project_stat(payload: ProjectItem):
    """统计项目根路径下的文件夹大小等相关信息"""
    return payload.dir_stat()


@app.post('/data_operations')
async def sfs2obs(payload: Sfs2OBSItem):
    return payload.to()
