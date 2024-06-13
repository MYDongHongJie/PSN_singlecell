import logging
from subprocess import check_output

logger = logging.getLogger(__name__)


def check_sge_job():
    # 检查SGE的任务状态
    try:
        qstat_output = check_output(['/software/gridengine/bin/lx-amd64/qstat', '-u', '*'], text=True, encoding='utf8')
    except FileNotFoundError:
        return None
    # qstat开始两行和最后一行均是无效行
    for line in qstat_output.split('\n')[2:-1]:
        job_id, prior, name, user, state, *_ = line.split()
        if state == 'Eqw':
            # 获取任务信息
            job_detail = check_output(['qstat', '-j', job_id], text=True)
            file = [n.split(':', 1)[1] for n in job_detail.split('\n') if n.startswith('script_file')][0].strip()
            logger.info('删除JobID:%s', job_id)
            check_output(['qdel', job_id])
            cmd = open(file).readlines()[1][12:].split()
            logger.info('重启任务：%s', cmd)
            check_output(cmd, text=True)
