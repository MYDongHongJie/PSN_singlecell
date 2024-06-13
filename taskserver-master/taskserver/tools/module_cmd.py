import os, sys
import subprocess
import contextlib

from oebio.utils.log import getLogger
from taskserver.tools.api_upload import api_upload
logger = getLogger('oe.cloud.sc.qsub')
# cmd = 'cat /proc/1/cgroup | grep -qi docker     && echo "Docker" || echo "Not Docker"'
# output = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
# print(output.stdout.decode('utf-8'))
if  os.path.exists("/software/gridengine/bin/lx-amd64/qstat") :
    file="/data/software/modules/modules-v4.2.1/init/python.py"
    modulepath="/data/software/modulefiles/"
    exec(open(file).read())
    os.environ['MODULEPATH'] = os.path.abspath(modulepath)
else:
    logger.info("using software in docker images")

# =======================================================================================================================
@contextlib.contextmanager
def module_cmd(*softwares):

    def run_cmd(command, projectid, taskid):
        logger.info(f"运行命令：{command}")
        ret = subprocess.run(command,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             encoding="utf-8")
        logger.info(ret.stdout)
        if ret.returncode != 0:
            logger.info(f"{command}:error[{ret}]")
            logger.info(f"{ret.stderr}")
            state="fail"
#            api_upload(taskid,state,5)
            sys.exit()
        return ret.returncode

    for s in softwares:
        if  os.path.exists("/software/gridengine/bin/lx-amd64/qstat") :
            module('load', s)
        else:
            continue
    try:
        yield run_cmd
    finally:
        print('任务执行完毕')
        for s in softwares:
            if os.path.exists("/software/gridengine/bin/lx-amd64/qstat"):
                module('unload', s)
            else:
                continue
