import getpass
import logging
import os
import zipfile
from datetime import datetime
from pathlib import Path
from shutil import rmtree
from textwrap import dedent

from pydantic import BaseModel, Field
from subprocess import check_output
from subprocess import run, CalledProcessError
from .settings import SC_MODULE
from .utils.obs import format_key, obs_method

logger = logging.getLogger(__name__)


class TaskAddItem(BaseModel):
    workroot: str
    projectId: str
    taskId: str
    priority: int = Field(0, ge=-1023, le=1024)
    queue: str = 'big'
    slot: int = 4
    cmd: str = 'python scRNA.py'
    image: str = None
    memory: int = Field(default=None, ge=1, le=512)
    callback: str = None
    rootDir: str = None

    def get_run_file(self, wd, i=1, suffix='.sh'):
        file = wd / f'{self.projectId}_{self.taskId}_{i}{suffix}'
        if file.exists():
            return self.get_run_file(wd, i=i+1, suffix=suffix)
        else:
            return file

    def build_cmd(self):
        pid, tid = self.projectId, self.taskId
        wd = Path(f'{self.workroot}20{pid[:2]}/{pid[2:4]}/{pid}/{tid}')
        if wd.exists():  # 已存在文件夹时，请删除对应的文件夹，避免遗留文件干扰
            from shutil import rmtree
            rmtree(wd)
        wd.mkdir(parents=True)
        if self.image:
            cmd = self.build_cci(wd)
        else:
            cmd = self.build_qsub(wd)
        return cmd

    def build_qsub(self, wd):
        run_file = self.get_run_file(wd)
        cmd = f'qsub -V -cwd -pe smp {self.slot} -sync no -N T{self.taskId} ' \
              f'-p {self.priority} -q {self.queue}'.split()
        cmd.extend(['-wd', f'{wd}', f'{run_file}'])
        raw_cmd = self.get_raw_cmd(wd)
        root_dir = self.rootDir or os.environ.get('SCC_OBS_ROOT', default=self.rootDir)
        callback = self.callback or os.environ.get('SCC_CALLBACK', default=self.callback)
        run_file.write_text(f'# 脚本自动生成时间：{datetime.now()}\n'
                            f'# Qsub CMD: {" ".join(cmd)}\n'
                            f'export SCC_CALLBACK={callback}\n'
                            f'export SCC_OBS_ROOT={root_dir}\n'
                            f'export SC_ROOT_DIR={self.workroot}\n'
                            'module purge\n'
                            f'module load {SC_MODULE}\n'
                            f'{raw_cmd}\n')
        return *cmd, cmd[9]

    def run(self):
        try:
            *cmd, output_id = self.build_cmd()
            logger.debug('CMD: %s', ' '.join(cmd))
            env = os.environ.copy()
            env.setdefault('SCC_CALLBACK', self.callback)
            env.setdefault('SCC_OBS_ROOT', self.rootDir)
            env.setdefault('SC_ROOT_DIR', self.workroot)
            output = run(cmd, text=True, capture_output=True, check=True, env=env)
            return {'success': 1, 'error': output.stderr, 'qsubId': output_id}
        except CalledProcessError as e:
            return {'success': 0, 'error': e.stderr + e.stdout, 'qsubId': ''}
        except Exception as e:
            return {'success': 0, 'error': str(e), 'qsubId': ''}

    def build_cci(self, wd):
        name = f'{datetime.now():%Y-%m-%d-%H-%M}-{self.projectId}-{self.taskId}'
        group_id = os.getgid()
        user_id = os.getuid()
        log_io = '{"collectionContainers": ["container-0"]}'
        cmd = self.get_raw_cmd(wd)
        root_dir = self.rootDir or os.environ.get('SCC_OBS_ROOT', default=self.rootDir)
        callback = self.callback or os.environ.get('SCC_CALLBACK', default=self.callback)
        yaml_content = f'''
        # auto created at {datetime.now()}
        apiVersion: batch/v1
        kind: Job
        metadata:
          name: {name}
          annotations:
            description: "单细胞云平台投递任务{self.projectId}-{self.taskId}"
          labels:
            owner: "{getpass.getuser()}"
            customprop: "{self.projectId}"
          namespace: snakemake-scrna
        spec:
          backoffLimit: 0
          completions: 1
          parallelism: 1
          template:
            metadata:
              annotations:
                cri.cci.io/container-type: secure-container
                log.stdoutcollection.kubernetes.io: '{log_io}'
              name: {name}
            spec:
              containers:
              - name: container-0
                image: {self.image}
                env:
                - name: SCC_CALLBACK
                  value: "{callback}"
                - name: SCC_OBS_ROOT
                  value: "{root_dir}"
                - name: SC_ROOT_DIR
                  value: "{self.workroot}"
                command:
                - bash
                - -c
                - '{cmd}'
                workingDir: {wd}
                # 镜像拉取策略，如果指定了:latest标签，则为Always, 否则为IfNotPresent
                imagePullPolicy: Always
                lifecycle: {{}}
                resources:
                  limits:
                    cpu: {self.slot}
                    memory: {self.memory}Gi
                  requests:
                    cpu: {self.slot}
                    memory: {self.memory}Gi
                volumeMounts:
                - mountPath: /data
                  name: cci-sfs-import-data
                  readOnly: true
                - mountPath: /home
                  name: cci-sfs-import-home
                - mountPath: /public/cloud_scRNA
                  name: cci-sfs-import-cloud-scrna
                - mountPath: /public/scRNA_works
                  name: cci-sfs-import-scrna-works
                  readOnly: true
                - mountPath: /public/dev_scRNA
                  name: cci-sfs-import-devscrna
                  readOnly: true

              securityContext: # POD的安全策略，使用与外部用户同样的id，被容器继承
                runAsUser: {user_id}
                runAsGroup: {group_id}
              imagePullSecrets:
              - name: imagepull-secret
              # 生命周期
              restartPolicy: Never
              # POD终止所需的时长（秒）
              terminationGracePeriodSeconds: 30
              activeDeadlineSeconds: {336 * 60 * 60}
              # DNS策略
              dnsPolicy: ClusterFirst
              # 调度
              schedulerName: default-scheduler
              priority: {self.priority}
              volumes:
              - name: cci-sfs-import-data
                persistentVolumeClaim:
                  claimName: cci-sfs-import-data
              - name: cci-sfs-import-home
                persistentVolumeClaim:
                  claimName: cci-sfs-import-home
              - name: cci-sfs-import-cloud-scrna
                persistentVolumeClaim:
                  claimName: cci-sfs-import-cloud-scrna
              - name: cci-sfs-import-scrna-works
                persistentVolumeClaim:
                  claimName: cci-sfs-import-scrna-works
              - name: cci-sfs-import-devscrna
                persistentVolumeClaim:
                  claimName: cci-sfs-import-devscrna
        '''
        runfile = self.get_run_file(wd, suffix='.yaml')
        runfile.write_text(dedent(yaml_content))
        cmd = ['kubectl', 'create', '-f', f'{runfile}']
        return *cmd, name

    def get_raw_cmd(self, wd):
        cmd = f'{self.cmd} -t {self.taskId} -p {self.projectId} -w {wd}'
        return f'{cmd}  >{wd}/run_task.log  2>&1 ' if self.image else cmd


class TaskRemoveItem(BaseModel):
    qsubId: str


def scan_dir(d: Path):
    if d.exists():
        total, *_ = check_output(f'du --summarize {d.resolve()} -b', shell=True, text=True).split()
        mtime = d.stat().st_mtime
    else:
        total, mtime = 0, 0
    return int(total), mtime


class ProjectItem(BaseModel):
    root: Path
    code: str
    task_code: str = None
    sub: bool = False

    @property
    def workdir(self):
        workdir = self.root / f'20{self.code[:2]}/{self.code[2:4]}/{self.code}'
        if self.task_code:
            workdir = workdir / self.task_code
        return workdir

    def dir_stat(self):
        workdir = self.workdir
        total, mtime = scan_dir(workdir)
        sub_dirs = {}
        if self.sub:
            for d in workdir.glob('*'):
                if d.is_dir():
                    t, m = scan_dir(d)
                    sub_dirs[d.stem] = {'name': d.stem, 'path': str(d), 'size': t, 'mtime': m}
        logger.info(f'路径{workdir}中统计信息如下{total}, {mtime}')
        return {'name': self.task_code or self.code, 'path': workdir,
                'size': total, 'mtime': mtime, 'sub_dirs': sub_dirs}


class Sfs2OBSItem(BaseModel):
    source: list[str]
    target: list[str] = None
    op: str = 'sfs2obs'
    obs_root: str

    def to(self):
        if func := getattr(self, self.op, None):
            return func()
        raise NotImplementedError()

    def sfs2obs(self):
        for s, t in zip(self.source, self.target):
            source_path = Path(s)
            if not source_path.exists():
                continue
            zip_file_name = source_path / f'{source_path.name}.zip'
            zip_file_name.unlink(missing_ok=True)
            files = list(source_path.rglob('*'))
            with zipfile.ZipFile(zip_file_name, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zip_handle:
                for file in files:
                    relative_path = file.relative_to(source_path)
                    if str(relative_path).startswith(('input', 'output')):
                        continue
                    zip_handle.write(file, relative_path)
            key = format_key(f'{self.obs_root}{t}/{zip_file_name.name}')
            if obs_method('putFile', objectKey=key, file_path=str(zip_file_name)):
                return {'status': 1, 'msg': f'文件{key}上传失败'}
        # 仅在所有任务执行完成转存OBS后再进行SFS文件删除
        for s in self.source:
            rmtree(s, ignore_errors=True)
        return {'status': 0, 'msg': '已成功由SFS转存OBS并删除SFS对应目录'}

    def obs2sfs(self):
        for s, t in zip(self.source, self.target):
            obs_path = Path(f'{self.obs_root}{s}')
            object_key = format_key(str(obs_path / f'{obs_path.name}.zip'))
            if obs_method('getObjectMetadata', object_key):
                continue
            download_path = Path(t) / f'{obs_path.name}.zip'
            if obs_method('getObject', object_key, download_path):
                continue
            zipfile.ZipFile(download_path).extractall(t)
            download_path.unlink()
            obs_method('deleteObject', object_key)
        return {'status': 0, 'msg': '已全部下载成功'}
