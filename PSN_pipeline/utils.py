import os
import subprocess
import datetime
def slurm(ffile,jobname,cmd,coren,mem,batch):
    fp = open(ffile,"w")
    fp.writelines("#!/bin/bash\n")
    fp.writelines("#SBATCH   -p {}\n".format(batch))
    fp.writelines("#SBATCH -N 1\n")
    fp.writelines("#SBATCH -n 1\n")
    fp.writelines("#SBATCH -c {}\n".format(coren))
    fp.writelines("#SBATCH --mem={}\n".format(mem))
    fp.writelines("#SBATCH   --job-name={}\n".format(jobname))
    fp.writelines("#SBATCH   --output={}.%j.out\n".format(jobname))
    fp.writelines("#SBATCH   --error={}.%j.err\n".format(jobname))
    fp.writelines("#SBATCH   --parsable\n")
    fp.writelines(cmd)
    fp.close()

def slurmdep(ffile,jobname,cmd,pid,coren,mem,batch):
    fp = open(ffile,"w")
    fp.writelines("#!/bin/bash\n")
    fp.writelines("#SBATCH   -p {}\n".format(batch))
    fp.writelines("#SBATCH -N 1\n")
    fp.writelines("#SBATCH -n 1\n")
    fp.writelines("#SBATCH -c {}\n".format(coren))
    fp.writelines("#SBATCH --mem={}\n".format(mem))
    fp.writelines("#SBATCH   --job-name={}\n".format(jobname))
    fp.writelines("#SBATCH   --output={}.%j.out\n".format(jobname))
    fp.writelines("#SBATCH   --error={}.%j.err\n".format(jobname))
    fp.writelines("#SBATCH   --parsable\n")
    fp.writelines("#SBATCH   --dependency=afterok:{}:+5\n".format(pid))
    fp.writelines(cmd)
    fp.close()


def execCmd(cmd):
    try:
        print("命令{}开始运行{}".format(cmd, datetime.datetime.now()))
        pid=int(str(subprocess.check_output(cmd, shell=True),'utf-8').replace("\n",""))
        print(pid)
        print("命令{}结束运行{}" .format(cmd, datetime.datetime.now()))
        return pid
    except:
        print("{}\t运行失败".format(cmd))

def CatFastqPerl(datamanage,contract,workdir,sampleinfo,batch):
    #sample_type={}
    downdir=datamanage + '*/' + contract+ '/'+'*'
    #cmd = ''
    cmd = 'sh /PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/FastqPrepare.sh {} {}'.format(downdir,sampleinfo)
    
    # scmd1 = '/PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/stdpipe/bin/perl \
    #        /PERSONALBIO/work/singlecell/s00/software/script/1.source/stdpipe/stdpipeV3/ClassifyRawdataToSampleID.pl \
    #        -m {} -r "{}" -o {} >rawfastq.sh;'.format(sampleinfo,str(downdir),workdir)
    # cmd = cmd + scmd1
    # for line in open(sampleinfo):
    #     if not line.startswith("sample"):
    #         line=line.strip().split()
    #         sample_type[line[0]]=line[1]
    # for sample in sample_type.keys():
    #     os.system('mkdir -p {}/1.raw_data/{}'.format(workdir,sample))

    # nn = r''' awk -F "\t" '{print "cat " $3 ">>1.raw_data/"$1"/"$1"_S1_L001_R1_001.fastq.gz";print "cat " $4 ">>1.raw_data/"$1"/"$1"_S1_L001_R2_001.fastq.gz" }' '''
    # scmd2 = '{}  raw_fastq.txt > merge.sh && sh merge.sh'.format(nn)
    # cmd = cmd + scmd2

    fqslurm = '{}/PrepareFastq.slurm'.format(workdir)
    slurm(fqslurm,'PrepareFastq',cmd,2,"100M",batch)
    sbatch ='sbatch {}'.format(fqslurm)
    pid = execCmd(sbatch)
    return pid

def runQC(workdir,batch,pid):
    data_dir=os.path.abspath(workdir)
    cmdrun="/PERSONALBIO/work/singlecell/s00/software/FastQC/fastqc -t 6 -o QC {}/1.raw_data/*/*gz".format(data_dir)
    qcslurm = "{}/qc.slurm".format(workdir)
    slurmdep(qcslurm,"runQC",cmdrun,pid,6,"6G",batch)
    cmd="mkdir QC && sbatch {}".format(qcslurm)
    execCmd(cmd)
