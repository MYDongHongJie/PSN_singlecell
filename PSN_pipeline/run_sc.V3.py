#coding:utf-8
# A scripts to run multi sample vdj
# Author: Personalbio idsz.Dai
import argparse
import os
from pickle import FALSE
import threading
import shutil
import datetime
import json
import time
import sys
sys.path.append('/PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA')
from utils import *
qcslurm="qc.slurm"
stslurm="st.slurm"
formatslurm = "format.slurm"
#crslurm="count.slurm"
cpslurm="cellrangercp.slurm"



def run_cellranger(sample_info,data_dir,species,nointron,fpid):
    genome_config="/PERSONALBIO/work/singlecell/s00/software/script/1.source/ref.json"
    genome_dic=json.load(open(genome_config))
    db_choose=genome_dic[species]['dir']
    #parse sample info
    sample_file=os.path.abspath(sample_info)
    data_dir=os.path.abspath(data_dir)+"/1.raw_data"
    sample_info_dic={}
    cmds = []
    pids = []
    for line in open(sample_file):
        if not line.startswith("sample"):
            line=line.strip().split()
            sample_info_dic[line[0]]=line[1]
            cmdrun="/PERSONALBIO/work/singlecell/s00/software/script/1.source/cellranger-7.1.0/bin/cellranger count \
                --no-bam --localcores 12 --id={} --fastqs={} --sample={} --transcriptome={} --include-introns=true"\
                .format(line[0],data_dir,line[0],db_choose)
            if nointron:
                cmdrun="/PERSONALBIO/work/singlecell/s00/software/script/1.source/cellranger-7.1.0/bin/cellranger count \
                --no-bam --localcores 12 --id={} --fastqs={} --sample={} --transcriptome={} --include-introns=false"\
                .format(line[0],data_dir,line[0],db_choose)
            crslurm="cr."+line[0]+".slurm"
            slurmdep(crslurm,"crcount",cmdrun,fpid,12,"80G",batch)
            cmd="sbatch {}".format(crslurm)
            pidm = execCmd(cmd)
            pids.append(pidm)
    return(pids,sample_info_dic)

#将需要反馈的cellranger结果拷贝到summary文件夹
def PrePareCellranger(pid,sample_info_dic):
     cmdrun=""
     if os.path.exists("CellrangerOut"):
         pass
     else:
         srun = "mkdir -p CellrangerOut;"
         cmdrun=cmdrun+srun
         for sample in sample_info_dic.keys():
             srun = "mv {} CellrangerOut;".format(sample)
             cmdrun=cmdrun+srun       

     if os.path.exists("summary/02_cellranger"):
        pass
     else:
        srun = "mkdir -p summary/02_cellranger;"
        cmdrun=cmdrun+srun
 
     for sample in sample_info_dic.keys():
         srun = "mkdir -p summary/02_cellranger/{};".format(sample)
         cmdrun=cmdrun+srun
         cmdrun=cmdrun+srun
         srun="cp -r  \
               CellrangerOut/{}/outs/cloupe.cloupe  CellrangerOut/{}/outs/web_summary.html  \
               CellrangerOut/{}/outs/filtered_feature_bc_matrix  CellrangerOut/{}/outs/molecule_info.h5 \
               CellrangerOut/{}/outs/filtered_feature_bc_matrix.h5 summary/02_cellranger/{};".format(\
               sample,sample,sample,sample,sample,sample)
         cmdrun=cmdrun+srun
     srun = "cp sample_info.txt summary/sample_info.txt;"
     cmdrun=cmdrun+srun
     slurmdep(cpslurm,"copy",cmdrun,":".join(list(map(str,pid))),1,"50M",batch)
     cmd1="sbatch {}".format(cpslurm)
     outid = execCmd(cmd1)
     return(outid)

#跑seurat
def RunSeurat(sample_type,mt_cutoff,pid,gather,workdir,double):
    cmdrun = ''
    nsample = len(open(workdir+"/sample_info.txt",'r').readlines())
    srun = "cp /PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/seurat.R summary/seurat.R;cd summary;"
    srun2 = "Rscript seurat.R -t {}  -f {} -g {} ".format(sample_type,mt_cutoff,gather)
    cmdrun = srun + srun2
    if double:
        cmdrun = cmdrun + "-d;"
    memd = nsample*20
    if memd >400:
        memd = 400
    slurmdep(stslurm,"seurat",cmdrun,pid,1,str(memd)+"G",batch)
    run = "sbatch {};".format(stslurm)
    outid = execCmd(run)
    return outid

#准备标准分析报告和结果反馈目录
def format_file(pid,species,nsample,partner):
    cmdrun = "cp /PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/format.sh . && bash format.sh "+species+" "+str(nsample)+" "+partner
    slurmdep(formatslurm,"format",cmdrun,pid,1,"200M",batch)
    cmd = "sbatch {}".format(formatslurm)
    outid = execCmd(cmd)

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('--sample_info', type=str,default="sample_info.txt",help='sample info')
    parser.add_argument('--sample', type=str,default="sample.txt",help='sample txt')
    parser.add_argument('--datamanage', required=True,default='/PERSONALBIO/work/datamanage/Datum/singlecell/2023/',type=str,
                        help='library data path')
    parser.add_argument('--species',required=True,type=str,
                        default='hsa',
                        help='choose species,related to genome and database to use ')
    parser.add_argument('--nointron',action='store_true',default=False,
                        help='whether only run cellranger analysis')
    parser.add_argument('--double',action='store_true',default=False,
                        help='whether only run cellranger analysis')
    parser.add_argument('--batch',required=True,type=str,
                       default='Batch3',
                       help='batch node')
    ##parser.add_argument('--mt',required=True,type=str,
                        #default='^MT',
                        #help='Single cell mt gene pattern')
    parser.add_argument('--mt_cut',required=True,type=str,
                        default='20',
                        help='Single cell mt gene cut off value')
    parser.add_argument('--gather',type=str,
                        default='harmony',
                        help='remove batch method')
    parser.add_argument('--contract',required=True,type=str,
                       help='batch node')
    parser.add_argument('--workdir',required=True,type=str,
                       help='work directory')
    parser.add_argument('--partner',required=True,type=str,
                       default='xxxxxx',
                       help='委托单位')
    parser.add_argument('--step', type=str,help = '1:catfastq,2:cellranger,3:seurat,4:result/report',default="1,2,3,4")

    args = parser.parse_args()
    global batch
    fpid = 1
    cpid = [2]
    ppid = 3
    spid = 4
    batch = args.batch
    nsample =  len(open(args.workdir+"/sample.txt",'r').readlines()) 

    sample_info_dic={}   
    for line in open(args.sample_info):
        if not line.startswith("sample"):
            line=line.strip().split()
            sample_info_dic[line[0]]=line[1]
 
    #cat library fastq data
    if '1' in args.step:
        fpid = CatFastqPerl(args.datamanage,args.contract,args.workdir,args.sample,batch)
        #run QC
        runQC(args.workdir,batch,fpid)
    if '2' in args.step:
        cpid,sample_info_dic = run_cellranger(args.sample_info,args.workdir,args.species,args.nointron,fpid)
        ppid = PrePareCellranger(cpid,sample_info_dic)
    if '3' in args.step:
        if os.path.exists("summary/02_cellranger"):
            pass
        else:
            ppid = PrePareCellranger(cpid,sample_info_dic)
        spid = RunSeurat(args.species,args.mt_cut,ppid,args.gather,args.workdir,args.double)
    if '4' in args.step:
        format_file(":".join(list(map(str,[spid,ppid,fpid]))),args.species,nsample,args.partner)
            

if __name__ == '__main__':
    main()
