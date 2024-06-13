#!/usr/bin/env python3
# encoding: utf-8
# Author liuling
# Modified luyao
# Desc   : This script is used to push SOR recird
# coding: utf-8

import os,time, glob ,shutil,yaml, argparse
import requests
import getpass
import json
import numpy as np
import pandas as pd
import re
from oebio.tools.oe_sor_upload import oe_sor_upload
from collections import OrderedDict

#=======================================================================================================================
def project(compact,yaml_file,user,pro_type,info_type,logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    logfile=os.path.join(logdir,"project_upload_success.txt")
    if not os.path.exists(logfile):
        desc = yaml_file["report"].copy()
        desc.pop("Reference")
        desc["Rawdata_Path"] = ','.join(map(str, yaml_file["cellranger_params"]["raw_data_obs_address"]))
        oe_sor_upload(compact,
                      pro_type,
                      info_type,
                      desc,
                      sample_name="",
                      operator=user)
        os.mknod(logfile)
    else:
        print("Project information has been upload !Skipping!")
#=======================================================================================================================
def sample(compact,metadata,yaml_file,user,pro_type,info_type,logdir):
    '''
    :param compact:
    :param metadata:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    samples = pd.read_csv(metadata, dtype=str).set_index("sampleid", drop=False)
    samples.index.names = ["sample_id"]
    for (sample) in samples.index:
        logfile=os.path.join(logdir,sample+"sample_upload_success.txt")
        if not os.path.exists(logfile):
            desc=dict()
            desc["Sample_Type"]=yaml_file["report"]["Species"]
            desc["Group1"]=samples.loc[(sample), ["group"]][0]
            desc["Batch_ID"]=samples.loc[(sample), ["batch"]][0]
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name=sample,
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Sample information for {sample} has been upload !Skipping!")
#=======================================================================================================================
def QC(compact,metadata,yaml_file,user,pro_type,info_type,logdir):
    '''
    :param compact:
    :param metadata:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    ## get sample names
    if type(yaml_file["cellranger_params"]["samples"]) == list:
        if all(isinstance(i, str) for i in yaml_file["cellranger_params"]["samples"]):
            sample_list = yaml_file["cellranger_params"]["samples"]
        elif all(isinstance(i, dict) for i in yaml_file["cellranger_params"]["samples"]):
            sample_list = [v for i in yaml_file["cellranger_params"]["samples"] for v,k in i.items()]
        samples = pd.Series(sample_list,index = sample_list)
    elif type(yaml_file["cellranger_params"]["samples"]) == str and yaml_file["cellranger_params"]["samples"] == "all" :
        samples = pd.read_csv(metadata,dtype=str).set_index("sampleid",drop=False)
        samples.index.names = ["sample_id"]
    else :  sys.exit("cellranger_params - samples can only be string 'all' or a list of sample names.")
    ## main
    for (sample) in samples.index:
        logfile=os.path.join(logdir,sample+"_QC_upload_success.txt")
        if not os.path.exists(logfile):
            os.system(f"cp result/cellranger/{sample}/outs/web_summary.html /public/dev_scRNA/documents_scRNA_dev/5.项目信息/10x_mutiltome/6.常见问题/cellranger_arc_results/%s_{sample}_web_summary.html"  %(compact))
            QC=pd.read_csv(f"result/cellranger/{sample}/outs/summary.csv", dtype=str)
            desc=pd.DataFrame(QC.values.T, index= [i.title() for i in list(QC.columns.str.replace(" ", "_"))], columns=QC.index).to_dict()[0]
            desc["Gex_Q30_Bases_In_Read_2"] = desc.pop("Gex_Q30_Bases_In_Read_2")
            desc["Gex_Reads_Mapped_Confidently_To_Transcriptome(%)"] = desc.pop("Gex_Reads_Mapped_Confidently_To_Transcriptome")
            desc["Gex_Fraction_Of_Transcriptomic_Reads_In_Cells(%)"] = desc.pop("Gex_Fraction_Of_Transcriptomic_Reads_In_Cells")
            desc["Gex_Reads_Mapped_Confidently_To_Genome(%)"] = desc.pop("Gex_Reads_Mapped_Confidently_To_Genome")
            desc["Gex_Clean_Base"] = str(format(int(desc["Gex_Sequenced_Read_Pairs"].replace(",",""))*300 / 1000000000, '.2f')) + "G"
            desc["Atac_Confidently_Mapped_Read_Pairs%)"] = desc.pop("Atac_Confidently_Mapped_Read_Pairs")
            desc["Atac_Fraction_Of_High-Quality_Fragments_Overlapping_Peaks(%)"] = desc.pop("Atac_Fraction_Of_High-Quality_Fragments_Overlapping_Peaks")
            desc["Atac_Tss_Enrichment_Score"] = desc.pop("Atac_Tss_Enrichment_Score")
            desc=OrderedDict(desc)
            desc.update({'QC__Type':"multimodal"})
            desc.move_to_end('QC__Type', last=False)
            desc=dict(desc)
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name=sample,
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"QC  information for {sample} has been upload !Skipping!")
#=======================================================================================================================
def database(compact,yaml_file,user,pro_type,info_type,logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    ## reference genome DB
    for i in yaml_file["database_url"].keys():
        logfile=os.path.join(logdir,"Genome_Database_upload_success.txt")
        if not os.path.exists(logfile):
            desc = dict()
            desc["Database"] = i
            desc["Linkage"] = yaml_file["database_url"][i]["Linkage"]
            desc["Version"] = yaml_file["database_url"][i]["Version"]
            oe_sor_upload(
                compact,
                pro_type,
                info_type,
                desc,
                sample_name="",
                operator=user)
            os.mknod(logfile)
        else:
            print(f"Database information for reference genome has been upload !Skipping!")

    ## report-related DB
    # if yaml_file["module"]["report"] == True:
    #     reportDB_file = "scripts/src/database.txt"
    #     reportDB = pd.read_table(reportDB_file, dtype=str).set_index("使用数据库", drop=False)
    #     reportDB.to_dict()
    #     for i in reportDB.index:
    #         logfile=os.path.join(logdir,i+"_Database_upload_success.txt")
    #         if not os.path.exists(logfile):
    #             desc=dict()
    #             desc["Database"]=i
    #             desc["Linkage"]=reportDB.loc[i,"网页链接"]
    #             desc["Version"]=reportDB.loc[i,"版本"]
    #             oe_sor_upload(compact,
    #                           pro_type,
    #                           info_type,
    #                           desc,
    #                           sample_name="",
    #                           operator=user)
    #             os.mknod(logfile)
    #         else:
    #             print(f"Database information for {i} has been upload !Skipping!")
#=======================================================================================================================
def software(compact,yaml_file,user,pro_type,info_type,logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    if yaml_file["module"]["cellranger"] == True:
        logfile=os.path.join(logdir,"cellranger_software_upload_success.txt")
        if not os.path.exists(logfile):
            desc=dict()
            desc["Software"] = yaml_file['cellranger_params']['envmodules'].split("/")[0]
            desc["Version"] = yaml_file['cellranger_params']['envmodules'].split("/")[1]
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name="",
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Software information for Cellranger has been upload !Skipping!")
    if yaml_file["module"]["report"] == True:
        software_table = pd.read_csv("scripts/src/software.txt", dtype=str).set_index("Software", drop=False)
        software_version = software_table.to_dict()['Version']
        software_version.pop('Cell Ranger ARC')
        for i in software_version.keys():
            logfile=os.path.join(logdir,i+"_software_upload_success.txt")
            if not os.path.exists(logfile):
                desc=dict()
                desc["Software"]=i
                desc["Version"]=software_version[i]
                oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name="",
                          operator=user)
                os.mknod(logfile)
            else:
                print(f"Software information for {i} has been upload !Skipping!")

#=======================================================================================================================
def benchmark(compact,yaml_file,user,pro_type,info_type,logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
        #: 1. Running time in seconds  s
        #: 2. Running time in h:m:s   h:m:s
        #: 3. Maximal RSS in MB   max_rss
        #: 4. Maximal VMS in MB   max_vms
        #: 5. Maximal USS in MB   max_uss
        #: 6. Maximal PSS in MB   max_pss
        #: 7. I/O read in bytes   io_in
        #: 8. I/O written in bytes  io_out
        #: 9. mean_load python.psutil.cpu_percent
            #: Count of CPU seconds, divide by running time to get mean load estimate
            self.cpu_seconds = cpu_seconds or 0
            #: First time when we measured CPU load, for estimating total running time
            self.first_time = None
            #: Previous point when measured CPU load, for estimating total running time
            self.prev_time = None
     '''
    ## extract rules
    rule_files=glob.glob("rules/*smk")
    rules={}
    rule_name=""
    for i in rule_files:
        with open(i,'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("rule "):
                    ll = re.split('\s+|:',line)
                    rule_name = ll[1]
                if "qsub_p=" in line:
                    qsub_p = re.sub('(qsub_p=|,)',"",line)
                    if str(qsub_p).isdigit():
                        qsub_p = str(qsub_p)
                        rules[rule_name] = qsub_p
                    else:
                        qsub_p = qsub_p.replace('config','yaml_file')
                        rules[rule_name] = qsub_p
    ## benchmarks
    list=glob.glob("benchmarks/*/*txt")
    for i in list:
        filename=i.split('/')[-1].replace('.benchmark.txt','').split('.')
        logfile=os.path.join(logdir,"-".join(filename)+"_cpu_mem_upload_success.txt")
        if not os.path.exists(logfile):
            data =pd.read_table(i,sep='\s+', header=0)
            Run_Time = round(float(data['s']) / 3600, 2)
            Memory = round(float(data['max_rss']) / 1024, 2)
            CPU_Load = float(data['mean_load'])
            if len(filename)==1:
                sample_name="project_level"
                Analysis_Module=filename[0]
            if len(filename)==2:
               sample_name=filename[0]
               Analysis_Module = filename[1]
            if Analysis_Module not in rules.keys():
                print("Couldn't find cpu info for rule " + Analysis_Module )
                continue
            CPU=eval(rules[Analysis_Module])
            desc={"Analysis__Module": Analysis_Module,
                "Run_Time": Run_Time ,
                "CPU": CPU,
                "CPU_Load": CPU_Load,
                "Memory": Memory }
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name=sample_name,
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Cpu_Mem information for {'-'.join(filename)} has been upload !Skipping!")
#=======================================================================================================================
parser = argparse.ArgumentParser(description="单细胞多模态项目SOR自动上传")
parser.add_argument("-i", "--config",default="config/config.yaml")
parser.add_argument("-m", "--metadata",default="config/metadata.csv")
parser.add_argument("-l", "--logdir",default="logs/upload_sor", help="sor upload successed log dir")
parser.add_argument("-s", "--step",default="Project,Sample,QC,Ref_Database,Software,Cpu_Mem", help="sor upload modules")
args = parser.parse_args()
config = os.path.abspath(args.config)

if __name__ == '__main__':
    yaml_file = yaml.load(open(config),Loader=yaml.FullLoader)
    compact = yaml_file["report"]["Project_Num"]
    pro_type = "单细胞多模态"
    user = getpass.getuser()
    if not os.path.exists(args.logdir):
        os.makedirs(args.logdir)
    if "Project" in args.step:
        project(compact,yaml_file ,user,pro_type,'Project',args.logdir)
    if "Sample" in args.step:
        sample(compact,args.metadata,yaml_file,user,pro_type,"Sample",args.logdir)
    if "QC" in args.step:
        QC(compact,args.metadata,yaml_file,user,pro_type,"QC",args.logdir)
    if "Ref_database" in args.step:
        database(compact,yaml_file ,user,pro_type,'Ref_Database',args.logdir)
    if "Software" in args.step:
        software(compact,yaml_file ,user,pro_type,'Softwares',args.logdir)
    if "Cpu_Mem" in args.step:
        benchmark(compact,yaml_file ,user,pro_type,'Cpu_Mem',args.logdir)
