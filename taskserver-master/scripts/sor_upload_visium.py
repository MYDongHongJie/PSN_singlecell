#!/usr/bin/env python3
# encoding: utf-8
# Author liuling
# Modified xfyang
# Date: 2021/01/22 13:31
# Desc   : This script is used to push SOR recird
# coding: utf-8

import os, time, glob, shutil, yaml, argparse
import requests
import getpass
import json
import numpy as np
import pandas as pd
from oebio.tools.oe_sor_upload import oe_sor_upload
from collections import OrderedDict


# =======================================================================================================================
def project(compact, yaml_file, user, pro_type, info_type, logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    logfile = os.path.join(logdir, "project_upload_success.txt")
    if not os.path.exists(logfile):
        desc = yaml_file["report"]
        desc["Rawdata_Path"] = ','.join(map(str, yaml_file["report"]["Rawdata_Path"]["raw_data_obs_address"])) \
                               + "||" \
                               + ','.join(map(str, yaml_file["report"]["Rawdata_Path"]["image_obs_address"]))
        oe_sor_upload(compact,
                      pro_type,
                      info_type,
                      desc,
                      sample_name="",
                      operator=user)
        os.mknod(logfile)
    else:
        print("Project information has been upload !Skipping!")


# =======================================================================================================================
def sample(compact, metadata, user, pro_type, info_type, logdir):
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
        logfile = os.path.join(logdir, sample + "sample_upload_success.txt")
        if not os.path.exists(logfile):
            desc = dict()
            desc["Sample_Type"] = samples.loc[(sample), ["species"]][0]
            desc["Group1"] = samples.loc[(sample), ["group"]][0]
            desc["Batch_ID"] = samples.loc[(sample), ["batchid"]][0]
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name=sample,
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Sample information for {sample} has been upload !Skipping!")


# =======================================================================================================================
def QC(compact, metadata, user, pro_type, info_type, logdir):
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
        logfile = os.path.join(logdir, sample + "QC_upload_success.txt")
        if not os.path.exists(logfile):
            QC = pd.read_csv(f"{logdir}/spaceranger/{sample}/outs/metrics_summary.csv")
            QC = QC.drop(['Sample ID'], axis=1)
            Clean_Base=(QC[["Number of Reads"]].values[0]*91/2+QC[["Number of Reads"]].values[0]*28/2)/1e9
            Clean_Base=round(float(Clean_Base),2)
            QC.insert(0, 'Clean_Base', f'{Clean_Base}G')
            desc = \
            pd.DataFrame(QC.values.T, index=QC.columns.str.replace(" ", "_").str.title(), columns=QC.index).to_dict()[0]
            desc = OrderedDict(desc)
            desc.update({'QC__Type': "spaceranger"})
            desc.move_to_end('QC__Type', last=False)
            desc = dict(desc)
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name=sample,
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"QC  information for {sample} has been upload !Skipping!")


# =======================================================================================================================
def database(compact, yaml_file, user, pro_type, info_type, logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    for i in yaml_file["database_url"].keys():
        logfile = os.path.join(logdir, i + "Ref_Database_upload_success.txt")
        if not os.path.exists(logfile):
            desc = dict()
            desc["Database"] = i
            desc["Linkage"] = yaml_file["database_url"][i]["Linkage"]
            desc["Version"] = yaml_file["database_url"][i]["Version"]
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name="",
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Database information for {i} has been upload !Skipping!")


# =======================================================================================================================
def software(compact, yaml_file, user, pro_type, info_type, logdir):
    '''
    :param compact:
    :param yaml_file:
    :param user:
    :param pro_type:
    :param info_type:
    :param logdir:
    :return:
    '''
    for i in yaml_file["softwares"].keys():
        logfile = os.path.join(logdir, i + "software_upload_success.txt")
        if not os.path.exists(logfile):
            desc = dict()
            desc["Software"] = i
            desc["Version"] = yaml_file["softwares"][i]["Version"]
            # desc["Parameter"]=yaml_file["softwares"][i]["Parameter"]
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name="",
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Software information for {i} has been upload !Skipping!")


# =======================================================================================================================
def benchmark(compact, yaml_file, user, pro_type, info_type, logdir):
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
    list = glob.glob("benchmarks/*txt")
    for i in list:
        filename = i.replace(".benchmark.txt", "").replace("benchmarks/", "").split(sep="_for_")
        if len(filename) == 1:
            logfile = os.path.join(logdir, filename[0] + "cpu_mem_upload_success.txt")
        if len(filename) == 2:
            logfile = os.path.join(logdir, filename[0] + "-" + filename[1] + "cpu_mem_upload_success.txt")
        if not os.path.exists(logfile):
            data = pd.read_table(i, sep='\s+', header=0)
            Run_Time = round(float(data['s']) / 3600, 2)
            Memory = round(float(data['max_rss']) / 1024, 2)
            CPU_Load = float(data['mean_load'])
            if len(filename) == 1:
                sample_name = "project_level"
                CPU = yaml_file["cpu_mem"][filename[0]]["cpu"]
            if len(filename) == 2:
                sample_name = filename[1]
                CPU = yaml_file["cpu_mem"][filename[0]]["cpu"]

            desc = {"Analysis__Module": filename[0],
                    "Run_Time": Run_Time,
                    "CPU": CPU,
                    "CPU_Load": CPU_Load,
                    "Memory": Memory}
            oe_sor_upload(compact,
                          pro_type,
                          info_type,
                          desc,
                          sample_name=sample_name,
                          operator=user)
            os.mknod(logfile)
        else:
            print(f"Cpu_Mem information for {filename[0]} has been upload !Skipping!")


# =======================================================================================================================
parser = argparse.ArgumentParser(description="空间转录组项目SOR自动上传")
parser.add_argument("-i", "--config", default="config/config.yaml")
parser.add_argument("-m", "--metadata", default="config/samples.csv")
parser.add_argument("-l", "--logdir", default="logs/upload_sor", help="sor upload successed log dir")
parser.add_argument("-s", "--step", default="Project,Sample,QC,Ref_Database,Software,Cpu_Mem",
                    help="sor upload modules")
args = parser.parse_args()
config = os.path.abspath(args.config)

if __name__ == '__main__':
    yaml_file = yaml.load(open(config), Loader=yaml.FullLoader)
    compact = yaml_file["report"]["Project_Num"]
    pro_type = yaml_file["report"]["Project_Type"]
    user = getpass.getuser()
    if not os.path.exists(args.logdir):
        os.makedirs(args.logdir)
    if "Project" in args.step:
        project(compact, yaml_file, user, pro_type, 'Project', args.logdir)
    if "Sample" in args.step:
        sample(compact, args.metadata, user, pro_type, "Sample", args.logdir)
    if "QC" in args.step:
        QC(compact, args.metadata, user, pro_type, "QC", args.logdir)
    if "Ref_database" in args.step:
        database(compact, yaml_file, user, pro_type, 'Ref_Database', args.logdir)
    if "Software" in args.step:
        software(compact, yaml_file, user, pro_type, 'Softwares', args.logdir)
    if "Cpu_Mem" in args.step:
        benchmark(compact, yaml_file, user, pro_type, 'Cpu_Mem', args.logdir)
