#!/usr/bin/env python3
# encoding: utf-8
# Author liuling
# Modified luyao
# Desc   : This script is used to push SOR recird
# coding: utf-8

import os,time, glob ,shutil,yaml, argparse,re
# import requests
import getpass
# import json
# import numpy as np
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
    # logfile=os.path.join(logdir,"project_upload_success.txt")
    # if not os.path.exists(logfile):
    desc = yaml_file["report"].copy()
    desc.pop("Task_Num")
    desc.pop("Reference")
    desc.pop("Library")
    desc["Rawdata_Path"] = ','.join(map(str, yaml_file["mobivision_params"]["raw_data_obs_address"]))
    #if yaml_file["mobivision_params"]["intron"]=="included":
    #    desc["Others"] = "single-nucleus"
    oe_sor_upload(compact,
                  pro_type,
                  info_type,
                  desc,
                  sample_name="",
                  operator=user)
    # os.mknod(logfile)
    # else:
    #     print("Project information has been upload !Skipping!")
#=======================================================================================================================
def sample(compact,metadata,user,pro_type,info_type,logdir):
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
        # logfile=os.path.join(logdir,sample+"sample_upload_success.txt")
        # if not os.path.exists(logfile):
        desc=dict()
        desc["Sample_Type"]=samples.loc[(sample), ["species"]][0]
        desc["Group1"]=samples.loc[(sample), ["group"]][0]
        desc["Batch_ID"]=samples.loc[(sample), ["batchid"]][0]
        oe_sor_upload(compact,
                    pro_type,
                    info_type,
                    desc,
                    sample_name=sample,
                    operator=user)
        # os.mknod(logfile)
        # else:
        #     print(f"Sample information for {sample} has been upload !Skipping!")
#=======================================================================================================================
#def open_html(file):
#    from selenium import webdriver
#    options = webdriver.ChromeOptions()
#    options.add_argument("headless")
#    options.add_argument('no-sandbox')
#    options.add_argument('dns-prefetch-disable')
#    options.add_argument('disable-gpu')
#    options.add_argument('disable-dev-shm-usage')
#    options.add_argument('disable-features=VizDisplayCompositor')
#    options.add_argument('disable-features=NetworkService')
#    options.add_argument('window-size=1920x945')
#    options.add_experimental_option("excludeSwitches",["ignore-certificate-errors"])
#    browser = webdriver.Chrome('/home/luyao/chromedriver', chrome_options=options) 
#    browser.get(f"file://{file}")
    #this step need to be fixed###########################################################
#    chemistry = browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(3) > td.tableMetric-value-chopped")
#    reference_path =  browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(5) > td.tableMetric-value")
    # include-introns
#    if browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(7) > td.tableMetric-value-chopped").text == "MobiVision v1.1":
#        include_introns = browser.find_element_by_css_selector("#sample_table > table > tbody > tr:nth-child(4) > td.tableMetric-value").text
#    else: 
#        include_introns = "unknown in v3"
#    return chemistry.text ,reference_path.text , include_introns
###############################################################################################
def QC(compact,yaml_file,user,pro_type,info_type,logdir):
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
    # if type(yaml_file["mobivision_params"]["samples"]) == list:
        # if all(isinstance(i, str) for i in yaml_file["mobivision_params"]["samples"]):
            # sample_list = yaml_file["mobivision_params"]["samples"]
        # elif all(isinstance(i, dict) for i in yaml_file["mobivision_params"]["samples"]):
            # sample_list = [k for i in yaml_file["mobivision_params"]["samples"] for k,v in i.items()]
        # samples = pd.Series(sample_list,index = sample_list)
    # elif type(yaml_file["mobivision_params"]["samples"]) == str and yaml_file["mobivision_params"]["samples"] == "all" :
            # samples = pd.Series()
            # for addr in yaml_file["mobivision_params"]["raw_data_obs_address"]:
                # folders = os.popen("obsutil ls -d -s %s | grep obs  | grep -v '\(md5\|bam_files\)' | sed 's#%s##'"%(addr,addr)).readlines()
                # for i in folders:
                    # if i.strip() != "":
                        # samples[i.strip().replace("/","")] = i.strip().replace("/","")
    # else :  sys.exit("mobivision_params - samples can only be string 'all' or a list of sample names.")
    samples = pd.read_csv("config/samples.csv",sep=",",dtype=str)
    samples.index = samples["sampleid"]
    # report - Count_QC this is weird because we haven't reach this step yet
    if os.path.exists("result/Count_QC/cell_statitics_before_after_QC.xls"):
        Count_QC = pd.read_csv("result/Count_QC/cell_statitics_before_after_QC.xls", dtype=str,sep="\t",index_col=0)
        if "Total_cells_after_rmdoublets" in Count_QC.columns:
            Count_QC = Count_QC[["Total_cells_after_rmdoublets","percent.mito_higher"]]
        else:
            Count_QC = Count_QC[["Total_cells_afterQC","percent.mito_higher"]]

    ## main
    for (sample) in samples.index:
        # logfile=os.path.join(logdir,sample+"_QC_upload_success.txt")
        # if not os.path.exists(logfile):
        # os.system(f"cp result/cellranger/{sample}/outs/web_summary.html /public/scRNA/works/Documents/backup/%s_{sample}_web_summary.html"  %(compact))
        QC=pd.read_csv(f"result/mobivision/{sample}/outs/{sample}_summary.csv", dtype=str)
        # rename and reorder 
        desc=pd.DataFrame(QC.values.T, index= QC.columns.str.replace(" ", "_"), columns=QC.index).to_dict()[0]
        if 'Q30_Bases_in_RNA_Read_2' in desc.keys():
            desc["Q30_Bases_in_RNA_Read"] = desc["Q30_Bases_in_RNA_Read"]+"; "+ desc["Q30_Bases_in_RNA_Read_2"]
            del desc["Q30_Bases_in_RNA_Read_2"]
        desc["Include_introns"] = "Intron excluded"
        if  len(re.split(r'[()]',desc["version"])) > 1:
            desc["Include_introns"] = re.split(r'[()]',desc["version"])[1].strip()
        desc.pop("Q30_Bases_in_Barcode+UMI")
        desc.pop("version")
        desc.pop("Fraction_of_Unique_Reads_in_Cells")
        desc.pop("id")
        desc.pop("Number_of_Valid_Reads")

        desc["Reads_Mapped_to_Transcriptome(%)"] = desc.pop("Reads_Mapped_Confidently_to_Transcriptome")
        desc["Fraction_Reads_in_Cells(%)"] = desc.pop("Fraction_of_Valid_Reads")
        desc["Chemistry"] = desc.pop("kit")
        desc["Reference_path"] =  desc.pop("reference")
        desc["Number_of_Reads"] = desc.pop("Number_of_Raw_Reads")
        desc["Q30_Bases_in_RNA_Read(%)"] = desc.pop("Q30_Bases_in_RNA_Read")
        desc["Reads_Mapped_to_Genome(%)"] = desc.pop("Reads_Mapped_to_Genome")
        desc["Clean_Base"] = str(format(int(desc["Number_of_Reads"].replace(",",""))*300 / 1000000000, '.2f')) + "G"
        desc=OrderedDict(desc)
        desc.update({'QC__Type':"mobivision"})
        desc.move_to_end('Number_of_Reads', last=False)
        desc.move_to_end('Reads_Mapped_to_Transcriptome(%)', last=False)
        desc.move_to_end('Fraction_Reads_in_Cells(%)', last=False)
        desc.move_to_end('Q30_Bases_in_RNA_Read(%)', last=False)
        desc.move_to_end('Median_Genes_per_Cell', last=False)
        desc.move_to_end('Estimated_Number_of_Cells', last=False)
        desc.move_to_end('QC__Type', last=False)
        desc=dict(desc)
        
        
        #need to change information extraction method
        #desc["Chemistry"], desc["Reference_path"], desc["intron"] = open_html(os.path.abspath(f"result/cellranger/{sample}/outs/web_summary.html"))
        # report-related
        if os.path.exists("result/Count_QC/cell_statitics_before_after_QC.xls"):
            desc["Total_cells_afterQC"] = Count_QC.loc[sample].iloc[0]
            desc["mito_cutoff"] = Count_QC.loc[sample].iloc[1]
        # upload
        oe_sor_upload(compact,
                      pro_type,
                      info_type,
                      desc,
                      sample_name=sample,
                      operator=user)
        # os.mknod(logfile)
        # else:
        #     print(f"QC  information for {sample} has been upload !Skipping!")
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
    logfile=os.path.join(logdir,"Genome_Database_upload_success.txt")
    if not os.path.exists(logfile):
        referenceDB_file = "%s/config.yaml" %(yaml_file["report"]["Reference"])
        referenceDB = yaml.load(open(referenceDB_file),Loader=yaml.FullLoader)
        desc=dict()
        desc["Database"]="genome"
        desc["Linkage"]=referenceDB['genome_linkage']
        desc["Version"]=referenceDB['genome_version']
#        desc["Path"]=referenceDB_file
        os.mknod(logfile)
        oe_sor_upload(compact,
                              pro_type,
                              info_type,
                              desc,
                              sample_name="",
                              operator=user)
        print(f"upload successfully: Database information for reference genome")
    else:
        print(f"Database information for reference genome has been upload !Skipping!")

    ## report-related DB
    if yaml_file["report_params"]["module"]["report"] == True:
        reportDB_file = "scripts/report/pic/tables/scRNA_database.txt"
        reportDB = pd.read_table(reportDB_file, dtype=str).set_index("使用数据库", drop=False)
        reportDB.to_dict()
        for i in reportDB.index:
            logfile=os.path.join(logdir,i+"_Database_upload_success.txt")
            if not os.path.exists(logfile):
                desc=dict()
                desc["Database"]=i
                desc["Linkage"]=reportDB.loc[i,"网页链接"]
                desc["Version"]=reportDB.loc[i,"版本"]
                oe_sor_upload(compact,
                              pro_type,
                              info_type,
                              desc,
                              sample_name="",
                              operator=user)
                os.mknod(logfile)
            else:
                print(f"Database information for {i} has been upload !Skipping!")
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
    logfile=os.path.join(logdir,"mobivision_software_upload_success.txt")
    if not os.path.exists(logfile):
        desc=dict()
        desc["Version"] = yaml_file['mobivision_params']['version']
        oe_sor_upload(compact,
                      pro_type,
                      info_type,
                      desc,
                      sample_name="",
                      operator=user)
        os.mknod(logfile)
    else:
        print(f"Software information for MobiVision has been upload !Skipping!")
    if yaml_file["report_params"]["module"]["report"] == True:
        software_table = pd.read_table("scripts/report/pic/tables/scRNA_software.txt", dtype=str).set_index("Software", drop=False)
        software_version = software_table.to_dict()['Version']
        software_version.pop('MobiVision')
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
    threads=1
    #for i in rule_files:
    #    with open (i,'r') as f:
    #        for line in f :
    #            line=line.strip()
    #            if line.startswith("rule "):
     #               ll= re.split('\s+|:',line)
    #                rule_name=ll[1]
    #                threads=yaml_file['cpu']['mobivision']
    #                rules[rule_name]=threads
    rules=yaml_file["cpu"].copy()
    ## benchmarks
    list=glob.glob("benchmarks/*txt")
    cpu_hours_all=0
    for i in list:
        filename=i.split(sep="/")[-1].replace(".benchmark.txt","").replace("benchmarks/","").split(sep=".")
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
            CPU=rules[Analysis_Module]
            cpu_hours = CPU*Run_Time
            cpu_hours_all=cpu_hours_all+cpu_hours
            desc={"Analysis__Module": Analysis_Module,
                "Run_Time": Run_Time ,
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
            print(f"Cpu_Mem information for {'-'.join(filename)} has been upload !Skipping!")
            
    #upload all cpu_hours
    logfile=os.path.join(logdir,"all_cpu_hours_upload_success.txt")
    desc={"Analysis__Module": 'all',
          "Run_Time": cpu_hours_all }
    oe_sor_upload(compact,
                  pro_type,
                  info_type,
                  desc,
                  sample_name="project_level",
                  operator=user)
    os.mknod(logfile)
#=======================================================================================================================
parser = argparse.ArgumentParser(description="单细胞转录组项目SOR自动上传")
parser.add_argument("-i", "--config",default="config/config.yaml")
parser.add_argument("-m", "--metadata",default="config/samples.csv")
parser.add_argument("-l", "--logdir",default="logs/upload_sor", help="sor upload successed log dir")
parser.add_argument("-s", "--step",default="Project,Sample,QC,Ref_Database,Software,Cpu_Mem", help="sor upload modules")
args = parser.parse_args()
config = os.path.abspath(args.config)

if __name__ == '__main__':
    yaml_file = yaml.load(open(config),Loader=yaml.FullLoader)
    compact = yaml_file["report"]["Task_Num"].split('-')[0]
    pro_type = "单细胞转录组"
    user = getpass.getuser()
    if not os.path.exists(args.logdir):
        os.makedirs(args.logdir)
    if "Project" in args.step:
        project(compact,yaml_file ,user,pro_type,'Project',args.logdir)
    if "Sample" in args.step:
        sample(compact,args.metadata,user,pro_type,"Sample",args.logdir)
    if "QC" in args.step:
        QC(compact,yaml_file,user,pro_type,"QC",args.logdir)
    if "Ref_database" in args.step:
        database(compact,yaml_file ,user,pro_type,'Ref_Database',args.logdir)
    if "Software" in args.step:
        software(compact,yaml_file ,user,pro_type,'Softwares',args.logdir)
    if "Cpu_Mem" in args.step:
        benchmark(compact,yaml_file ,user,pro_type,'Cpu_Mem',args.logdir)
