#!/usr/bin/env python3
# Author lip
# Date: 2023.05.17

# Import Required Libraries
import argparse
import glob
import subprocess
import os
#import configparser
import sys
import pandas as pd
import re
import shutil

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates script for BD pipline v2.0")
parser.add_argument("-i", "--input_path", help="Input directory with fastq.gz",
                    default="./raw_data")
parser.add_argument("-s", "--sampleid", help="sampleid for analysis.",
                    default="OE1")
parser.add_argument("-d", "--genome", help="Input the file of BD genome which build by STAR.(eg 'Human/genome.tar.gz')",
                    default="/data/database/BD-scRNA-refdata/genome_BD_v2/Human/genome.tar.gz")
parser.add_argument("-t", "--Basic_Algo_Only", help="whether if Sequencing file mixed samples (default:false) ",
                    choices=["false", "true"], default="false")
parser.add_argument("-v", "--vdj", help="Specify if VDJ run: human, mouse, humanBCR, humanTCR, mouseBCR, mouseTCR",
                    default= "mRNA")
parser.add_argument("-f", "--force_cell", help="force cells to expect number (eg, 10000 or 8000) ",
                    default="false")
parser.add_argument("--AbSeq", help="AbSeq reference file in FASTA format ",
                    default=None)
parser.add_argument("--get_bam", help="whether to skip creation of the Bam Output ",
                    default="false")
parser.add_argument("-n", "--nthread", help="thread (default:8) ",
                    default="16")
args = parser.parse_args()

#python BD_get_qsub.py  -i ./raw_data/ -s G03841_CD_huichang -d /data/database/BD-scRNA-refdata/Human/GRCh38.tar.gz  -g /data/database/BD-scRNA-refdata/Human/gencode.v32.annotation.gtf -t false -v human -n 6

#1.get parame

inputDirectory = os.path.abspath(args.input_path)
sampleid = str(args.sampleid)
genome_fa = os.path.abspath(args.genome)
BAO_type = str(args.Basic_Algo_Only)
vdj_type = str(args.vdj)
force_cell = str(args.force_cell)
nthread = int(args.nthread)
get_bam = str(args.get_bam)
AbSeq_fa = str(args.AbSeq)

if nthread < 4:
    nthread = 8
if vdj_type !="mRNA":
    nthread = 16

cwl_name = "../../envs/BD_workflow/2.0.cwl"
abs_path = os.path.abspath(cwl_name)

#2.Create script file and sampleid_yml file
pwd = os.getcwd()
scriptName = ("%s.BD_V2.0.sh" % sampleid)
if os.path.isfile(scriptName):
    os.remove(scriptName)
##get script
script = open(scriptName, 'a')
script.write("module purge && source /home/lipeng/miniconda3/bin/activate CWLtool && cwl-runner --parallel --rm-tmpdir" + " ")
script.write("--tmp-outdir-prefix " + sampleid + "temp/tempfile" + " ")
script.write("--tmpdir-prefix " + sampleid + "tmp/tmpfile" + " ")
script.write("--outdir " + sampleid + " ")
script.write( abs_path + " ")
script.write(sampleid + "_wta_2.0_h.yml" + "\n")
script.close()


##get yml file
ymlName = ("%s_wta_2.0_h.yml" % sampleid)
if os.path.isfile(ymlName):
    os.remove(ymlName)

yamlfile = open(ymlName, 'a')
yamlfile.write("cwl:tool: rhapsody" + "\n")
yamlfile.write("\n")
yamlfile.write("Reads:"+ "\n")

if vdj_type =="mRNA":
    term_gz = glob.glob(r'%s/%s/*/*.gz'% (inputDirectory,sampleid))
else:
    term_gz = glob.glob(r'%s/%s/*_fastqs/*.gz'% (inputDirectory,sampleid)) + glob.glob(r'%s/%s_TCR/*_fastqs/*.gz'% (inputDirectory,sampleid)) + glob.glob(r'%s/%s_BCR/*_fastqs/*.gz'% (inputDirectory,sampleid)) +  glob.glob(r'%s/%s_dan/*_fastqs/*.gz'% (inputDirectory,sampleid))
term_gz.sort()
for gzfile in term_gz:
    print(gzfile)
    yamlfile.write("\n")
    yamlfile.write(" - class: File" + "\n")
    yamlfile.write(("   location: \"%s\"" % gzfile) + "\n")
yamlfile.write("\n")
yamlfile.write("Reference_Archive:" + "\n")
yamlfile.write("   class: File" + "\n")
yamlfile.write(("   location: \"%s\"" % genome_fa) + "\n")
#yamlfile.write("Transcriptome_Annotation:" + "\n")
#yamlfile.write("   class: File" + "\n")
#yamlfile.write(("   location: \"%s\"" % genome_gtf) + "\n")
yamlfile.write(("Maximum_Threads: %s" % nthread) + "\n")
yamlfile.write("\n")

if AbSeq_fa !="None":
    AbSeq_fa_param = AbSeq_fa
    yamlfile.write("\n")
    yamlfile.write("AbSeq_Reference:" + "\n")
    yamlfile.write(" - class: File" + "\n")
    yamlfile.write(("   location: \"%s\"" % AbSeq_fa) + "\n")

if BAO_type =="false":
    yamlfile.write("#Basic_Algo_Only: true" + "\n")
else:
    yamlfile.write("Basic_Algo_Only: %s" % BAO_type + "\n")
if vdj_type =="mRNA" or vdj_type == "dan":
    yamlfile.write("#VDJ_Version: human" + "\n")
else:
    yamlfile.write("VDJ_Version: %s" % vdj_type + "\n")
if force_cell =="false":
    yamlfile.write("#Exact_Cell_Count: true" + "\n")
else:
    yamlfile.write("Exact_Cell_Count: %s" % force_cell + "\n")
if get_bam =="false":
    yamlfile.write("Generate_Bam: false" + "\n")
else:
    yamlfile.write("Generate_Bam: true" + "\n")

yamlfile.close()


#3.print scripts
#scripts = glob.glob(os.path.join(pwd + "/" + scriptName))
#scripts.sort()
#print(scripts)
#for i in range(0, len(scripts), 1):
#    BD_run_task = 'cd %s && qsub -V -cwd -pe smp %s -q big@hcu-0007,big@hcu-0008,big@hcu-0009,big@hcu-0011,big@hcu-0016  %s' % (pwd, nthread, scripts[i])
#    print("step:"+ BD_run_task)
#    subprocess.run(BD_run_task, shell=True, check=True)
#
BD_run_task = 'cd %s && sh %s.BD_V2.0.sh' % (pwd, sampleid)
print("step:"+ BD_run_task)
subprocess.run(BD_run_task, shell=True, check=True)
