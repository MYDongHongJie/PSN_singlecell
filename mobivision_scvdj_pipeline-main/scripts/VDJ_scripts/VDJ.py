# -*- coding: utf-8 -*-
# test for vdj
# Author: donghongjie
# Date: 07/21/2023
import argparse
import configparser
import subprocess
import os,sys
import re
import csv
import shutil
#import pandas as pd
#import numpy as np


args  = argparse.ArgumentParser(description='VDJ analysis')
args.add_argument('-v','--input', type=str, default = "result/mobivision/",help="Input dir")
args.add_argument('-o','--output', type=str, default = 'result/VDJ_aggr/')
args.add_argument('-q','--groupby', type=str, default ='Sampleid')
args.add_argument('-t','--VDJ', type=str, default ='TCR,BCR',help="VDJ type,")
args.add_argument('-m','--metadata', type=str, default ='config/samples.csv',help="config/samples.csv")
args.add_argument('-p','--thread', type=int, default  = '4')
args  = args.parse_args()

input = os.path.abspath(args.input)
group = args.groupby
metadata = args.metadata
thread = args.thread
output = os.path.abspath(args.output)
VDJ = re.split(',',args.VDJ)

with open(metadata, "r", ) as file:
    reader = csv.DictReader(file)
    sample_name = [row['sampleid'] for row in reader]

#1.get TCR/BCR info from raw matrix
#paths_py = os.path.dirname(os.path.realpath(__file__))
vdjdata_r = "scripts/VDJ_scripts/vdj_merge.R"
immuarch_plot = "scripts/VDJ_scripts/testvdjtk.R"
mlo_immuarch = "module purge && module load immunoinformatics/1.0.0" 
for subtype in VDJ:
    outputdir  = "{}/{}".format(output,subtype)
    cmd1 = 'module purge && module load OESingleCell/3.0.d && Rscript %s -i %s -o %s  -m %s --type %s \n' % (vdjdata_r, input, outputdir, metadata ,subtype )
    print("step1 get TCR/BCR info from raw matrix :"+cmd1)
    subprocess.run(cmd1, shell=True, check=True)

    print("Start to analysis ",subtype," in ",outputdir)
    #CalcSpectratype
    cmd2 = '%s && vdjtools CalcSpectratype -m %s/Clonotypes/vdj_metadata.xls %s/Spectratyping/sample  \n' % ( mlo_immuarch,outputdir,outputdir)
    print("step2:"+cmd2)
    subprocess.run(cmd2, shell=True, check=True)
    #PlotFancySpectratype
    for sample in sample_name:
        sampleVDJ=sample+"_"+subtype
        cmd3 = '%s && vdjtools PlotFancySpectratype -t 10 %s/Clonotypes/%s.xls %s/Spectratyping/%s' % (mlo_immuarch,outputdir,sampleVDJ,outputdir,sampleVDJ)
        print("step3:"+cmd3)
        subprocess.run(cmd3, shell=True, check=True)
        cmd4 = '%s && vdjtools PlotSpectratypeV  -t 12 %s/Clonotypes/%s.xls %s/Spectratyping/%s && \
                rename  txt xls %s/Spectratyping/*.txt' % (mlo_immuarch,outputdir,sampleVDJ,outputdir,sampleVDJ,outputdir)
        print("step4:"+cmd4)
        subprocess.run(cmd4, shell=True, check=True)
    #Gene_usage
    if len(sample_name) == 1:
        print("样本数小于2，不提供热图结果")
    else:
        cmd5 = '%s && vdjtools CalcSegmentUsage -m %s/Clonotypes/vdj_metadata.xls -p  %s/Gene_usage/gene' % (mlo_immuarch,outputdir,outputdir)
        print("step5:"+cmd5)
        subprocess.run(cmd5, shell=True, check=True)
    for sample in sample_name:
        sampleVDJ=sample+"_"+subtype
        cmd6 = '%s && vdjtools PlotFancyVJUsage %s/Clonotypes/%s.xls %s/Gene_usage/%s && \
                rename  txt xls %s/Gene_usage/*.txt' % (mlo_immuarch,outputdir,sampleVDJ,outputdir,sampleVDJ,outputdir)
        print("step6:"+cmd6)
        subprocess.run(cmd6, shell=True, check=True)
        cmd7 = '%s && vdjtools PlotQuantileStats %s/Clonotypes/%s.xls %s/Clonity/%s && \
                rename  txt xls %s/Clonity/*.txt' % (mlo_immuarch,outputdir,sampleVDJ,outputdir,sampleVDJ,outputdir)
        print("step7:"+cmd7)
        subprocess.run(cmd7, shell=True, check=True)
   #DiversityStats
    if len(sample_name) == 1:   
        print("样本数小于2，不提供Diversity结果")
    else:
        cmd8 = '%s && vdjtools CalcDiversityStats -m  %s/Clonotypes/vdj_metadata.xls %s/Diversity/sample' % (mlo_immuarch,outputdir,outputdir)
        cmd9 = '%s && vdjtools RarefactionPlot -m  %s/Clonotypes/vdj_metadata.xls -f group  -l sample.id %s/Diversity/sample' % (mlo_immuarch,outputdir,outputdir)
        print("step8:"+cmd8)
        subprocess.run(cmd8, shell=True, check=True)
        print("step9:"+cmd9)
        subprocess.run(cmd9, shell=True, check=True)
        #VENN plot
        if len(sample_name) >5:
            print("样本数大于5, 默认不进行多样性 venn 分析。")
        else:
            cmd10 = '%s && vdjtools JoinSamples  -m %s/Clonotypes/vdj_metadata.xls -x 1 -p %s/Diversity/sample && \
                    rename  txt xls %s/Diversity/*.txt' % (mlo_immuarch,outputdir,outputdir,outputdir)
            print("step10:"+cmd10)
            subprocess.run(cmd10, shell=True, check=True)
    #copy files
    #pic_dir = "/public/scRNA_works/works/liuhongyan/Script/scVDJ/scvdj_report_script/pictures/"
    #cmd11 = 'cp %s/cd3_indel.png                %s/Spectratyping && \
    #         cp %s/vj_gene_usage_frequency.png  %s/Gene_usage  && \
    #         cp %s/vj_gene_usage_frequency_heatmap.png   %s/Gene_usage   && \
    #         cp %s/clonotypes_abundances.png    %s/Clonotypes    \n' % (pic_dir,outputdir,pic_dir,outputdir,pic_dir,outputdir,pic_dir,outputdir)
    #print("step11:"+cmd11)
    #subprocess.run(cmd11, shell=True, check=True)
    #immunarch plot
    cmd12 = '%s && Rscript %s -i %s -o %s -t %s -m %s ' % (mlo_immuarch,immuarch_plot,input,outputdir,subtype,metadata)
    print("step12:"+cmd12)
    subprocess.run(cmd12, shell=True, check=True)
    cmd13 = 'for pic in %s/*/*.pdf;do n=$(echo $pic|sed \'s/.pdf//\');/data/software/ImageMagick/ImageMagick-v7.0.8-14/bin/convert  -quality 500 -antialias -density 300 -transparent white -trim -flatten -sharpen 0x1.0 $pic ${n}.png;done' % (outputdir)
    subprocess.run(cmd13, shell=True, check=True)
