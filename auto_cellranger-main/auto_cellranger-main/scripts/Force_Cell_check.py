import os
#import re
import sys
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="本脚本用于在出具报告时确认是否需要进行force cell")
parser.add_argument("-i", "--samplescsv",help="sample.csv")
args = parser.parse_args()

samples = args.samplescsv

#测试文件
#samples = "/public/scRNA_works/works/liuxuan/Test/snakemake/DOE20224404_count/config/samples.csv"
#output = "/public/scRNA_works/works/liuxuan/Test/snakemake/DOE20224404_count/Force_cell.check"

samples = pd.read_csv(samples, dtype=str).set_index("sampleid", drop=False)

data = pd.DataFrame(columns=["sampleid","Median Genes per Cell","Median UMI Counts per Cell","Fraction Reads in Cells"])
for index in range(len(samples.index)):
   sample = samples.index[index]
   sample_stat = "./result/cellranger/%s/outs/metrics_summary.csv"%(sample)
   metrix_sum = pd.read_csv(sample_stat, dtype=str)
   metrix_sum['sampleid'] = sample
   df = metrix_sum.filter(items=["sampleid","Median Genes per Cell","Median UMI Counts per Cell","Fraction Reads in Cells"])
   data = data.append(df)


data.index = data.sampleid

for id in range(0,data.shape[0]): 
    sampleid=data.index[id]
    Median_Genes_per_Cell = int("".join(data['Median Genes per Cell'][id].split(",")))
    Median_UMI_Counts_per_Cell = int("".join(data['Median UMI Counts per Cell'][id].split(",")))
    Fraction_Reads_in_Cells = float(data['Fraction Reads in Cells'][id].strip('%'))
    #print(sampleid)
    if Median_Genes_per_Cell < 700 or Median_UMI_Counts_per_Cell < 1000:
        if Median_Genes_per_Cell >= 700 and Median_UMI_Counts_per_Cell < 1000:
            info = "Error: " + str(sampleid) +" 样本的Median_UMI_Counts_per_Cell 为 " +str(Median_UMI_Counts_per_Cell) + ", 小于 1000，请停止出具报告, 与技术支持确认cellranger结果是否需要调整。"
            print(f"\033[0;31m{info}\033[0m")
            sys.stderr.write(f"\033[0;31m{info}\033[0m\n")
        elif Median_Genes_per_Cell < 700 and Median_UMI_Counts_per_Cell >= 1000:
            info = "Error: " + str(sampleid) +" 样本Median_Genes_per_Cell 为 " +str(Median_Genes_per_Cell) + ", 小于 700，请停止出具报告, 与技术支持确认cellranger结果是否需要调整。"
            print(f"\033[0;31m{info}\033[0m")
            sys.stderr.write(f"\033[0;31m{info}\033[0m\n")
        elif Median_Genes_per_Cell < 700 and Median_UMI_Counts_per_Cell < 1000:
            info = "Error: " + str(sampleid) +" 样本Median_Genes_per_Cell 为 " +str(Median_Genes_per_Cell) + ", 小于 700; 且 Median_UMI_Counts_per_Cell 为 " +str(Median_UMI_Counts_per_Cell) + ", 小于 1000，请停止出具报告, 与技术支持确认cellranger结果是否需要调整。"
            print(f"\033[0;31m{info}\033[0m")
            sys.stderr.write(f"\033[0;31m{info}\033[0m\n")
    elif Fraction_Reads_in_Cells < 40:
        info = "Warning: " + str(sampleid) + " 样本存在Fraction_Reads_in_Cells 小于 40。"
        print(f"\033[0;33m{info}\033[0m") 
