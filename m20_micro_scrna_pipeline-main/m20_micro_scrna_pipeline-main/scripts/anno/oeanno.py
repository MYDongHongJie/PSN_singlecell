#!/usr/bin/env python3
# Version 1.01

import argparse
import glob
import os
import subprocess
import shutil
import sys
import string
import pandas as pd
from pathlib import Path

# Path(__file__).parent.absolute()

parser = argparse.ArgumentParser(description='Generates annotation scripts ')

parser.add_argument('-n', '--name', default="Unigene", help='prefix name, default=Unigene')
parser.add_argument('-o', '--outdir', type=str, default='anno_result', help='Output directory for annotation results,  default=anno_result')
parser.add_argument('-m', '--mode', type=str,help='plant animal or micro')
parser.add_argument('-i', '--input', help='input file')

args = parser.parse_args()

###### Write the script(s)
mode = args.mode
name = args.name
input = args.input
outdir= args.outdir
scriptpath = sys.path[0]

df = pd.read_table(input)
df = df[['#Query_id','Subject_id','E_value','Identity']]
anno = pd.read_table('/data/database/denovo_database_2020/kegg/KEGG_anno.txt')
out = pd.merge(df, anno, left_on='Subject_id', right_on='id',how='left')
out = out.drop('id', 1)
out = out[out['KO'].notnull()]
outfile = outdir + name + ".KEGG.gene.anno.xls_tmp"
out.to_csv(outfile, index=False,sep='\t')

if mode == "plant":
	cmd =' awk -F"\\t" -v OFS="\\t" \'{split($8,f,"[,]");split($9,b,"[|]");for(i=1; i<=length(b); i++) print $1,$2,$3,$4,$5,$6,$7,f[i]"\\t"b[i]}\' %s |awk -F"\\t" -v OFS="\\t" \'NR==FNR{a[$1]=$1;b[$1]=$5;next}{if($8==a[$8]) print $0,b[$8]}\' %s - |awk -F"\\t" -v OFS="\\t" \'$10~"plant"\' |cut -f1-9 |awk -F"\\t" -v OFS="\\t" \'{a[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]=a[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]"|"$8;b[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]=b[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]"|"$9}END{for(i in a) print i"\\t"a[i]"\\t"b[i]}\'|sed "s:\\t|:\\t:g"|sed "1i #Gene_id\\tDatabase_GeneID\\te_value\\tidentity\\tKO\\tgene_name\\tDescription\\tPathway\\tPathway_definition" > %s ' % (outdir + name + ".KEGG.gene.anno.xls_tmp", str(Path(__file__).parent.absolute()) + "/KEGGpathway_three_levels_v2_191022.xls",outdir + name + ".KEGG.gene.anno.xls")
else:
	cmd =' awk -F"\\t" -v OFS="\\t" \'{split($8,f,"[,]");split($9,b,"[|]");for(i=1; i<=length(b); i++) print $1,$2,$3,$4,$5,$6,$7,f[i]"\\t"b[i]}\' %s |awk -F"\\t" -v OFS="\\t" \'NR==FNR{a[$1]=$1;b[$1]=$5;next}{if($8==a[$8]) print $0,b[$8]}\' %s - |awk -F"\\t" -v OFS="\\t" \'$10~"animal"\' |cut -f1-9 |awk -F"\\t" -v OFS="\\t" \'{a[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]=a[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]"|"$8;b[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]=b[$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7]"|"$9}END{for(i in a) print i"\\t"a[i]"\\t"b[i]}\'|sed "s:\\t|:\\t:g"|sed "1i #Gene_id\\tDatabase_GeneID\\te_value\\tidentity\\tKO\\tgene_name\\tDescription\\tPathway\\tPathway_definition" > %s ' % (outdir + name + ".KEGG.gene.anno.xls_tmp", str(Path(__file__).parent.absolute()) + "/KEGGpathway_three_levels_v2_191022.xls",outdir + name + ".KEGG.gene.anno.xls")
cmd += '&&awk -F"\\t" -v OFS="\\t"  \'{if($8!="--") print $1,$5,$8,$9}\'  %s|awk -F"\\t" -v OFS="\\t" \'NR>1{split($3,f,"[|]");split($4,b,"[|]");for(i=1; i<=length(b); i++) print b[i],f[i],$1,$2}\'  |awk -F"\\t" -v OFS="\\t" \'{a[$1"\\t"$2]=a[$1"\\t"$2]";"$3;b[$1"\\t"$2]=b[$1"\\t"$2]";"$4}END{for(i in a) print i"\\t"a[i]"\\t"b[i]}\' |awk -F"\\t" -v OFS="\\t" \'{print $1,$2,$3,$4,$4}\'| awk -F"\\t" -v OFS="\\t" \'{a=gsub(";","\\t",$5); print $0}\' |awk -F"\\t" -v OFS="\\t" \'{print $1,$2"\\t"(NF-5)"\\t"$3"\\t"$4}\'|awk -F"\\t" -v OFS="\\t" \'{sub(";","",$4);sub(";","",$5); print $0}\' |awk -F"\\t" -v OFS="\\t" \'NR==FNR{a[$1]=$1;b[$1]=$2"--"$3;next}{if($2==a[$2]) print b[$2],$0}\' %s - |sort|sed "1i #class\\tPathway_definition\\tPathway\\tGene_number\\tGene_id\\tKOs"  >%s' % (outdir + name + ".KEGG.gene.anno.xls", str(Path(__file__).parent.absolute()) + "/KEGGpathway_three_levels_v2_191022.xls",outdir + name + ".KEGG.pathway.xls")
cmd += '&&totalnum=`cat %s|awk -F"\\t" -v OFS="\\t" \'$8!="--"&&$8!="Pathway"\' |wc -l`&&sed "s/--/\\t/" %s |awk -F"\\t" -v OFS="\\t" \'NR>1{print $2,$1,$6}\' |awk -F"\\t" -v OFS="\\t" \'{split($3,b,"[;]");for(i=1; i<=length(b); i++) print $1,$2,b[i]}\' |awk -F"\\t" -v OFS="\\t" \'++a[$1,$2,$3]==1\' |awk -F"\\t" -v OFS="\\t" \'{a[$1"\\t"$2]=a[$1"\\t"$2]";"$3}END{for(i in a) print i"\\t"a[i]}\' |awk -F"\\t" -v OFS="\\t"  \'{print $1,$2,$3,$3}\'| awk -F"\\t" -v OFS="\\t"  \'{a=gsub(";","\\t",$4); print $0}\' |awk -F"\\t" -v OFS="\\t"  \'{print $1,$2"\\t"(NF-4)"\\t"$3}\'|awk -F"\\t" -v OFS="\\t" \'{sub(";","",$4); print $0}\'|sort -t$\'\\t\' -k2 -r |awk -F"\\t" -v OFS="\\t" \'{print $1,$2,$3,$3/"\'$totalnum\'",$4}\' |sed "1i Classification_level2\\tClassification_level1\\tgene_number\\tpercentage\\tGenes">%s' % (outdir + name + ".KEGG.gene.anno.xls",outdir + name + ".KEGG.pathway.xls",outdir + name + ".KEGG.classification.xls")
print(outdir + name + ".KEGG.classification.xls")
cmd += '&&perl %s  -tab %s -od %s -key %s' % (scriptpath + "/kegg_classification.pl",outdir + name + ".KEGG.classification.xls", outdir,  name)
subprocess.call(cmd,shell=1)
