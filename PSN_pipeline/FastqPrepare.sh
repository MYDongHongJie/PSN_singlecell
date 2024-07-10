#!/bin/bash

# 判断传入的参数数量
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <path> <sampleinfo>"
    exit 1
fi

# 获取所有传入的参数
params=("$@")

sampleinfo="${params[$#-1]}"
paths="${params[@]:0:$#-1}"
path=$(IFS=,; echo "${paths[*]}")

format_path=`echo ${path}|sed 's|,|\t|g'`
if [ -f mergefileR1.sh ] || [ -f mergefileR2.sh ]; then
	rm mergefileR1.sh  && rm mergefileR2.sh
fi

cat /PERSONALBIO/work/singlecell/s03/proj_doing/Micro/M4/youhua202404/merge_head.info > mergefileR1.sh 
cat /PERSONALBIO/work/singlecell/s03/proj_doing/Micro/M4/youhua202404/merge_head.info > mergefileR2.sh 
cat ${sampleinfo}|cut -f 1,2 |while read a b;do 
	mkdir -p 1.raw_data/${a}/
	n=$(find $format_path -name *.gz |grep -Ei "raw_1.fq|R1_001.fastq|R1.fastq"|grep ${b});
	j=$(find $format_path -name *.gz |grep -Ei "raw_2.fq|R2_001.fastq|R2.fastq"|grep ${b});
	echo "cat " $n " > 1.raw_data/${a}/${a}_S1_L001_R1_001.fastq.gz ">> mergefileR1.sh
	echo "cat " $j " > 1.raw_data/${a}/${a}_S1_L001_R2_001.fastq.gz ">> mergefileR2.sh
	# echo "cat " $n " > 1.raw_data/${a}/${a}_S1_L001_R1_001.fastq.gz ">> mergefileR1.sh
	# echo "cat " $j " > 1.raw_data/${a}/${a}_S1_L001_R2_001.fastq.gz ">> mergefileR2.sh
done

sh  mergefileR1.sh && sh  mergefileR2.sh

