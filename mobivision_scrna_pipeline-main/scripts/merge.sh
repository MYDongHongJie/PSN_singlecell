#!/bin/bash
# list probable extend names
module purge & module load pigz
ls -d raw_data/*/ | sed 's/raw_data\///g' | sed 's/.$//g' > sample_name.txt
for i in $(cat sample_name.txt);
do
    mkdir ${PWD}/raw_data/${i}/merged
    zcat ${PWD}/raw_data/${i}/*R1_001.fastq.gz  | pigz -p 10  > ${PWD}/raw_data/${i}/merged/sample${i}.R1.fastq.gz
    zcat ${PWD}/raw_data/${i}/*R2_001.fastq.gz | pigz -p 10  > ${PWD}/raw_data/${i}/merged/sample${i}.R2.fastq.gz
done
