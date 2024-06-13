#!/bin/bash
# 0. 配置免密
# ssh-keygen -t rsa
# ssh-copy-id scrna@10.100.10.69 #发送公钥  # 密码：oesinglecell
# 1. 指定目录路径和输出文件名
complete_dir=$1
# project_id=$(echo $dir_path | cut -d'/' -f$(expr $(echo $dir_path | tr -dc '/' | wc -c) - 1))
project_id=$2

# project_id=$(echo $dir_path | cut -d'/' -f$(expr $(echo $dir_path | tr -dc '/' | wc -c) + 1))
remote_host="scrna@10.100.10.69"
# bam1=$dir_path/outs/possorted_genome_bam.bam
# bam2=$dir_path/outs/possorted_genome_bam.bam.bai
# output_file="$dir_path/outs/${project_id}_md5.xls"
remote_dir=/cds/obs-scrnabak/auto_cellranger
if ! ssh -t ${remote_host} "[ -d ${remote_dir}/${project_id} ] && echo ok || mkdir -p ${remote_dir}/${project_id} " ; then
  # 读取要添加的参数
  # new_row="${complete_dir} ${project_id}"
  # 在文件末尾追加新行
  # echo "${new_row}" >> /public/cloud_scRNA/scrna_auto_task/file/fail_cellranger.txt
  echo "本地服务器连接异常，上传失败"
  # echo $new_row
  exit 1
fi 
echo 开始上传$complete_dir
scp -r ${complete_dir} ${remote_host}:${remote_dir}/${project_id}/ && echo '上传成功，工作目录可以删除'

