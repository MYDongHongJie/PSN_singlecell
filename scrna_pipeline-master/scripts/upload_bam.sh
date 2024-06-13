#!/bin/bash
# 0. 配置免密
# cat ~/.ssh/id_rsa.pub | ssh scrna@10.100.10.69 'cat >> ~/.ssh/authorized_keys'
# cat ~/.ssh/id_dsa.pub | ssh ${remote_host} 'cat >> ~/.ssh/authorized_keys'
# 1. 指定目录路径和输出文件名
sample=$1
project_id=$2
#sample_id=$(echo $dir_path | cut -d'/' -f$(expr $(echo $dir_path | tr -dc '/' | wc -c) + 1))
remote_host="scrna@10.100.10.69"
dir_path=${sample}/outs/
bam1=${sample}/outs/possorted_genome_bam.bam
bam2=${sample}/outs/possorted_genome_bam.bam.bai
if test -f "$bam1" &&  test -f "$bam2"; then
  echo $project_id
  echo $sample
  echo 开始生成md5
  # 2. 遍历目录中的.bam文件
  # 3. 生成MD5校验码并将其追加到输出文件中
  cd $dir_path
  output_file="$dir_path/${project_id}_md5.xls"
  if test -f ${sample}/possorted_genome_bam.bam && test -f ${sample}/possorted_genome_bam.bam.bai && test -f ${output_file}; then
    echo "已生成md5，直接上传"
  else
    mkdir $sample
    mv possorted_genome_bam.bam ${sample}
    mv possorted_genome_bam.bam.bai ${sample}
    md5sum ${sample}/* >> ${project_id}_md5.xls
  fi
  cd -
  # 创建目标目录
  ssh -t ${remote_host} "[ -d /cds/scRNA-bambak/${project_id} ] && echo ok || mkdir -p /cds/scRNA-bambak/${project_id} "
    # 4. 将.bam文件上传到远程服务器
  echo 开始上传bam
  # 如覆盖bam，则先删除原先bam的md5
  remote_file=/cds/scRNA-bambak/${project_id}/${project_id}_md5.xls
  if ssh $remote_host "[ -d /cds/scRNA-bambak/${project_id}/${sample} ]"; then
    ssh ${remote_host} "grep -v ${sample} ${remote_file} | sed '/^$/d' > temp  && mv temp ${remote_file}"
  fi
  # scp -r ${dir_path}/${sample}/ ${remote_host}:/cds/scRNA-bambak/${project_id}/ && rm -r ${dir_path}/${sample}/
  scp -r ${dir_path}/${sample}/ ${remote_host}:/cds/scRNA-bambak/${project_id}/
  # Check if the upload was successful
  if [ $? -eq 0 ]; then
      # Check if remote file size is non-zero
      # remote_size=$(ssh "${remote_host}" "du -s /cds/scRNA-bambak/${project_id}/${sample}/" | cut -f1)
      zero_kb_files=$(ssh "${remote_host}" "find /cds/scRNA-bambak/${project_id}/${sample}/ -type f -size 0")
      # if [ "${remote_size}" -ne 0 ]; then
      if [ -z "${zero_kb_files}" ]; then
          # Modification time and size are non-zero, delete local files
          rm -r "${dir_path}/outs/${sample}/"
          # echo "rm -r ${dir_path}/outs/${sample}/"
          echo "Files uploaded successfully and remote modification time and size are valid. Local files deleted."
      else
          # Remote size is zero, exit with an error
          echo "Error: Remote file size is zero. Exiting with status 1."
          exit 1
      fi
  else
      # Upload failed, exit with an error
      echo "Error: Upload failed. Exiting with status 1."
      exit 1
  fi
  # 5. 将输出文件上传到远程服务器
  echo 开始上传$output_file
  # chmod 777 $output_file
  if ssh $remote_host "[ -f $remote_file ]"; then
    cat "$output_file" | ssh ${remote_host} "cat >> ${remote_file}" && rm $output_file 
  else
    scp "$output_file" ${remote_host}:${remote_file} && rm $output_file 
  fi
else
  echo 目录不存在，跳过bam上传
fi

