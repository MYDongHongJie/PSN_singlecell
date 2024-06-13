#!/bin/bash
# 0. 配置免密
# cat ~/.ssh/id_rsa.pub | ssh ${remote_host} 'cat >> ~/.ssh/authorized_keys'
# cat ~/.ssh/id_dsa.pub | ssh ${remote_host} 'cat >> ~/.ssh/authorized_keys'
# 1. 指定目录路径和输出文件名
dir_path=$1
# 判断最后一个字符是否是/
if [ "${dir_path: -1}" = "/" ]; then
    # 如果是/，去掉最后一个字符
    dir_path="${dir_path%/}"
    echo "去掉最后一个/后的变量值为: $dir_path"
else
    echo "变量值不以/结尾"
fi
# project_id=$(echo $dir_path | cut -d'/' -f$(expr $(echo $dir_path | tr -dc '/' | wc -c) - 1))
project_id=$2
sample_id=$(echo $dir_path | cut -d'/' -f$(expr $(echo $dir_path | tr -dc '/' | wc -c) + 1))
remote_host="scrna@10.100.10.69"
a=$dir_path/*.BAM
echo $a
base_name=`basename ${a} .BAM`
bam1=$dir_path/${base_name}.BAM
bam2=$dir_path/${base_name}.BAM.bai
bam3=$dir_path/${base_name}.BAM.csi
bam1file=`basename $bam1`
bam2file=`basename $bam2`

output_file="$dir_path/${project_id}_md5.xls"

if test -f "$bam3" ; then
  echo 'csi exist'
  bam2=`echo $bam3`
fi

bam_files2=$(find "${dir_path}/${sample_id}/" -type f -name "*.BAM")

if { test -f "$bam1" && test -f "$bam2"; } || { test -f "${bam_files2}" && test -f "${bam_files2}.bai" && test -f "${dir_path}/${project_id}_md5.xls"; }; then
# if test -f "$bam1" &&  test -f "$bam2" ; then
  echo $project_id
  echo $sample_id
  echo 开始生成md5
  # 2. 遍历目录中的.bam文件
  # 3. 生成MD5校验码并将其追加到输出文件中
  cd ${dir_path}/
  output_file2=${project_id}_md5.xls
  if  test -f "${bam_files2}" && test -f "${bam_files2}.bai" && test -f "${project_id}_md5.xls"; then
    echo "已生成md5，直接上传"
  else
    mkdir $sample_id
    mv $bam1file ${sample_id}/${base_name}.BAM
    mv $bam2file ${sample_id}/${base_name}.BAM.bai
    md5sum ${sample_id}/* > ${project_id}_md5.xls
  fi
  # 创建目标目录
  if ! ssh -t ${remote_host} "[ -d /cds/scRNA-bambak/${project_id} ] && echo ok || mkdir -p /cds/scRNA-bambak/${project_id} " ; then
    echo "本地服务器连接异常，上传失败"
    #exit 1
  fi 
    # 4. 将.bam文件上传到远程服务器
  echo 开始上传bam
  # 如覆盖bam，则先删除原先bam的md5
  remote_file=/cds/scRNA-bambak/${project_id}/${project_id}_md5.xls
  if ssh $remote_host "[ -d /cds/scRNA-bambak/${project_id}/${sample_id} ]"; then
    ssh ${remote_host} "grep -v ${sample_id} ${remote_file} | sed '/^$/d' > temp  && mv temp ${remote_file}"
  fi
  # scp -r ${dir_path}/${sample_id}/ ${remote_host}:/cds/scRNA-bambak/${project_id}/ && rm -r ${dir_path}/${sample_id}/
  scp -r ${sample_id}/ ${remote_host}:/cds/scRNA-bambak/${project_id}/
  # Check if the upload was successful
  if [ $? -eq 0 ]; then
      # Check if remote file size is non-zero
      # remote_size=$(ssh "${remote_host}" "du -s /cds/scRNA-bambak/${project_id}/${sample_id}/" | cut -f1)
      zero_kb_files=$(ssh "${remote_host}" "find /cds/scRNA-bambak/${project_id}/${sample_id}/ -type f -size 0")
      # if [ "${remote_size}" -ne 0 ]; then
      if [ -z "${zero_kb_files}" ]; then
          # Modification time and size are non-zero, delete local files
          rm -r "${sample_id}/"
          # echo "rm -r ${dir_path}/${sample_id}/"
          echo "Files uploaded successfully and remote modification time and size are valid. Local files deleted."
      else
          # Remote size is zero, exit with an error
          echo "Error: Remote file size is zero. Exiting with status 1."
          #exit 1
      fi
  else
      # Upload failed, exit with an error
      echo "Error: Upload failed. Exiting with status 1."
      #exit 1
  fi
  # 5. 将输出文件上传到远程服务器
  echo 开始上传$output_file2
  # chmod 777 $output_file
  if ssh $remote_host "[ -f $remote_file ]"; then
    cat "$output_file2" | ssh ${remote_host} "cat >> ${remote_file}" && rm $output_file2 
  else
    scp "$output_file2" ${remote_host}:${remote_file} && rm $output_file2 
  fi
else
  echo 目录不存在，跳过bam上传
fi

