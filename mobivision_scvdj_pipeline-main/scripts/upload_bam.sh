# 0. 配置免密
# cat ~/.ssh/id_rsa.pub | ssh ${remote_host} 'cat >> ~/.ssh/authorized_keys'
# cat ~/.ssh/id_dsa.pub | ssh ${remote_host} 'cat >> ~/.ssh/authorized_keys'
# 1. 指定目录路径和输出文件名
sample=$1
project_id=$2
#sample_id=$(echo $dir_path | cut -d'/' -f$(expr $(echo $dir_path | tr -dc '/' | wc -c) + 1))
remote_host="scrna@10.100.10.69"
dir_path=result/cellranger/${sample}/outs/per_sample_outs/${sample}/count/
bam1=result/cellranger/${sample}/outs/per_sample_outs/${sample}/count/sample_alignments.bam
bam2=result/cellranger/${sample}/outs/per_sample_outs/${sample}/count/sample_alignments.bam.bai
if test -f "$bam1" &&  test -f "$bam2"; then
  echo $project_id
  echo $sample
  echo 开始生成md5
  # 2. 遍历目录中的.bam文件
  # 3. 生成MD5校验码并将其追加到输出文件中
  cd $dir_path
  output_file="$dir_path/${project_id}_md5.xls"
  if test -f ${sample}/sample_alignments.bam && test -f ${sample}/sample_alignments.bam.bai && test -f ${output_file}; then
    echo "已生成md5，直接上传"
  else
    mkdir $sample
    mv sample_alignments.bam ${sample}
    mv sample_alignments.bam.bai ${sample}
    md5sum ${sample}/* > ${project_id}_md5.xls
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
  scp -r ${dir_path}/${sample}/ ${remote_host}:/cds/scRNA-bambak/${project_id}/ && rm -r ${dir_path}/${sample}/
  # 5. 将输出文件上传到远程服务器
  echo 开始上传$output_file
  # chmod 777 $output_file
  if ssh $remote_host "[ -f $remote_file ]"; then
    cat "$output_file" | ssh ${remote_host} "cat >> ${remote_file}" && rm $output_file 
  else
    scp "$output_file" ${remote_host}:${remote_file}  && rm $output_file
  fi
else
  echo 目录不存在，跳过bam上传
fi

