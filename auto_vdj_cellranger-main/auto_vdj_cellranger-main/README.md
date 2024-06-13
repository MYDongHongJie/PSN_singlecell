# 根据归并完成产出表解析参数自动完成cellranger count分析——弹性云版

## 一、流程概要
cellranger
1.单细胞自动化归并流程会在该目录下生成 cellranger_config.csv 文件。  
  yaml_prepare脚本会自动根据该文件解析生成 config.yaml 文件。  
2.运行data_prepare脚本下载并整理raw_data。  
3.运行cellranger count 分析.  
  
email
1.OA流程下来会自动抓取生成 email_config.csv 文件在该目录里。  
  yaml_prepare脚本会自动根据该文件解析重新添加 config.yaml 文件内容。  
2.运行 cellranger_sample_stat 整理cellranger结果。  
3.运行 cellranger_summary 统计cellranger结果。  
4.运行 cellranger_check 根据具体情况发送质控邮件。  
  同时运行 cellranger_sor 上传登记质控信息。  
5.判断upload_bam是否运行bam_upload，上传bam文件。  
6.判断 cellranger_report 是否生成QC报告，是则生成并上传云报告，否则上传cellranger结果到 obs://oe-scrna/works/scRNA/PROJECT/ 目录。  
7.判断此前任务结果，运行 update_status 修改任务状态。  

## 二、snakemake 流程提交

### 1. 加载环境并自动生成config.yaml文件 

```bash
source scripts/envs/.bashrc
/data/software/conda_envs/snakemake/bin/python scripts/yaml_prepare.py -c config/config.yaml -i cellranger_config.csv -e email_config.csv 
```
### 2. 检查Snakefile配置是否正确 

```bash 
smt cellranger.smk 
```
  
### 3.提交任务

```bash
smq cellranger.smk

```

### 4.自动化代码
```
source scripts/envs/.bashrc && 
/data/software/conda_envs/snakemake/bin/python scripts/yaml_prepare.py -c config/config.yaml -i cellranger_config.csv -e email_config.csv && 
smq cellranger.smk
```

