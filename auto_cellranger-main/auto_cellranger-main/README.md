# 根据归并完成产出表解析参数自动完成cellranger count分析

## 一、流程概要
1.单细胞自动化归并流程会在该目录下生成 cellranger_config.csv 文件。  
  yaml_prepare脚本会自动根据该文件解析生成 config.yaml 文件。  
2.运行data_prepare脚本下载并整理raw_data。  
3.运行cellranger count 分析.  
  
## 二、snakemake 流程提交

### 1. 加载环境并自动生成config.yaml文件 

```bash
source scripts/envs/.bashrc
python scripts/yaml_prepare.py -c config/config.yaml -i cellranger_config.csv

```
### 2. 检查Snakefile配置是否正确 

```bash 
smt cellranger.smk 
```
  
### 3.提交任务

```bash
smq cellranger.smk & 

```

### 4.自动化代码
```
source scripts/envs/.bashrc && 
python scripts/yaml_prepare.py -c config/config.yaml -i cellranger_config.csv && 
smq cellranger.smk & 
```
