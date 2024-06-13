# mobivision_scrna_pipeline

## 一、运行环境配置
### 1.1、脚本下载

构建项目文件夹，并从从gitlab下载mobivision_scrna_pipeline对应脚本（内网clone,无需联网），分析请在联网节点运行。

```
project=HT2021-XXXX
species=Human
mkdir ${project}-${species}
cd ${project}-${species}
git clone   http://gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/mobivision_scrna_pipeline.git
# 如git未配置用户名、密码:
# git clone  http://<your user name>:<password>@gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/mobivision_scrna_pipeline.git
```
### 1.2 配置运行环境
```
source scrna_pipeline/envs/.bashrc   #加载环境

mlns  
```
## 二、填写项目相关文件

### A. cellranger流程

**填写好config/config.yaml 配置文件。** 配置文件包括以下4个模块：< 1.report 项目相关信息 2. email 跟邮件反馈相关参数 3. cellranger_params 跟cellranger质控相关参数 4. report_params 跟出报告流程相关参数 5.亚群分析相关参数 6.各步骤投递cpu数.>

在四个模块中第一二模块填写的是项目层面的参数，不论是质控还是出报告都需要逐行确认以及填写的；三四模块分别为质控和出报告两个流程各自需要用到的参数。


```bash
$ cat ./config/config.yaml 
## ========================= Project information =========================
report:
  Project_Num: HT2021-test  #项目编号
  Task_Num: DOExxx-bx      #任务单号
  Customer:                 #客户姓名
  Sales:                    #销售姓名
  Library:"墨卓3'转录组"     #墨卓的文库类型
  Executor: OExxxx          #执行人员编号
  Species: 小鼠-心肌组织          #实验物种
  Sample_Num: 1              #样本数目
  Others: 3' scRNA          #如果是5’建库项目，请填写 "5' scRNA"。Please fill in "5' scRNA" if the project is 5′-end scRNA-seq.
  Reference:                #参考基因组路径，以下两行为人小鼠的最新版本路径
  #Reference: /public/scRNA_works/works/guochy/Mobivision_demo/GRCh38_10X_reference/mobivision-GRCh38-2020-A # Lastest human ref
  #Reference: /public/scRNA_works/works/guochy/Mobivision_demo/mm10_10X_reference/mobivision-mm10-2020-A # Lastest mouse ref

  # Project_Path:
  # Project_Start_Time:
  # Project_End_Time:
  # Data_Cleaning_Time:

## ========================= Email feedback =========================
email:
  emailto:
    - "xxx@oebiotech.com"    # 质控反馈邮件收件人，填写自己的邮箱。
  Project_management: "朱卉" # 项目部责任人
  Laboratory: "毕远芝经理" # 实验部门负责人

## ========================= Mobivision parameters =========================
mobivision_params:
  module:
    mobivision: true # change cellranger to mobivision
    upload_bam: False # 是否把bam上传到OBS，默认False (默认不生成bam文件)
    mobivision_report: False # 是否生成mobivision报告，默认不出，只质控不分析项目选择True
  version: "v2.1"  #选择mobivision的版本，默认是2.1
  merge_samples: False # 是否merge samples,默认False
  intron: included  #默认included
  raw_data_obs_address: #可以列表形式填写多行，路径需要以'/'结尾。没有obs路径 填写 - "" 即可,默认不用填。
    - ""
  samples: all #用法1. 填写all 自动从raw_data_obs_address 获取所有样本名。
  #samples:    #用法2. 指定样本名。
  # - S1
  # - S2
  #samples:    #用法3. 指定样本名以及force-cells数量。注：方法2和方法3不可混用。
  # - S1: 5000
  Number_of_Reads: "350M"

```


### B. 出报告流程：
**除了第一二部分以外，填写config/config.yaml 第四部分内容** 

```bash
$ tail config/config.yaml 
## ========================= Report parameters =========================
report_params:
  module:
    cellranger_aggr: false
    create: false
    qc: false
    clustering: false
    marker: false
    celltyping: false
    diffexp: false
    report: false
    report_upload: false
  envmodules: "OESingleCell/3.0.d"
  samples_file: "config/samples.csv"
  libraries_file: "config/libraries.csv"
  diffexp_file: "config/diffexp.tsv"
  #ngene_numi_fold2sd: 2
  cut1: "median"  #可选median(中位数)，mean(平均值)
  cut2: "mad"     #可选median(中位数)、sd(标准差)或mad(绝对中位差)
  log10GenesPerUMI: 0.7  #最小值
  percent_mito: NULL # NULL: 自动判断
  percent_HB: 0.05
  rmdoublets: TRUE
  rmdoublets_methods: DoubletFinder
  roboust_linear_model: NULL #根据拟合线性模型处理离域值。例如nCount_RNA:nFeature_RNA;默认值：NULL
  reduct1: "pca"
  reduct2: "umap"
  resolution: 0.4
  pointsize: 0.5 # NULL：自动判断
  marker_test: "presto"
  diffexp_test: "presto"
  celltypingdb: "hpca" # Human: hpca/blueprint_encode/schcl # Mouse: immgen/mouse.rnaseq/scmca

  ```
**填写samples.csv 样本信息文件**
只需填写 sampleid，group，batchid，其余部分会自动填充。
```bash
$ cat ./config/samples.csv 
sampleid,species,group,batchid
S1,,S1,1
S2,,S2,2
```

首次使用需要配置~/.obsutilconfig文件，请咨询相关人员。

## 三、snakemake 流程提交 

### 1. 检查Snakefile配置是否正确 

```bash 
smt Snakefile 
```

### 2.绘制流程运行图（可略）

```bash 
smpl Snakefile pdf run_dag.pdf 
```
### 3.提交任务
请在能联网的节点运行,自动反馈邮件脚本、PPI需要联网才能运行。

```bash
smq Snakefile 
```





