## 10x Xenium 组织原位分析标准分析流程

## 一、运行环境配置

### 1.1 脚本下载

构建项目文件夹，并从从gitlab下载对应脚本m20_scrna_pipeline下载运行脚本

```
mkdir ${project}-${species}_${tissue}
cd ${project}-${species}_${tissue}
module load git 
git clone http://gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/xenium_pipeline.git
```

### 1.2 配置运行环境

```
cd ${project}-${species}_${tissue}
source xenium_pipeline/envs/.bashrc
```

## 二、数据准备

### 2.1 填写样本信息文件samples.csv

samples.csv文件共有4列，分别为sampleid,species,tissue,group，分别表示样本ID，物种名，组织名，组名，每个样本一行，填写完整即可。

### 2.2 cpu以及内存配置文件

cluster.yaml文件为cpu和内存配置文件，分模块配置相应cpu和内存。

### 2.3 参数配置文件

config.yaml文件为参数配置文件，共有8个模块，分别如下

- **Xenium** xenium运行模块 
- **QCreport** 质控报告模块，创建质控分析报告
- **Create** seurat对象创建模块，读取xenium结果矩阵，创建seurat对象
- **QC**  seurat质控模块，对数据进一步质控统计
- **Clustering**  降维聚类模块，对数据进行降维聚类
- **Celltype**  细胞类型鉴定模块，进行细胞类型的鉴定
- **Marker**  标记基因查找模块，寻找每个cluster的显著性标记基因
- **Report** 标准分析报告生成模块，生成相应报告


配置文件具体内容如下：

```yaml
## 不同模块部分
module:
  Xenium: True
  QCreport: True
  Create: True
  QC: True
  Clustering: True
  Celltype: True
  Marker: True
  Report: True

## 不同模块参数
Xenium:
  method: "relabel" ##可选“relabel”，“resegment” ，“Import_segment”，若为空则直接利用下机结果
  raw_data_obs_address: #可以列表形式填写多行，路径需要以'/'结尾。填写 - "" 则根据任务单号自动生成obs链接。
    - "obs://oe-scrna/works/scRNA/sunkun/xenium_demo/"
  raw_data_local_address: "" ##本地文件路径，与obs路径2选1，需要以“/”结尾
  samples: all #用法1. 填写all 自动从raw_data_obs_address 获取所有样本名。
  #samples:    #用法2. 指定样本名。挑选部分样本进行质控分析，仅限质控分析
  # - S1
  # - S2
  expansion_distance: 5 #细胞边界外阔大小,“resegment”参数
  resegment_nuclei: true #参数为"true" "false",是否重新定位细胞核位置，“resegment”参数
  coordinate_transform: "" #转录本划分数据文件，“Import_segment”参数
  nuclei: "" # 细胞核分割数据文件，“Import_segment”参数
  cells: "" # 细胞边界数据文件，“Import_segment”参数
  units: "" #"microns" or "pixels" 分割边界数值单位微米或者像素，“Import_segment”参数
QC_report:
  emailto: "xxx@oebiotech.com"

QC:
  vars2regress: "NULL"
  filters: "nCount_xenium,nFeature_xenium"
  lower: "NULL,NULL"
  upper: "NULL,NULL"
  normmeth: "sctransform"
  rmdoublets: "FALSE"
  crop: "TRUE" #图形是否合并，默认合并
  location: "v" # 图形方向，"v" 表示x < y 图形竖着，“h”表示x > y 图形横着

Clustering:
  reduct1: "pca" ##去批次 “pca,harmony”
  reduct2: "umap"
  clusteringuse: "snn"
  resolution: "0.2"
  palette: "customecol2"
  component: "30"
  batchid: "batchid"
  location: "v"
  splitby: "sampleid"
  groupby: "clusters"
  groups: "sampleid"
  propby: "clusters"
  ptsize: "0.5"

Marker:
  test_method: "wilcox"
  topn: "10"
  topby: "gene_diff"
  vismethod: "featureplot,vlnplot,heatmap"
  location: "v"
  topn_vis: "5"  # cluster过多时，可选择展示部分marker

Celltyping:
  location: "v"
  clusterby: "clusters"
  palette: "customecol2"
  ptsize: "0.5"
    
```



## 三、snakemake流程提交

输入文件及配置文件准备完毕后即可运行snakemake流程

```bash
source envs/.bashrc  ##初始化环境
smt Snakefile >logs/smt_$(date +%Y%m%d%H%M%S).log ##检查整个流程运行是否顺利，并打印执行命令，如有问题，请咨询相关人员。
smq Snakefile ##运行流程
```

## 四、Demo数据snakemake流程示范

```shell
# step1：下载流程
module load git
git clone http://gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/xenium_pipeline.git

# step2：下载Demo数据（可直接在config文件中填写路径地址）
module load obsutil/5.2.12
mkdir raw_data
obsutil cp obs://oe-scrna/works/scRNA/sunkun/xenium_demo/  raw_data/ -r -f

#step3: 配置参数文件
#对于Demo数据，参数部分不用变
  
#step4：流程运行
source xenium_pipeline/envs/.bashrc
smt Snakefile
smq Snakefile
```

## 五、报告核查

**脚本检查完毕后，请参考以下要点再人工核查一遍：**

- 文件生成是否正常
- 图片生成是否有误
- 参数是否合适
- 结果是否正常
  ……

**如果项目样本数超过2个，由于报告比较大不便上传邮件附件，请使用 `tar hcvf - HT2020-XXX*_Report/  | pigz -6 -p 10 -k >HT2020-XXX*_Report.tar.gz`  命令压缩报告，将压缩包上传至 `obs://oe-scrna/Analysis_Report/` ， 发报告时邮件提供obs链接即可。**

```
python3 -m http.server 8080
project_report=`ls -d1 *Report*`
project_name=${project}
tar hcvf - $project_report  | pigz -6 -p 10 -k > ${project_report}.tar.gz
module load obsutil/5.2.12
obsutil cp ${project_report}.tar.gz   obs://oe-scrna/Analysis_Report//${project_name}/ -r -f -flat -vlength -vmd5 
echo obs://oe-scrna/Analysis_Report/${project_name}/${project_report}.tar.gz
```

执行人:
监督人:
