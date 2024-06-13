## M20单细胞转录组标准分析流程

## 一、运行环境配置

### 1.1 脚本下载

构建项目文件夹，并从从gitlab下载对应脚本m20_scrna_pipeline下载运行脚本

```
mkdir ${project}-${species}_${tissue}
cd ${project}-${species}_${tissue}
module load git 
git clone http://gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/m20_scrna_pipeline.git -b master
```

### 1.2 配置运行环境

```
cd m20_scrna_pipeline
source envs/.bashrc  
mlns  
```

## 二、数据准备

### 2.1 填写样本信息文件samples.csv

samples.csv文件共有4列，分别为sampleid,species,group,batchid，分别表示样本ID，物种名，组别，批次，每个样本一行，填写完整即可。

### 2.2 填写差异分组信息文件diff_group.tsv

diff_group.tsv文件共有5列，分别为treatment,treatment_name,control,control_name,type，分别表示处理组，处理组名，对照组，对照组名，组所属类型，同单细胞转录组差异分析输入表格。

### 2.3 cpu以及内存配置文件

cluster.yaml文件为cpu和内存配置文件，分模块配置相应cpu和内存。

### 2.4 参数配置文件

config.yaml文件为参数配置文件，共有8个模块，分别如下

- **genome** 参考基因组制作
- **STARsolo** STARsolo模块，运行STARsolo进行比对注释
- **QCreport**  质控报告模块，创建质控分析报告
- **Create**  seurat对象创建模块，读取STARsolo结果矩阵，创建seurat对象
- **QC** seurat质控模块，对数据进一步质控
- **Clustering** 降维聚类模块，对数据进行降维聚类
- **Celltype**  细胞类型鉴定模块，进行细胞类型的鉴定
- **Marker**  标记基因查找模块，寻找每个cluster的显著性标记基因
- **Diffexp**  差异基因查找模块，寻找组间差异基因
- **Report**  标准分析报告生成模块，生成相应报告

配置文件具体内容如下：

```yaml
## 不同模块部分
module:
  genome: False
  STARsolo: True
  QCreport: True
  Create: True
  QC: True
  Clustering: True
  Celltype: True
  Marker: True
  Diffexp: True
  Report: True

## 不同模块参数
## STARsolo部分参数
STARsolo:
    empty_drops: "" #"5000"填入具体的细胞数
    raw_data_obs_address: #可以列表形式填写多行，路径需要以'/'结尾。填写 - "" 则根据任务单号自动生成obs链接。
      - ""
    samples: all #用法1. 填写all 自动从raw_data_obs_address 获取所有样本名。
    #samples:    #用法2. 指定样本名。
    # - S1
    # - S2
    raw_data: False # 所下载的数据是否为raw_data，若raw_data为True则自动转为clean_data，不同批次类型（rawdata/cleandata）相同可直接合并
    ## 数据自动下载部分可以自动合并试测和补测raw_data/clean_data数据，若已经合并了，请手动下载设置
    platte: "ditto" #绘制基因的reads覆盖度图的调色板

## QC report的回复邮箱
QC_report:
  emailto: "xxx@oebiotech.com"
## QC参数
QC:
  vars2regress: "nCount_RNA"
  filters: "nFeature_RNA,nCount_RNA,log10GenesPerUMI,percent.mito,percent.HB,percent.lncRNA,percent.rRNA"
  lower: "200,1000,0.7,0,0,0,0"
  upper: "Inf,Inf,Inf,0.15,0.05,Inf,Inf"
  features2filter: True
  mincell4gene: 0.01
  normmeth: "sctransform"
  rmdoublets: "TRUE"
  method: "doubletfinder"
## 降维参数
Clustering:
  RNA_reduct1: "mnn"
  RNA_reduct2: "umap"
  RNA_clusteringuse: "slm"
  RNA_resolution: "0.4"
  RNA_palette: "customecol2"
  RNA_component: "10"
  RNA_groupby: "sampleid,group"
  RNA_batchid: "batchid"
## marker基因查找参数
Marker:
  RNA_topn_marker: "10"
  RNA_test_method: "presto"
  RNA_topn: "10"
  RNA_topby: "gene_diff"
  RNA_vismethod: "featureplot,vlnplot,heatmap"
## 细胞类型鉴定参数
Celltyping:
  refbuiltin: "hpca"
  annolevel: "main"
  demethod: "wilcox"
## 差异基因参数
Diffexp:
  FC: "1.5"
  pvalue: NULL
  padj: "0.05"
  test: "presto"
  topn: 25
```

### 2.5 参考基因组制作

```shell
#!/usr/bin/bash
module load STRS/1.0.0

OUTDIR="STAR_index"
mkdir -p $PWD/${OUTDIR}

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles genome.fa  \ ##下载的基因组fasta文件
--sjdbGTFfile  genomic.gtf  \ ## 下载的gtf文件
--genomeSAindexNbases 10 \
--sjdbGTFfeatureExon gene \
--sjdbGTFtagExonParentTranscript gene_id \
--sjdbOverhang 122
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
git clone http://gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/m20_scrna_pipeline.git -b master

# step2：下载Demo数据
module load obsutil/5.2.12
cd m20_scrna_pipeline
mkdir result
obsutil cp obs://oe-scrna/works/scRNA/sunkun/M20_scRNA/clean_data/  result/ -r -f -vlength
mv result/clean_data result/0.clean_data
## 注意：Demo数据的fasta格式已经调整好，新的数据需要将文件夹名命名为样本名，fasta格式为sampleid_1.fq.gz, sampleid_2.fq.gz

#step3: 配置参数文件
#对于Demo数据，参数部分不用变，仅需修改以下参数
QC_report:
    emailto: "xxx@oebiotech.com" #此处换成自己的邮箱
  
#step4：流程运行
source envs/.bashrc
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
