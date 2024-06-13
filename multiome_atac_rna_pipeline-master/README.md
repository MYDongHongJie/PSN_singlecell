# 单细胞多模态流程使用
## 1.下载分析流程并载入环境
```sh
git clone http://gitlab.oebiotech.cn/sunkun/multiome_atac_rna_pipeline.git
source envs/.bashrc
```
## 2.修改config.yaml文件
### 2.1质控分析修改内容
```yaml
##===================== 1.Basee information =============================

## ==================== 1.1 Project information ============================
report:
  #project name id ,eg:HT2020-10440
  Project_Num: HT2021-Demo
  #任务单号
  Task_Num:
  #customer name, eg:汤耀辉
  Customer: test
  #sales name,eg:祖启东
  Sales: test
  #deal by who, eg:OE094
  Executor: OE0926
  #specieas, eg: Human 或 Mouse,首字母大写
  Species: Human
  #sample numbers ,eg: 4
  Sample_Num: 2
  Project_Type: '单细胞多模态'
  #project path ,will be filled automaticly
  Project_Path: ''
  #project start time ,will be filled automaticly
  Project_Start_Time: ''
  #project end time, for report time, eg:2021-02-05
  Project_End_Time: ''
  Data_Cleaning_Time: ''
  Others: ''
  Reference: /data/database/cellranger-refdata/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 #human
  #/data/database/cellranger-refdata/refdata-cellranger-arc-mm10-2020-A-2.0.0 #mouse
## ===================== 1.2 Email feedback ==============================
email:
  emailto:
    - "kun.sun@oebiotech.com"    # 质控反馈邮件收件人，填写自己的邮箱。

##====================== 2.Mdule control===============================
module:
  cellranger: True
  create: false
  qc: false
  clustering: false
  marker: false
  celltyping: false
  wnn: false
  diffexp: false
  report: false

## ====================== 3.Parameters control ===========================

## ======================= 3.1 Cellranger params =========================
cellranger_params:
  aggr: false #是否对多个样本进行合并
  envmodules: "cellranger-arc/2.0.0"
  raw_data_obs_address: #可以列表形式填写多行，路径需要以'/'结尾。没有obs路径 填写 - "" 即可。
    - "obs://oe-scrna/works/scRNA/sunkun/multimodal/HT2021-multimodal-pipeline/"
```
质控分析完成后，根据质控结果，无异常发送【质控反馈】报告，有异常发送【质控异常】主题报告。
## 2.2 标准分析
标准分析部分，根据项目情况修改运行参数
```yaml
##=======================3.2 Create params =============================
create_params:
  metrics: "percent.mito,percent.fragment.peaks"

##=======================3.3 QC params ===============================
qc_params:
  vars2regress: "nCount_RNA|nCount_ATAC"
  filters: "nFeature_RNA,nCount_RNA,log10GenesPerUMI,atac_fragments,TSS.enrichment,nucleosome_signal"
  lower: "200,1000,0.7,1000,2,0"
  upper: "Inf,Inf,Inf,100000,100,2"
  features2filter: "NULL"
  mincell4gene: 1
  rmdoublets: "TRUE"
  method: "doubletfinder"
  normmeth: "sctransform|logtf-idf"
##=========================3.4 clustering ==============================
cluster_params:
  #RNA
  RNA_reduct1: "mnn"
  RNA_reduct2: "umap"
  RNA_clusteringuse: "snn"
  RNA_resolution: "0.4"
  RNA_palette: "blindless"
  RNA_component: "10"
  RNA_groupby: "sampleid,batch"
  #ATAC
  ATAC_reduct1: "harmony"
  ATAC_reduct2: "umap"
  ATAC_clusteringuse: "slm"
  ATAC_resolution: "0.4"
  ATAC_palette: "blindless"
  ATAC_component: "30"
  ATAC_groupby: "sampleid,batch"
##========================3.5 marker================================
marker_params:
  #RNA
  RNA_topn_marker: "10"
  RNA_test_method: "presto"
  RNA_topn: "10"
  RNA_topby: "gene_diff"
  RNA_vismethod: "featureplot,vlnplot,heatmap"
  #ATAC
  ATAC_topn_marker: "10"
  ATAC_test_method: "presto"
  ATAC_topn: "10"
  ATAC_topby: "gene_diff"
  ATAC_vismethod: "featureplot,vlnplot,heatmap"
##==============================3.6 celltyping ========================
celltyping_params:
  ##人数据集:<hpca,blueprint_encode,schcl> 小鼠数据集：<immgen,mouse.rnaseq,scmca>
  refbuiltin: "hpca"
  annolevel: "main"
  demethod: "wilcox"

##===============================3.7 wnn ===========================
wnn_params:
  clusteringuse: "slm"
  resolution: "0.4"
  groupby: "sampleid,batch"
  wnn_topn_marker: "10"
  wnn_test_method: "presto"
  wnn_topn: "10"
  wnn_topby: "gene_diff"
  wnn_vismethod: "featureplot,vlnplot,heatmap"
  refbuiltin: "hpca"
  annolevel: "main"
  demethod: "wilcox"
##=============================3.8 diffexp ===========================
diffexp_params:
  FC: "1.5"
  pvalue: "0.05"
  test: "presto"
```
标准分析完成后，将数据上传至OBS并发送【项目报告】邮件