# 10X空间转录组分析手册

## 一、运行环境配置
### 1.1、脚本下载
构建项目文件夹，并从从gitlab下载对应脚本snake_visium_workflow下载运行脚本
```
mkdir ${project}-${species}_${tissue}
cd ${project}-${species}_${tissue}
module load git 
git  clone  http://gitlab.oebiotech.cn/SingleCell/singlecell-snakemake-pipeline/st_pipeline -b produce_hwy
```

### 1.2 配置运行环境
```
source st_pipeline/envs/.bashrc   
mlns  
```

## 二、数据准备

### 2.1 填写项目配置文件
详细说明情况参考config/config.yaml文件

### 2.2 填写样本信息及差异筛选分组信息
#### 样本信息（config/samples.csv）
```bash
sampleid,slide,slide_area,fastq,image,cytaimage,species,group,batchid,image_json_file
human_prostate_acinar,V11J26-002,B1,raw_data/human_prostate_acinar/human_prostate_acinar_fastqs,raw_data/images/human_prostate_acinar.tif,,人-前列腺癌,human_prostate_acinar,1,
human_prostate_cancer,V11J26-003,B1,raw_data/human_prostate_cancer/human_prostate_cancer_fastqs,raw_data/images/human_prostate_cancer.tif,,人-前列腺癌,human_prostate_cancer,2,
```
参数说明：

  - **sampleid** 样本名称，必填项； 
  - **slide** Visium芯片序列号信息（如：V19L29-035，如有多个芯片序列号，会用括号标注样本序号），根据obs://oe-scrna/空间转录组HE染色图片/下对应项目质检报告填写，必填项； 
  - **slide_area** Visium芯片捕获区ID（如A1,B1,C1,D1）与样本的对应关系，根据obs://oe-scrna/空间转录组HE染色图片/下对应项目质检报告填写，必填项； 
  - **fastq** fastq路径，如：raw_data/12-5d_1/12-5d_1_fastqs/，可不填，将在2.3中自动填充
  - **image** H&E图像路径（在cytassist FFPE项目中为图像校正所需显微镜图像路径），如：raw_data/images/12-5d_1.jpg，可不填，将在2.3中自动填充
  - **cytaimage** cytassist FFPE项目图像路径，如：raw_data/images/cytaimage/12-5d_1.jpg，可不填，如果config.ymal中library_type为cytassist将在2.3中自动填充
  - **species** 物种信息，可不填，将在2.3中根据config.ymal中species参数自动填充 
  - **group** 分组信息，如有，需填写 
  - **batchid** 批次效应，如有，需填写
  - **image_json_file** 图像校正json文件，如有，需填写
#### 差异筛选分组信息（config/diff_group.tsv）

```bash  
treatment       treatment_name  control control_name    type
human_prostate_acinar   human_prostate_acinar   human_prostate_cancer   human_prostate_cancer   sampleid
```
参数说明：

  - **treatment** 实验组样本名称，需填写
  - **treatment_name** 实验组样本名称，需填写
  - **control** 对照组样本名称，需填写
  - **control_name** 对照组样本名称，需填写
  - **type** 分组类别，如sampleid或者group 

### 2.3.数据量确认

每次的下机安排邮件中会附有数据产出表格，请查看对应项目号的 **Raw Base（G）/ Reads（M）** 列，如果合计下机数据量不足合同要求，请及时反馈项目部。
单细胞下机数据产出不足的可容许范围是合同要求*10%，具体规则如下： 

    - 合同400M，最低要求360M。
    - 合同250M，最低要求225M。

Illumina/MGI 可以根据数据路径分辨：

    - obs://oebiotech/Illumina-scRNA/
    - obs://oebiotech/MGI2000-scRNA/          
                            
如果下机数据量低于“最低要求”，则反馈项目部安排补测。

**空间转录组数据分析存在一定探索性，分析过程中请注意查看各项分析结果。如存在质控或其他分析结果异常请及时调整参数或方法，切勿只按固定参数运行！**

### 2.5.分析注意事项 
- 1. 从今天起下单的单细胞项目，一个项目号一律仅出具一份常规报告。所有样本质控完成，提供确认单后方可启动常规分析。事业部已给全体销售发布公告，如有销售咨询，可直接回复是单细胞规定。确有特殊情况，需走OA评估。

- 2. 所有项目报告一律上传obs，包括完整质控报告，obs存放路径规则如下图， obs://oe-scrna/Analysis_Report//HT2021-XXXX/报告压缩包 。
发送常规报告时邮件主题务必统一关键字【项目报告】，格式参照：【项目报告】HT2021-xxxx xx老师 单细胞转录组/免疫组/ATAC/空间转录组项目报告。对于只质控不分析项目，质控反馈与完整质控报告同时发送，邮件主题统一关键字：【质控反馈和项目报告】，格式参照：【质控反馈和项目报告】HT2021-xxxx xx老师 单细胞转录组项目质控反馈及报告。

- 3. 今天起下单的单细胞项目，质控反馈邮件发送给顾婕，抄送朱卉、陆瑶、张倩和单细胞邮箱（空间转录组除外，质控反馈仍发给朱卉）。涉及项目报告（包括完整质控报告），发给朱卉，抄送陆瑶和单细胞邮箱。

- 4. 今天起下单的所有项目，质控反馈完就提交OA流程（包括空间转录组），项目部会以邮件形式发送确认单，然后执行分析，10个工作日内出具报告。（仅签订1个样本的项目，质控反馈后直接分析，出具报告后提交OA，


## 三、snakemake流程提交 
直接提交run_pipeline.sh ,会自动下载数据并执行分析
```bash
source envs/.bashrc  ##初始化环境
smt pre_Snakefile >logs/pre_smt.log ##检查数据下载流程
smq pre_Snakefile  ##执行数据下载，如果为首次使用，需要配置~/.obsutilconfig文件，请咨询相关人员。
smt Snakefile >logs/smt_$(date +%Y%m%d%H%M%S).log ##检查整个流程运行是否顺利，并打印执行命令
smq Snakefile ##运行流程
```

## 四、报告核查

**脚本检查完毕后，请参考以下要点再人工核查一遍：** 

  - 质控过滤后spot数目是否正常
  - 样本分析名与下机名称是否对应
  - Visium芯片序列号和Visium芯片捕获区ID是否对应
  - 数据量是否足够
  - 基因组比对率是否合格
  - 差异检验实验对照是否正确
  - 备注有无额外分析要求
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

## 六、上传和数据清理

1: 出具完整报告后，将`Spaceranger/`结果上传至  `obs://oe-store/work/scRNA/个人目录` ，上传完成后即可删除华为云上`Spaceranger结果`、`rawdata、filtered_seurat.rds`、`seurat_seurat.rds`，保留分析报告和最终生成的`singlecell_object.clustering_resolution*.rds`。

2: 出具报告后OBS上的Spaceranger结果保留4个月，过期后请定期删除。华为云上报告和RDS文件保留6个月，如半年以上无后续要求即可上传OBS存档。

**由于空间转录组没有cleandata， raw data保存在OBS， 所以分析完成后可直接删除华为云上的raw data。**

执行人:   
监督人:   

