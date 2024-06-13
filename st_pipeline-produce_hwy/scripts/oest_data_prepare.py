#!/public/dev_scRNA/software/miniconda3/envs/OESingleCell_3.00_visium_produce/bin/python
# encoding: utf-8
"""
Project : snake_visium_pipeline
Author  : Xiufeng Yang
Contact : xiufeng.yang@oebiotech.com
File   : data_prepare.py.py
IDE    : PyCharm
Time   : 2020-12-01 14:16:04
Desc   :
"""
import os, datetime, re, string
import yaml,time
import click
import subprocess
import pandas as pd
from glob import glob
from collections import OrderedDict
from docx import Document

exec(open('/data/software/modules/modules-v4.2.1/init/python.py').read())
os.environ['MODULEPATH'] =os.path.abspath("./envs/")

def file_filter(f):
    if f == 'images' or f[-4:] in [".gpr", ".xls"]:
        return False
    else:
        return True


def formated_sampleids(text):
    # 创建一个转换表，将标点符号和空格,换行符都映射为空字符
    translator = str.maketrans('', '', string.punctuation + ' ' + '\n')
    # 使用转换表去除标点符号和空格
    cleaned_text = text.translate(translator)
    if(cleaned_text[0].isdigit()):
        print(cleaned_text,"样品名字是数字开头，默认开头加S，改为","S"+cleaned_text)
        cleaned_text="S"+cleaned_text    
    return cleaned_text


def ordered_yaml_load(stream, Loader=yaml.SafeLoader,
                      object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def _construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        _construct_mapping)
    return yaml.load(stream, OrderedLoader)


def ordered_yaml_dump(data, stream=None, Dumper=yaml.SafeDumper,
                      object_pairs_hook=OrderedDict, **kwds):
    class OrderedDumper(Dumper):
        pass

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())
    OrderedDumper.add_representer(object_pairs_hook, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwds)


def raw_data_download(raw_data_obs_address, image_obs_address, outdir):
    """
    下载原始数据，返回raw_data里的sampleid
    """
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始下载raw_data")
    if not os.path.exists(f'{outdir}/raw_data/'):
        module('load', "obsutil/5.2.12")
        raw_data_batch = len(raw_data_obs_address)
        for i in range(raw_data_batch):
            if raw_data_batch == 1:
                subprocess.run(
                    f"obsutil cp  {raw_data_obs_address[i]}  {outdir}/raw_data/ -r -f -flat -vlength -vmd5", shell=True,
                    check=True)
            elif raw_data_batch > 1:
                subprocess.run(
                    f"obsutil cp {raw_data_obs_address[i]}  {outdir}/raw_data -r -f -flat -vlength -vmd5", shell=True,
                    check=True)
                old_names = [item for sublist in
                             [glob(f'{outdir}/raw_data/*' + ext) for ext in [f"/*.gz", "/*.md5"]]
                             for item in sublist]
                batch_names = [item for sublist in
                               [glob(f'{outdir}/raw_data/*' + ext) for ext in [f"/batch_*.gz"]]
                               for item in sublist]
                filter_names = [x for x in old_names if
                                all(y not in x for y in batch_names)]
                for old_name in filter_names:
                    file_path, file_name = os.path.split(old_name)
                    os.rename(old_name, os.path.join(file_path, "batch_" + str(i + 1) + "_" + file_name))
    else:
        print("\033[0;31m%s\033[0m" % f"Warning: {outdir}/raw_data/样本 不为空，如需重新下载请删除")

    if image_obs_address:
        if not os.path.exists(f'{outdir}/raw_data/images'):
            ###  download  slide images  file
            image_obs_batch = len(image_obs_address)
            for i in range(image_obs_batch):
                module('load', "obsutil/5.2.12")
                subprocess.run(
                    f"obsutil cp  {image_obs_address[i]}  {outdir}/raw_data/images -r -f -flat -vlength -vmd5",
                    shell=True, check=True)
        else:
            print("\033[0;31m%s\033[0m" % f"Warning: {outdir}/raw_data/images 不为空，如需重新下载请删除")

    # 使用file_filter函数过滤，获取实际样本名列表
    raw_sampleid_list = os.listdir(f"{outdir}/raw_data/")
    raw_sampleid_list = list(filter(file_filter, raw_sampleid_list))
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> raw_data下载完毕\n\n")

    return raw_sampleid_list


def raw_data_formart(samplesfiles, outdir):
    """
    创建rawdata软链接
    """
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始创建软链接")

    samples = pd.read_csv(samplesfiles, dtype=str).set_index("sampleid", drop=False)
    library = pd.DataFrame(columns=['library_id', 'molecule_h5', 'cloupe_file', 'spatial_folder'])
    samples_update = samples
    raw_samples = os.listdir(f"{outdir}/raw_data/")
    raw_samples = list(filter(file_filter, raw_samples))
    for sample in raw_samples:
        if sample in list(samples.index):
            if not os.path.exists(f'{outdir}/raw_data/{sample}/{sample}_fastqs'):
                os.mkdir(f'{outdir}/raw_data/{sample}/{sample}_fastqs')
            else:
                print("\033[0;31m%s\033[0m" % f"Warning: {sample}的fastqs文件夹已存在，现已删除并重新创建fastqs文件夹")
                os.system(f"rm -rf {outdir}/raw_data/{sample}/{sample}_fastqs")
                os.mkdir(f'{outdir}/raw_data/{sample}/{sample}_fastqs')

            ######################## # fq1_list = glob(f'{outdir}/raw_data/{sample}/*R{i}.fastq,fq.gz')
            fq1_list = [item for sublist in
                        [glob(f'{outdir}/raw_data/{sample}' + ext) for ext in [f"/*R1.fastq.gz", "/*_1.fq.gz","/*R1_001.fastq.gz", "/*R1.fq.gz"]]
                        for item in sublist]
            fq1_list.sort(key=os.path.getmtime)
            j = 1
            for fastq1 in fq1_list:
                ##R1
                des1 = f'{outdir}/raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_R1_001.fastq.gz'
                fastq_show1 = fastq1.replace(f"{outdir}", "")
                des_show1 = f"raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_R1_001.fastq.gz"
                if not os.path.exists(des1):
                    os.symlink(os.path.abspath(fastq1), des1)
                    print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show1}  {des_show1} ")
                else:
                    print("\033[0;31m%s\033[0m" % f'Warning: {des_show1} has existed!!!')

                ##R2
                fastq2 = fastq1.replace("R1.fastq.gz", "R2.fastq.gz").replace("_1.fq.gz", "_2.fq.gz").replace("R1_001.fastq.gz", "R2_001.fastq.gz").replace("R1.fq.gz", "R2.fq.gz")
                des2 = f'{outdir}/raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_R2_001.fastq.gz'
                fastq_show2 = fastq2.replace(f"{outdir}", "")
                des_show2 = f"raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_R2_001.fastq.gz"
                if not os.path.exists(des2):
                    os.symlink(os.path.abspath(fastq2), des2)
                    print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show2}  {des_show2} ")
                else:
                    print("\033[0;31m%s\033[0m" % f'Warning: {des_show2} has existed!!!')

                ##I1
                I1= fastq1.replace("R1.fastq.gz", "I1.fastq.gz").replace("R1_001.fastq.gz", "I1_001.fastq.gz").replace("R1.fq.gz", "I1.fq.gz")
                des3= f'{outdir}/raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_I1_001.fastq.gz'
                fastq_show3 = I1.replace(f"{outdir}", "")
                des_show3 = f"raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_I1_001.fastq.gz"
                if os.path.exists(I1):
                    try:
                        os.symlink(os.path.abspath(I1), des3)
                    except FileExistsError as e:
                        print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show3}  {des_show3} ")

                ##I2
                I2= fastq1.replace("R1.fastq.gz", "I2.fastq.gz").replace("R1_001.fastq.gz", "I2_001.fastq.gz").replace("R1.fq.gz", "I2.fq.gz")
                des4= f'{outdir}/raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_I2_001.fastq.gz'
                fastq_show4 = I2.replace(f"{outdir}", "")
                des_show4= f"raw_data/{sample}/{sample}_fastqs/{sample}_S2_L00{j}_I2_001.fastq.gz"
                if os.path.exists(I2):
                    try:
                        os.symlink(os.path.abspath(I2), des4)
                    except FileExistsError as e:
                        print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show4}  {des_show4} ")

                j = j + 1

            samples_update.loc[sample, "fastq"] = f'raw_data/{sample}/{sample}_fastqs'
            library.loc[sample, 'library_id'] = sample
            library.loc[sample, 'molecule_h5'] = f"{sample}/outs/molecule_info.h5"
            library.loc[sample, 'cloupe_file'] = f"{sample}/outs/cloupe.cloupe"
            library.loc[sample, 'spatial_folder'] = f"{sample}/outs/spatial"

        else:
            print("\033[0;31m%s\033[0m" % f"Warning: fastq: {sample} 和报告中的样本名不一致，请核查")
    samples_update.to_csv(samplesfiles, encoding='utf-8', index=False)
    library.to_csv(os.path. dirname(samplesfiles)+"/"+"library.csv", encoding='utf-8', index=False)
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 软链接创建完毕\n\n")

def check_obs_path(obs_address):
    """
    检验obs路径是否可用
    """
    module('load', "obsutil/5.2.12")
    print_obs_text = subprocess.run(f'obsutil ls {obs_address}', shell=True, capture_output=True, text=True).stdout.strip()
    print_obs_text_lines = print_obs_text.split('\n')

    # 查找包含 'Object list:' 的行，如果有，说明路径下有文件
    for line in print_obs_text_lines:
        if 'Object list:' in line:
            return True
    return False

def get_tasknum_obspath(config):
    """
    根据taskid更新config.yaml中的projectid和raw_data下载地址等
    """
    with open(config) as f:
        config_info = ordered_yaml_load(f)

    task_num = config_info["report"]["Task_Num"]
    # 解析 Task_Num，获取前面部分
    task_num_parts = task_num.split('-')
    project_prefix = task_num_parts[0]
    config_info["report"]["Project_Num"] = project_prefix

    # 更新样本数、路径和时间等信息
    if  config_info["report"]["Project_Path"] =="":
        config_info["report"]["Project_Path"] = os.path.abspath('.')
    if  config_info["report"]["Project_Start_Time"] =="":
        config_info["report"]["Project_Start_Time"] = time.strftime("%Y_%m_%d")

    # 获取合同年份
    match = re.search(r'[A-Za-z](\d{4})', task_num)
    if match:
        year = match.group(1)  # 提取匹配到的四位数
    else:
        print("\033[0;31m%s\033[0m" % f"Warning: 未匹配到合同的年份，默认设置为当前年份，请核查")
        year = datetime.datetime.now().year

    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始校验obs路径")
    # 如果obs路径为空，则自动填充，如果不为空，说明是手动填入，不修改obs路径
    if not config_info["report"]["Rawdata_Path"]["raw_data_obs_address"]:
        # 更新 config_info 中的字段
        raw_data_obs_address = f"obs://scrna-cloud/samples/{year}/{task_num}/"
        image_obs_address = f"obs://scrna-cloud/samples/空间转录组HE染色图片/{task_num}/"
        # 如果路径正确则填充，不正确则报错
        if check_obs_path(raw_data_obs_address):
            config_info["report"]["Rawdata_Path"]["raw_data_obs_address"] = [f"{raw_data_obs_address}"]
        else:
            print(f"\033[0;31m%s\033[0m" % f"Error: raw_data {raw_data_obs_address} 路径错误，请核查")
        if check_obs_path(image_obs_address):
            config_info["report"]["Rawdata_Path"]["image_obs_address"] = [f"{image_obs_address}"]
        else:
            print(f"\033[0;31m%s\033[0m" % f"Error: image {image_obs_address} 路径错误，请核查")
    else:
        # 检验手动填入的路径是否正确
        exs_rawdata_path = config_info["report"]["Rawdata_Path"]["raw_data_obs_address"][0]
        exs_img_path = config_info["report"]["Rawdata_Path"]["image_obs_address"][0]
        if check_obs_path(exs_rawdata_path):
            print("\033[0;31m%s\033[0m" % f"您手动填入的rawdata路径经检验有效")
        else:
            print("\033[0;31m%s\033[0m" % f"Error: 您填入的rawdata路径错误，请核查: \n{exs_rawdata_path}")
        if check_obs_path(exs_img_path):
            print("\033[0;31m%s\033[0m" % f"您手动填入的image路径经检验有效")
        else:
            print("\033[0;31m%s\033[0m" % f"Error: 您填入的image路径错误，请核查: \n{exs_img_path}")
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> obs路径校验完毕\n\n")
    with open(config, 'w', encoding='utf-8') as fp:
        ordered_yaml_dump(config_info, fp, default_flow_style=False,encoding='utf-8',allow_unicode=True)

    return config_info


def capture_text_from_docx(raw_sampleid_list):
    '''
    从质检报告的docx文档获取需要的信息
    '''
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始读取质检报告")

    file_path = subprocess.run('ls ./raw_data/images/*docx*', shell=True, capture_output=True, text=True).stdout.strip()
    document = Document(file_path)

    filename = os.path.basename(file_path)

    # 从docx文档名获取销售姓名
    match = re.search(r'[(（]([^)）]*)[)）]', filename) # 使用正则表达式匹配括号内的内容
    if match:
        sales_name = match.group(1)
    else:
        print("无法从docx文件名中提取销售姓名，请自行添加")

    # 获取文档表格信息
    tables = document.tables

    # 获取第一个表格
    first_table = tables[0]

    #获取表格行数
    table_row=len(first_table.rows)
    # 第一个表格去重并保留原来的顺序
    all_rows = []
    for row in first_table.rows:
        row_cells = [cell.text for cell in row.cells]
        unique_row_cells = list(OrderedDict.fromkeys(row_cells))   #这里有bug
        all_rows.extend(unique_row_cells)

    # all_rows总量之前的字符串再去重并保留原来的顺序
    unique_rows_before_zl = []
    zl_index = next((index for index, item in enumerate(all_rows) if item.startswith('总量')), None)

    for i, row in enumerate(all_rows):
        if i < zl_index:
            if row not in unique_rows_before_zl:
                unique_rows_before_zl.append(row)

    # 重新拼接字符串
    all_rows = unique_rows_before_zl + all_rows[zl_index:]

    # 从项目名称获取文库类型
    library_index = all_rows.index('项目名称')
    library_str = all_rows[library_index + 1]
    if "CytAssist" in library_str:
        library = "cytassist"
    else:
        if "FFPE" in library_str:
            library = "ffpe"
        else:
            library = "fresh"

    # 获取客户姓名
    customer_index = all_rows.index('客户姓名')
    customer = all_rows[customer_index + 1]
    # 替换逗号和下划线等为连字符
    customer = customer.replace(',', '-').replace('_', '-').replace('，', '-').replace('、', '-').replace('/', '-')

    # 获取样品类型后面的一个字符串,物种
    sample_type_index = all_rows.index('样品类型')
    species = all_rows[sample_type_index + 1]

    # 获取样品数量后面的一个字符串，docx中的样品数量
    sample_num_index = all_rows.index('样品数量')
    docx_sample_num = int(all_rows[sample_num_index + 1])

    # 获取实际样本数量
    real_sample_num = len(raw_sampleid_list)

    # 获取样本名和对应的slide
    sampleid_slide_area = {}  # 定义存储结果的字典
    sampleid_index = {}
    # 提取raw_data路径下的样本名，进行格式化，和格式化后的质控报告里的样本名进行校验
    raw_sampleid_list = raw_sampleid_list
    formated_raw_sampleid_list = [formated_sampleids(item) for item in raw_sampleid_list]

    # 创建一个字典，"格式化后的sampleid"："raw_data里的sampleid"
    sampleid_dict = {}
    for raw_sampleid in raw_sampleid_list:
        formatted_sampleid = formated_sampleids(raw_sampleid)
        sampleid_dict[formatted_sampleid] = raw_sampleid

    # 判断实际样本数量和docx样本数量是否相等
    if docx_sample_num == real_sample_num:
        sample_num = real_sample_num
    else:
        print("\033[0;31m%s\033[0m" % f"Warning: 报告中样本数为{docx_sample_num}，rawdata文件夹中的实际样本数为{real_sample_num}，样本名以rawdata中的样本名为准，请核查")
        sample_num = real_sample_num

    # 提取质控报告里的sampleid和slide_area
    first_sampleid_row=table_row-docx_sample_num

    baogao_sampleid_title={}
    for index, item in enumerate(first_table.rows[first_sampleid_row-1].cells):
        baogao_sampleid_title[item.text]=index
    # 循环提取数据
    for i in range(first_sampleid_row,table_row):
        baogao_sampleid = first_table.rows[i].cells[baogao_sampleid_title['样本名称']].text
        slide_area = first_table.rows[i].cells[baogao_sampleid_title['捕获区ID']].text
        baogao_sampleid_index = first_table.rows[i].cells[baogao_sampleid_title['编号']].text
        print(baogao_sampleid_index+" : "+baogao_sampleid+" : "+slide_area)
        # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
        if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
            sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
            sampleid_slide_area[sampleid] = slide_area
            sampleid_index[baogao_sampleid_index] = sampleid
        else:
            print("\033[0;31m%s\033[0m" % f"Warning: 报告中的sampleid: '{baogao_sampleid}' 和raw_data中的sampleid: {raw_sampleid_list} 不一致，请核查")
    # 获取空间表达玻片SN的index
    sn_index = all_rows.index('空间表达玻片SN')

    # 获取空间表达玻片SN
    sn = all_rows[sn_index + 1]
    sampleid_sn = {}
    # 如果不存在括号，就是一个sn
    if re.findall(r"[(（]",sn) == []:
        for sampleid in sampleid_dict.values():
            if formated_sampleids(sampleid) in formated_raw_sampleid_list:
            # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
                sampleid_sn[sampleid] = sn
    else:
        for i in sn.split("\n"):
            sn_name, bar = re.split(r"[（(]",i)
            seq = re.split('\D| ', bar)
            try:
                while True:
                    seq.remove("")
            except ValueError:
                pass
            for j in seq:
                if j in sampleid_index.keys() and  formated_sampleids(sampleid_index[j]) in formated_raw_sampleid_list:
                    sampleid_sn[sampleid_index[j]]=sn_name.strip()
    result = {
        "sales_name": sales_name,
        "library": library,
        "customer": customer,
        "species": species,
        "sample_num": sample_num,
        "sampleid_sn": sampleid_sn,
        "sampleid_slide_area": sampleid_slide_area
    }
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 质检报告读取完毕\n\n")
    return result


def updata_config_yaml(config, config_info_dict):
    """
    修改config.yaml中的Customer，Sales，Species，Sample_Num，library等
    """
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始更新config.yaml")

    with open(config) as f:
        config_info = ordered_yaml_load(f)
    config_info["report"]["Customer"] = config_info_dict["customer"]
    config_info["report"]["Sales"] = config_info_dict["sales_name"]
    config_info["report"]["Species"] = config_info_dict["species"]
    config_info["report"]["Sample_Num"] = config_info_dict["sample_num"]
    config_info["params"]["library_type"] = config_info_dict["library"]

    # 填充人和小鼠参考基因组以及探针序列
    species = config_info_dict["species"].split('-')[0]
    if species == "人":
        print("\033[0;31m%s\033[0m" % f"如果是人的自建基因组，需在该脚本完成后，自行填写参考基因组的绝对路径，最后不加斜杠")
        config_info["database"]["reference"] = "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A"
        if config_info_dict["library"] == "ffpe":
            config_info["database"]["probe_set"] = "/public/dev_scRNA/software/spaceranger-2.0.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv"
        elif config_info_dict["library"] == "cytassist":
            config_info["database"]["probe_set"] = "/public/dev_scRNA/software/spaceranger-2.0.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
    elif species == "小鼠":
        print("\033[0;31m%s\033[0m" % f"如果是小鼠的自建基因组，需在该脚本完成后，自行填写参考基因组的绝对路径，最后不加斜杠")
        config_info["database"]["reference"] = "/data/database/cellranger-refdata/refdata-gex-mm10-2020-A"
        if config_info_dict["library"] == "ffpe":
            config_info["database"]["probe_set"] = "/public/dev_scRNA/software/spaceranger-2.0.0/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv"
        elif config_info_dict["library"] == "cytassist":
            config_info["database"]["probe_set"] = "/public/dev_scRNA/software/spaceranger-2.0.0/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv"
    else:
        print(f"\033[0;31m%s\033[0m" % f"该项目的物种为：{species}，您需要手动填写参考基因组，注意是绝对路径，最后不加斜杠")

    with open(config, 'w', encoding='utf-8') as fp:
        ordered_yaml_dump(config_info, fp, default_flow_style=False,encoding='utf-8',allow_unicode=True)
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> config.yaml更新完毕\n\n")


def updata_sample_csv(samples, config_info_dict, outdir):
    '''
        修改sample.csv中的内容
    '''
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始更新samples.csv")

    library = config_info_dict["library"]
    species = config_info_dict["species"]
    sampleid_sn = config_info_dict["sampleid_sn"]
    sampleid_slide_area = config_info_dict["sampleid_slide_area"]
    # 读取CSV文件
    samples_update = pd.DataFrame(columns=["sampleid","slide","slide_area","fastq","image","cytaimage","species","group","batchid","image_json_file"
])

    formated_sampleid_list = []  # 创建一个空列表来存储sampleid
    sampleid_dict = {}    # 创建一个字典，"格式化后的sampleid"："raw_data里的sampleid"
    for sampleid in sampleid_sn:
        samples_update.loc[sampleid, "sampleid"] = sampleid
        samples_update.loc[sampleid, "slide"] = sampleid_sn[sampleid]
        samples_update.loc[sampleid, "slide_area"] = sampleid_slide_area[sampleid]
        samples_update.loc[sampleid, "species"] = species
        formatted_sampleid = formated_sampleids(sampleid)
        formated_sampleid_list.append(formatted_sampleid)  # 将sampleid添加到列表中
        sampleid_dict[formatted_sampleid] = sampleid    # "格式化后的sampleid"："raw_data里的sampleid"
        slide = sampleid_sn[sampleid]

        # 下载gpr文件
        outfile = f'{outdir}/raw_data/{slide}.gpr'
        slide_group = str(slide)[0:6]
        if not os.path.exists(outfile):
            module('load', "axel")
            subprocess.run(
                f"axel -n 2  -o {outfile}  https://s3-us-west-2.amazonaws.com/10x.spatial-slides/gpr/{slide_group}/{slide}.gpr",
                shell=True, check=True)
        else:
            print("\033[0;31m%s\033[0m" % f"{outfile} 已存在")

    # 读取images文件夹的内容
    image_path = "raw_data/images"
    images = [file for file in os.listdir(f"{outdir}/raw_data/images") if file.endswith(".tif")]
    if library == "cytassist":
        for image in images:
            filename, suffix = os.path.splitext(image)
            if "-扫描仪原图" in image:
                images_sampleid = image.replace("-扫描仪原图.jpg", "").replace("-扫描仪原图.tif", "").replace("-扫描仪原图.png", "")
                if(images_sampleid[0].isdigit()):
                    images_sampleid="S"+images_sampleid
                formated_images_sampleid = formated_sampleids(images_sampleid)
                raw_data_sampleid = sampleid_dict[formated_images_sampleid]
                # 如果images里的格式化后的sampleid和raw_data里格式化后的不一样，则报错
                if formated_images_sampleid not in formated_sampleid_list:
                    print("\033[0;31m%s\033[0m" % f"Warning: image:{image} 和报告中的样本名不一致，请核查")
                else:
                    # 如果sampleid和raw_data的sampleid相等，则使用该image名字，如果不同则修改成和raw_data一样的sapleid名
                    if images_sampleid == raw_data_sampleid:
                        samples_update.loc[raw_data_sampleid, "image"] = f'raw_data/images/{image}'
                    else:
                        # 重命名文件
                        os.rename(f"{image_path}/{image}", f"{image_path}/{raw_data_sampleid}-扫描仪原图{suffix}")
                        samples_update.loc[raw_data_sampleid, "image"] = f'{image_path}/{raw_data_sampleid}-扫描仪原图{suffix}'
            else:
                images_sampleid = image.replace(".jpg", "").replace(".tif", "").replace(".png", "")
                if(images_sampleid[0].isdigit()):
                    images_sampleid="S"+images_sampleid
                formated_images_sampleid = formated_sampleids(images_sampleid)
                raw_data_sampleid = sampleid_dict[formated_images_sampleid]
                # 如果images里的格式化后的sampleid和raw_data里格式化后的不一样，则报错
                if formated_images_sampleid not in formated_sampleid_list:
                    print("\033[0;31m%s\033[0m" % f"Warning: cytaimage image:{image} 和报告中的样本名不一致，请核查")
                else:
                    # 如果sampleid和raw_data的sampleid相等，则使用该image名字，如果不同则修改成和raw_data一样的sapleid名
                    if images_sampleid == raw_data_sampleid:
                        samples_update.loc[raw_data_sampleid, "cytaimage"] = f'raw_data/images/{image}'
                    else:
                        # 重命名文件
                        os.rename(f"{image_path}/{image}", f"{image_path}/{raw_data_sampleid}{suffix}")
                        samples_update.loc[raw_data_sampleid, "cytaimage"] = f'{image_path}/{raw_data_sampleid}{suffix}'
    else:
        for image in images:
            filename, suffix = os.path.splitext(image)
            if suffix in [".jpg", ".tif", ".png"]:
                images_sampleid = image.replace(".jpg", "").replace(".tif", "").replace(".png", "")
                if(images_sampleid[0].isdigit()):
                    images_sampleid="S"+images_sampleid
                formated_images_sampleid = formated_sampleids(images_sampleid)
                raw_data_sampleid = sampleid_dict[formated_images_sampleid]
                # 如果images里的格式化后的sampleid和raw_data里格式化后的不一样，则报错
                if formated_images_sampleid not in formated_sampleid_list:
                    print("\033[0;31m%s\033[0m" % f"Warning: cytaimage image:{image} 和报告中的样本名不一致，请核查")
                else:
                    # 如果sampleid和raw_data的sampleid相等，则使用该image名字，如果不同则修改成和raw_data一样的sapleid名
                    if images_sampleid == raw_data_sampleid:
                        samples_update.loc[raw_data_sampleid, "image"] = f'raw_data/images/{image}'
                    else:
                        # 重命名文件
                        os.rename(f"{image_path}/{image}", f"{image_path}/{raw_data_sampleid}{suffix}")
                        samples_update.loc[raw_data_sampleid, "image"] = f'{image_path}/{raw_data_sampleid}{suffix}'
    # 保存
    samples_update.to_csv(samples, encoding='utf-8', index=False)
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> samples.csv更新完毕\n\n")


@click.command("Data preparing for cellranger/spaceranger pipeline")
@click.option('-c', '--config',
              default="config/config.yaml",
              help='config yaml files with raw_data_obs_address/image_obs_address information.default:config/config.yaml')
@click.option('-s', '--samples',
              default="config/samples.csv",
              help='sample metadata information csv file.default:config/samples.csv')
@click.option('-o', '--outdir',
              default="./",
              help='Output directory,default:.')

def data_prepare(config, samples, outdir):
    """
    主函数
    """
    # 1.先根据config.yaml文件中的Task_Num修改config.yaml的projectid和下载地址。
    task_num_info = get_tasknum_obspath(config)

    # 2.获取原始数据和切片图片的下载地址等。
    raw_data_obs_address = task_num_info["report"]["Rawdata_Path"]["raw_data_obs_address"]
    image_obs_address = task_num_info["report"]["Rawdata_Path"]["image_obs_address"]

    # 3.执行下载
    raw_sampleid_list = raw_data_download(raw_data_obs_address, image_obs_address, outdir)

    # 4.根据image里的docx文档，再次更新config.yaml的其他内容
    config_info_dict = capture_text_from_docx(raw_sampleid_list)
    updata_config_yaml(config, config_info_dict)

    # 4.修改sample.csv文件，生成sampleid,slid,slide_area等
    updata_sample_csv(samples, config_info_dict, outdir)

    # 5.生成fastq软链接
    raw_data_formart(samples, outdir)
    print("\033[0;35m%s\033[0m" % f">>>>>>>>>>>>>>> data_prepare完成")

if __name__ == "__main__":
    data_prepare()
