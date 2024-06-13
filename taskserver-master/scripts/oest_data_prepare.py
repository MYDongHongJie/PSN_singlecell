#!/public/scRNA_works/works/guochy/conda/miniconda3/envs/taskserver_2.0/bin/python
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
import os, sys, re, string,datetime
import time
import click
import subprocess
import pandas as pd
from glob import glob
from collections import OrderedDict
from docx import Document
import yaml


exec(open('/data/software/modules/modules-v4.2.1/init/python.py').read())#this may cause problem if we turn to docker image
os.environ['MODULEPATH'] =os.path.abspath("./envs/")


##===========================================Step1: download all input files including image===========================================
def raw_data_download(fastqs, outdir):
    suffix = ["jpg", "tif", "tiff","png"]
    for batch, fastq in fastqs.items():
        path = None
        for file in fastq:
            print(file)
            split_path = file.split("/")
            file_name = split_path.pop(-1)
            for suf in suffix:
                if suf in file:
                    path = "/".join(split_path)           
            module('load', "obsutil/5.2.12")
            subprocess.run(
                f"obsutil cp  {file}  {outdir}/raw_data/download_from_obs/{file_name}  -f -flat -vlength -vmd5",
                shell=True, check=True)
        if path:
            module('load', "obsutil/5.2.12")
            subprocess.run(
                f"obsutil cp  {path}  {outdir}/raw_data/download_from_obs/  -flat -r -include=*.docx  -f -vlength",
                shell=True, check=True)

##==========================================Step2: format fastq names=================================================================
def raw_fastq_formart(sample, outdir):
        if not os.path.exists(f'{outdir}/raw_data/{sample}_fastqs'):
            os.mkdir(f'{outdir}/raw_data/{sample}_fastqs')
            fq1_list = [item for sublist in
                        [glob(f'{outdir}/raw_data/download_from_obs' + ext) for ext in [f"/*R1.fastq.gz", "/*_1.fq.gz","/*R1_001.fastq.gz","/*R1.fq.gz"]]
                        for item in sublist]
            fq1_list.sort(key=os.path.getmtime)
            j = 1
            for fastq1 in fq1_list:
                ##R1
                des1 = f'{outdir}/raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R1_001.fastq.gz'
                fastq_show1 = fastq1.replace(f"{outdir}", "")
                des_show1 = f"raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R1_001.fastq.gz"
                print(des1)
                print(des_show1)
                if not os.path.exists(des1):
                    os.symlink(os.path.abspath(fastq1), des1)
                    print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show1}  {des_show1} ")
                else:
                    print("\033[0;31m%s\033[0m" % f'Warning: {des_show1} has existed!!!')

                ##R2
                fastq2 = fastq1.replace("R1.fastq.gz", "R2.fastq.gz").replace("_1.fq.gz", "_2.fq.gz").replace("R1_001.fastq.gz", "R2_001.fastq.gz").replace("R1.fq.gz", "R2.fq.gz")
                des2 = f'{outdir}/raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R2_001.fastq.gz'
                fastq_show2 = fastq2.replace(f"{outdir}", "")
                des_show2 = f"raw_data/{sample}_fastqs/{sample}_S2_L00{j}_R2_001.fastq.gz"
                if not os.path.exists(des2):
                    os.symlink(os.path.abspath(fastq2), des2)
                    print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show2}  {des_show2} ")
                else:
                    print("\033[0;31m%s\033[0m" % f'Warning: {des_show2} has existed!!!')

                ##I1
                I1= fastq1.replace("R1.fastq.gz", "I1.fastq.gz").replace("R1_001.fastq.gz", "I1_001.fastq.gz")
                des3= f'{outdir}/raw_data/{sample}_fastqs/{sample}_S2_L00{j}_I1_001.fastq.gz'
                fastq_show3 = I1.replace(f"{outdir}", "")
                des_show3 = f"raw_data/{sample}_fastqs/{sample}_S2_L00{j}_I1_001.fastq.gz"
                if   os.path.exists(I1):
                    os.symlink(os.path.abspath(I1), des3)
                    print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show3}  {des_show3} ")
                ##I2
                I2= fastq1.replace("R1.fastq.gz", "I2.fastq.gz").replace("R1_001.fastq.gz", "I2_001.fastq.gz")
                des4= f'{outdir}/raw_data/{sample}_fastqs/{sample}_S2_L00{j}_I2_001.fastq.gz'
                fastq_show4 = I2.replace(f"{outdir}", "")
                des_show4= f"raw_data/{sample}_fastqs/{sample}_S2_L00{j}_I2_001.fastq.gz"
                if   os.path.exists(I2):
                    os.symlink(os.path.abspath(I2), des4)
                    print("\033[0;32m%s\033[0m" % f"Making symlink: {sample} {j} {fastq_show4}  {des_show4} ")
                j = j + 1

##==========================================Step3: format image names=================================================================
def raw_image_formart(sample, outdir):
    if not os.path.exists(f'{outdir}/raw_data/images'):
        os.mkdir(f'{outdir}/raw_data/images')
        image_list = [item for sublist in
                    [glob(f'{outdir}/raw_data/download_from_obs' + ext) for ext in [f"/*.jpg", "/*.tif", "/*.tiff","/*.png"]]
                    for item in sublist]
        image_list.sort(key=os.path.getmtime)
        j = 1
        for image in image_list:
            image_name = image.split("/")
            des1 = f'{outdir}/raw_data/images/{image_name[-1]}'
            image_show1 = image.replace(f"{outdir}", "")
            des_show1 = f"raw_data/images/{image_name[-1]}"
            print(des1)
            print(des_show1)
            if not os.path.exists(des1):
                os.system(f"cp {image} {des1}")
                print("\033[0;32m%s\033[0m" % f"Copying image: {sample} {image_show1}  {des_show1} ")
            else:
                print("\033[0;31m%s\033[0m" % f'Warning: {des_show1} has existed!!!')
            j = j + 1
##============================================Step4: generate sample.csv===================================================================
def formated_sampleids(text):
    # 创建一个转换表，将标点符号和空格都映射为空字符
    translator = str.maketrans('', '', string.punctuation + ' ')
    # 使用转换表去除标点符号和空格
    cleaned_text = text.translate(translator)
    return cleaned_text
def capture_text_from_docx(raw_sampleid_list,outdir):
    '''
    从质检报告的docx文档获取需要的信息
    '''
    print("\033[0;36m%s\033[0m" % f">>>>>>>>>>>>>>> 开始读取质检报告")
    file_path = subprocess.run(f'ls {outdir}/raw_data/download_from_obs/*docx*', shell=True, capture_output=True, text=True).stdout.strip()
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
    # 第一个表格去重并保留原来的顺序
    all_rows = []
    for row in first_table.rows:
        row_cells = [cell.text for cell in row.cells]
        unique_row_cells = list(OrderedDict.fromkeys(row_cells))
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
     # 提取质控报告里的sampleid和slide_area，分为表格里包含接头和不包含接头的。
    if "接头" not in all_rows:
        ref_index_zl = next((index for index, item in enumerate(all_rows) if item.startswith('总量')), None)
        # 循环提取数据
        for i in range(docx_sample_num):
            baogao_sampleid = all_rows[ref_index_zl + 3 + i * 9]
            slide_area = all_rows[ref_index_zl + 2 + i * 9]
            # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
            if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
                sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
                sampleid_slide_area[sampleid] = slide_area
                print("\033[0;31m%s\033[0m" % f"报告中的sampleid: '{baogao_sampleid}' 和raw_data中的sampleid: {raw_sampleid_list} 一致，将进行后续分析！")
            else:
                print("\033[0;31m%s\033[0m" % f"Warning: 报告中的sampleid: '{baogao_sampleid}' 和raw_data中的sampleid: {raw_sampleid_list} 不一致，可能此次分析不涉及这个样本，请注意！")
    else:
        ref_index_jt = next((index for index, item in enumerate(all_rows) if item.startswith('接头')), None)
        # 循环提取数据
        for i in range(docx_sample_num):
            baogao_sampleid = all_rows[ref_index_jt + 3 + i * 8]
            slide_area = all_rows[ref_index_jt + 2 + i * 8]
            # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
            if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
                sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
                sampleid_slide_area[sampleid] = slide_area
                print("\033[0;31m%s\033[0m" % f"报告中的sampleid: '{baogao_sampleid}' 和raw_data中的sampleid: {raw_sampleid_list} 一致，将进行后续分析！")
            else:
                print("\033[0;31m%s\033[0m" % f"Warning: 报告中的sampleid: '{baogao_sampleid}' 和raw_data中的sampleid: {raw_sampleid_list} 不一致，可能此次分析不涉及这个样本，请注意！")
    # 获取空间表达玻片SN的index
    sn_index = all_rows.index('空间表达玻片SN')
    # 获取空间表达玻片SN
    sn = all_rows[sn_index + 1]
    sampleid_sn = {}
    res = []
    sn_list = []
    # 如果不存在括号，就是一个sn
    if sn.find("（") == -1:
        if "接头" not in all_rows:
            ref_index_zl = next((index for index, item in enumerate(all_rows) if item.startswith('总量')), None)
            # 循环提取数据
            for i in range(docx_sample_num):
                baogao_sampleid = all_rows[ref_index_zl + 3 + i * 9]
                # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
                if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
                    sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
                    sampleid_sn[sampleid] = sn
                else:
                    pass
        else:
            ref_index_jt = next((index for index, item in enumerate(all_rows) if item.startswith('接头')), None)
            # 循环提取数据
            for i in range(docx_sample_num):
                baogao_sampleid = all_rows[ref_index_jt + 3 + i * 8]
                # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
                if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
                    sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
                    sampleid_sn[sampleid] = sn
                else:
                    pass
    else:
        for i in sn.split("\n"):
            sn_name, bar = i.split("（")
            seq = re.split(r'[-,，]', bar.split("）")[0])
            if len(seq) == 1:
                res.append(sn_name)
            else:
                res.append((sn_name + ",") * (int(seq[1]) - int(seq[0]) + 1))
        split_list = [item.split(',') for item in res]
        for sublist in split_list:
            for element in sublist:
                if element != '':
                    sn_list.append(element)
        if "接头" not in all_rows:
            ref_index_zl = next((index for index, item in enumerate(all_rows) if item.startswith('总量')), None)
            # 循环提取数据
            for j in range(docx_sample_num):
                baogao_sampleid = all_rows[ref_index_zl + 3 + j * 9]
                # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
                if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
                    sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
                    sampleid_sn[sampleid] = sn_list[j]
                else:
                    pass
        else:
            ref_index_jt = next((index for index, item in enumerate(all_rows) if item.startswith('接头')), None)
            # 循环提取数据
            for j in range(docx_sample_num):
                baogao_sampleid = all_rows[ref_index_jt + 3 + j * 8]
                # 如果格式化后的报告里的sampleid在格式化后的raw_data的sampleid列表里，则将sampleid命名为raw_data的sampleid名
                if formated_sampleids(baogao_sampleid) in formated_raw_sampleid_list:
                    sampleid = sampleid_dict[formated_sampleids(baogao_sampleid)]
                    sampleid_sn[sampleid] = sn_list[j]
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

##==========================================Step5: download  slide gpr file=================================================================
def slide_download(slide,outdir):
    outfile = f'{outdir}/raw_data/{slide}.gpr'
    slide_group = str(slide)[0:6]
    if not os.path.exists(outfile):
        module('load', "axel")
        subprocess.run(
            f"axel -n 2  -o {outfile}  https://s3-us-west-2.amazonaws.com/10x.spatial-slides/gpr/{slide_group}/{slide}.gpr",
            shell=True, check=True)
    else:
        print("\033[0;31m%s\033[0m" % f"Warning: {outfile} has existed!!!")
##==========================================args============================================================================================
@click.command("Data preparing for cellranger/spaceranger pipeline")
@click.option('-i', '--input',type=click.Path("r"),
              default="input/input.json",
              help='config yaml files with raw_data_obs_address/image_obs_address information.default:config/config.yaml')
@click.option('-o', '--outdir', type=click.Path(exists=False),
              default="./",
              help='Output directory,default:.')
def data_prepare(input, outdir):
    input_info = pd.read_json(f'{input}', orient="index")
    fastqs = input_info.loc['base', "sample"]['dataFiles']
    select_sample = input_info.loc["base","sample"]["name"]
    raw_data_download(fastqs, outdir)
    info_dict = capture_text_from_docx([select_sample],outdir)
    library_type = info_dict['library']
    if input_info.loc["parameters","image_json_file"] != "": #调整为去指定位置获取
        image_json = "/public/scRNA_works/works/guochy/ST_taskserver/json_files/" + input_info.loc["parameters","image_json_file"]
    else:
        image_json = None
    sample_info = pd.DataFrame(columns=["sampleid","slide","slide_area","fastq","image","cytaimage","species","group","batchid","image_json_file"]) #应该可以用left_join函数把根据字典生成的dataframe替代这个新的dataframe中指定的列
    sample_info.loc[select_sample,"sampleid"] = select_sample
    sample_info.loc[select_sample,"slide"] = info_dict['sampleid_sn'][select_sample]
    sample_info.loc[select_sample,"slide_area"] = info_dict['sampleid_slide_area'][select_sample]
    sample_info.loc[select_sample,"fastq"] = ""
    sample_info.loc[select_sample,"image"] = ""
    sample_info.loc[select_sample,"cytaimage"] = ""
    sample_info.loc[select_sample,"species"] = info_dict['species']
    sample_info.loc[select_sample,"group"] = ""
    sample_info.loc[select_sample,"batchid"] = ""
    if image_json:
        sample_info.loc[select_sample,"image_json_file"] = image_json
    else:
        sample_info.loc[select_sample,"image_json_file"] = ""
    raw_fastq_formart(select_sample, outdir)
    raw_image_formart(select_sample, outdir)
    slide = sample_info.loc[select_sample,"slide"]
    slide_download(slide,outdir) #this must be the last step
    fastq_dir = f'{outdir}/raw_data/{select_sample}_fastqs/'
    sample_info.loc[select_sample,"fastq"] = fastq_dir
    #slidefile = f'{outdir}/raw_data/{slide}.gpr'
    if library_type == "cytassist":
        image_list = [item for sublist in
                    [glob(f'{outdir}/raw_data/images' + ext) for ext in [f"/*.jpg", "/*.tif", "/*.tiff","/*.png"]]
                    for item in sublist]
        origin_list = [item for sublist in
                    [glob(f'{outdir}/raw_data/images' + ext) for ext in [f"/*原图.jpg", "/*原图.tif", "/*原图.tiff","/*原图.png"]]
                    for item in sublist]
        sample_info.loc[select_sample,"image"] = origin_list[0]
        extra_list = [x for x in image_list if x not in origin_list]
        sample_info.loc[select_sample,"cytaimage"] = extra_list[0]
    else:
        image_list = [item for sublist in
                    [glob(f'{outdir}/raw_data/images' + ext) for ext in [f"/*.jpg", "/*.tif", "/*.tiff","/*.png"]]
                    for item in sublist]
        sample_info.loc[select_sample,"image"] = image_list[0]

    ## write samples.csv
    sample_info.to_csv(f'{outdir}/config/samples.csv', index=False)

    ##update config yaml
    subprocess.run(f'cp -r /public/scRNA_works/works/guochy/ST_taskserver/ST_snakemake_template/st_pipeline/config/config.yaml {outdir}/config/',shell=True, check=True)
    with open(f'{outdir}/config/config.yaml','r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
        config['report']['customer'] = info_dict['customer']
        config['report']['Sales'] = info_dict['sales_name']
        config['report']['Species'] = info_dict['species']
        config['params']['library_type'] = info_dict['library']
        config['params']['Sample_Num'] = info_dict['sample_num']
    with open(f'{outdir}/config/config.yaml','w',encoding='utf-8') as yw:
        yaml.dump(config,yw,allow_unicode=True,encoding='utf-8',default_flow_style=False)


if __name__ == "__main__":
    data_prepare()
