#!/usr/bin/env python3
# encoding: utf-8
import yaml,re,click
import pandas as pd

def update_local_file( sample ,cellranger_config, status = "wait"):
    df_local = pd.read_csv(cellranger_config, sep=",")
    df_local.loc[df_local['Sample'].isin(sample), 'status'] = status
    df_local.to_csv(cellranger_config, index=False)
    if status == "running":
        print(f"更新状态成功,样本 "+ str(sample) +" 进行running。" , flush=True)
    if status == "finish":
        print(f"更新状态成功,样本 "+ str(sample) + " 完成分析", flush=True)

def remove_dup_items(d):
    """
    去除字典对象内的重复元素。
    """
    new_d = {}
    for k, v in d.items():
        new_d[k] = list(set(v))
    return new_d

def replace_chars(s):
    """
    如果字符串中包含要替换的reference字典中的任何一个key，则将该字符替换为匹配的key。
    """
    for key in reference.keys():
        if key in s:
            s = s.replace(s, key)
    return s

def extract_digits(project):
    """
    从字符串中提取数字的前四位。
    """
    match = re.search(r'\d{4}', project)
    if match:
        return match.group()

reference = {"人": "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A", 
             "小鼠":"/data/database/cellranger-refdata/refdata-gex-mm10-2020-A",
             "大鼠":"/data/database/cellranger-refdata/refdata-cellranger-Rnor_6.0"}

@click.command()
@click.option('-i', 'cellranger_config', type=str,default = 'cellranger_config.csv', 
    help='传入的样本产出表')
@click.option('-c', 'configfile', type=str, default = 'config/config.yaml',
    help='config yaml files with raw_data_obs_address, samples information. df:config/config.yaml')

def yaml_prepare(cellranger_config,configfile):
    """example: python yaml_prepare.py -c config/config.yaml -i cellranger_config.csv """
    # 读取表格数据
    df = pd.read_csv( cellranger_config, sep=",")
    if 'wait' in df['status'].values.tolist():
        df = df.loc[df['status'].isin(["wait"])] 
    else:
        print("The status with no sample is wait, exit! ")
        exit()
    # 删除所有包含 NaN 值的列
    df = df.dropna(axis=1)
    df["物种"] = [replace_chars(sp) for sp in df["物种"]]
    # 将 DataFrame 转换为字典
    data_dict = df.to_dict(orient='list')
    new_dict = remove_dup_items(data_dict)
    # 加载 config.yaml
    with open(configfile,'r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
    # 填写 config.yaml
    config["report"]["Project_Num"] = new_dict['任务单号2.0'][0]
    config["report"]["Species"] = new_dict['物种'][0]
    config["report"]["Sample_Num"] = len(new_dict['Sample'])
    config["report"]["Reference"] = reference[new_dict['物种'][0]]
    config["cellranger_params"]["raw_data_obs_address"] = [f"obs://scrna-cloud/samples/"+ extract_digits(new_dict['任务单号2.0'][0]) + f"/" + new_dict['任务单号2.0'][0] + f"/"]
    config["cellranger_params"]["samples"] = new_dict['Sample']
    # 输出 config.yaml
    with open(configfile,'w',encoding='utf-8') as yw:
        yaml.dump(config,yw)
    update_local_file( new_dict['Sample'], cellranger_config,"running")
    
if __name__ == '__main__':
    yaml_prepare()
