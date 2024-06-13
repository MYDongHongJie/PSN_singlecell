#!/usr/bin/env python3
# encoding: utf-8
import yaml,re,click,os,time,shutil
import pandas as pd

def update_local_file( sample ,config_table, column = "Sample", status = "wait"):
    df_local = pd.read_csv(config_table, sep=",")
    df_local.loc[df_local[column].isin(sample), 'status'] = status
    df_local.to_csv(config_table, index=False)
    if status == "running":
        print(f"更新状态成功,"+ str(column)+ " " + str(sample) +" 进行running。" , flush=True)
    if status == "finish":
        print(f"更新状态成功,"+ str(column)+ " " + str(sample) + " 完成分析", flush=True)

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
    for key in Reads.keys():
        if key in s:
            s = s.replace(s, key)
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
             "小鼠":"/data/database/cellranger-refdata/refdata-gex-mm10-2020-A"}

vdj_Reference = {"人": "/data/database/cellranger-refdata/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0", 
             "小鼠":"/data/database/cellranger-refdata/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0"}

Reads = {"100G": "350M", "150G":"500M", "10G":"35M"}

def yaml_prepare_cellranger(config_table,configfile):
    # 读取表格数据
    df = pd.read_csv( config_table, sep=",")
    # 数字开头的样本添加 S
    df['Sample'] = ["S"+sp if re.match(r"^[0-9]", sp) else sp for sp in df['Sample']]
    sample_sc = list(set([s for s in df['Sample'] if not (s.endswith("_TCR") or s.endswith("_BCR"))]))
    df.to_csv(config_table, index=False)
    # 挑选status为wait的样本进行分析
    if 'wait' in df['status'].values.tolist() or 'error' in df['status'].values.tolist():
        df = df.loc[df['status'].isin(["wait","error"])] 
    else:
        print("The status with no sample is wait, exit! ")
        exit(1)
    # 删除所有包含 NaN 值的列
    df = df.dropna(axis=1)
    sample_test = df[(df["目标数据量"] < 20 )  & ~(df["Sample"].str.endswith("_TCR") | df["Sample"].str.endswith("_BCR")) ]["Sample"].tolist()
    # df["物种"] = [replace_chars(sp) for sp in df["物种"]]
    # 将 DataFrame 转换为字典
    data_dict = df.to_dict(orient='list')
    new_dict = remove_dup_items(data_dict)
    # 加载 config.yaml
    with open(configfile,'r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
    # 填写 config.yaml
    config["report"]["Project_Num"] = new_dict['任务单号2.0'][0]
    config["report"]["Species"] = new_dict['物种'][0]
    config["report"]["Sample_Num"] = len(sample_sc)
    new_dict["物种"] = [replace_chars(sp) for sp in new_dict["物种"]]
    config["report"]["scrna_Reference"] = reference[new_dict['物种'][0]]
    config["report"]["vdj_Reference"] = vdj_Reference[new_dict['物种'][0]]
    config["cellranger_params"]["raw_data_obs_address"] = [f"obs://scrna-cloud/samples/"+ extract_digits(new_dict['任务单号2.0'][0]) + f"/" + new_dict['任务单号2.0'][0] + f"/"]
    config["cellranger_params"]["samples"] = new_dict['Sample']
    config["cellranger_params"]["test_samples"] = sample_test
    # 输出 config.yaml
    with open(configfile,'w',encoding='utf-8') as yw:
        yaml.dump(config,yw,allow_unicode=True)
    update_local_file( new_dict['Sample'], config_table, "Sample","running")
    # 去除此前运行的重复样本结果
    for i in new_dict['Sample']:
        dir_path=f"result/cellranger/" + i
        if os.path.exists(dir_path):
            print(f"样本"+ str(i) +"在之前进行过cellranger，即将删除该样本cellranger结果并重新运行！")
            shutil.rmtree(dir_path)
        if os.path.exists(i + "_cellranger.check"):
            os.remove(i + "_cellranger.check")
    # 修改 cluster.yaml
    with open("config/cluster_txy.yaml",'r',encoding='utf-8') as y:
        cluster_txy=yaml.load(y,yaml.Loader)
    # cluster_txy['__default__']['custom_prop'] = time.strftime("%Y-%m-%d-%H-%M-")+ new_dict['任务单号2.0'][0] + "-" + "-".join(new_dict['Sample'])
    cluster_txy['__default__']['custom_prop'] = new_dict['任务单号2.0'][0].lower() + f"-auto" #+ "-" + "-".join([s.lower() for s in new_dict['Sample']])
    with open("config/cluster_txy.yaml",'w',encoding='utf-8') as yw:
            yaml.dump(cluster_txy,yw)
####
def yaml_prepare_email(config_table,configfile):
    # 读取表格数据
    df = pd.read_csv( config_table, sep=",")
    # 挑选status为wait的任务单号
    if 'ready' in df['status'].values.tolist() or 'error' in df['status'].values.tolist():
        df = df.loc[df['status'].isin(["ready","error"])] 
    else:
        print("没有需要生成邮件的任务单号! ")
        exit(1)
    # 删除所有包含 NaN 值的列
    df = df.dropna(axis=1)
    # 将 DataFrame 转换为字典
    data_dict = df.to_dict(orient='list')
    new_dict = remove_dup_items(data_dict)
    new_dict['report'][0] = "F" if new_dict['reads'][0] == "10G" else new_dict['report'][0]
    new_dict["reads"] = [replace_chars(sp) for sp in new_dict["reads"]]
    df_cellranger = pd.read_csv( "cellranger_config.csv", sep=",")
    if new_dict['rwdh'][0] not in df_cellranger['任务单号2.0'].tolist():
        print("项目任务单号与本目录项目不符，请核查！")
        exit(1)
    #df_cellranger = df_cellranger.loc[df_cellranger['status'].isin(["finish"])] 
    #sample_list = list(set(df_cellranger['Sample']))
    # 加载 config.yaml
    with open(configfile,'r',encoding='utf-8') as y:
        config=yaml.load(y,yaml.Loader)
    # if new_dict['rwdh'][0] == config["report"]["Project_Num"]:
    # 填写 config.yaml
    config["report"]["Executor"] = "OE0988"
    config["report"]["Customer"] = new_dict['khxm'][0] + f"_" + new_dict['lxrxm'][0] if new_dict['khxm'][0] != new_dict['lxrxm'][0] else new_dict['khxm'][0]
    config["report"]["Project_Num"] = new_dict['rwdh'][0] 
    config["report"]["Task_Num"] = new_dict['rwdh'][0] 
    config["report"]["Sales"] = new_dict['lastname'][0]
    #config["report"]["Sample_Num"] = len(sample_list)
    config['email'] = {}
    config["email"]["emailto"] = ['qingzhen.zi@oebiotech.com','kaiqi.guo@oebiotech.com']
    config["email"]["remark"] = new_dict['fxbz'][0] if "fxbz" in new_dict else ""
    config["email"]["Project_management"] = new_dict['xmzy'][0]
    config["email"]["Laboratory"] = f"王翠亭博士"
    config["cellranger_params"]["module"]["cellranger_report"] = new_dict['report'][0] == "T"
    config["cellranger_params"]["module"]["upload_bam"] = True
    config["cellranger_params"]["Number_of_Reads_for_mutil"] = Reads[new_dict["reads"][0]] + f',35M' #+ Reads[str(int(new_dict['目标数据量'][1])) + f"G"]
    #config["cellranger_params"]["samples"] = sample_list
    # else:
        # print("项目任务单号与本目录项目不符，请核查！")
        # exit()
    # 输出 config.yaml
    with open(configfile,'w',encoding='utf-8') as yw:
        yaml.dump(config,yw,allow_unicode=True)
    update_local_file( new_dict['rwdh'], config_table, "rwdh", "running")
    # 修改 cluster.yaml
    with open("config/cluster_txy.yaml",'r',encoding='utf-8') as y:
        cluster_txy=yaml.load(y,yaml.Loader)
    cluster_txy['__default__']['custom_prop'] = new_dict['rwdh'][0].lower() + f"-auto-email" 
    with open("config/cluster_txy.yaml",'w',encoding='utf-8') as yw:
            yaml.dump(cluster_txy,yw)

#######################
@click.command()
@click.option('-i', 'config_table', type=str,default = 'cellranger_config.csv', 
    help='传入的样本产出表/OA信息表。')
@click.option('-e', 'email_config_table', type=str,default = 'email_config.csv', 
    help='传入的样本产出表/OA信息表。')
@click.option('-c', 'configfile', type=str, default = 'config/config.yaml',
    help='config yaml files with raw_data_obs_address, samples information. df:config/config.yaml')
#######################
def run(config_table, email_config_table, configfile):
    """example: /data/software/conda_envs/snakemake/bin/python scripts/yaml_prepare.py -c config/config.yaml -i cellranger_config.csv -e email_config.csv """
    yaml_prepare_cellranger(config_table, configfile)
    yaml_prepare_email(email_config_table, configfile)

if __name__ == '__main__':
    run()