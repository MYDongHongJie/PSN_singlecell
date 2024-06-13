import requests
import pandas as pd
import json
from configparser import ConfigParser
import click
@click.command()
@click.option('-i','--projectnum', help='the project number. ')

def main(projectnum): 
    data = oeconfirm_api(projectnum)
    with open("%s.json" %projectnum,'w') as original_json:
        json.dump(data, original_json,ensure_ascii=False, indent=4)
    
    # config.ini
    project=pd.DataFrame(data['分析基本信息'])
    project.index=list(project['key'])
    config=ConfigParser()
    config.add_section("par")
    config.set("par","项目编号", project["value"]["项目编号"])
    config.set("par","客户姓名", project["value"]["客户名称"])
    config.set("par","实验物种", project["value"]["样本物种"])
    config.set("par","样本数目", str(len(data['样本信息']['样本分析名称'])))
    config.set("par","执行编号", "")
    config.add_section("db")
    config.set("db","genome", "")
    config.add_section("parameter")
    config.set("parameter","nGene_nUMI_fold2sd", str(2))
    config.set("parameter","percent_mito", "")
    config.set("parameter","marker_test", "bimod")
    config.set("parameter","diffexp", "MAST")
    config.write(open("config.ini", "w"))
    # metadata.csv
    samples=pd.DataFrame((data['样本信息']))
    if not all(samples['样本下机名称'] == samples['样本分析名称'] ):
        samples.iloc[:,0:2].to_csv("change_name.tsv",sep='\t',index=False,header=False)
        print("请注意：样本分析名有所变化，下机名与分析名对应关系请见 change_name.tsv")

    samples.index=list(samples['样本分析名称'])
    samples = samples.rename(columns={'批次':'batchid','样本分析名称':'sampleid','样本下机名称':'rawsampleid'}) 
    samples['group']=''
    samples['species']=project["value"]["样本物种"]
    order = ['sampleid','species', 'group','batchid']
    samples=samples[order]
    
    groups=pd.DataFrame((data['样本分组信息']))
    groups.index=list(groups['group'])
    groups1=groups[groups["type"]=="常规"]
    for line in groups1.index: 
        for  sample in groups1.loc[line,"samples"].split(","):
            samples.loc[sample,'group'] = groups1.loc[line,"group"]
    # other groups
    if any(groups["type"]!="常规"):
        groups2=groups[groups["type"]!="常规"]
        sample_groups=dict.fromkeys(samples.index)
        for line in groups2.index:
            for  sample in groups2.loc[line,"samples"].split(","):
                if  sample_groups[sample] == None :
                    sample_groups[sample] =  [groups2.loc[line,"group"]]
                else:
                    sample_groups[sample].append(groups2.loc[line,"group"])
    
        sample_groups_df= pd.DataFrame.from_dict(sample_groups,orient="index",columns=[f"group{i}" for i in range(2, max([len(i) for i in sample_groups.values() ])+2)])
        samples=pd.concat([samples,sample_groups_df],axis=1)

    samples.to_csv("metadata.csv",sep=',',index=False) #   na_rep='NA'
    # diff_group.tsv
    diff=pd.DataFrame((data['差异比较信息']))
    diff.insert(0,'type','sampleid')
    diff.insert(2,'case_sample','')
    diff.insert(4,'control_sample','')
    
    for line in diff.index:
        diff.loc[line,'case_sample'] = groups.loc[diff.loc[line,'case'],'samples']
        diff.loc[line,'control_sample'] = groups.loc[diff.loc[line,'control'],'samples']
        for i in samples.columns:
            if diff.loc[line,'case'] in samples.loc[:,i].tolist() and diff.loc[line,'control'] in samples.loc[:,i].tolist():
                diff.loc[line,'type'] = i
                continue
    diff.to_csv("diffexp.tsv",sep='\t',index=False)

def oeconfirm_api(analysis_id):
    # 在线确认单秘钥并不是不变的，可能会一个月更改一次
    oeconfirm_token = "oebiotech.confrim.T34u"

    # 在线确认单api链接固定
    confirm_api_url = "https://cloud.oebiotech.cn/oeconfirm/get_record"

     # 请求数据
     # analysis_id：分析编号
    # token: 在线确认单秘钥
    post_data = {'analysis_id': analysis_id, 'token': oeconfirm_token}
    response = requests.post(confirm_api_url, data=post_data, timeout=5)
    result_json = None
    if response.status_code == 200:
        result_json = response.json()
    elif response.status_code == 500:
        print("未查询到记录信息")
        exit()
    else:
        print("访问出错")
        exit()
    return result_json

if __name__ == "__main__":
    main()
