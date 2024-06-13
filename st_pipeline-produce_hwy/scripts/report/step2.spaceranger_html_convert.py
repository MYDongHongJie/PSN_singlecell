#  -*- coding: UTF-8 -*
"""
Project :convert yaml to html
Author  :Xiufeng.Yang
Contact : xiufeng.yang@oebiotech.com
File   : report_oebio.py
IDE    : PyCharm
Time   : 2020-04-20
Desc   : produce 10x visium html report by using oebio
"""

import os
import yaml
import click
from oebio.report import Report
from oebio.app import *
from collections import OrderedDict

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

@click.command()
@click.option('-d', '--outdir', prompt='the report directory',
              help='the directory of your report files ')
@click.option('-c', '--configfile', prompt='config yaml file. [default:config/config.yaml]',
              default='config/config.yaml',
              help='the config file which contain your program information.')
def convert_html(outdir, configfile):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    # 项目摘要统计数据====================================================================================================
    with open(configfile) as f:
        config = ordered_yaml_load(f)
    # file = open(configfile, 'r', encoding='utf-8')
    os.chdir(outdir)
    # cfg = file.read()
    # config = yaml.load(cfg)
    library_type = config['params']['library_type']
    
    summ = pd.read_csv('../../../config/samples.csv', index_col=0, sep=',')
    ## 质控时切片数据和样本数量相同
    sample_num = len(summ)  # 样本数量
    slices_num = len(summ)  # 切片数量
    # Define header info
    header_info = {key: config["report"][key] for key in config["report"] if
                   key in ["Project_Num", "Customer", "Species","Task_Num"]}

    header_info["订单编号"] = header_info.pop("Project_Num")
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["物种信息"] = header_info.pop("Species")
    header_info["样本数量"] = sample_num
    if library_type == "fresh":
        header_info["文库类型"] = "Fresh"
    elif library_type == "ffpe":
        header_info["文库类型"] = "FFPE"
    elif library_type == "cytassist":
        header_info["文库类型"] = "CytAssist"
    header_info["任务单号"] = header_info.pop("Task_Num")
    title = "10x Visium空间转录组项目质控报告"
    report = Report(title, title=title, header_info=header_info, oe_welcome='''
    ####  感谢您选择欧易生物！
    
    **关于网页报告使用有以下几点提示：**
    
    1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节，便于您快速查看相应内容。
    2. **报告正文**中的**图片**均可以点击后进行**放大查看**，且图片点击后可以包含更多的细节信息，比如左上角会显示具体的**图片数量**，右上角有相关的**图片工具**，**键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
    3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张，更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
    4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**，当表格多于20行的时候，表格会嵌入到网页之中，可以利用滚动栏进行调整设置，同时会在表格上方增加搜索设置，便于快速查询表格信息。
    5. 在报告中浏览，需要返回顶部的时候，可以使用网页**右下角的白色箭头标识**快速返回顶部。
    6. 本提示可以点击右上方的不再显示，则后续打开网页均不会显示本提示。
    ''')
    if library_type == "fresh":
        report.add_yaml_config(os.path.join(scriptdir, 'visium_report.yaml'))
    elif library_type == "ffpe":
        report.add_yaml_config(os.path.join(scriptdir, 'visium_report_FFPE.yaml'))
    elif library_type == "cytassist":
        report.add_yaml_config(os.path.join(scriptdir, 'visium_report_cytassist.yaml'))

    os.environ['OEBIO'] = os.path.join(scriptdir,"pic")

    ################################################ 项目概况#################################################################
    summary = report.add_section('项目概况')

    summary_1 = summary.add_section('项目摘要_仅质控',
                                    sample_num=sample_num,
                                    slices_num=slices_num)

    ################################################ 分析流程 #################################################################
    experi_module = report.add_section('技术简介') 
    experi1 = experi_module.add_section('背景简介')
    experi2 = experi_module.add_section('技术原理')
    experi3 = experi_module.add_section('实验流程')

    ################################################ 项目分析结果 #############################################################
    result_module = report.add_section('项目质控结果')

    ################################################ SpaceRanger分析结果 #############################################################
    #section1 = result_module.add_section('SpaceRanger分析结果')
    sample_name = summ.index[0] 
    cloupe_file_link = ""
    HE_file = list(glob(f"../../spaceranger/{sample_name}/outs/spatial/*{sample_name}*"))[0]
    SpaceRanger_folder = list(glob(r"[0-9]*.SpaceRanger"))[0]
    section1 = result_module.add_section('基因定量质控',HE_format = os.path.splitext(HE_file)[1], cloupe_file_link = cloupe_file_link)
    #section1 = result_module.add_section('SpaceRanger分析')
    # samplefile = pd.read_csv('../../../config/samples.csv', index_col=0, sep=',')
    # sample_name = samplefile.index[0]

    ######################################### 附录 #############################################################################
    folder = list(glob(r"[0-9]*Supplemental_material"))[0]
    section8 = report.add_section('附录')
    section8_1 = section8.add_section("Loupe Browser使用教程",folder=folder )
    section8_2 = section8.add_section("实验技术方法说明",folder=folder)
    #section8_3 = section8.add_section('AI修图使用说明')
    section8_4 = section8.add_section("参考数据库链接",folder=folder)

    section9 = report.add_section('申明')
    report.write_to('质控说明.html',zip_report_name=f'{outdir}.zip')

if __name__ == "__main__":
    convert_html()
