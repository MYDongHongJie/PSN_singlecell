# encoding:utf-8
import shutil
import os
from oebio.report import Report
from oebio.app import *
import click
import subprocess
import sys
import time
import yaml
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

@click.command()
@click.option('-i', '--input', prompt='the program path',
              help='the directory of your program.[default:./]',default="./")
@click.option('-c', '--configfile', prompt='config yaml file. [default:config/config.yaml]', default='config/config.yaml',
              help='the config file which contain your program information.')
def pathseq_report(input, configfile):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    #===================================================================================================================
    ### step1 整理文件夹
    #===================================================================================================================
    with open(configfile) as f:
        config = ordered_yaml_load(f)
    program_num = config['report']['Project_Num']
    program_path = os.path.abspath(input)
    config_ref_file = config['database']['reference']
    with open(f'{config_ref_file}/config.yaml') as f:
        config_ref = ordered_yaml_load(f)

    os.chdir(program_path)
    log_file = open("%s/logs/pathseq_report_log.txt" % (program_path), "w")
    if os.path.exists("%s/result/report/%s_PathSeq_Report" % (program_path, program_num)):
        shutil.rmtree("%s/result/report/%s_PathSeq_Report" % (program_path, program_num ), ignore_errors=True)
        os.makedirs("%s/result/report/%s_PathSeq_Report" % (program_path, program_num ))
    else:
        os.makedirs("%s/result/report/%s_PathSeq_Report" % (program_path, program_num ))
        print("Create Pathseq Report...")
        log_file.write("Create Pathseq  Report..." + "\n")

    output_dir = "%s/result/report/%s_PathSeq_Report" % (program_path, program_num)

    ###################################### 1.Bacteria_reads_Bacteria_UMIs_distribution
    if os.path.exists("%s/result/pathseq/pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution" % (program_path)):
        subprocess.call('ln -s   %s/result/pathseq/pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution  %s/1.Bacteria_reads_Bacteria_UMIs_distribution  ' % (program_path, output_dir),                        shell=True)
    else:
        print("Can not find pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find pathseq_visualization/1.Bacteria_reads_Bacteria_UMIs_distribution results!!!!!!!!!!!!!!!!" + "\n")

    ###################################### 2.Top10_most_dominant_Bacterial_genera
    if os.path.exists("%s/result/pathseq/pathseq_visualization/2.Top10_most_dominant_Bacterial_genera" % (program_path)):
        subprocess.call('ln -s  %s/result/pathseq/pathseq_visualization/2.Top10_most_dominant_Bacterial_genera   %s/2.Top10_most_dominant_Bacterial_genera ' % (program_path, output_dir),                        shell=True)
    else:
        print("Can not find pathseq_visualization/2.Top10_most_dominant_Bacterial_genera results!!!!!!!!!!!!!!!!")
        log_file.write("Can not find pathseq_visualization/2.Top10_most_dominant_Bacterial_genera results!!!!!!!!!!!!!!!!" + "\n")
    subprocess.call('cp -r %s/src/ %s/src' % (scriptdir, output_dir), shell=True)
    #===================================================================================================================
    ### step2 生成网页报告
    #===================================================================================================================
    os.chdir(output_dir)
    # Define header info
    header_info = {key: config["report"][key] for key in config["report"] if
                   key in ["Project_Num", "Customer", "Species", "Executor"]}
    header_info["项目编号"] = header_info.pop("Project_Num")
    header_info["客户姓名"] = header_info.pop("Customer")
    header_info["物种信息"] = header_info.pop("Species")
    header_info["项目执行人员"] = header_info.pop("Executor")
    title = "空间转录组肿瘤微生物鉴定分析结题报告"
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
    report.add_yaml_config(os.path.join(scriptdir, 'pathseq_report.yaml'))
    background_intro = report.add_section('背景概述')
    principle_intro = report.add_section('原理简介')
    apply_intro = report.add_section('应用示例')
    #### 项目分析结果
    result_module = report.add_section('项目分析结果')
    section1 = result_module.add_section('微生物组成概况')#===============================================================
    section1.add_table('1.Bacteria_reads_Bacteria_UMIs_distribution/summary_statistic.xls',
                         caption='各个样本中微生物鉴定情况概览',
                         headerdict={'sampleid': '样本名称',
                                     'Total Spaceranger Detect Spots': '样本中检出的Spots数目',
                                     'Number of Positive Genus Capture Spots ': '样本中检出的包含微生物的Spots数目',
                                     'Total  Number of Genus Reads': '样本中微生物相关的Reads数目',
                                     'Total  Number of Genus UMIs': '样本中微生物相关的UMIs数目'
                                     })
    section1_1 = section1.add_section('微生物UMIs数目')
    section1_2 = section1.add_section('微生物Reads数目')
    section1_3 = section1.add_section('小提琴图展示')
    section2 = result_module.add_section('Top10微生物组成情况')#==========================================================
    section2_1 = section2.add_section('饼状统计图')
    section2_2 = section2.add_section('柱状统计图')
    section2_3 = section2.add_section('空间分布图')
    section3= report.add_section("附录")#===============================================================================
    section3_1 = section3.add_section("分析软件列表")
    #section3_2= section3.add_section("参考数据库链接")
    section4 = report.add_section('申明')
    report.write_to('report.html')
if __name__ == "__main__":
    pathseq_report()
