#!/usr/bin/env python3
# encoding: utf-8
# Author Liuhongyan
# Modified luyao
# Date: 2020/1/24 13:31
# Desc   : This script is used to get cellranger summary information 
check_yaml="""
Estimated_Number_of_Cells: 
  - "2000"
  - "12000"
Q30_Base_in_RNA_Read: "80%"
Fraction_Reads_in_Cells: 
  - "30%"
  - "70%"
Median_Genes_per_Cell:
  - "800"
  - "900"
Reads_Mapped_Confidently_to_Transcriptome: 
  - "25%"
  - "30%"
"""
import os
import re
import yaml
import argparse
import smtplib
import pandas as pd
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.header import Header
from email.mime.multipart import MIMEMultipart

parser = argparse.ArgumentParser(description="质控自动反馈脚本")
parser.add_argument("-i", "--config",help="config.yaml")
parser.add_argument("-s", "--stat",help="all samples cellranger summary")
parser.add_argument("-m", "--html", help="all samples cellranger summary html, 用 , 分隔")
parser.add_argument("-p", "--png", help="all samples cellranger summary png, 用 , 分隔")
# parser.add_argument("-step", "--step",type=str,help="eg：cellranger")
parser.add_argument("-o", "--output",help="output statistic file")
args = parser.parse_args()

stat = args.stat
# step = args.step
inconfig = args.config
check = args.output
html = args.html
png = args.png

def cellarnger_check(stat,html,cellranger_only,Estimated_Number_of_Cells,Q30_Base_in_RNA_Read,Fraction_Reads_in_Cells,Median_Genes_per_Cell,Number_of_Reads,Reads_Mapped_Confidently_to_Transcriptome,Project_Num,Species,Customer,Project_management,Laboratory,raw_data_obs_address,include_introns):
    """
    cellarnger check result 
    
    """
    ### 各比较指标 ### 
    Med_Genes_Cell_stand = int(Median_Genes_per_Cell[1])
    Med_Genes_Cell_min = int(Median_Genes_per_Cell[0])
    Est_Num_Cells_max_stand = int(Estimated_Number_of_Cells[1])
    Est_Num_Cells_min_stand = int(Estimated_Number_of_Cells[0])
    Frac_Reads_Cells_stand = int(Fraction_Reads_in_Cells[1].replace("%",""))
    Frac_Reads_Cells_min = int(Fraction_Reads_in_Cells[0].replace("%",""))
    Confi_Trome_stand = int(Reads_Mapped_Confidently_to_Transcriptome[1].replace("%",""))
    Confi_Trome_min = int(Reads_Mapped_Confidently_to_Transcriptome[0].replace("%",""))
    Antisense_stand = 10 
    if int(Number_of_Reads.replace("M","")) ==350:
        data_min=300
    elif int(Number_of_Reads.replace("M","")) ==500:
        data_min=450
    else:
        data_min = int(Number_of_Reads.replace("M",""))*0.9
    ### 变量设置 ###
    cell_big_error_sampleid = []
    cellranger_only_header = ''
    info = ''
    all_sample_text = ''
    global_error_flag = False
    global_sth_flag  = False
    df = pd.read_csv(stat,index_col=0)
    for id in range(0,df.shape[0]): 
        sampleid=df.index[id]
        sample_error = ''
        Q30_message = '' 
        sample_name = "    " + sampleid + " 样本"
        sample_warning = ''
        error_flag = False
        warning_flag = False

        Est_Num_Cells_Sam=int(str(df.loc[sampleid]["Estimated Number of Cells"]).replace(',',''))
        Med_Genes_Cell_Sam=int(str(df.loc[sampleid]["Median Genes per Cell"]).replace(',',''))
        Q30_Bases_Sam=float(df.loc[sampleid]["Q30 Bases in RNA Read"].replace("%",""))
        Frac_Reads_Cells_Sam=float(df.loc[sampleid]["Fraction Reads in Cells"].replace("%",""))
        Number_Reads_Sam=int(str(df.loc[sampleid]["Number of Reads"]).replace(",",""))/1000000
        Confi_Trome_Sam = float(df.loc[sampleid]["Reads Mapped Confidently to Transcriptome"].replace("%",""))
        Antisense_Sam = float(df.loc[sampleid]["Reads Mapped Antisense to Gene"].replace("%",""))
        Seq_Saturation = float( df.loc[sampleid]["Sequencing Saturation"].replace("%", ""))
        if int(Number_of_Reads.replace("M", "")) != 33:
            ### Number of Reads ###
            if ( Number_Reads_Sam < data_min ):
                sample_error += "，下机数据不足，低于最低要求 " + str(int(data_min)) + "M，请核实该项目下机数据"
                error_flag = True 
            ### Reads Mapped Confidently to Transcriptome ###
            if ( Confi_Trome_Sam < Confi_Trome_min): 
                sample_error += "，Reads Mapped Confidently to Transcriptome 低于 " + str(Confi_Trome_stand) +  "%"
                error_flag = True
            elif ( Confi_Trome_min <= Confi_Trome_Sam < Confi_Trome_stand ):
                sample_warning += "，Reads Mapped Confidently to Transcriptome 稍低于 " + str(Confi_Trome_stand) +  "%"
                warning_flag = True
            ### Estimated Number of Cells ###
            if ( Est_Num_Cells_Sam > Est_Num_Cells_max_stand ):
                cell_big_error_sampleid.append(sampleid)
            elif ( Est_Num_Cells_Sam < Est_Num_Cells_min_stand ):
                sample_error += "，Estimated Number of Cells 低于 " + str(Est_Num_Cells_min_stand)
                error_flag = True 

            ### Median Genes per Cell ### 
            if ( Med_Genes_Cell_Sam < Med_Genes_Cell_min ):
                sample_error += "，Median Genes per Cell 低于 " + str(Med_Genes_Cell_stand)
                error_flag = True
            elif ( Med_Genes_Cell_min <= Med_Genes_Cell_Sam < Med_Genes_Cell_stand ):
                sample_warning += "，Median Genes per Cell 稍低于 " + str(Med_Genes_Cell_stand)
                warning_flag = True

            ### Fraction Reads in Cells ### 
            if ( Frac_Reads_Cells_Sam < Frac_Reads_Cells_min): 
                sample_error += "，Fraction Reads in Cells 低于 " + str(Frac_Reads_Cells_stand)+  "%"
                error_flag = True
            elif ( Frac_Reads_Cells_min <= Frac_Reads_Cells_Sam < Frac_Reads_Cells_stand ):
                sample_warning += "，Fraction Reads in Cells 稍低于 " + str(Frac_Reads_Cells_stand) +  "%"
                warning_flag = True
            # ### Antisense_Sam
            # if (include_introns == True and Antisense_Sam > Antisense_stand ) :
            #     sample_warning += "，Reads Mapped Antisense to Gene 较高，在细胞核转录组中属正常现象"
            #     warning_flag = True
            ### Q30 Bases in RNA Read ###
            if ( Q30_Bases_Sam < float(Q30_Base_in_RNA_Read.replace("%","")) ):
                Q30_message += "，Q30 Bases in RNA Read 低于" + Q30_Base_in_RNA_Read +  "，反馈时请附上 Fastqc 结果。"

            ### 单样本质控结果汇总 ### 
            if error_flag:
                global_error_flag = True
                all_sample_text += sample_name + sample_error + "。" + Q30_message + "\n"
            else:
                if warning_flag:
                    global_sth_flag = True
                    all_sample_text += sample_name + sample_warning + "，可以继续后续分析。" + Q30_message + "\n"
                else:
                    if Q30_message:
                        global_sth_flag = True
                        all_sample_text += sample_name + Q30_message + "其余各指标正常。"
                    else:
                        all_sample_text += sample_name + "各质控指标正常。\n"
        else:
            ### 试测 10G 项目 ###
            ### Estimated Number of Cells ###
            if ( Est_Num_Cells_Sam < 4000 ):
                sample_error += "，Estimated Number of Cells 低于 4000"
                error_flag = True 
            ### Median Genes per Cell ### 
            if ( Med_Genes_Cell_Sam < 600 ):
                sample_error += "，Median Genes per Cell 低于 600"
                error_flag = True
            ### Seq_Saturation ###
            if ( Seq_Saturation > 40 ):
                sample_error += "，Sequencing Saturation 高于 40%"
                error_flag = True
            ### 单样本质控结果汇总 ###
            if error_flag:
                global_error_flag = True
                all_sample_text += sample_name + sample_error + "，请进一步确认是否加测。" + "\n"
            else:
                all_sample_text += sample_name + "各质控指标正常，建议加测。\n"

    ### 所有样本质控结果汇总 ### 
    if global_error_flag :
        header_flag="质控异常"
    elif global_sth_flag :
        header_flag="质控反馈"
    else:
        header_flag="质控反馈"
        if int(Number_of_Reads.replace("M", "")) != 33:
            all_sample_text = "    样本各质控指标正常。\n"
        else:
            all_sample_text = "    样本各质控指标正常，建议加测。\n"
    ## paste header and text
    if cellranger_only:
        import yaml,time,sys
        header = "【" + header_flag +"和项目报告】"+ Task_Num + " " + Customer + "老师 " + "单细胞转录组项目质控反馈和项目报告"
        all_samples = ", ".join(list(df.index))
        info = Project_management + '：\n    您好。\n    ' + Task_Num + " " + Customer + "老师 单细胞转录组项目 " + \
               str(len(df.index)) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n" + \
               "\n    该项目为只质控不分析项目，已上传云报告并将完整QC报告上传至OBS，\n    路径：obs://oe-scrna/Analysis_Report/" + \
               Project_Num + "/" + Task_Num + "_QC_Report" + time.strftime('%Y_%m_%d') +".tar.gz \n "
    else:
        header = "【" + header_flag + "】"+ Project_Num + " " + Customer + "老师 " + "单细胞转录组项目质控反馈"
        all_samples = ", ".join(list(df.index))
        info = Project_management + '：\n    您好。\n    ' + Project_Num + " " + Customer + "老师 单细胞转录组项目 " + \
               str(len(df.index)) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n"

    # if cell_big_error_sampleid:
    #    big_samples = "、 ".join(cell_big_error_sampleid)
    #    cell_big_error = Laboratory + '：\n    您好。\n    如下图所示：' + \
    #        big_samples + " 样本" + "细胞数高于 " + str(Est_Num_Cells_max_stand) + ", 请知悉! \n"
    #    info += cell_big_error

    if int(Number_of_Reads.replace("M", "")) == 33:
        test_proj_text = "\n该项目为试测10G项目。\n"
        info += test_proj_text

    ### 返回 邮件相关信息 ###
    return (info,header)

# def final(header, email):
#     error_info=header
#     output = "NA"
#     send_email(email,header,output,error_info)

def send_email(emailto,header,context_html,png,info):
    '''
        邮件发送
    '''
    emailfrom="singlecell@oebiotech.com"
    password="9LqYhXGFfkwgVhh5"
    message = MIMEMultipart('related')
    # message = MIMEText(out, 'plain', 'utf-8')
    message['From'] = Header(emailfrom, 'utf-8')
    message['To'] =  ",".join(emailto)
    #邮件正文内容
    # message.attach(MIMEText(info, 'html', 'utf-8'))
    png_list = png.split(",")
    pic_inline = ''
    for index,png_file in enumerate(png_list):
        with open(png_file, 'rb') as image:
            image_info = MIMEImage(image.read())
            image_info.add_header('Content-Id', f'<image{index + 1}>')
            message.attach(image_info)
            tmp_pic_inline = f'''
                                <br><img src="cid:image{index + 1}" width="600"></br>
                            '''
            pic_inline += tmp_pic_inline
    info_html = info.replace('\n','<br> &emsp; ')
    message.attach(MIMEText( info_html + pic_inline, 'html', 'utf-8'))

    # 这里的filename可以任意写，写什么名字，邮件中显示什么名字
    # att1 = MIMEText(open(context_stat, 'rb').read(), 'plain', 'utf-8')
    # context_stat_name=os.path.abspath(context_stat).split("/")[-1]
    # att1["Content-Type"] = 'application/octet-stream'
    # att1["Content-Disposition"] = 'attachment; filename="' + context_stat_name + '"'
    # message.attach(att1)
    array = context_html.split(",")
    for id in range(len(array)):
        s_html = array[id]
        att = MIMEText(open(s_html, 'rb').read(), 'plain', 'utf-8')
        html_name=os.path.abspath(s_html).split("/")[-1]
        att["Content-Type"] = 'application/octet-stream'
        att["Content-Disposition"] = 'attachment; filename="' + html_name + '"'
        message.attach(att)

    subject = header
    message['Subject'] = Header(subject, 'utf-8')
    smtp=smtplib.SMTP_SSL('smtp.exmail.qq.com',port=465)
    smtp.connect('smtp.exmail.qq.com')
    smtp.login(emailfrom,password)
    smtp.sendmail(emailfrom,emailto, message.as_string())

if __name__ == '__main__':
    check_dict = yaml.load(check_yaml, Loader=yaml.FullLoader)
    Estimated_Number_of_Cells = check_dict["Estimated_Number_of_Cells"]
    Q30_Base_in_RNA_Read = check_dict["Q30_Base_in_RNA_Read"]
    Fraction_Reads_in_Cells = check_dict["Fraction_Reads_in_Cells"]
    Median_Genes_per_Cell = check_dict["Median_Genes_per_Cell"]
    Reads_Mapped_Confidently_to_Transcriptome = check_dict["Reads_Mapped_Confidently_to_Transcriptome"]
    
    config = {}
    with open(inconfig,'r',encoding='utf-8') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    Project_Num = config["report"]["Project_Num"]
    Task_Num = config["report"]["Task_Num"]
    Customer = config["report"]["Customer"]
    Species = config["report"]["Species"]
    Project_management =  config["email"]["Project_management"]
    Laboratory = config["email"]["Laboratory"]
    emailto = config["email"]["emailto"]
    raw_data_obs_address = config["cellranger_params"]["raw_data_obs_address"][0]
    include_introns =  config["cellranger_params"]["include_introns"]
    Number_of_Reads = config["cellranger_params"]["Number_of_Reads"]
    cellranger_only = config["cellranger_params"]["module"]["cellranger_report"]

    #if "cellranger" in step:
    (info,header)=cellarnger_check(stat,html,cellranger_only,Estimated_Number_of_Cells,Q30_Base_in_RNA_Read,Fraction_Reads_in_Cells,Median_Genes_per_Cell,Number_of_Reads,Reads_Mapped_Confidently_to_Transcriptome,Project_Num,Species,Customer,Project_management,Laboratory,raw_data_obs_address,include_introns)

    send_email(emailto,header,html,png,info)

    with open(check, "w",encoding='utf-8-sig') as log:
        log.write(info)

    # if "finish" in step:
    #     header = Project_Num + '  项目分析完成' +  " ( " + Customer + "老师  " + title  + " )"
    #     final(header,single)
