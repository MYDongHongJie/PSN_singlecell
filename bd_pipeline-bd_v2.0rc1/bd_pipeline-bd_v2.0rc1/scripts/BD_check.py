#!/usr/bin/env python3
# encoding: utf-8
# Author Liuhongyan
# Modified luyao
# Date: 2020/1/24 13:31
# Desc   : This script is used to get bd summary information 
check_yaml="""
Putative_Cell_Count: 
  - "2000"
  - "20000"
Pct_Q30_Bases_in_Filtered_R2: "80%"
Pct_Cellular_Reads_Aligned_Uniquely: 
  - "30%"
  - "70%"
Median_Genes_per_Cell:
  - "800"
  - "900"
Annotated_Transcriptome_Pct: 
  - "25%"
  - "30%"
Pct_Read_Pair_Overlap: 
  - "10%"
  - "20%"
Pct_Reads_Filtered_Out: 
  - "10%"
  - "20%"
Pct_Assigned_to_Cell_Labels: 
  - "70%"
  - "80%"
Pct_TCR: "80%"
Pct_BCR: "80%"
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
parser.add_argument("-s", "--stat",help="all samples bd summary")
parser.add_argument("-m", "--html", help="all samples bd summary html, 用 , 分隔")
#parser.add_argument("-p", "--png", help="all samples bd summary png, 用 , 分隔")
# parser.add_argument("-step", "--step",type=str,help="eg：bd")
parser.add_argument("-o", "--output",help="output statistic file")
args = parser.parse_args()

stat = args.stat
# step = args.step
inconfig = args.config
check = args.output
html = args.html
#png = args.png

def bd_check(stat,html,Putative_Cell_Count,Pct_Q30_Bases_in_Filtered_R2,Pct_Cellular_Reads_Aligned_Uniquely,Median_Genes_per_Cell,Total_Reads_in_FASTQ,Annotated_Transcriptome_Pct,Pct_Read_Pair_Overlap,Pct_Reads_Filtered_Out,Pct_Assigned_to_Cell_Labels,Pct_TCR,Pct_BCR,Project_Num,Task_Num,Species,Customer,Project_management,Laboratory,QC_report,stat_type):
    """
    BD check result 
    
    """
    ### 各比较指标 ### 
    Med_Genes_Cell_stand = int(Median_Genes_per_Cell[1])
    Med_Genes_Cell_min = int(Median_Genes_per_Cell[0])
    Est_Num_Cells_max_stand = int(Putative_Cell_Count[1])
    Est_Num_Cells_min_stand = int(Putative_Cell_Count[0])
    Frac_Reads_Cells_stand = int(Pct_Cellular_Reads_Aligned_Uniquely[1].replace("%",""))
    Frac_Reads_Cells_min = int(Pct_Cellular_Reads_Aligned_Uniquely[0].replace("%",""))
    Confi_Trome_stand = int(Annotated_Transcriptome_Pct[1].replace("%",""))
    Confi_Trome_min = int(Annotated_Transcriptome_Pct[0].replace("%",""))
    #add VDJ checkpoint
    Pct_Read_Pair_Overlap_stand = int(Pct_Read_Pair_Overlap[1].replace("%",""))
    Pct_Read_Pair_Overlap_min = int(Pct_Read_Pair_Overlap[0].replace("%",""))
    Pct_Reads_Filtered_Out_stand = int(Pct_Reads_Filtered_Out[1].replace("%",""))
    Pct_Reads_Filtered_Out_min = int(Pct_Reads_Filtered_Out[0].replace("%",""))
    Pct_Assigned_to_Cell_Labels_stand = int(Pct_Assigned_to_Cell_Labels[1].replace("%",""))
    Pct_Assigned_to_Cell_Labels_min = int(Pct_Assigned_to_Cell_Labels[0].replace("%",""))
    Pct_TCR_stand = float(Pct_TCR.replace("%",""))
    Pct_BCR_stand = float(Pct_BCR.replace("%",""))

    Antisense_stand = 10 
    if int(Total_Reads_in_FASTQ.replace("M","")) ==350:
        data_min=300
    elif int(Total_Reads_in_FASTQ.replace("M","")) ==500:
        data_min=450
    else:
        data_min = int(Total_Reads_in_FASTQ.replace("M",""))*0.9
    ### 变量设置 ###
    cell_big_error_sampleid = []
    info = ''
    all_sample_text = ''
    global_error_flag = False
    global_sth_flag  = False
    df = pd.read_csv(stat,index_col=0)
    if stat_type =="dan":
        df = df[~df.index.str.contains("_" + stat_type)]
    for id in range(0,df.shape[0]): 
        sampleid=df.index[id]
        sample_error = ''
        Q30_message = '' 
        sample_name = "    " + sampleid + " 样本"
        sample_warning = ''
        error_flag = False
        warning_flag = False
        Est_Num_Cells_Sam=int(df.loc[sampleid]["Putative_Cell_Count"])
        Med_Genes_Cell_Sam=int(df.loc[sampleid]["Median_Genes_per_Cell"])
        Q30_Bases_Sam=float(df.loc[sampleid]["Pct_Q30_Bases_in_Filtered_R2"])
        Frac_Reads_Cells_Sam=float(df.loc[sampleid]["Pct_Cellular_Reads_Aligned_Uniquely"])
        Number_Reads_Sam=int(df.loc[sampleid]["Total_Reads_in_FASTQ"])/1000000
        Confi_Trome_Sam = float(df.loc[sampleid]["Annotated_Transcriptome_Pct"])
        Pct_Read_Pair_Overlap_Sam=float(df.loc[sampleid]["Pct_Read_Pair_Overlap"])
        Pct_Reads_Filtered_Out_Sam=float(df.loc[sampleid]["Pct_Reads_Filtered_Out"])
        Pct_Assigned_to_Cell_Labels_Sam=float(df.loc[sampleid]["Pct_Assigned_to_Cell_Labels"])
        if  'MeanMCP_BCR' in  df.columns :
            Pct_TCR_Sam=float(df.loc[sampleid]["Pct_TCR"])
            Pct_BCR_Sam=float(df.loc[sampleid]["Pct_BCR"])

        #Antisense_Sam = float(df.loc[sampleid]["Reads Mapped Antisense to Gene"].replace("%",""))
        #Seq_Saturation = float( df.loc[sampleid]["Sequencing Saturation"].replace("%", ""))
        if int(Total_Reads_in_FASTQ.replace("M", "")) != 33:
            ### Total_Reads_in_FASTQ ###
            if ( Number_Reads_Sam < data_min ):
                sample_error += "，下机数据不足，低于最低要求 " + str(int(data_min)) + "M，请核实该项目下机数据"
                error_flag = True 
            ### Annotated_Transcriptome_Pct ###
            if ( Confi_Trome_Sam < Confi_Trome_min): 
                sample_error += "，Annotated_Transcriptome_Pct 低于 " + str(Confi_Trome_stand) +  "%"
                error_flag = True
            elif ( Confi_Trome_min <= Confi_Trome_Sam < Confi_Trome_stand ):
                sample_warning += "，Annotated_Transcriptome_Pct 稍低于 " + str(Confi_Trome_stand) +  "%"
                warning_flag = True
            ### Putative_Cell_Count ###
            if ( Est_Num_Cells_Sam > Est_Num_Cells_max_stand ):
                cell_big_error_sampleid.append(sampleid)
            elif ( Est_Num_Cells_Sam < Est_Num_Cells_min_stand ):
                sample_error += "，Putative_Cell_Count 低于 " + str(Est_Num_Cells_min_stand)
                error_flag = True 

            ### Median_Genes_per_Cell ### 
            if ( Med_Genes_Cell_Sam < Med_Genes_Cell_min ):
                sample_error += "，Median_Genes_per_Cell 低于 " + str(Med_Genes_Cell_stand)
                error_flag = True
            elif ( Med_Genes_Cell_min <= Med_Genes_Cell_Sam < Med_Genes_Cell_stand ):
                sample_warning += "，Median_Genes_per_Cell 稍低于 " + str(Med_Genes_Cell_stand)
                warning_flag = True

            ### Pct_Cellular_Reads_Aligned_Uniquely ### 
            if ( Frac_Reads_Cells_Sam < Frac_Reads_Cells_min): 
                sample_error += "，Pct_Cellular_Reads_Aligned_Uniquely 低于 " + str(Frac_Reads_Cells_stand)+  "%"
                error_flag = True
            elif ( Frac_Reads_Cells_min <= Frac_Reads_Cells_Sam < Frac_Reads_Cells_stand ):
                sample_warning += "，Pct_Cellular_Reads_Aligned_Uniquely 稍低于 " + str(Frac_Reads_Cells_stand) +  "%"
                warning_flag = True
            # ### Antisense_Sam
            # if (include_introns == True and Antisense_Sam > Antisense_stand ) :
            #     sample_warning += "，Reads Mapped Antisense to Gene 较高，在细胞核转录组中属正常现象"
            #     warning_flag = True
            ### Pct_Q30_Bases_in_Filtered_R2 ###
            if ( Q30_Bases_Sam < float(Pct_Q30_Bases_in_Filtered_R2.replace("%","")) ):
                Q30_message += "，Pct_Q30_Bases_in_Filtered_R2 低于" + Pct_Q30_Bases_in_Filtered_R2 +  "，反馈时请附上 Fastqc 结果。"
        #add  VDJ check 20230706
            #Pct_Read_Pair_Overlap
            if ( Pct_Read_Pair_Overlap_Sam > Pct_Read_Pair_Overlap_stand): 
                sample_error += "，Pct_Read_Pair_Overlap 高于 " + str(Pct_Read_Pair_Overlap_stand)+  "%"
                error_flag = True
            elif ( Pct_Read_Pair_Overlap_min <= Pct_Read_Pair_Overlap_Sam < Pct_Read_Pair_Overlap_stand ):
                sample_warning += "，Pct_Read_Pair_Overlap 稍高于 " + str(Pct_Read_Pair_Overlap_min) +  "%"
                warning_flag = True
            #Pct_Reads_Filtered_Out
            if ( Pct_Reads_Filtered_Out_Sam > Pct_Reads_Filtered_Out_stand): 
                sample_error += "，Pct_Reads_Filtered_Out 高于 " + str(Pct_Reads_Filtered_Out_stand)+  "%"
                error_flag = True
            elif ( Pct_Reads_Filtered_Out_min <= Pct_Reads_Filtered_Out_Sam < Pct_Reads_Filtered_Out_stand ):
                sample_warning += "，Pct_Reads_Filtered_Out 稍高于 " + str(Pct_Reads_Filtered_Out_min) +  "%"
                warning_flag = True
            #Pct_Assigned_to_Cell_Labels
            if ( Pct_Assigned_to_Cell_Labels_Sam < Pct_Assigned_to_Cell_Labels_min): 
                sample_error += "，Pct_Assigned_to_Cell_Labels 低于 " + str(Pct_Assigned_to_Cell_Labels_min)+  "%"
                error_flag = True
            elif ( Pct_Assigned_to_Cell_Labels_min <= Pct_Assigned_to_Cell_Labels_Sam < Pct_Assigned_to_Cell_Labels_stand ):
                sample_warning += "，Pct_Assigned_to_Cell_Labels 稍低于 " + str(Pct_Assigned_to_Cell_Labels_stand) +  "%"
                warning_flag = True
            #Pct_TCR BCR
            if  'MeanMCP_BCR' in  df.columns :
                if ( Pct_TCR_Sam < 70): 
                    sample_error += "，Pct_TCR 低于 " + str(70)+  "%"
                    error_flag = True
                elif ( 70 <= Pct_TCR_Sam < Pct_TCR_stand ):
                    sample_warning += "，Pct_TCR 稍高于 " + str(70) +  "%"
                    warning_flag = True
                if ( Pct_BCR_Sam < 70): 
                    sample_error += "，Pct_BCR 低于 " + str(70)+  "%"
                    error_flag = True
                elif ( 70 <= Pct_BCR_Sam < Pct_TCR_stand ):
                    sample_warning += "，Pct_BCR 稍高于 " + str(70) +  "%"
                    warning_flag = True
            print(error_flag,sample_error)
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
            ### Putative_Cell_Count ###
            if ( Est_Num_Cells_Sam < 4000 ):
                sample_error += "，Putative_Cell_Count 低于 4000"
                error_flag = True 
            ### Median_Genes_per_Cell ### 
            if ( Med_Genes_Cell_Sam < 600 ):
                sample_error += "，Median_Genes_per_Cell 低于 600"
                error_flag = True
            ### Seq_Saturation ###
            #if ( Seq_Saturation > 40 ):
            #    sample_error += "，Sequencing Saturation 高于 40%"
            #    error_flag = True
            ### 单样本质控结果汇总 ###
            if error_flag:
                global_error_flag = True
                all_sample_text += sample_name + sample_error + "，请进一步确认是否加测。" + "\n"
            else:
                all_sample_text += sample_name + "各质控指标正常，建议加测。\n"
            ### 判断数据量是否是加测后数据
            if ( Number_Reads_Sam > 100 ):
                all_sample_text += '<font color="#FF0000"> 请注意！！！样本数据量大于100M，请确认是否为试测数据。</font> \n'

    ### 所有样本质控结果汇总 ### 
    if global_error_flag :
        header_flag="质控异常"
    elif global_sth_flag :
        header_flag="质控反馈"
    else:
        header_flag="质控反馈"
        if int(Total_Reads_in_FASTQ.replace("M", "")) != 33:
            all_sample_text = "    样本各质控指标正常。\n"
        else:
            all_sample_text = "    样本各质控指标正常，建议加测。\n"
    ## paste header and text
    all_samples = ", ".join(list(df.index))
    
    if QC_report:
        import yaml,time,sys
        header_flag = header_flag + "和项目报告"

    if 'MeanMCP_BCR'  in df.columns.to_list():
        Type_vdj = "True"
        header = "【" + header_flag + "】"+ Task_Num + " " + Customer + "老师 " + "BD单细胞转录组和免疫组项目质控反馈"
        info = Project_management + '：\n    您好。\n    ' + Task_Num + " " + Customer + "老师 单细胞转录组和免疫组项目 " + \
            str(len(df.index)) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n"
    else :
        if stat_type =="mRNA":
            header = "【" + header_flag + "】"+ Task_Num + " " + Customer + "老师 " + "BD单细胞转录组项目质控反馈"
            info = Project_management + '：\n    您好。\n    ' + Task_Num + " " + Customer + "老师 单细胞转录组项目 " + \
                str(len(df.index)) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n"
        else:
            header = "【" + header_flag + "】"+ Task_Num + " " + Customer + "老师 " + "BD单细胞转录组和表面蛋白项目项目质控反馈"
            info = Project_management + '：\n    您好。\n    ' + Task_Num + " " + Customer + "老师 单细胞转录组和表面蛋白项目 " + \
                str(len(df.index)) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n"

    if cell_big_error_sampleid:
       big_samples = "、 ".join(cell_big_error_sampleid)
       cell_big_error = Laboratory + '：\n    您好。\n    如下图所示：' + \
           big_samples + " 样本" + "细胞数高于 " + str(Est_Num_Cells_max_stand) + ", 请知悉! \n"
       info += cell_big_error

    if int(Total_Reads_in_FASTQ.replace("M", "")) == 33:
        test_proj_text = "\n该项目为试测10G项目。\n"
        info += test_proj_text

    if QC_report:
        import yaml,time,sys
        header = header + "和项目报告"
        #info = info + \
        #       "\n    该项目为只质控不分析项目，已上传云报告并将完整QC报告上传至OBS，\n    路径：obs://oe-scrna/Analysis_Report/" + \
        #       Project_Num + "/" + Task_Num + "_QC_Report" + time.strftime('%Y_%m_%d') +".zip \n "
        info = info + \
               "\n    该项目为只质控不分析项目，完整QC报告已上传云报告。\n "
    ### 返回 邮件相关信息 ###
    return (info,header)

# def final(header, email):
#     error_info=header
#     output = "NA"
#     send_email(email,header,output,error_info)

def send_email(emailto,header,context_html,info,stat):
    '''
        邮件发送
    '''
    emailfrom="sc_autosend@oebiotech.com"
    password="OEbiosc1616"
    message = MIMEMultipart('related')
    # message = MIMEText(out, 'plain', 'utf-8')
    message['From'] = Header(emailfrom, 'utf-8')
    message['To'] =  ",".join(emailto)
    #邮件正文内容
    # message.attach(MIMEText(info, 'html', 'utf-8'))
    #png_list = png.split(",")
    #pic_inline = ''
    #for index,png_file in enumerate(png_list):
    #    with open(png_file, 'rb') as image:
    #        image_info = MIMEImage(image.read())
    #        image_info.add_header('Content-Id', f'<image{index + 1}>')
    #        message.attach(image_info)
    #        tmp_pic_inline = f'''
    #                            <br><img src="cid:image{index + 1}" width="600"></br>
    #                        '''
    #        pic_inline += tmp_pic_inline
    info_html = info.replace('\n','<br> &emsp; ')
    info_html = info_html.replace('<br>&emsp;<br>&emsp;', '<br><br>')
    # table # 
    df = pd.read_csv(stat,index_col=0,dtype = 'str')
    bd_html = df.T.to_html(index=True, border=1, justify="left",bold_rows=False)
    
    if 'MeanMCP_BCR'  in df.columns.to_list():
        Type_vdj = "True"
        bd_html = bd_html.replace('<table border="1" class="dataframe">','<h3>BD转录组和免疫组汇总结果:</h3> <table bordercolor="grey" border="1" cellspacing="0" cellpadding="5" class="dataframe">' )
    else :
        Type_vdj = "True"
        bd_html = bd_html.replace('<table border="1" class="dataframe">','<h3>BD转录组汇总结果:</h3> <table bordercolor="grey" border="1" cellspacing="0" cellpadding="5" class="dataframe">' )
    info_html = info_html + bd_html
    #message.attach(MIMEText( info_html + pic_inline, 'html', 'utf-8'))
    message.attach(MIMEText( info_html, 'html', 'utf-8'))
    
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
    #add csv file
    att2 = MIMEText(open(stat, 'rb').read(), 'base64', 'utf-8')
    att2["Content-Type"] = 'application/octet-stream'
    att2["Content-Disposition"] = 'attachment; filename="metrics_summary_stat.csv"'
    message.attach(att2)

    subject = header
    message['Subject'] = Header(subject, 'utf-8')
    smtp=smtplib.SMTP_SSL('smtp.exmail.qq.com',port=465)
    smtp.connect('smtp.exmail.qq.com')
    smtp.login(emailfrom,password)
    smtp.sendmail(emailfrom,emailto, message.as_string())

if __name__ == '__main__':
    check_dict = yaml.load(check_yaml, Loader=yaml.FullLoader)
    Putative_Cell_Count = check_dict["Putative_Cell_Count"]
    Pct_Q30_Bases_in_Filtered_R2 = check_dict["Pct_Q30_Bases_in_Filtered_R2"]
    Pct_Cellular_Reads_Aligned_Uniquely = check_dict["Pct_Cellular_Reads_Aligned_Uniquely"]
    Median_Genes_per_Cell = check_dict["Median_Genes_per_Cell"]
    Annotated_Transcriptome_Pct = check_dict["Annotated_Transcriptome_Pct"]
    Pct_Read_Pair_Overlap = check_dict["Pct_Read_Pair_Overlap"]
    Pct_Reads_Filtered_Out = check_dict["Pct_Reads_Filtered_Out"]
    Pct_Assigned_to_Cell_Labels = check_dict["Pct_Assigned_to_Cell_Labels"]
    Pct_TCR = check_dict["Pct_TCR"]
    Pct_BCR = check_dict["Pct_BCR"]
    
    config = {}
    with open(inconfig,'r',encoding='utf-8') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    Project_Num = config["report"]["Task_Num"].split('-')[0]
    Task_Num = config["report"]["Task_Num"]
    Customer = config["report"]["Customer"]
    Species = config["report"]["Species"]
    Project_management =  config["email"]["Project_management"]
    Laboratory = config["email"]["Laboratory"]
    emailto = config["email"]["emailto"]
    stat_type = config['BD_params']['VDJ']
    #raw_data_obs_address = config["BD_params"]["raw_data_obs_address"][0]
    #include_introns =  config["bd_params"]["include_introns"]
    Total_Reads_in_FASTQ = config["BD_params"]["Number_of_Reads"]
    QC_report = config["BD_params"]["QC_report"]

    #if "bd" in step:
    (info,header)=bd_check(stat,html,Putative_Cell_Count,Pct_Q30_Bases_in_Filtered_R2,Pct_Cellular_Reads_Aligned_Uniquely,Median_Genes_per_Cell,Total_Reads_in_FASTQ,Annotated_Transcriptome_Pct,Pct_Read_Pair_Overlap,Pct_Reads_Filtered_Out,Pct_Assigned_to_Cell_Labels,Pct_TCR,Pct_BCR,Project_Num,Task_Num,Species,Customer,Project_management,Laboratory,QC_report,stat_type)
    #print(emailto,header,html,info,stat)
    send_email(emailto,header,html,info,stat)

    with open(check, "w",encoding='utf-8-sig') as log:
        log.write(info)

    # if "finish" in step:
    #     header = Project_Num + '  项目分析完成' +  " ( " + Customer + "老师  " + title  + " )"
    #     final(header,single)

#QC_report 
#tar -zcvf DOE202212817-b12_QC_20230306.tar.gz --exclude=G12037_CD_huichang/Logs --exclude=G2969757_CD_jiechang/Logs  G12037_CD_huichang/ G2969757_CD_jiechang/
