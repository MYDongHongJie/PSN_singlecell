#!/usr/bin/env python3
# encoding: utf-8
# Author Liuhongyan
# Date: 2022/11/4 
# Desc   : This script is used to get cellranger mutil summary information 

check_yaml="""
Estimated_Number_of_Cells: 
  - "2000"
  - "25000"
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
parser.add_argument("-s", "--stat",help="all samples cellranger summary, 用 , 分隔")
parser.add_argument("-m", "--html", help="all samples cellranger summary html, 用 , 分隔")
parser.add_argument("-o", "--outdir",help="outdir")
args = parser.parse_args()

stat = args.stat
# step = args.step
inconfig = args.config
outdir = args.outdir
html = args.html


def cellranger_summary(stat, outdir):
    index=0
    tcr="TCR"
    bcr="BCR"
    stat_list = stat.split(",")
    for stat in stat_list: 
        index=index+1
        sampleid = os.path.split(stat)[-1].split(".")[0].replace("_metrics_summary","")
        print(sampleid)
        df = pd.read_csv(stat,index_col=0)
        if index == 1:
            scrna_df = df[df['Library Type']=="Gene Expression"]
            scrna = scrna_df[['Metric Name','Metric Value']]
            scrna.columns = ['SampleID',sampleid]  
            if "VDJ T" in df['Library Type'].tolist():
                tcr_df = df[df['Library Type']=="VDJ T"]
                tcr = tcr_df[['Metric Name','Metric Value']]
                tcr.columns = ['SampleID',f'{sampleid}_TCR'] 
                tcr.set_index('SampleID', inplace=True)
                tcr = tcr[~tcr.index.duplicated()] 
            if "VDJ B" in df['Library Type'].tolist(): 
                bcr_df = df[df['Library Type']=="VDJ B"]
                bcr = bcr_df[['Metric Name','Metric Value']]
                bcr.columns = ['SampleID',f'{sampleid}_BCR'] 
                bcr.set_index('SampleID', inplace=True)
                bcr = bcr[~bcr.index.duplicated()]
        else:
            scrna_df = df[df['Library Type']=="Gene Expression"]
            scrna_tamp = scrna_df[['Metric Value']]
            scrna = pd.concat([scrna,scrna_tamp],axis=1)
            scrna = scrna.rename(columns={'Metric Value':sampleid})
            if "VDJ T" in df['Library Type'].tolist():
                tcr_df = df[df['Library Type']=="VDJ T"]
                tcr_df = tcr_df[['Metric Name','Metric Value']]
                tcr_df.columns = ['SampleID',f'{sampleid}_TCR'] 
                tcr_df.set_index('SampleID', inplace=True)
                tcr_df= tcr_df[~tcr_df.index.duplicated()]
                tcr = pd.concat([tcr,tcr_df],axis=1)
                #tcr = tcr.rename(columns={'Metric Value':f'{sampleid}_TCR'})
            if "VDJ B" in df['Library Type'].tolist(): 
                bcr_df = df[df['Library Type']=="VDJ B"]
                bcr_df = bcr_df[['Metric Name','Metric Value']]
                bcr_df.columns = ['SampleID',f'{sampleid}_BCR'] 
                bcr_df.set_index('SampleID', inplace=True)
                bcr_df= bcr_df[~bcr_df.index.duplicated()]
                bcr = pd.concat([bcr,bcr_df],axis=1)
    scrna = scrna.drop_duplicates()
    scrna.to_csv(f'{outdir}/scrna_cellranger_stat.xls', sep="\t", encoding='utf-8', index=False)
    if type(tcr) == pd.DataFrame:
        #tcr = tcr.drop_duplicates()
        ### 如果存在缺失值，将缺失值替换成 0 ###
        if tcr.isnull().values.any():
            tcr.fillna(0, inplace=True)
        tcr = tcr.reset_index()
        if "index" in tcr.columns:
            tcr = tcr.rename(columns={'index':'SampleID'})
        tcr.to_csv(f'{outdir}/TCR_cellranger_stat.xls',sep="\t", encoding='utf-8', index=False)
    if type(bcr) == pd.DataFrame:
        #bcr = bcr.drop_duplicates()
        ### 如果存在缺失值，将缺失值替换成 0 ###
        if bcr.isnull().values.any():
            bcr.fillna(0, inplace=True)
        ## 将行名设置成列
        bcr = bcr.reset_index()
        if "index" in bcr.columns:
            bcr = bcr.rename(columns={'index':'SampleID'})
        bcr.to_csv(f'{outdir}/BCR_cellranger_stat.xls', sep="\t", encoding='utf-8', index=False)
    ### 返回 邮件相关信息 ###
    return (scrna,tcr,bcr)


def cellranger_check(scrna,cellranger_only,tcr,bcr,Estimated_Number_of_Cells,Q30_Base_in_RNA_Read,Fraction_Reads_in_Cells,Median_Genes_per_Cell,Number_of_Reads,Reads_Mapped_Confidently_to_Transcriptome,Project_Num,Species,Customer,Project_management,Laboratory,remark,include_introns):
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

    Number_of_Reads_list =  Number_of_Reads.split(",")
    if int(Number_of_Reads_list[0].replace("M","")) ==350:
        scrna_data_min=300
    elif int(Number_of_Reads_list[0].replace("M","")) ==500:
        scrna_data_min=450
    else:
        scrna_data_min = int(Number_of_Reads_list[0].replace("M",""))*0.9

    if int(Number_of_Reads_list[1].replace("M","")) ==35:
        vdj_data_min=30
    else:
        vdj_data_min = int(Number_of_Reads_list[1].replace("M",""))*0.9

    ### 变量设置 ###
    cell_big_error_sampleid = []
    info = ''
    all_samples =list()
    all_sample_text = ''
    global_error_flag = False
    global_sth_flag  = False
    sample_normal = []
    sample = []
    sample_message=" "
    error_flag = False
    sample_abnormity = []
    data_error=" "
    # if include_introns == True :
    #     introns = "核"
    # else: 
    #     introns = ""
    introns = ""

    ### 单细胞转录组 #### 
    all_sample_text = "    转录组：\n"
    for id in range(1,scrna.shape[1]):
        #id=1
        df = scrna 
        sampleid=df.columns[id]
        sample_error = ''
        Q30_message = '' 
        sample_name = "    " + sampleid + " 样本"
        sample_warning = ''
        error_flag = False
        warning_flag = False
        Est_Num_Cells_Sam=int(str( df[df['SampleID']=="Estimated number of cells"][sampleid]['Library'] ).replace(',',''))
        Med_Genes_Cell_Sam=int(str( df[df['SampleID']=="Median genes per cell"][sampleid]['Cells'] ).replace(',',''))
        Q30_Bases_Sam = float( df[df['SampleID']=="Q30 RNA read"][sampleid]['Library'] .replace("%",""))
        Frac_Reads_Cells_Sam= float( df[df['SampleID']=="Confidently mapped reads in cells"][sampleid]['Cells'] .replace("%",""))
        Number_Reads_Sam=int(str( df[df['SampleID']=="Number of reads"][sampleid]['Library'] ).replace(",",""))/1000000
        Confi_Trome_Sam = float( df[df['SampleID']=="Confidently mapped to transcriptome"][sampleid]['Library'] .replace("%",""))
        Antisense_Sam = float( df[df['SampleID']=="Confidently mapped antisense"][sampleid]['Library'] .replace("%",""))
        Seq_Saturation = float( df[df['SampleID']=="Sequencing saturation"][sampleid]['Library'].replace("%", ""))
        all_samples.append(sampleid)
        if int(Number_of_Reads_list[0].replace("M","")) !=33:
            ### Number of Reads ###
            if ( Number_Reads_Sam < scrna_data_min ):
                sample_error += "，下机数据不足，低于最低要求 " + str(int(scrna_data_min)) + "M，请核实该项目下机数据"
                error_flag = True 
            ### Confidently mapped to transcriptome ###
            if ( Confi_Trome_Sam < Confi_Trome_min): 
                sample_error += "，Confidently mapped to transcriptome 低于 " + str(Confi_Trome_stand) +  "%"
                error_flag = True
            elif ( Confi_Trome_min <= Confi_Trome_Sam < Confi_Trome_stand ):
                sample_warning += "，Confidently mapped to transcriptome 稍低于 " + str(Confi_Trome_stand) +  "%"
                warning_flag = True
            ### Estimated number of cells ###
            if ( Est_Num_Cells_Sam > Est_Num_Cells_max_stand ):
                cell_big_error_sampleid.append(sampleid)
            elif ( Est_Num_Cells_Sam < Est_Num_Cells_min_stand ):
                sample_error += "，Estimated number of cells 低于 " + str(Est_Num_Cells_min_stand)
                error_flag = True 
            ### Median genes per cell ### 
            if ( Med_Genes_Cell_Sam < Med_Genes_Cell_min ):
                sample_error += "，Median genes per cell 低于 " + str(Med_Genes_Cell_stand)
                error_flag = True
            elif ( Med_Genes_Cell_min <= Med_Genes_Cell_Sam < Med_Genes_Cell_stand ):
                sample_warning += "，Median genes per cell 稍低于 " + str(Med_Genes_Cell_stand)
                warning_flag = True
            ### Fraction Reads in Cells ### 
            if ( Frac_Reads_Cells_Sam < Frac_Reads_Cells_min): 
                sample_error += "，Confidently mapped reads in cells 低于 " + str(Frac_Reads_Cells_stand)+  "%"
                error_flag = True
            elif ( Frac_Reads_Cells_min <= Frac_Reads_Cells_Sam < Frac_Reads_Cells_stand ):
                sample_warning += "，Confidently mapped reads in cells 稍低于 " + str(Frac_Reads_Cells_stand) +  "%"
                warning_flag = True
            # ### Antisense_Sam
            # if (include_introns == True and Antisense_Sam > Antisense_stand ) :
            #     sample_warning += "，Confidently mapped antisense 较高，在细胞核转录组中属正常现象"
            #     warning_flag = True
            ### Q30 RNA read ###
            if ( Q30_Bases_Sam < float(Q30_Base_in_RNA_Read.replace("%","")) ):
                Q30_message += "，Q30 RNA read 低于" + Q30_Base_in_RNA_Read +  "，反馈时请附上 Fastqc 结果。"
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
    ### 免疫组库
    data_error=" " 
    if type(tcr) == pd.DataFrame: 
        for id in range(1,tcr.shape[1]):
            sampleid=tcr.columns[id]
            all_samples.append(sampleid)
            ### Number of Reads ###
            Number_Reads_Sam=int(str( tcr[tcr['SampleID']=="Number of reads"][sampleid].iloc[0] ).replace(",",""))/1000000
            if ( Number_Reads_Sam < vdj_data_min ):
                data_error +=f"    {sampleid} 样本下机数据不足，低于最低要求 {str(int(vdj_data_min))} M，请核实该项目下机数据。\n "

    if type(bcr) == pd.DataFrame: 
        for id in range(1,bcr.shape[1]):
            sampleid=bcr.columns[id]
            all_samples.append(sampleid)
            ### Number of Reads ###
            Number_Reads_Sam=int(str( bcr[bcr['SampleID']=="Number of reads"][sampleid].iloc[0] ).replace(",",""))/1000000
            if ( Number_Reads_Sam < vdj_data_min ):
                data_error +=f"    {sampleid} 样本下机数据不足，低于最低要求 {str(int(vdj_data_min))} M，请核实该项目下机数据。\n "
    ###免疫组库###
    htmls  = html.split(",")
    for i in htmls:
        with open(i, encoding='UTF-8') as f:
            fr = f.readlines()
            line = fr[12]
            Sample_ID = re.findall(r'"Run ID":"(.*?)"', line)
            if (re.findall(r'T Cell Expression', line)) and (re.findall(r'B Cell Expression', line)):
                Num_cells_T = tcr.loc[tcr['SampleID'] == 'Estimated number of cells', Sample_ID[0] + "_TCR"].values[0]
                Num_cells_B = bcr.loc[bcr['SampleID'] == 'Estimated number of cells', Sample_ID[0] + "_BCR"].values[0]
                Num_cells_T = Num_cells_T.replace(",", "")
                Num_cells_B = Num_cells_B.replace(",", "")
                if (int(Num_cells_T) ==0) and (int(Num_cells_B) !=0):
                    Num_cells = "B Cell Expression"
                    error_flag_alerts_bool = True
                    sample_abnormity.append(Sample_ID[0] + "_TCR")
                elif (int(Num_cells_T) !=0) and (int(Num_cells_B) ==0):
                    Num_cells = "T Cell Expression"
                    error_flag_alerts_bool = True
                    sample_abnormity.append(Sample_ID[0] + "_BCR")
                elif (int(Num_cells_T) ==0) and (int(Num_cells_B) ==0):
                    Num_cells = "T Cell Expression"
                    error_flag_alerts_bool = True
                    sample_abnormity.append(Sample_ID[0] + "_BCR")
                    sample_abnormity.append(Sample_ID[0] + "_TCR")
                else:
                    Num_cells = line
            if (re.findall(r'T Cell Expression', line)) and (not re.findall(r'B Cell Expression', line)):
                Num_cells_T = tcr.loc[tcr['SampleID'] == 'Estimated number of cells', Sample_ID[0] + "_TCR"].values[0]
                Num_cells_T = Num_cells_T.replace(",", "")
                if (int(Num_cells_T) ==0):
                    Num_cells = "T Cell Expression"
                    error_flag_alerts_bool = True
                    sample_abnormity.append(Sample_ID[0] + "_TCR")
                else:
                    Num_cells = line
            if (not re.findall(r'T Cell Expression', line)) and (re.findall(r'B Cell Expression', line)):
                Num_cells_B = bcr.loc[bcr['SampleID'] == 'Estimated number of cells', Sample_ID[0] + "_BCR"].values[0]
                Num_cells_B = Num_cells_B.replace(",", "")
                if (int(Num_cells_B) ==0):
                    Num_cells = "B Cell Expression"
                    error_flag_alerts_bool = True
                    sample_abnormity.append(Sample_ID[0] + "_BCR")
                else:
                    Num_cells = line
            if (re.findall(r'T Cell Expression', Num_cells)) and (not re.findall(r'B Cell Expression', Num_cells)):
                alerts = re.findall('("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\{"level":"WARN)|("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\{"level":"ERROR")', line)
                if len(alerts) != 0:
                    error_flag_alerts_bool = True
                    sample_abnormity.append(Sample_ID[0] + "_TCR")
                else:
                    error_flag_alerts_str = " "
                    sample_normal.append(Sample_ID[0])
            if (not re.findall(r'T Cell Expression', Num_cells)) and (re.findall(r'B Cell Expression', Num_cells)):
                alerts = re.findall('("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\{"level":"WARN")|("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\{"level":"ERROR")', line)
                if len(alerts) != 0 :
                    error_flag_alerts = True
                    sample_abnormity.append(Sample_ID[0] + "_BCR")
                else:
                    error_flag_alerts =" "
                    sample_normal.append(Sample_ID[0])
            if (re.findall(r'T Cell Expression', Num_cells)) and (re.findall(r'B Cell Expression', Num_cells)):
                alerts = re.findall('("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\{"level":"WARN")|("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\{"level":"ERROR")|("title":"V\(D\)J Barcode Rank Plot"}}},"alerts":\[\]})', line)
                if (alerts[0][0] == '' and alerts[0][1] == '"title":"V(D)J Barcode Rank Plot"}}},"alerts":[{"level":"ERROR"') or (alerts[0][0] == '"title":"V(D)J Barcode Rank Plot"}}},"alerts":[{"level":"WARN"'):
                    if (alerts[1][0] == '' and alerts[1][1] == '"title":"V(D)J Barcode Rank Plot"}}},"alerts":[{"level":"ERROR"') or (alerts[1][0] == '"title":"V(D)J Barcode Rank Plot"}}},"alerts":[{"level":"WARN"'):
                        error_flag_alerts = True
                        sample_abnormity.append(Sample_ID[0] + "_BCR")
                    else:
                        error_flag_alerts =" "
                        sample_normal.append(Sample_ID[0])
                    error_flag_alerts = True
                    sample_abnormity.append(Sample_ID[0] + "_TCR")
                else:
                    if (alerts[1][0] == '' and alerts[1][1] == '"title":"V(D)J Barcode Rank Plot"}}},"alerts":[{"level":"ERROR"') or (alerts[1][0] == '"title":"V(D)J Barcode Rank Plot"}}},"alerts":[{"level":"WARN"'):
                        error_flag_alerts = True
                        sample_abnormity.append(Sample_ID[0] + "_BCR")
                    else:
                        error_flag_alerts =" "
                        sample_normal.append(Sample_ID[0])
            samples = ", ".join(sample)
    samples = ", ".join(sample)
    samples_abnormity = ", ".join(sample_abnormity)
    samples_normal = ", ".join(sample_normal)
    ### 所有样本质控结果汇总 ### 
    if global_error_flag or data_error !=" " :
        header_flag="质控异常"
        if len(samples_abnormity) ==0 :
            all_sample_text += "    免疫组库：\n    样本各质控指标正常。\n"
        if len(samples_abnormity) !=0 :
            all_sample_text += "    免疫组库：\n"
    elif global_sth_flag :
        header_flag="质控反馈"
        if len(samples_abnormity) ==0 :
            all_sample_text += "    免疫组库：\n    样本各质控指标正常。\n"
    else:
        header_flag="质控反馈"
        if int(Number_of_Reads_list[0].replace("M","")) !=33:
            all_sample_text = "    转录组：\n    样本各质控指标正常。\n    免疫组库：\n"
            if len(samples_abnormity) ==0 :
                all_sample_text += " 样本各质控指标正常。\n"
        else:
            all_sample_text = "    转录组：\n    样本各质控指标正常，建议加测。\n    免疫组库：\n"
            if len(samples_abnormity) ==0 :
                all_sample_text += " 样本各质控指标正常。\n"
    if len(samples_abnormity) !=0 :
        header_flag="质控异常"
        if len(samples_abnormity) !=0 :
            all_sample_text += f" {samples_abnormity} 样本质控异常，异常指标详见网页提示。\n"

    ## paste header and text
    header = "【" + header_flag + "】"+ Project_Num + " " + Customer + "老师 " + "单细胞" + introns + "转录组和免疫组库项目质控反馈"

    all_samples = ", ".join(all_samples )

    if cellranger_only:
        import yaml,time,sys
        header = "【" + header_flag +"和项目报告】"+ Task_Num + " " + Customer + "老师 " + "单细胞转录组及免疫组库项目质控反馈和项目报告"
        Project_Num = re.sub("-.*","",Project_Num)
        info = Project_management + '：\n    您好。\n    ' + Task_Num + " " + Customer + "老师 单细胞转录组及免疫组库项目 " + \
               str(len(list(scrna.columns[1:]))) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n" + \
               "\n    该项目为只质控不分析项目，已上传 QC 云报告。\n "
    else:
        info = Project_management + '：\n    您好。\n    ' + Project_Num + " " + Customer + "老师 单细胞" + introns + "转录组和免疫组库项目 " + str(len(list(scrna.columns[1:]))) + " 个样本 (" +  all_samples + ") 质控, " + "物种：" + Species + "，质控结果如下：\n" + all_sample_text + "    详情见附件。\n\n" 

    if remark !="":
        remark_text = "\n注：" + remark + "\n"
        info += remark_text

    # if cell_big_error_sampleid:
    #     big_samples = ",".join(cell_big_error_sampleid)
    #     cell_big_error = Laboratory + '：\n    您好。\n    如下所示：' + \
    #         big_samples + " 样本" + "细胞数高于 " + str(Est_Num_Cells_max_stand) + ", 请知悉! \n"
    #     info += cell_big_error
    if cell_big_error_sampleid:
       big_samples = "、 ".join(cell_big_error_sampleid)
       cell_big_error =   '    如下图所示：' + \
           big_samples + " 样本" + "细胞数高于 " + str(Est_Num_Cells_max_stand) + ", 还请先判断是否有force cell调整空间,再反馈项目部! \n"
       info += cell_big_error
    # if int(Number_of_Reads.replace("M", "")) == 33:
    #     test_proj_text = "\n该项目为试测10G项目。\n"
    #     info += test_proj_text

    if int(Number_of_Reads_list[0].replace("M","")) ==33:
        test_proj_text = "\n该项目为试测10G项目。\n"
        info += test_proj_text

    ### 返回 邮件相关信息 ###
    return (info,header)


# def final(header, email):
#     error_info=header
#     output = "NA"
#     send_email(email,header,output,error_info)

def send_email(emailto,header,context_html,info,scrna,tcr,bcr,outdir):
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
    info_html = info.replace('\n','<br>&emsp;')
    info_html = info_html.replace('<br>&emsp;<br>&emsp;', '<br><br>')
    # table # 
    scrna_html = scrna.to_html(index=False, border=1, justify="left",bold_rows=False)
    scrna_html = scrna_html.replace('<table border="1" class="dataframe">','<h3>转录组汇总结果:</h3> <table bordercolor="grey" border="1" cellspacing="0" cellpadding="5" class="dataframe">' )

    # message.attach(MIMEText( info_html + scrna_html  , 'html', 'utf-8'))
    email_html = info_html + scrna_html
    if type(tcr) == pd.DataFrame: 
        tcr_html= tcr.to_html(index=False, border=1, justify="left",bold_rows=False)
        tcr_html = tcr_html.replace('<table border="1" class="dataframe">','<h3>免疫组库 TCR 汇总结果:</h3> <table bordercolor="grey" border="1" cellspacing="0" cellpadding="5" class="dataframe">' )
        email_html = email_html + tcr_html

    if type(bcr) == pd.DataFrame: 
        bcr_html = bcr.to_html(index=False, border=1, justify="left",bold_rows=False)
        bcr_html = bcr_html.replace('<table border="1" class="dataframe">','<h3>免疫组库 BCR 汇总结果:</h3> <table bordercolor="grey" border="1" cellspacing="0" cellpadding="5" class="dataframe">' )
        email_html = email_html + bcr_html

    message.attach(MIMEText( email_html, 'html', 'utf-8'))
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

    # Write HTML String to file.html
    with open(f'{outdir}/cellranger_summary_email.html', "w") as file:
        file.write( email_html)

if __name__ == '__main__':
    check_dict = yaml.load(check_yaml, Loader=yaml.FullLoader)
    Estimated_Number_of_Cells = check_dict["Estimated_Number_of_Cells"]
    Q30_Base_in_RNA_Read = check_dict["Q30_Base_in_RNA_Read"]
    Fraction_Reads_in_Cells = check_dict["Fraction_Reads_in_Cells"]
    Median_Genes_per_Cell = check_dict["Median_Genes_per_Cell"]
    Reads_Mapped_Confidently_to_Transcriptome = check_dict["Reads_Mapped_Confidently_to_Transcriptome"]
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    config = {}
    with open(inconfig,'r',encoding='utf-8') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    Project_Num = config["report"]["Project_Num"]
    Task_Num = config["report"]["Task_Num"]
    Customer = config["report"]["Customer"]
    Species = config["report"]["Species"]
    Project_management =  config["email"]["Project_management"]
    Laboratory = config["email"]["Laboratory"]
    remark = config["email"]["remark"]
    ## 自动获取邮件收件人的邮箱信息
    if config["email"]["emailto"] is None or config["email"]["emailto"][0] == "":
        if 'CLOUDUSER' in os.environ: 
            emailto = list() 
            emailto.append(os.environ['CLOUDUSER'])
            print(f"The value of CLOUDUSER is {emailto[0]}")
        else:
            print("There's no CLOUDUSER variable in the current environment, please set CLOUDUSER on .bashrc or email to on config.yaml ")
    else:
        emailto = config["email"]["emailto"]
    include_introns =  config["cellranger_params"]["include_introns"]
    Number_of_Reads = config["cellranger_params"]["Number_of_Reads_for_mutil"]
    cellranger_only = config["cellranger_params"]["module"]["cellranger_report"]
    
    ## cellranger summary clean 
    (scrna,tcr,bcr) = cellranger_summary(stat, outdir) 

    ## cellranger mutil qc result check 
    (info,header)=cellranger_check(scrna,cellranger_only,tcr,bcr,Estimated_Number_of_Cells,Q30_Base_in_RNA_Read,Fraction_Reads_in_Cells,Median_Genes_per_Cell,Number_of_Reads,Reads_Mapped_Confidently_to_Transcriptome,Project_Num,Species,Customer,Project_management,Laboratory,remark,include_introns)

    ## email send 
    send_email(emailto,header,html,info,scrna,tcr,bcr,outdir)


