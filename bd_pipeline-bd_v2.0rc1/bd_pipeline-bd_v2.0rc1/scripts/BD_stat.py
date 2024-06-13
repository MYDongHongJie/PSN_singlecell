import argparse
import pandas as pd
import os

# Assume your data is in test.txt in the current working directory 
parser = argparse.ArgumentParser(description="Sample stats combine script for BD pipline")
parser.add_argument("-i", "--input_path", help="Input with Metrics_Summary.csv",
                    default="sample_Metrics_Summary.csv")
parser.add_argument("-o", "--output", help="sampleid for results.",
                    default="outputdir")
parser.add_argument("-t", "--type", help="Customized for AbSeq",
                    default=None)
args = parser.parse_args()


input = os.path.abspath(args.input_path)
output = str(args.output)
data_type = str(args.type)
df = pd.read_csv(input, names=range(20),comment='#',sep=",")
df = df.dropna(axis=1,how='all')

if data_type !="None":
    df_vdj = None
    df_FQ = df.iloc[[0,1,2]].dropna(axis=1,how='all')[[0,1,5]]
    df_QC = df.iloc[[4,5,6]].dropna(axis=1,how='all')[[0,1,2,3]]
    df_QC.index = [0,1,2]
    df_align = df.iloc[[8,9,10]].dropna(axis=1,how='all')[df.columns[0:4]]
    df_align.index = [0,1,2]
    df_cell = df.iloc[[15,16,17]].dropna(axis=1,how='all')[[0,7]]
    df_cell.index = [0,1,2]
    merge_df = pd.concat([df_cell,df_vdj,df_FQ,df_QC,df_align],axis=1)
    merge_df.iloc[0] = ("Putative_Cell_Count","Median_Genes_per_Cell","Total_Reads_in_FASTQ","Pct_Read_Pair_Overlap","Pct_Reads_Filtered_Out","Total_Filtered_Reads","Pct_Q30_Bases_in_Filtered_R2","Pct_Assigned_to_Cell_Labels","Pct_Cellular_Reads_Aligned_Uniquely","Cellular_Reads","Annotated_Transcriptome_Pct","Introns_Pct","Intergenic_Regions_Pct")
else :
    if len(df) > 18:
        vdj_Mean = df.iloc[19:22].dropna(axis=1,how='all')[[13,7]].T[[20,21]]
        vdj_Mean = vdj_Mean.replace('BCR',"MeanMCP_BCR").replace("TCR","MeanMCP_TCR")
        vdj_Mean.index = [0,1]
        vdj_pct_info = df.iloc[12:14].dropna(axis=1,how='all')[[11,10]]
        vdj_pct = pd.DataFrame({"Pct_BCR":["Pct_BCR",vdj_pct_info.iloc[0,0]],"Pct_TCR":["Pct_TCR",vdj_pct_info.iloc[1,1]]})
        vdj_pct.index = [0,1]
        df_vdj = pd.concat((vdj_Mean,vdj_pct),axis=1)
        df_FQ = df.iloc[[0,4]].dropna(axis=1,how='all')[[0,1,5]]
        df_FQ.index = [0,1]
        df_QC = df.iloc[[5,9]].dropna(axis=1,how='all')[[0,1,2,3]]
        df_QC.index = [0,1]
        df_align = df.iloc[[10,14]].dropna(axis=1,how='all')[df.columns[0:4]]
        df_align.index = [0,1]
        df_cell = df.iloc[[17,18]].dropna(axis=1,how='all')[[0,7]]
        df_cell.index = [0,1]
        merge_df = pd.concat([df_cell,df_vdj,df_FQ,df_QC,df_align],axis=1)
        merge_df.iloc[0] = ("Putative_Cell_Count","Median_Genes_per_Cell","MeanMCP_BCR","MeanMCP_TCR","Pct_BCR","Pct_TCR","Total_Reads_in_FASTQ","Pct_Read_Pair_Overlap","Pct_Reads_Filtered_Out","Total_Filtered_Reads","Pct_Q30_Bases_in_Filtered_R2","Pct_Assigned_to_Cell_Labels","Pct_Cellular_Reads_Aligned_Uniquely","Cellular_Reads","Annotated_Transcriptome_Pct","Introns_Pct","Intergenic_Regions_Pct")
    elif len(df) == 18:
        VDJ_type = df.iloc[2][7].split('_')[-1]
        vdj_Mean = df.iloc[17:18].dropna(axis=1,how='all')[[7]]
        vdj_Mean.columns = ["MeanMCP_" + VDJ_type]
        if VDJ_type =="TCR":
            vdj_Mean.insert(column="MeanMCP_BCR",value="0",loc = 1)
        elif VDJ_type =="BCR":
            vdj_Mean.insert(column="MeanMCP_TCR",value="0",loc = 1)
        # vdj_Mean = vdj_Mean.replace('BCR',"MeanMCP_BCR").replace("TCR","MeanMCP_TCR")
        vdj_Mean = vdj_Mean[['MeanMCP_BCR','MeanMCP_TCR']]
        new_term = vdj_Mean.columns.tolist()
        vdj_Mean.loc[-1] = new_term
        vdj_Mean.index = vdj_Mean.index+1
        vdj_Mean = vdj_Mean.sort_index()
        vdj_Mean.index = [0,1]
        vdj_pct_info = df.iloc[8:11].dropna(axis=1,how='all')[[11,10]].drop(9)
        vdj_pct_info.index = [0,1]
        vdj_pct_info = vdj_pct_info.replace('VDJ_BCR_Pct',"Pct_BCR").replace("VDJ_TCR_Pct","Pct_TCR")
        vdj_pct_info.columns = vdj_pct_info.loc[0].tolist()
        #vdj_pct = pd.DataFrame({"Pct_BCR":["Pct_BCR",vdj_pct_info.iloc[0,0]],"Pct_TCR":["Pct_TCR",vdj_pct_info.iloc[1,1]]})
        vdj_pct = vdj_pct_info
        df_vdj = pd.concat((vdj_Mean,vdj_pct),axis=1)
        df_FQ = df.iloc[[0,3]].dropna(axis=1,how='all')[[0,1,5]]
        df_FQ.index = [0,1]
        df_QC = df.iloc[[4,7]].dropna(axis=1,how='all')[[0,1,2,3]]
        df_QC.index = [0,1]
        df_align = df.iloc[[8,11]].dropna(axis=1,how='all')[df.columns[0:4]]
        df_align.index = [0,1]
        df_cell = df.iloc[[14,15]].dropna(axis=1,how='all')[[0,7]]
        df_cell.index = [0,1]
        merge_df = pd.concat([df_cell,df_vdj,df_FQ,df_QC,df_align],axis=1)
        merge_df.iloc[0] = ("Putative_Cell_Count","Median_Genes_per_Cell","MeanMCP_BCR","MeanMCP_TCR","Pct_BCR","Pct_TCR","Total_Reads_in_FASTQ","Pct_Read_Pair_Overlap","Pct_Reads_Filtered_Out","Total_Filtered_Reads","Pct_Q30_Bases_in_Filtered_R2","Pct_Assigned_to_Cell_Labels","Pct_Cellular_Reads_Aligned_Uniquely","Cellular_Reads","Annotated_Transcriptome_Pct","Introns_Pct","Intergenic_Regions_Pct")
    else:
        df_vdj = None
        df_FQ = df.iloc[[0,1]].dropna(axis=1,how='all')[[0,1,5]]
        df_QC = df.iloc[[2,3]].dropna(axis=1,how='all')[[0,1,2,3]]
        df_QC.index = [0,1]
        df_align = df.iloc[[4,5]].dropna(axis=1,how='all')[df.columns[0:4]]
        df_align.index = [0,1]
        df_cell = df.iloc[[8,9]].dropna(axis=1,how='all')[[0,7]]
        df_cell.index = [0,1]
        merge_df = pd.concat([df_cell,df_vdj,df_FQ,df_QC,df_align],axis=1)
        merge_df.iloc[0] = ("Putative_Cell_Count","Median_Genes_per_Cell","Total_Reads_in_FASTQ","Pct_Read_Pair_Overlap","Pct_Reads_Filtered_Out","Total_Filtered_Reads","Pct_Q30_Bases_in_Filtered_R2","Pct_Assigned_to_Cell_Labels","Pct_Cellular_Reads_Aligned_Uniquely","Cellular_Reads","Annotated_Transcriptome_Pct","Introns_Pct","Intergenic_Regions_Pct")
merge_df.to_csv(output,header=0,index = 0)
