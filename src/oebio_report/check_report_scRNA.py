# -*- coding: utf-8 -*-
import click
import os
import shutil
import re
import glob
import pandas as pd

@click.command()
@click.option('-i','--input', prompt='Report directory',
              help='The directory of report')
@click.option('-n','--samplenum', prompt='Sample number',type= int,
              help='The number of samples')
@click.option('-s','--species', prompt='species',
              help='The specie where the samples come from. e.g. human, mouse, other')
@click.option('-d','--diffgroup', prompt='Diff group number',type= int,
              help='The number of diff group')
@click.option('-o','--out', prompt='Output directory',
             help='The output directory of "Warning.txt"')

def check(input,samplenum,species,diffgroup,out):

    """
    python check_report_scRNA.py -i HTXXXX_report/ -n 6 -s human -d 2 -o ./

    """

    warn = open('%s/Warning.txt' %(out),'w')
    w = "该报告可能存在如下问题 T^T ，请进一步核查！！\n"
    n = 0
    ### cellranger
    f = glob.glob(input + '/*.CellRanger')
    if f:
        if samplenum == 1:
            if len(os.listdir(f[0])) == 2 :
                print("CellRanger √")
            else:
                print("× CellRanger 目录下 %s 个文件, 应有 2 个文件，请核查!!!" %(len(os.listdir(f[0]))) )
                n += 1
                w += str(n) + ". CellRanger 目录下 %s 个文件, 应有 2 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])))
        else:
            if len(os.listdir(f[0])) == samplenum*2 + 1 :
                print("CellRanger √")
            else:
                print("× CellRanger 目录下 %s 个文件, 应有 %s 个文件，请核查!!!" %(len(os.listdir(f[0])), samplenum*2 + 1) )
                n += 1
                w += str(n) + ". CellRanger 目录下 %s 个文件, 应有 %s 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])), samplenum*2 + 1)
    else:
        print("× CellRanger 目录不存在，请核查!!!")
        n += 1
        w += str(n) + ". CellRanger 目录不存在！\n"
    ### Count_QC
    f = glob.glob(input + '/*.Count_QC')
    if f:
        count_qc_tmp = 0
        if len(os.listdir(f[0])) == 15 :
            qc_stat = pd.read_csv('%s/cell_statitics_before_after_QC.xls' %(f[0]),index_col=0,sep='\t')
            sample_num = len(qc_stat)   #样本数量
            if sample_num != samplenum:
                print("× Count_QC 中有 %s 个样本，应有 %s 个样本，请核查!!!" %(sample_num, samplenum))
                count_qc_tmp += 1
                n += 1
                w += str(n) + ". Count_QC 中有 %s 个样本，应有 %s 个样本，数目错误请核查！\n" %(sample_num, samplenum)
            else:
                if min(qc_stat['Total_cells_afterQC'])/max(qc_stat['Total_cells_afterQC']) <= 0.5:
                    print("× Count_QC 中存在样本间细胞数相差>=2倍，请注意。")
                    count_qc_tmp += 1
                    n += 1
                    w += str(n) + ". Count_QC 中存在样本间细胞数相差>=2倍，请注意。\n"
                for index,row in qc_stat.iterrows():
                    if (row["Total_cells_beforeQC"] - row["Total_cells_afterQC"]) / row["Total_cells_beforeQC"] > 0.4 :
                        print("× Count_QC 中 %s 样本过滤掉的细胞数目超过 40%%，请注意。" %(index))
                        count_qc_tmp += 1
                        n += 1
                        w += str(n) + ". Count_QC 中 %s 样本过滤掉的细胞数目超过 40%%，请注意。\n" %(index)
        else:
            print("× Count_QC 目录下 %s 个文件, 应有 15 个文件，请核查!!!" %(len(os.listdir(f[0]))) )
            n += 1
            w += str(n) + ". Count_QC 目录下 %s 个文件, 应有 15 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])))
            count_qc_tmp += 1
        if count_qc_tmp == 0:
            print("Count_QC √")
    else:
        print("× Count_QC 目录不存在，请核查!!!")
        n += 1
        w += str(n) + ". Count_QC 目录不存在！\n"
    ### Clustering
    f = glob.glob(input + '/*.Clustering')
    if f:
        clustering_correct = 4
        if samplenum > 1:
            clustering_correct += 1
        if species != "human" and species != "mouse": 
            clustering_correct += 1
        if len(os.listdir(f[0])) == clustering_correct:
            print("Clustering √")
        else:
            print("× Clustering 目录下 %s 个文件, 应有 %s 个文件，请核查!!!" %(len(os.listdir(f[0])),clustering_correct) )
            n += 1
            w += str(n) + ". Clustering 目录下 %s 个文件, 应有 %s 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])),clustering_correct)
    else:
        print("× Clustering 目录不存在，请核查!!!")
        n += 1
        w += str(n) + ". Clustering 目录不存在！\n"
    ### Clusters_correlation
    f = glob.glob(input + '/*.Clusters_correlation')
    if f:
        correlation_results = pd.read_csv('%s/normalized_data_groupby_clusters.xls' %(f[0]), sep='\t')
        correlation_cluster_num = correlation_results.columns.size-1
        if len(os.listdir(f[0])) == 4 :
            print("Clusters_correlation √")
        else:
            print("× Clusters_correlation 目录下 %s 个文件, 应有 4 个文件，请核查!!!" %(len(os.listdir(f[0]))) )
            n += 1
            w += str(n) + ". Clusters_correlation 目录下 %s 个文件, 应有 4 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])))
    else:
        print("× Clusters_correlation 目录不存在，请核查!!!")
        n += 1
        w += str(n) + ". Clusters_correlation 目录不存在！\n"
    ### Marker
    f = glob.glob(input + '/*.Marker')
    if f:
        Marker_results = pd.read_csv('%s/top10_markers_for_each_cluster_anno.xls' %(f[0]), sep='\t')
        marker_cluster_num = max(Marker_results['cluster'])
        if len(os.listdir(f[0])) == 5 :
            if 'correlation_cluster_num' in locals() :
                if correlation_cluster_num != marker_cluster_num:
                    print("× Marker 中的细胞群数与 Clusters_correlation 中的细胞群数不一致，请核查!!!")
                    n += 1
                    w += str(n) + ". Marker 中的细胞群数与 Clusters_correlation 中的细胞群数不一致，数目错误请核查！\n"
                else:
                    print("Marker √")
            else:
                print("? Marker 目录下文件数一致，但未知与 Clusters_correlation 中的细胞群数关系，请先确定 Clusters_correlation 群数后再次核查！")
                n += 1
                w += str(n) + ". Marker 目录下文件数一致，但未知与 Clusters_correlation 中的细胞群数关系，请先确定 Clusters_correlation 群数后再次核查！\n"
        else:
            print("× Marker 目录下 %s 个文件, 应有 5 个文件，请核查!!!" %(len(os.listdir(f[0]))) )
            n += 1
            w += str(n) + ". Marker 目录下 %s 个文件, 应有 5 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])))
    else:
        print("× Marker 目录不存在，请核查!!!")
        n += 1
        w += str(n) + ". Marker 目录不存在！\n"

    ### Reference_celltype
    f = glob.glob(input + '/*.Reference_celltype')
    if species == "human" or species == "mouse":
        if f:
            if len(os.listdir(f[0])) == 9 :
                celltyping_file = glob.glob(f[0] +'/*_resolution*simplified_celltyping_results.csv')
                celltyping = pd.read_csv('%s' %(celltyping_file[0]), sep=',')
                celltyping_cluster_num = max(celltyping['clusters'])
                if 'marker_cluster_num' in locals():
                    if celltyping_cluster_num != marker_cluster_num:
                        print("× Reference_celltype 中的细胞群数与 Marker 中的细胞群数不一致，请核查!!!")
                        n += 1
                        w += str(n) + ". Reference_celltype 中的细胞群数与 Marker 中的细胞群数不一致，数目错误请核查！\n"
                    else:
                        if glob.glob(f[0] +'/*_top*_celltyping_on_*_resolution*_plot.pdf')[0].split('/')[-1].startswith(species.capitalize()):
                            print("Reference_celltype √")
                        else:
                            print("× Reference_celltype 中物种选用错误，应为 %s ，请核查!!!" %(species.capitalize()))
                            n += 1
                            w += str(n) + ". Reference_celltype 中物种选用错误，应为 %s ，请核查！\n" %(species.capitalize())
                else:
                    print("? Reference_celltype 目录下文件数一致，但未知与 Marker 中的细胞群数关系，请先确定 Marker 群数后再次核查!!！")
                    n += 1
                    w += str(n) + ". Reference_celltype 目录下文件数一致，但未知与 Marker 中的细胞群数关系，请先确定 Marker 群数后再次核查!\n"
            else:
                print("× Reference_celltype 目录下 %s 个文件, 应有 9 个文件，请核查!!!" %(len(os.listdir(f[0]))) )
                n += 1
                w += str(n) + ". Reference_celltype 目录下 %s 个文件, 应有 9 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])) )
        else:
            print("× Reference_celltype 目录不存在，请核查!!!")
            n += 1
            w += str(n) + ". Reference_celltype 目录不存在！\n"
    else:
        print("特殊物种无细胞类型鉴定结果。")
    ### Diffexp
    f = glob.glob(input + '/*.Diffexp')
    if diffgroup != 0:
        if f:
            if len(os.listdir(f[0])) == diffgroup*5 +1 :
                print("Diffexp √")
            else:
                print("× Diffexp 目录下 %s 个文件, 应有 %s 个文件，请核查!!!" %(len(os.listdir(f[0])), diffgroup*5 +1) )
                n += 1
                w += str(n) + ". Diffexp 目录下 %s 个文件, 应有 %s 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])), diffgroup*5 +1)
        else:
            print("× Diffexp 目录不存在，请核查!!!")
            n += 1
            w += str(n) + ". Diffexp 目录不存在！\n"
    else:
        print("无差异比较组。")
    ### enrichment
    f = glob.glob(input + '/*.enrichment')
    if diffgroup != 0:
        f = glob.glob(input + '/*.enrichment')
        if f:
            if "GO_enrichment" in os.listdir(f[0]):
                tmp_go = 0
                if len(os.listdir(f[0] + "/GO_enrichment")) == diffgroup + 1:
                    for i in glob.glob(f[0]+ "/GO_enrichment/*/"):
                        if len(os.listdir(i)) != 45:
                            tmp_go +=1
                            print("× GO_enrichment 目录下 %s 分组有 %s 个文件，应有 45 个文件! 请核查!!!" %(i.split("GO_enrichment/")[-1].replace("/",""), len(os.listdir(i))))
                            n += 1
                            w += str(n) + ". GO_enrichment 目录下 %s 分组有 %s 个文件，应有 45 个文件，数目错误核查！\n" %(i.split("GO_enrichment/")[-1].replace("/",""), len(os.listdir(i)))
                else:
                    tmp_go +=1
                    print("× GO_enrichment 目录下 %s 个文件，应有 %s 个文件! 请核查!!!" %(len(os.listdir(f[0] + "/GO_enrichment")),diffgroup + 1))
                    n += 1
                    w += str(n) + ". GO_enrichment 目录下 %s 个文件，应有 %s 个文件，数目错误核查！\n" %(len(os.listdir(f[0] + "/GO_enrichment")),diffgroup + 1)
            else:
                tmp_go +=1
                print("× enrichment 目录下缺少 GO_enrichment! 请核查!!!")
                n += 1
                w += str(n) + ". enrichment 目录下缺少 GO_enrichment！\n"
            if tmp_go == 0 :
                print("GO_enrichment √")

            if "KEGG_enrichment" in os.listdir(f[0]):
                tmp_kegg = 0
                if len(os.listdir(f[0] + "/KEGG_enrichment")) == diffgroup + 1:
                    for i in glob.glob(f[0] + "/KEGG_enrichment/*-vs-*/"):
                        if len(os.listdir(i)) != 27:
                            tmp_kegg += 1
                            print("× Warning: KEGG_enrichment 目录下 %s 分组有 %s 个文件，应有 27 个文件! 请核查!!!" % (
                            i.split("KEGG_enrichment/")[-1].replace("/", ""), len(os.listdir(i))))
                            n += 1
                            w += str(n) + ". KEGG_enrichment 目录下 %s 分组有 %s 个文件，应有 27 个文件，数目错误核查！\n" % (
                            i.split("KEGG_enrichment/")[-1].replace("/", ""), len(os.listdir(i)))
                else:
                    tmp_kegg +=1
                    print("× KEGG_enrichment 目录下 %s 个文件，应有 %s 个文件! 请核查!!!" %(len(os.listdir(f[0] + "/KEGG_enrichment")),diffgroup + 1))
                    n += 1
                    w += str(n) + ". KEGG_enrichment 目录下 %s 个文件，应有 %s 个文件，数目错误核查！\n" %(len(os.listdir(f[0] + "/KEGG_enrichment")),diffgroup + 1)
            else:
                tmp_kegg += 1
                print("× enrichment 目录下缺少 KEGG_enrichment! 请核查!!!")
                n += 1
                w += str(n) + ". enrichment 目录下缺少 KEGG_enrichment！\n"
            if tmp_kegg == 0:
                print("KEGG_enrichment √")
        else:
            print("× enrichment 目录不存在，请核查!!!")
            n += 1
            w += str(n) + ". enrichment 目录不存在！\n"
    else:
        print("无差异富集分析。")
    ### PPI
    f = glob.glob(input + '/*.PPI')
    if diffgroup != 0:
        if f:
            if len(os.listdir(f[0])) == diffgroup*7 :
                print("PPI √")
            else:
                print("× PPI 目录下 %s 个文件, 应有 %s 个文件，请核查!!!" %(len(os.listdir(f[0])), diffgroup*7) )
                n += 1
                w += str(n) + ". PPI 目录下 %s 个文件, 应有 %s 个文件，数目错误请核查！\n" %(len(os.listdir(f[0])), diffgroup*7)
        else:
            print("× PPI 目录不存在，请核查!!!")
            n += 1
            w += str(n) + ". PPI 目录不存在！\n"
    else:
        print("无差异基因PPI分析。")
        
    warn.write(w)
    warn.close()
    if n == 0:
        os.remove('%s/Warning.txt' %(out))
        print("Good Luck! 报告检查通过，请人工检查确认无误后发送报告！")
    else:
        print("\n ！！！该报告可能存在 %s 个问题，请查看 Warning.txt，进一步核查补充完整后再发送！！！" %(n))

if __name__ == '__main__':
    check()
