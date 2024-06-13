import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description="split_violin plot")
parser.add_argument("-i", "--input", help="group_exp file")
parser.add_argument("-g", '--gene', help="file of the  genes and terms")
parser.add_argument("-o", '--output', help="outputdir",
                    default="split_violin")
args = parser.parse_args()

input = os.path.abspath(args.input)
gene_group = os.path.abspath(args.gene)
output = os.path.abspath(args.output)

df = pd.read_csv(input,sep="\t")
gene_term = pd.read_csv(gene_group,sep="\t")

#gene_names = df.columns.tolist()[3:119]
term_list = gene_term.columns.tolist()

for term in term_list:
    if os.path.exists(f'{output}/{term}/'):
            print(f"\033[0;33mWaring: Folder {output}/{term} exists. Skipping...\033[0m")
    else:
        os.mkdir(f'{output}/{term}/')
    genelist = gene_term[term]
    genelist = genelist.dropna(axis=0,how='all').tolist()
    for gene in genelist:
        print(gene)
        df_1 = df[[gene,'group','clusters']]
        df_1.sort_values('group', ascending=False,inplace=True)
        #g = sns.catplot(x="clusters",y=gene,hue="group",kind="violin", split=True, data=df_1)
        g = sns.catplot(x="clusters",y=gene,hue="group",kind="violin", split=True, data=df_1,
                    palette={'WT_6':'#E6C013', 'KO_15': '#0672BB'},inner = "quartile")
        plt.xticks(rotation=30,horizontalalignment='right',fontsize=13)
        g.fig.set_size_inches(7,5)
        g.set_axis_labels('clusters',gene+" data")
        plt.savefig(output + "/" + term + '/' + gene + '.data.split_volin.pdf', orientation = 'horizontal')