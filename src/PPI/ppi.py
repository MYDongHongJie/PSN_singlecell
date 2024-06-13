#!/usr/bin/env python3
# encoding: utf-8
"""
Project : Visium
Author  : Xiufeng Yang
Contact : xiufeng.yang@oebiotech.com
File   : ppi.py
IDE    : PyCharm
Time   : 2020-05-08 15:17:42
Desc   :
"""
# =======================================================================================================================
# =======================================================================================================================
import requests
import pandas as pd
from pandas.core.frame import DataFrame
from time import sleep
import click
import os
import subprocess
import cairosvg
import math 
import numpy as np

def string_switch(x, y, z, s=1):
    with open(x, "r", encoding="utf-8") as f:
        # readlines以列表的形式将文件读出
        lines = f.readlines()
    with open(x, "w", encoding="utf-8") as f_w:
        # 定义一个数字，用来记录在读取文件时在列表中的位置
        n = 0
        # 默认选项，只替换第一次匹配到的行中的字符串
        if s == 1:
            for line in lines:
                if y in line:
                    line = line.replace(y, z)
                    f_w.write(line)
                    n += 1
                    break
                f_w.write(line)
                n += 1
            # 将剩余的文本内容继续输出
            for i in range(n, len(lines)):
                f_w.write(lines[i])
        # 全局匹配替换
        elif s == 'g':
            for line in lines:
                if y in line:
                    line = line.replace(y, z)
                f_w.write(line)


def mapping(mygene, species, outname):
    """
    ## For a given list of proteins the script resolves them (if possible) to the best matching STRING identifier
    ## and prints out the mapping on screen in the TSV format
    :param mygene: input gene list
    :param species: species: NCBI taxon identifiers
    :param outname: outputfile name
    :return: mapping results
    """
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    if len(mygene) > 2000:
        mapping_results=DataFrame()
        for index in range(1,math.ceil(len(mygene)/2000)+1):
            if index==1:
                mygene_subset=mygene[index:index*2000+1]
            elif index==math.ceil(len(mygene)/2000):
                mygene_subset=mygene[2000*(index-1):len(mygene)]
            else:
                mygene_subset=mygene[2000*(index-1):index*2000+1]
            params = {
                "identifiers": "\r".join(mygene_subset),  # your protein list
                "species": species,  # species NCBI identifier
                "limit": 1,  # only one (best) identifier per input protein
                "echo_query": 1,  # see your input identifiers in the output
                "caller_identity": "www.awesome_app.org"  # your app name
            }
            ## Construct URL
            request_url = "/".join([string_api_url, output_format, method])
            ## Call STRING
            results_subset = requests.post(request_url, data=params)
            mapping_list_subset = [x.split("\t") for x in results_subset.text.strip().split("\n")]
            mapping_results_subset = DataFrame(mapping_list_subset)
            mapping_results=mapping_results.append(mapping_results_subset)
    else:
        params = {
                "identifiers": "\r".join(mygene),  # your protein list
                "species": species,  # species NCBI identifier
                "limit": 1,  # only one (best) identifier per input protein
                "echo_query": 1,  # see your input identifiers in the output
                "caller_identity": "www.awesome_app.org"  # your app name
        }
        ## Construct URL
        request_url = "/".join([string_api_url, output_format, method])
        ## Call STRING
        results = requests.post(request_url, data=params)
        
        ## Read and parse the results
        mapping_list = [x.split("\t") for x in results.text.strip().split("\n")]
        mapping_results = DataFrame(mapping_list)
    
    mapping_results.to_csv(outname, sep='\t', header=False, index=False)
    return mapping_results


def get_network(mapping_results, species, outname):
    """

    :param mapping_results:
    :param species:
    :param outname:
    :return:
    """
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "network"
    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    mygene = list(mapping_results[0])
    ## For all gene call STRING
    if len(mygene) > 2000:
        results_data=DataFrame()
        for index in range(1,math.ceil(len(mygene)/2000)+1):
            if index==1:
                mygene_subset=mygene[index:index*2000+1]
            elif index==math.ceil(len(mygene)/2000):
                mygene_subset=mygene[2000*(index-1):len(mygene)]
            else:
                mygene_subset=mygene[2000*(index-1):index*2000+1]
            
            params = {
                "identifiers": "%0d".join(mygene_subset),  # your protein
                "species": species,  # species NCBI identifier
                "caller_identity": "www.awesome_app.org"  # your app name
            }
            ## Construct URL
            results = requests.post(request_url, data=params)
            ## Save the network to file
            results_list_subset = [x.split("\t") for x in results.text.strip().split("\n")]
            results_data_subset = DataFrame(results_list_subset)
            results_data = results_data.append(results_data_subset)
    else:
        params = {
            "identifiers": "%0d".join(mygene),  # your protein
            "species": species,  # species NCBI identifier
            "caller_identity": "www.awesome_app.org"  # your app name
        }
        ## Call STRING
        results = requests.post(request_url, data=params)
        ## Save the network to file
        results_list = [x.split("\t") for x in results.text.strip().split("\n")]
        results_data = DataFrame(results_list)
    
    results_data.drop_duplicates(inplace=True)
    results_data.to_csv(outname, sep='\t', header=False, index=False)
    return results_data

def get_network_image(mapping_results, species, outname):
    """
    ## For each protein in a list save the PNG image of
    ## STRING network of its interaction partners.
    :param mapping_results:
    :param species:
    :param outname:
    :return:
    """
    string_api_url = "https://string-db.org/api"
    output_format = "svg"
    method = "network"
    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    mygene = list(mapping_results['gene'])
    ## For all gene call STRING
    params = {
        "identifiers": "%0d".join(mygene),  # your protein
        "species": species,  # species NCBI identifier
        "hide_disconnected_nodes": 1,
        "network_flavor": "confidence",  # show confidence links
        "caller_identity": "www.awesome_app.org"  # your app name
    }
    ## Call STRING
    response = requests.post(request_url, data=params)
    ## Save the network to file
    file_name = outname
    print("Saving interaction network to %s" % file_name)

    with open(file_name, 'wb') as fh:
        fh.write(response.content)
    read_list = []
    with open(file_name, mode="r") as f:
        read_s = f.readlines()
    txt = """
    <rect width="100%" height="100%"
        style="fill:rgb(255,255,255);stroke-width:1;stroke:rgb(0,0,0)"
    />
    """
    read_s.insert(1, txt)
    s = "".join(read_s)
    with open(file_name, mode="w") as f:
        f.write(s)
    sleep(1)


@click.command()
@click.option("--inputfile", "-i", type=click.Path("r"), help="input files with gene list")
@click.option("--species", "-s", type=click.INT,
              help="NCBI taxon identifiers (e.g. Human is 9606, Mouse is 10090, Rat is 10116),See: STRING organisms(https://string-db.org/cgi/input.pl?input_page_active_form=organisms).")
@click.option("--prefix", "-p", type=click.STRING, help="prefix names for output file.")
@click.option("--noft", "-n", type=click.INT, default = 25, help="no. of tops.")
@click.option("--outputdir", "-o", type=click.Path(exists=False), help="output directory.")
def main(inputfile, species, prefix, noft, outputdir):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    inputfile = os.path.realpath(inputfile)
    outputdir = os.path.realpath(outputdir)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    os.chdir(str(outputdir))

    diffgene = pd.read_table(inputfile)
    if "gene_name" in diffgene.columns:
        mygene = list(diffgene["gene_name"])
    else:
        mygene = list(diffgene["gene"])
    mapping_result = mapping(mygene, species, prefix + ".mapping_results.tsv")
    results_data = get_network(mapping_result, species, prefix + ".string_protein-protein-interaction.tsv")
    
    if results_data.shape[0]>1:
        if "up_down" in diffgene.columns:
            uptop_diff = diffgene[diffgene['up_down']=="Up"].sort_values(by="FoldChange",ascending=False)[0:noft]
            downtop_diff = diffgene[diffgene['up_down']=="Down"].sort_values(by="FoldChange",ascending=True)[0:noft]
            top_diff = uptop_diff.append(downtop_diff)
            get_network_image(top_diff, species, prefix + ".string_protein-protein-interaction.svg")
        else:
            get_network_image(diffgene, species, prefix + ".string_protein-protein-interaction.svg")
        cairosvg.svg2png(url=prefix + ".string_protein-protein-interaction.svg",
                    write_to=prefix + ".string_protein-protein-interaction.png",dpi=1000)
        cairosvg.svg2pdf(url=prefix + ".string_protein-protein-interaction.svg",
                    write_to=prefix + ".string_protein-protein-interaction.pdf")
        ###################################################################################################################
        # change networks color by up_down information
        if "up_down" in diffgene.columns:
            svg_files = prefix + ".string_protein-protein-interaction.svg"
            cmd = f"python {scriptdir}/change_STRING_colors.py -s {svg_files}"
            subprocess.run(cmd, shell=True, check=True)
            # color_gene = pd.read_table("color_table.tsv",header=None)
            # color_gene.columns = ["gene", "color"]
            # new_color_gene = pd.merge(diffgene, color_gene, on=["gene"], how='outer')
            genelist = pd.DataFrame(mapping_result[[0, 5]])
            genelist.columns = ["gene", "gene_color"]
            new_color_gene = pd.merge(diffgene, genelist, on=["gene"], how='outer')

            def new_color(row):
                if row["up_down"] == "Up":
                    color = "rgb(255,0,0)"
                elif row["up_down"] == "Down":
                    color = "rgb(0,255,0)"
                else:
                    color = "Na"
                return color

            new_color_gene["color"] = new_color_gene.apply(lambda x: new_color(x), axis=1)
            new_color = pd.DataFrame(new_color_gene[['gene_color', 'color']])
            new_color.to_csv("color_table.tsv", sep='\t', header=False, index=False)
            cmd = f"python {scriptdir}/change_STRING_colors.py -s {svg_files} -c color_table.tsv "
            subprocess.run(cmd, shell=True, check=True)
            ##covert image format into pdf and png
            cairosvg.svg2png(url=prefix + ".string_protein-protein-interaction.new_colors.svg",
                             write_to=prefix + ".string_protein-protein-interaction.new_colors.png",
                             dpi=1000)
            cairosvg.svg2pdf(url=prefix + ".string_protein-protein-interaction.new_colors.svg",
                             write_to=prefix + ".string_protein-protein-interaction.new_colors.pdf")
            # convert -background white 1.png 2.png


if __name__ == "__main__":
    main()
