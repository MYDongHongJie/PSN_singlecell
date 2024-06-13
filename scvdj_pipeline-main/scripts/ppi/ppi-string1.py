#!/usr/bin/env python3
# encoding: utf-8
"""
Project : Visium
Author  : Xiufeng Yang
Contact : xiufeng.yang@oebiotech.com
File   : ppi-string1.py
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


def mapping(mygene, species):
    """
    ## For a given list of proteins the script resolves them (if possible) to the best matching STRING identifier
    ## and prints out the mapping on screen in the TSV format
    :param mygene: input gene list
    :param species: species: NCBI taxon identifiers
    :return: mapping results
    """
    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
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
    results_data.to_csv(outname, sep='\t', header=False, index=True)


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
    mygene = list(mapping_results[0])
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
    sleep(1)


@click.command()
@click.option("--inputfile", "-i", type=click.Path(exists=True, resolve_path=True), help="input files with gene list")
@click.option("--species", "-s", type=click.INT,
        help="NCBI taxon identifiers (e.g. Human is 9606, see: STRING organisms(https://string-db.org/cgi/input.pl?input_page_active_form=organisms).")
@click.option("--prefix", "-p", type=click.STRING,help="prefix names for output file.")
@click.option("--outputdir", "-o", type=click.Path(exists=False, resolve_path=True), help="output directory.")
def main(inputfile, species, prefix, outputdir):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    os.chdir(str(outputdir))
    diffgene = pd.read_table(inputfile)
    mygene = list(diffgene["gene"])
    mapping_result = mapping(mygene, species)
    get_network_image(mapping_result, species, prefix + "string_protein-protein-interaction.png")
    get_network(mapping_result, species, prefix + "string_protein-protein-interaction.tsv")

if __name__ == "__main__":
    main()
