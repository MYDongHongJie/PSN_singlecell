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

# @click.command()
# @click.option('-i', '--input', prompt='the program path',
#               help='the directory of your program.[default:./]',default="./")
# @click.option('-c', '--configfile', prompt='config yaml file. [default:config/config.yaml]', default='config/config.yaml',
#               help='the config file which contain your program information.')
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
def st_report(input, configfile):
    """ python3 st_report.py  -i ./ -c ./config.txt"""
    with open(configfile) as f:
        config = ordered_yaml_load(f)

    # file = open(configfile, 'r', encoding='utf-8')
    # cfg = file.read()
    # config = yaml.load(cfg)
    program_num = config['report']['Project_Num']
    program_path = os.path.abspath(input)
    if config['report']['Project_End_Time'] == "":
        config["report"]['Project_End_Time'] = time.strftime("%Y_%m_%d")
        report_time = config["report"]['Project_End_Time']
    else:
        report_time = config['report']['Project_End_Time']
    config_ref_file = config['database']['reference']
    with open(f'{config_ref_file}/config.yaml') as f:
        config_ref = ordered_yaml_load(f)
    config["database_url"]["Ref_genome"]['Linkage'] = config_ref['genome_linkage']
    config["database_url"]["Ref_genome"]['Version'] = config_ref['genome_version']
    config["database_url"]["Ref_genome_annotation"]['Linkage'] = config_ref['annotation_linkage']
    config["database_url"]["Ref_genome_annotation"]['Version'] = config_ref['annotation_version']
    ###################################### Report ########################################################
    ######################################################################################################
    os.chdir(program_path)
    log_file = open("%s/logs/report_log.txt" % (program_path), "w")
    script_dir = os.path.dirname(os.path.realpath(__file__))

    num = 0
    if config["module"]["spaceranger_report"]:
        if os.path.exists("%s/result/report/%s_QC_Report_%s" % (program_path, program_num,report_time )):
            shutil.rmtree("%s/result/report/%s_QC_Report_%s" % (program_path, program_num, report_time), ignore_errors=True)
        else:
            os.makedirs("%s/result/report/%s_QC_Report_%s" % (program_path, program_num, report_time))
            print("Create spaceranger html...")
            log_file.write("Create spaceranger QC..." + "\n")
        outdir = "%s/result/report/%s_QC_Report_%s" % (program_path, program_num, report_time)

        ###################################### SpaceRanger
        if os.path.exists("%s/result/spaceranger" % (program_path)):
            num += 1
            name = [str(i).split('/')[-2] for i in glob("result/spaceranger/*/outs")]
            #name.remove('aggr')  
            # samplefile = pd.read_csv('../../../config/samples.csv', index_col=0, sep=',')
            # sample_name = samplefile.index[0]
            for j in name:
                os.makedirs("%s/%s.SpaceRanger/%s" % (outdir, num, j))
                subprocess.call('ln -s %s/result/spaceranger/%s/outs/* %s/%s.SpaceRanger/%s ' % (program_path, j, outdir, num, j),
                                shell=True)                  
                subprocess.call('unlink %s/%s.SpaceRanger/%s/spatial' % (outdir, num, j), shell=True)      
                subprocess.call('cp -r %s/result/spaceranger/%s/outs/spatial/ %s/%s.SpaceRanger/%s/spatial/' % (program_path, j, outdir, num, j), shell=True)  
                #subprocess.call('rm  %s/%s.SpaceRanger/%s/spatial/*.tiff' % (outdir, num, j), shell=True) 
                if os.path.exists(f"{outdir}/{num}.SpaceRanger/{j}/possorted_genome_bam.bam" ):
                        subprocess.call(f'rm {outdir}/{num}.SpaceRanger/{j}/possorted_genome_bam.bam*', shell=True) 
                if os.path.exists(f"{outdir}/{num}.SpaceRanger/aggr" ):
                        subprocess.call(f'rm -r {outdir}/{num}.SpaceRanger/aggr/', shell=True)                            

        else:
            print("Can not find SpaceRanger results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find SpaceRanger results!!!!!!!!!!!!!!!!" + "\n")

        num += 1
        subprocess.call('cp -r %s/supplemental_material/ %s/%s.Supplemental_material/' % (script_dir, outdir,num), shell=True)
        #subprocess.call('cp -r %s/src/ %s/' % (script_dir, outdir), shell=True)
    # ## HE images
    # if os.path.exists("%s/result/HEimages" % (program_path)):
    #     os.makedirs("%s/%s.HEimages/" % (outdir, num))
    #     subprocess.call('cp %s/result/HEimages/*tif %s/%s.HEimages/ ' % (program_path, outdir, num), shell=True)

    if config["module"]["report"]:
        if os.path.exists("%s/result/report/%s_Report_%s" % (program_path, program_num,report_time )):
            shutil.rmtree("%s/result/report/%s_Report_%s" % (program_path, program_num, report_time), ignore_errors=True)
        else:
            os.makedirs("%s/result/report/%s_Report_%s" % (program_path, program_num, report_time))
            print("Create Report...")
            log_file.write("Create Report..." + "\n")
        outdir = "%s/result/report/%s_Report_%s" % (program_path, program_num, report_time)

        ###################################### SpaceRanger
        if os.path.exists("%s/result/spaceranger" % (program_path)):
            num += 1
            name = [str(i).split('/')[-2] for i in glob("result/spaceranger/*/outs")]

            for j in name:
                os.makedirs("%s/%s.SpaceRanger/%s" % (outdir, num, j))
                subprocess.call('ln -s %s/result/spaceranger/%s/outs/* %s/%s.SpaceRanger/%s ' % (program_path, j, outdir, num, j),shell=True)
                if os.path.exists(f"{outdir}/{num}.SpaceRanger/{j}/possorted_genome_bam.bam" ):
                        subprocess.call(f'rm {outdir}/{num}.SpaceRanger/{j}/possorted_genome_bam.bam*', shell=True)
        
                #subprocess.call('rm %s/%s.SpaceRanger/%s/cloupe.cloupe' % (outdir, num, j), shell=True)
                #subprocess.call('cp %s/result/spaceranger/%s.png %s/%s.SpaceRanger/%s ' % (program_path,j, outdir, num,j), shell=True)
        else:
            print("Can not find SpaceRanger results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find SpaceRanger results!!!!!!!!!!!!!!!!" + "\n")

        ###################################### Count_QC
        if os.path.exists("%s/result/count_qc" % (program_path)):
            num += 1
            os.makedirs("%s/%s.Count_QC" % (outdir, num))
            subprocess.call('cp %s/result/count_qc/[!f]*  %s/%s.Count_QC/' % (program_path, outdir, num), shell=True)
        else:
            print("Can not find Count_QC results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find Count_QC results!!!!!!!!!!!!!!!!!" + "\n")

        ###################################### Clustering

        if os.path.exists("%s/result/cluster_seurat/umap_Dimension_Reduction" % (program_path)):
            num += 1
            os.makedirs("%s/%s.Clustering" % (outdir, num))
            subprocess.call('cp %s/result/cluster_seurat/umap_Dimension_Reduction/* %s/%s.Clustering/' % (program_path, outdir, num), shell=True)
            # if config['params']['reduct1_method'] == "mnn":
            #     subprocess.call('cp %s/result/cluster_seurat_pca/umap_Dimension_Reduction/* %s/%s.Clustering/' % (
            #     program_path, outdir, num), shell=True)
        else:
            print("Can not find umap_Dimension_Reduction results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find umap_Dimension_Reduction results!!!!!!!!!!!!!!!!" + "\n")

        if len(name) >= 1:
            num +=1 
            subprocess.call('cp -r %s/result/cluster_seurat/visualize_cluster_by_* %s/%s.Visualize_cluster_by_clusters/' % (program_path, outdir, num), shell=True)

        ###################################### Marker
        if os.path.exists("%s/result/marker" % (program_path)):
            num += 1
            os.makedirs("%s/%s.Marker" % (outdir, num))
            subprocess.call('cp -r %s/result/marker/markers_vis4cluster* %s/%s.Marker/' % (program_path, outdir, num), shell=True)
            subprocess.call('cp -r %s/result/marker/topmarker_gene_heatmap.p* %s/%s.Marker/' % (program_path, outdir, num),
                            shell=True)
            subprocess.call('cp -r  %s/result/marker/top10_markers_for_each_cluster_anno.xls %s/%s.Marker/' % (program_path, outdir, num),
                            shell=True)
            subprocess.call(
                'cp -r  %s/result/marker/all_markers_for_each_cluster_anno.xls %s/%s.Marker/' % (program_path, outdir, num),
                shell=True)
        else:
            print("Can not find Marker results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find Marker results!!!!!!!!!!!!!!!!" + "\n")
        ###################################### Celltype
        if os.path.exists("%s/result/celltype" % (program_path)):
            num += 1
            os.makedirs("%s/%s.Reference_CellType" % (outdir, num))
            os.makedirs("%s/%s.Reference_CellType/1.celltype_spatialplot" % (outdir, num))
            os.makedirs("%s/%s.Reference_CellType/2.celltype_interaction" % (outdir, num))
            os.makedirs("%s/%s.Reference_CellType/3.celltype_heatmap" % (outdir, num))
            os.makedirs("%s/%s.Reference_CellType/4.top_celltype_spatialplot" % (outdir, num))
            subprocess.call(
                'cp -r %s/result/celltype/1.celltype_spatialplot/*_spatial_transcriptomics_celltyping_plot.*  %s/%s.Reference_CellType/1.celltype_spatialplot' % (
                    program_path, outdir, num),
                shell=True)
            subprocess.call(
                'cp -r %s/result/celltype/1.celltype_spatialplot/split*  %s/%s.Reference_CellType/1.celltype_spatialplot' % (program_path, outdir, num),
                shell=True)
            subprocess.call('cp -r  %s/result/celltype/*celltype.xls  %s/%s.Reference_CellType/' % (program_path, outdir, num), shell=True)

            subprocess.call('cp -r  %s/result/celltype/2.celltype_interaction/*celltype_interaction*.p*  %s/%s.Reference_CellType/2.celltype_interaction' % (
                program_path, outdir, num), shell=True)
            subprocess.call(
                'cp -r  %s/result/celltype/2.celltype_interaction/*interactions.xls  %s/%s.Reference_CellType/2.celltype_interaction' % (program_path, outdir, num),
                shell=True)
            subprocess.call(
                'cp -r  %s/result/celltype/3.celltype_heatmap/*celltype_heatmap.p*  %s/%s.Reference_CellType/3.celltype_heatmap' % (program_path, outdir, num),
                shell=True)
            subprocess.call(
                'cp -r  %s/result/celltype/4.top_celltype_spatialplot/*celltype_spatial*  %s/%s.Reference_CellType/4.top_celltype_spatialplot' % (program_path, outdir, num),
                shell=True)
            subprocess.call(
                'cp -r  %s/result/celltype/celltype_infor.csv  %s/%s.Reference_CellType/' % (program_path, outdir, num),
                shell=True)
            # subprocess.call(
            #     'cp -r  %s/result/celltype/RCTD细胞类型鉴定结果说明.docx  %s/%s.Reference_CellType/' % (program_path, outdir, num),
            #     shell=True)
                
        else:
            print("Can not find CellTyping results!!!!!!!!!!!!!!!!")
            log_file.write("Can not find CellTyping results!!!!!!!!!!!!!!!!" + "\n")
        ###################################### Diffexp
        if os.path.exists("%s/result/diffexp" % (program_path)):
            num += 1
            os.makedirs("%s/%s.Diffexp" % (outdir, num))
            subprocess.call('cp %s/result/diffexp/diffexp_results_stat.xls %s/%s.Diffexp/' % (program_path, outdir, num),
                            shell=True)
            subprocess.call('cp %s/result/diffexp/*-vs-*-diff-*_anno.xls %s/%s.Diffexp/' % (program_path, outdir, num), shell=True)
            subprocess.call('cp %s/result/diffexp/*-vs-*-all_diffexp_genes_anno.xls %s/%s.Diffexp/' % (program_path, outdir, num),
                            shell=True)
        else:
            print("There is no Diffexp result.")
            log_file.write("There is no Diffexp result." + "\n")

        ###################################### enrichment
        if os.path.exists("%s/result/diffexp/enrichment" % (program_path)):
            num += 1
            os.makedirs("%s/%s.Enrichment/background_files" % (outdir, num))
            subprocess.call('cp -r %s/result/diffexp/enrichment/KEGG_enrichment  %s/%s.Enrichment/KEGG_enrichment' % (program_path, outdir, num), shell=True)
            subprocess.call('cp -r %s/result/diffexp/enrichment/GO_enrichment  %s/%s.Enrichment/GO_enrichment' % (program_path, outdir, num), shell=True)
            gene_annotation = os.path.join(config['database']['reference'],"annotation")
            subprocess.call('cp -r %s/gene_annotation.xls  %s/%s.Enrichment/background_files' % (gene_annotation, outdir, num), shell=True)
            subprocess.call('cp -r %s/gene_go.backgroud.xls  %s/%s.Enrichment/background_files' % (gene_annotation, outdir, num), shell=True)
            subprocess.call('cp -r %s/gene_kegg.backgroud.xls  %s/%s.Enrichment/background_files' % (gene_annotation, outdir, num), shell=True)

        else:
            print("There is no enrichment result.")
            log_file.write("There is no enrichment result." + "\n")

        ###################################### PPI
        if os.path.exists("%s/result/diffexp/ppi" % (program_path)):
            num += 1
            os.makedirs("%s/%s.PPI" % (outdir, num))
            subprocess.call(
                'cp %s/result/diffexp/ppi/*string_protein-protein-interaction* %s/%s.PPI/' % (program_path, outdir, num),
                shell=True)
            subprocess.call(
                'cp %s/result/diffexp/ppi/*_ppi_network* %s/%s.PPI/' % (program_path, outdir, num),
                shell=True)
        else:
            print("There is no PPI result.")
            log_file.write("There is no PPI result." + "\n")

        num += 1
        subprocess.call('cp -r %s/supplemental_material/ %s/%s.Supplemental_material/' % (script_dir, outdir,num), shell=True)
        subprocess.call('rm %s/%s.Supplemental_material/10x_visium_空间转录组报告树状图.pdf' % (outdir,num), shell=True)
        subprocess.call('cp -r %s/supplemental_material/10x_visium_空间转录组报告树状图.pdf %s/' % (script_dir, outdir), shell=True)
    #subprocess.call('cp -r %s/src/ %s/result/' % (script_dir, program_path), shell=True)
        #subprocess.call('cp -r %s/src/ %s/' % (script_dir, outdir), shell=True)
    ####get database url information
    database_url=config['database_url']
    database = [dict(database_url[i],**{'Database': i}) for i in database_url.keys()]
    database = pd.DataFrame(database)
    database = database[['Database','Linkage','Version']]
    if config["module"]["spaceranger_report"]:
        database = database[:2]
    database.to_csv(f'{program_path}/scripts/report/pic/tables/database.txt', encoding='utf-8', sep='\t', index=False)
    ###get software information
    softwares_info = config['softwares']
    softwares = [dict(softwares_info[i], **{'Softwares': i}) for i in softwares_info.keys()]
    softwares = pd.DataFrame(softwares)
    if config["params"]['celltype'] != "SPOTLight":
        mask = softwares['Softwares'].str.contains("SPOTLight")
    else:
        mask = softwares['Softwares'].str.contains("RCTD")
    softwares = softwares[~mask]
    softwares= softwares[['Softwares','Version','Parameter']]
    softwares.to_csv(f'{program_path}/scripts/report/pic/tables/softwares.txt', encoding='utf-8', sep='\t', index=False)
    with open(configfile, 'w', encoding='utf-8') as fp:
        ordered_yaml_dump(config, fp, default_flow_style=False, encoding='utf-8', allow_unicode=True)
    log_file.close()
        #subprocess.call('cd %s/ && python %s/convert_to_html.py' % (outdir, script_dir), shell=True)
    #subprocess.call('cd %s/ && python %s/convert_to_html.py' % (outdir, script_dir), shell=True)
if __name__ == "__main__":
    st_report()
