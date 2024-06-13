#===========================================================================================================
# 运行spaceranger
# #===========================================================================================================
rule spaceranger:
        """
        run spaceranger for each samples 
        """
        input:
             #"logs/download.check.sucess.txt",
             fastq=get_fastqs,
             slide_file=get_slide_file 
        output:
              "result/spaceranger/{sample}/outs/web_summary.html",
              "result/spaceranger/{sample}/outs/metrics_summary.csv",
              "result/spaceranger/{sample}/outs/cloupe.cloupe",
              "result/spaceranger/{sample}/outs/raw_feature_bc_matrix.h5",
              "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5",
              "result/spaceranger/{sample}/outs/spatial/tissue_lowres_image.png",
              "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
        params:
            ref=reference,
            project=project,
            slide=get_slide,
            slidearea=get_slidearea,
            sample=get_sample,
            sample_type=get_sample_type,
            origin_image = get_origin_image,
            image_file= get_slide_image,
            cytaimage_cmd = "--cytaimage="  if config["params"]["library_type"] == "cytassist" else " ",
            cytaimage_file = get_cytaimage_image if config["params"]["library_type"] == "cytassist" else " ",
            ffpe_probe=  " "  if config["params"]["library_type"] == "fresh"  else f"--probe-set={probe_set} ",
            loupealignment=  get_sample_image_json_file ,
            r2_length= 91 if  config["params"]["library_type"] == "fresh"  else 50,
            generate_bam =  " " if config["params"]["generate_bam"] else "--no-bam",
            qsub_mem = config['spaceranger']['mem'],
            qsub_p = config['spaceranger']['cpu']
        benchmark:
                 "benchmarks/spaceranger_for_{sample}.benchmark.txt"
        envmodules:
                envmodules['spaceranger']
        shell:
            """
                 cd  result/spaceranger/;
                 rm -rf {params.sample};
                 spaceranger count    \
                      --id={params.sample}    \
                      --transcriptome={params.ref}   \
                      --fastqs={input.fastq}    \
                      --sample={params.sample}   \
                      {params.image_file} \
                      {params.cytaimage_cmd}{params.cytaimage_file} \
                      --slide={params.slide}   \
                      --slidefile={input.slide_file}   \
                      --area={params.slidearea}    \
                      --r1-length=28    \
                      --r2-length={params.r2_length}   \
                      --localcores={params.qsub_p}    \
                      --localmem={params.qsub_mem}    \
                      --description={params.project}_{params.sample_type}\
                      {params.ffpe_probe} \
                      {params.generate_bam} \
                      {params.loupealignment}
                 cd {params.sample}
                 ls | grep -v outs | xargs rm -rf
                 image_name=$(echo {params.origin_image} |sed 's:/: :g' |awk '{{print $NF}}')
                 cp {params.origin_image}   ./outs/spatial/image_$image_name
                 if [ ! "{params.cytaimage_file}" = " " ];then
                     cytaimage_name=$(echo {params.cytaimage_file} |sed 's:/: :g' |awk '{{print $NF}}')
                     cp {params.cytaimage_file}   ./outs/spatial/cytaimage_$cytaimage_name
                 fi
            """
#===========================================================================================================
# spaceranger 网页截图
#===========================================================================================================
localrules: spaceranger_shoot_clear

rule spaceranger_shoot_clear:
    input:
        web_summary = "result/spaceranger/{sample}/outs/web_summary.html"
    output:
        sample_png = "result/spaceranger/{sample}/outs/{sample}.png"
    benchmark:
        "benchmarks/spaceranger_shoot_clear_for_{sample}.benchmark.txt"
    log:
        "logs/spaceranger/run.snakemake.screenshot.{sample}.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        '''
        oest_html_screenshot.py \
              -i {input.web_summary} \
              -n {wildcards.sample} \
              -o result/spaceranger/{wildcards.sample}/outs >&2 2>{log} \
        # && rm -r raw_data/{wildcards.sample}/*.fastq.gz
        '''
#=======================================================================================================================
# 生产spacerangerQC报告
#=======================================================================================================================
if config["module"]["spaceranger_email"]:
    #localrules: spaceranger_email
    rule spaceranger_email:
        """
        spaceranger_report summary 
        """
        input:
            metadata=config['samples'],
            html=expand(["result/spaceranger/{sample}/outs/web_summary.html",
                         "result/spaceranger/{sample}/outs/metrics_summary.csv"],sample=samples.index)
        output:
            "logs/spaceranger_report.success.log",
            f"result/report/{project}-QC_spaceranger_report/summary.csv",
            f"result/report/{project}-QC_spaceranger_report/summary.csv.email.html"
        benchmark:
            f'benchmarks/spaceranger_report_for_{project}.benchmark.txt'
        log:
            f"logs/report/run.snakemake.spaceranger_{project}_report.output"
        envmodules:
            envmodules["oesinglecell"]
        shell:
            """
            sctool   -o  ./result/report  \
                    st-email  \
                        -s  {input.metadata} \
                        -i  ./result/spaceranger \
                        -c  config/config.yaml  >&2  2>{log} && touch logs/spaceranger_report.success.log
            """

if config["module"]["spaceranger_report"]:
    #localrules: spaceranger_report
    import time,re,os
    if config['report']['Project_End_Time']=="":
        config['report']['Project_End_Time']=time.strftime("%Y_%m_%d")
    rule spaceranger_report:
        """
        spaceranger QC report 只质控不分析 
        """
        input:
            spaceranger_outfiles=expand(["result/spaceranger/{sample}/outs/web_summary.html",
                                         "result/spaceranger/{sample}/outs/metrics_summary.csv",
                                         "result/spaceranger/{sample}/outs/cloupe.cloupe",
                                         "result/spaceranger/{sample}/outs/raw_feature_bc_matrix.h5",
                                         "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5"],sample=samples.index),
            config="config/config.yaml"
        output:
            f"logs/spaceranger_QC_report.txt"
        params:
            outdir=f"result/report/{config['report']['Project_Num']}_QC_Report_{config['report']['Project_End_Time']}",
            inputdir=os.getcwd()
        benchmark:
            f"benchmarks/report_summary_for_{config['report']['Project_Num']}_QC_Report_{config['report']['Project_End_Time']}.benchmark.txt"
        log:
            f"logs/{config['report']['Project_Num']}_QC_Report_{config['report']['Project_End_Time']}.log"
        shell:
            """
            ## step1: 文件拷贝  and  ## step2: html文件生成
            module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/step1.report_files_copy.py  -i {params.inputdir} -c {input.config} && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/step2.spaceranger_html_convert.py -c {input.config}  -d {params.outdir} >&2  2>{log} && touch logs/spaceranger_QC_report.txt
            """

    # rule spaceranger_report_upload:
    #     """
    #     cellranger report upload
    #     """
    #     input:
    #         report = f"result/report/{config['report']['Project_Num']}_QC_Report_{config['report']['Project_End_Time']}/report.html",
    #         mark = f'logs/spaceranger_QC_report.txt'       
    #     output:
    #         f"result/report/{config['report']['Project_Num']}_QC_Report_{config['report']['Project_End_Time']}.zip",
    #         "logs/spaceranger_report_upload.check"
    #     params:
    #         report_num=re.sub("-.*","",config['report']['Project_Num']), 
    #         report = f"{config['report']['Project_Num']}_QC_Report_{config['report']['Project_End_Time']}"
    #     log:
    #        "logs/report_upload.log"
    #     benchmark:
    #         "benchmarks/report_upload.benchmark.txt"
    #     resources:
    #         qsub_mem=10,
    #         qsub_p=1
    #     shell:
    #         '''
    #         cd result/report/ ; 
    #         if [ -d {params.report}"/result" ] ; then 
    #             rm -r {params.report}/result 
    #         fi
    #         /public/scRNA_works/works/guokaiqi/software/zip -r {params.report}.zip {params.report} >&2  2> ../../{log};
    #         /data/software/obsutil/5.2.12/obsutil cp {params.report}.zip   obs://oe-scrna/Analysis_Report/{params.report_num}/ -f -vlength &&
    #         ( echo upload succeed on : obs://oe-scrna/Analysis_Report/{params.report_num}/{params.report}.zip  |& tee ../../logs/spaceranger_report_upload.check ) ||
    #         echo upload fails |& tee ../../{log} && rm -r {params.report}.zip
    #         '''

      
#=======================================================================================================================
#  spaceranger aggr
#=======================================================================================================================
if config["module"]["spaceranger_aggr"]:
    import os
    import pandas as pd
    rule spaceranger_aggr_library_meta:
        """
        getting spaceranger aggr library meta 
        """
        input:
            metadata=config['samples']
        output:
            meta=os.path.abspath("config/library.csv")
        run:
            samples = pd.read_csv(input.metadata,dtype=str).set_index("sampleid",drop=False)
            library = pd.DataFrame(columns=['library_id', 'molecule_h5', 'cloupe_file', 'spatial_folder'])
            for sample in list(samples.index):
                library.loc[sample, 'library_id'] = sample
                library.loc[sample, 'molecule_h5'] = f"{sample}/outs/molecule_info.h5"
                library.loc[sample, 'cloupe_file'] = f"{sample}/outs/cloupe.cloupe"
                library.loc[sample, 'spatial_folder'] = f"{sample}/outs/spatial"
                library.to_csv(output.meta,encoding='utf-8',index=False)
    rule spaceranger_aggr:
        """
        spaceranger aggr
        """
        input:
            spaceranger_outfiles = expand(["result/spaceranger/{sample}/outs/web_summary.html",
                                           "result/spaceranger/{sample}/outs/metrics_summary.csv",
                                           "result/spaceranger/{sample}/outs/cloupe.cloupe",
                                           "result/spaceranger/{sample}/outs/raw_feature_bc_matrix.h5",
                                           "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5"], sample=samples.index),
              meta = os.path.abspath("config/library.csv")
        output:
              "result/spaceranger/aggr/outs/web_summary.html",
              "result/spaceranger/aggr/outs/summary.json",
              "result/spaceranger/aggr/outs/cloupe.cloupe",
              # "result/spaceranger/aggr/outs/raw_feature_bc_matrix.h5",
              "result/spaceranger/aggr/outs/filtered_feature_bc_matrix.h5"
        params:
           outdir="result/report/spaceranger",
           normalize=config["params"]["aggr_normalize"]
        benchmark:
             "benchmarks/spaceranger_aggr.benchmark.txt" ##spaceranger 运行无法采用 >&2 2>{log}，需要注意
        envmodules:
            envmodules['spaceranger']
        shell:
           """
           cd  result/spaceranger/;
           rm -rf aggr; 
           spaceranger aggr \
               --id=aggr  \
               --csv={input.meta}  \
               --normalize=none  
           """
#=======================================================================================================================
#  如果需要生成bam文件，则后续将进行整理上传，如果未生成，则不进行该步骤
#=======================================================================================================================
if config["params"]["generate_bam"] and config["params"]["bam_upload_rm"]:
    #localrules: bam_upload_rm
    import os
    rule bam_upload_rm:
        """
        prepare bam files for each sample and upload to  OBS
        """
        input:
            web_summary="result/spaceranger/{sample}/outs/web_summary.html"
            # bam="result/spaceranger/{sample}/outs/possorted_genome_bam.bam",
            # bai="result/spaceranger/{sample}/outs/possorted_genome_bam.bam.bai",
        output:
            "logs/{sample}_upload_bam_check"
        #log: 
        #   "logs/{sample}_upload_bam.log"
        params:
            project=project,
            dir=os.getcwd()
        shell:
            '''
            bash scripts/upload_bam.sh {params.dir}/result/spaceranger/{wildcards.sample} {params.project}  && touch logs/{wildcards.sample}_upload_bam_check    
            '''

if config["module"]["spaceranger_upload"]:
    rule spaceranger_upload :
        input:
            expand([
            "result/spaceranger/{sample}/outs/web_summary.html",
            "result/spaceranger/{sample}/outs/metrics_summary.csv",
            "result/spaceranger/{sample}/outs/cloupe.cloupe",
            "result/spaceranger/{sample}/outs/raw_feature_bc_matrix.h5",
            "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix.h5",
            "result/spaceranger/{sample}/outs/spatial/tissue_lowres_image.png",
            "result/spaceranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
            "result/spaceranger/{sample}/outs/{sample}.png",
            "logs/{sample}_upload_bam_check"] ,sample=samples.index)
        output:
            "logs/spaceranger_upload.check"
        log:
            "logs/spaceranger_upload.success.log"
        benchmark:
            "benchmarks/spaceranger_upload.benchmark.txt"
        params:
            project = project
        resources:
            qsub_mem=10,
            qsub_n=1
        shell:
            '''
            bash scripts/upload_cellranger_m.sh $PWD/result/spaceranger {params.project} |& tee {log}  && touch {output}
            '''

