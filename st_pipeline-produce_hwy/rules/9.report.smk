#=======================================================================================================================
if config['report']['Project_End_Time']=="":
    config['report']['Project_End_Time']=time.strftime("%Y_%m_%d")

import os 
rule report_summary:
    """
    report summary
    """
    input:
        all_input=all_input,
        config = "config/config.yaml"
    output:
        f"result/report/{config['report']['Project_Num']}_Report_{config['report']['Project_End_Time']}/分析说明.html"
    params:
        outdir=f"result/report/{config['report']['Project_Num']}_Report_{config['report']['Project_End_Time']}",
        inputdir=os.getcwd()
    benchmark:
        f"benchmarks/report_summary_for_{config['report']['Project_Num']}_Report_{config['report']['Project_End_Time']}.benchmark.txt"
    log:
        f"logs/{config['report']['Project_Num']}_Report_{config['report']['Project_End_Time']}.log"
    shell:
        """
        ## step1: 文件拷贝  and  ## step2: html文件生成
        module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/step1.report_files_copy.py  -i {params.inputdir} -c {input.config} && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/step2.report_html_convert.py -c {input.config}  -d {params.outdir} 
        """
