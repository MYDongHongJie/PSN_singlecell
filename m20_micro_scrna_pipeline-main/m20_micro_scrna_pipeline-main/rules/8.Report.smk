##===================================================== report =======================================================
localrules: Report

rule Report:
    '''
    report
    '''
    input:
        all_files= all_input,
        config="config/config.yaml"
    output:
        html=f"result/report/{project}_Report_{time.strftime('%Y_%m_%d')}/Report.html"
    log:
        "logs/report.log"
    params:
        report=f"result/report/{project}_Report_{time.strftime('%Y_%m_%d')}/"
    shell:
        '''
        module purge &&  /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report.py -i result -c {input.config} |& tee {log}
        '''
