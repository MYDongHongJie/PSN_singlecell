##===================================================== report =======================================================
rule Report:
    '''
    report
    '''
    input:
        config = "config/config.yaml",
        donefiles = inputfiles_report
    output:
        html = f"result/report/{project}_Report/Report.html"
    log:
        "logs/report.log"
    shell:
        '''
        {mp}
        /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report.py -i result -c {input.config} >&2 2>{log}
        '''

