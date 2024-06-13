rule report:
    """
    report
    """
    input:
        config="config/config.yaml"
    output:
        html = f"result/report/{project}_Report/Report.html"
    log:
        "logs/report.log"
    shell:
        '''
        python scripts/report.py -i result -c {input.config} >&2 2>{log}
        '''