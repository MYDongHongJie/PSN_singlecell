##===================================================== report =======================================================
localrules: report
rule report:
	'''
	report
	'''
	input:
		config = "config/config.yaml"
	output:
		html = f"result/report/{project}_Report_{time.strftime('%Y_%m_%d')}/report.html"
	log:
		"logs/report.log"
	resources:
		qsub_mem=10,
		qsub_p=config['cpu']['report']
	params:
		report=f"result/report/{project}_Report_{time.strftime('%Y_%m_%d')}/"
	shell:
		'''
		rm -rf {params.report}
		/public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report.py -i result -c {input.config} |& tee {log}
		'''
