rule convert:
    input:
        cluster_rds = "result/Count_QC/filtered.h5seurat",
        #temp_file="clustering_h5seurat_touched.check" # to make sure running after Clustering
    output:
        v3rds = "result/rds/data_ob_v3.rds"
    benchmark:
        "benchmarks/convert.benchmark.txt"
    log:  "logs/convert.log"
    resources:
        qsub_mem=40,
        qsub_p=2,
    shell:
        """
module load OESingleCell/3.0.d &&  Rscript scripts/R/v4tov3.r -i {input.cluster_rds} -f h5seurat -o result/rds &&
module unload OESingleCell/3.0.d && module load OESingleCell/2.0.0 && Rscript scripts/R/v4tov3.r -o result/rds
"""

localrules: report
rule report:
    """
    report summary
    """
    input:
        v3rds = "result/rds/data_ob_v3.rds"
        #config = "config/config.yaml",
    output:
        html = f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}/分析说明.html"
    params:
        report = f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}/",
        config = "config/config.yaml",
    benchmark:
        "benchmarks/report.benchmark.txt"
    log:  "logs/report.log"
    resources:
        qsub_mem=10,
        qsub_p=1,
    #envmodules:
        #config['report_params']['envmodules']
        #"oebio/1.2.15"
    shell:
        """
rm -rf {params.report}
module purge && /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report/scRNA_report.py \
  -i result \
  -c {params.config} \
  -r {input.v3rds} |& tee {log};
        """

localrules: report_upload
rule report_upload:
    """
    report upload
    zip hcu-0006 only
    """
    input:
        report = f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}/分析说明.html"
    output:
        protected(f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}.zip"),
        protected(f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}.zip.md5"),
        "result/report/report_upload.check"
    params:
        report_num = config['report']['Project_Num'],
        report = f"{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}"
    benchmark:
        "benchmarks/report_upload.benchmark.txt"
    log:  "logs/report_upload.log"
    resources:
        qsub_mem=10,
        qsub_p=1,
    #envmodules:
    #    #config['report_params']['envmodules']
    #    "oebio/1.2.15"
    shell:
        """
cd result/report/ ; 
zip -r {params.report}.zip {params.report} ;
md5sum {params.report}.zip > {params.report}.zip.md5 ;
obsutil cp {params.report}.zip   obs://oe-scrna/Analysis_Report/{params.report_num}/ -f -vlength -vmd5;
obsutil cp {params.report}.zip.md5   obs://oe-scrna/Analysis_Report/{params.report_num}/ -f -vlength -vmd5 &&
    ( echo upload succeed on : obs://oe-scrna/Analysis_Report/{params.report_num}/{params.report}.zip  |& tee report_upload.check ) ||
    echo upload fails |& tee {log}
        """

localrules: report_sor
rule report_sor:
    """
    upload report SOR
    """
    input:
        report = f"result/report/{config['report']['Project_Num']}_Report_{time.strftime('%Y_%m_%d')}/分析说明.html"
        # metrics_summary_stat="result/cellranger_sum/metrics_summary_stat.csv",
        #metadata = "config/samples.csv",
    output:
        directory("logs/upload_report_sor")
    params:
        config = "config/config.yaml"
    benchmark:
        "benchmarks/report_sor.benchmark.txt"
    resources:
        qsub_mem=4,
        qsub_p=1
    envmodules:
        config["report_params"]['envmodules']
    shell:
        """
python  scripts/scrna_sor.py \
  -i {params.config} \
  -l logs/upload_report_sor \
  -s "Project,QC,Cpu_Mem,Ref_database,Software,Sample"
        """
