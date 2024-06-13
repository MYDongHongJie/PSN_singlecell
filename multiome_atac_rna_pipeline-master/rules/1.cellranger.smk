##=====================================================step1:  data prepare ========================================================
rule data_prepare:
    """
    prepare fastqs for each sample
    """
    input:
        metadata = 'config/metadata.csv',
        yaml = 'config/config.yaml'
    output:
       temp_file=  temp("raw_data.check"),
       libraries=expand("config/library/{sample}_libraries.csv",sample=samples.index),
       libraries_aggr=os.path.abspath("config/library/libraries_aggr.csv")
    benchmark:
        "benchmarks/cellranger/data_prepare.benchmark.txt"
    log:
        "logs/cellranger/data_prepare.log"
    resources:
        qsub_mem=20,
        qsub_p=config['cpu']['data_prepare'],
        qsub_n=1
    envmodules:
        config['envmodules']['oesinglecell']
    shell:
        '''
        python scripts/data_prepare.py -i {input.metadata} -c {input.yaml} -o raw_data  |& tee {log}  && touch raw_data.check
        '''

# =================================================== step2: run cellranger ====================================================
rule cellranger:
    """
    run cellranger for each samples
    """
    input:
        libraries = "config/library/{sample}_libraries.csv"
    output:
        "result/cellranger/{sample}/outs/web_summary.html",
        "result/cellranger/{sample}/outs/summary.csv",
        "result/cellranger/{sample}/outs/atac_fragments.tsv.gz",
        "result/cellranger/{sample}/outs/per_barcode_metrics.csv",
        "result/cellranger/{sample}/outs/gex_molecule_info.h5"
    params:
        project = project,
        species = Species,
        reference=reference,
        sample=get_sample,
    benchmark:
        "benchmarks/cellranger/{sample}.cellranger.benchmark.txt"
    resources:
        qsub_mem=120,
        qsub_p=config['cpu']['cellranger'],
        qsub_n=1
    envmodules:
        config["envmodules"]['cellranger']
    shell:
        '''
        cd result/cellranger/ ;
        rm -rf {params.sample};
        cellranger-arc count \
                    --id={params.sample}   \
                    --reference={params.reference} \
                    --libraries=../../{input.libraries} \
                    --localcores={resources.qsub_p}  \
                    --localmem={resources.qsub_mem}   \
                    --description={params.project}_{params.species} \
                    --no-bam
        '''
##=============================================== step 3: report prepare =====================================================
localrules: report_prepare
rule report_prepare:
    input:
        webfile="result/cellranger/{sample}/outs/web_summary.html",
        summary="result/cellranger/{sample}/outs/summary.csv"
    output:
        sample_webfile="result/cellranger/{sample}/outs/{sample}_summary.html",
        sample_summary="result/cellranger/{sample}/outs/{sample}_summary.csv",
        sample_png="result/cellranger/{sample}/{sample}.png"
    envmodules:
        config["envmodules"]['chromedriver'],
        config["envmodules"]['oebio']
    params:
        outputfile="result/cellranger/{sample}/",
        pngname="{sample}"
    resources:
        qsub_mem = 10,
        qsub_p = 1,
        qsub_n = 1
    shell:
        '''
        cp {input.webfile} {output.sample_webfile}
        cp {input.summary} {output.sample_summary}
        python scripts/selenium_screenshot.py  -i {input.webfile} -n {params.pngname} -o  {params.outputfile}
        '''
##================================================ step4: report and send e-mail ===================================================
localrules: cellranger_report
rule cellranger_report:
    input:
        samples_csv=expand("result/cellranger/{sample}/outs/{sample}_summary.csv",sample=samples.index),
        samples_summary=expand("result/cellranger/{sample}/outs/{sample}_summary.html",sample=samples.index),
    output:
        report="logs/summary_report.log",
        sor=directory( "logs/upload_sor"),
        summary_file="result/cellranger/summary.csv"
    params:
        config = "config/config.yaml"
    resources:
        qsub_mem=10,
        qsub_p=1,
        qsub_n=1
    shell:
        """
        ##qc summary and report
        module load /public/dev_scRNA/sunkun/pipeline/multimode_pipeline/envs/OESingleCell/beta
        Rscript scripts/qc_summary.R  \
                -s config/metadata.csv  \
                -i result/cellranger/ \
                -o result/report \
                -c {params.config} && touch logs/summary_report.log
        module purge
       ## sor upload
       module load /public/dev_scRNA/sunkun/pipeline/multimode_pipeline/envs/OESingleCell/dev_multimodel
       python  scripts/sor_upload.py \
                -i  {params.config} \
                -l  logs/upload_sor \
                -s "Project,QC,Cpu_Mem,Ref_database,Software"
       module purge
       if [ -d "raw_data/download" ]; then rm -r raw_data/download ;fi
        """
#=======================================================step 5: run cellranger_arc aggr  ==========================================
if config["cellranger_params"]["aggr"]:
    rule run_aggr:
        '''
        run cellranger_arc aggr
        '''
        input:
            expand( "result/cellranger/{sample}/outs/atac_fragments.tsv.gz",sample=samples.index),
            expand( "result/cellranger/{sample}/outs/per_barcode_metrics.csv",sample=samples.index),
            expand( "result/cellranger/{sample}/outs/gex_molecule_info.h5",sample=samples.index),
            libraries=os.path.abspath("config/library/libraries_aggr.csv")
        output:
            "result/cellranger/aggr/{params.project}_{params.species}/outs/web_summary.html"
        benchmark:
            "benchmarks/cellranger/run_aggr.benchmark.txt"
        log:
            "logs/cellranger/run_aggr.log"
        resources:
            qsub_mem=64,
            qsub_p=config[ 'cpu']['cellranger_aggr'],
            qsub_n=1
        envmodules:
            config['envmodules']['oesinglecell'],
            config['envmodules']['cellranger']
        params:
            project=project,
            species=Species,
            reference=reference
        shell:
            '''
            cd result/cellranger/aggr
            cellranger-arc aggr --id= {params.project}_{params.species} \
                                --csv= {input.libraries} \
                                --normalize=depth \
                                --reference={params.reference}
            '''
#=======================================================step 6: QC report =============================================
if config["module"]["QCreport"]:
    localrules: run_qcreport
    rule run_qcreport:
        '''
        cellrangle result report
        '''
        input:
            config = "config/config.yaml",
            summary_file="result/cellranger/summary.csv"
        output:
            f"result/report/{project}_QCreport_{time.strftime('%Y_%m_%d')}/report.html"
        log:
            "logs/qc_report.log"
        resources:
            qsub_mem=10,
            qsub_p=config['cpu']['report']
        params:
            report=f"result/report/{project}_QCreport_{time.strftime('%Y_%m_%d')}/"
        shell:
            '''
            rm -rf {params.report}
		    /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python scripts/report_QC.py -i result -c {input.config} |& tee {log}
            '''