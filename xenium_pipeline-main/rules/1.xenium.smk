##===========================================step 1: 数据下载 ===========================================================
if localPath =="":
    rule data_download:
        """
        data download
        """
        input:
            config="config/config.yaml"
        output:
            outfile = expand(['raw_data/{sample}/gene_panel.json'],sample=samples.index)
        benchmark:
            "benchmarks/data_prepare.benchmark.txt"
        log:
            "logs/data_praper/data_prepare.log"
        shell:
            '''
            python scripts/data_prepare.py -c {input.config}  -o raw_data/  >&2 2>{log}
            '''

if localPath != "":
    rule data_download:
        """
        local data download
        """
        input:
            config = "config/config.yaml"
        output:
            outfile = "raw_data/{sample}/gene_panel.json"
        benchmark:
            "benchmarks/data_prepare_{sample}.benchmark.txt"
        log:
            "logs/data_praper/data_prepare_{sample}.log"
        params:
            localPath = localPath,
            sample= get_sample
        shell:
            '''
            cd raw_data/{params.sample} && ln -s {params.localPath}/{params.sample}/* ./ >&2 2>{log}
            '''

##===========================================step 2： 使用板载系统默认分析结果 ============================================
if xenium_method == "":
    rule xenium_null:
        """
        use xenium onboard result
        """
        input:
            pannel = "raw_data/{sample}/gene_panel.json"
        output:
            "result/1.Xenium/{sample}/outs/cells.csv.gz",
            "result/1.Xenium/{sample}/outs/transcripts.csv.gz",
            "result/1.Xenium/{sample}/outs/cell_boundaries.csv.gz",
            "result/1.Xenium/{sample}/outs/metrics_summary.csv",
            "result/1.Xenium/{sample}/outs/analysis_summary.html"
        benchmark:
            "benchmarks/xenium/{sample}_xenium_null.txt"
        log:
            "logs/xenium/{sample}_xenium_null.log"
        params:
            sample = get_sample
        shell:
            '''(
            cd result/1.Xenium/{params.sample}/outs &&\
            ln -s ../../../../raw_data/mbrain1/* ./ ) >&2 2>{log}
            '''
##============================================step 3: 使用relabel参数重新定量转录本 =======================================
if xenium_method == "relabel":
    rule xenium_relabel:
        """
        xenium relabel 
        """
        input:
            panel_file = "raw_data/{sample}/gene_panel.json"
        output:
            "result/1.Xenium/{sample}/outs/cells.csv.gz",
            "result/1.Xenium/{sample}/outs/transcripts.csv.gz",
            "result/1.Xenium/{sample}/outs/cell_boundaries.csv.gz",
            "result/1.Xenium/{sample}/outs/metrics_summary.csv",
            "result/1.Xenium/{sample}/outs/analysis_summary.html"
        benchmark:
            "benchmarks/xenium/{sample}_xenium_relabel.txt"
        log:
            "logs/xenium/{sample}_xenium_relabel.log"
        params:
            sample = get_sample
        shell:
            '''(
            cd result/1.Xenium
            rm -rf {params.sample}
            xeniumranger relabel --id={params.sample} --xenium-bundle=../../raw_data/{params.sample}/ --panel=../../{input.panel_file}
            ) >&2 2>{log}
            '''

##============================================step 4: 使用resegment参数重新分割细胞 =======================================
if xenium_method == "resegment":
    rule xenium_resegment:
        """
        xenium resegment
        """
        input:
            panel_file = "raw_data/{sample}/gene_panel.json"
        output:
            "result/1.Xenium/{sample}/outs/cells.csv.gz",
            "result/1.Xenium/{sample}/outs/transcripts.csv.gz",
            "result/1.Xenium/{sample}/outs/cell_boundaries.csv.gz",
            "result/1.Xenium/{sample}/outs/metrics_summary.csv",
            "result/1.Xenium/{sample}/outs/analysis_summary.html"
        benchmark:
            "benchmarks/xenium/{sample}_xenium_resegement.txt"
        log:
            "logs/xenium/{sample}_xenium_resegement.log"
        params:
            sample = get_sample,
            expansion_distance = config['params']['Xenium']['expansion_distance'],
            resegment_nuclei = config['params']['Xenium']['resegment_nuclei']
        shell:
            '''(
            cd result/1.Xenium
            rm -rf {params.sample}
            xeniumranger resegment --id={params.sample} \
                       --xenium-bundle=../../raw_data/{params.sample}/ \
                       --expansion-distance={params.expansion_distance} \
                       --resegment-nuclei={params.resegment_nuclei} )>&2 2>{log}
            '''


##============================================step 5: 使用Import_segment参数重新分割细胞 ==================================
if xenium_method == "Import_segment":
    rule xenium_Importsegment:
        """
        xenium Import segment result
        """
        input:
            panel_file = "raw_data/{sample}/gene_panel.json"
        output:
            "result/1.Xenium/{sample}/outs/cells.csv.gz",
            "result/1.Xenium/{sample}/outs/transcripts.csv.gz",
            "result/1.Xenium/{sample}/outs/cell_boundaries.csv.gz",
            "result/1.Xenium/{sample}/outs/metrics_summary.csv",
            "result/1.Xenium/{sample}/outs/analysis_summary.html"
        benchmark:
            "benchmarks/xenium/{sample}_xenium_importsegment.txt"
        log:
            "logs/xenium/{sample}_xenium_importsegment.log"
        params:
            sample=get_sample,
            coordinate_transform=config['params']['Xenium']['coordinate_transform'],
            nuclei=config['params']['Xenium']['nuclei'],
            cells=config['params']['Xenium']['cells'],
            units=config['params']['Xenium']['units']
        shell:
            '''(
            cd result/1.Xenium
            rm -rf {params.sample}
            xeniumranger import-segmentation --id={params.sample} \
                                 --xenium-bundle=../../raw_data/{params.sample}/ \
                                 --coordinate-transform={params.coordinate_transform} \
                                 --nuclei={params.nuclei} \
                                 --cells={params.cells} \
                                 --units={params.units} )>&2 2>{log}
            '''
##=============================================step 6: 截图  ===========================================================
rule screenshot:
    """
    html screenshot
    """
    input:
        html_file= "result/1.Xenium/{sample}/outs/analysis_summary.html"
    output:
        "result/1.Xenium/{sample}/{sample}.png"
    benchmark:
        "benchmarks/xenium/{sample}_xenium_screenshot.txt"
    log:
        "logs/xenium/{sample}_xenium_screenshot.log"
    params:
        sample=get_sample
    shell:
        '''
        python scripts/oest_html_screenshot.py -w 1600 -l 3500  -i {input.html_file} -n {params.sample} -o result/1.Xenium/{params.sample} >&2 2>{log}
        '''
##=============================================step 7: 汇总 =========================================================
rule qc_summary:
    """
    xenium quality control summary
    """
    input:
        samples_csv=expand("result/1.Xenium/{sample}/outs/metrics_summary.csv",sample=samples.index)
    output:
        "result/1.Xenium/summary_all.csv"
    benchmark:
        "benchmarks/xenium/qc_summary.txt"
    log:
        "logs/xenium/qc_summary.log"
    params:
        metadataFile = config['database']['metadataFile'],
        config = "config/config.yaml"
    shell:
        '''
         Rscript scripts/qc_summary.R  \
                -s {params.metadataFile}  \
                -i result/1.Xenium/ \
                -o result/report \
                -c {params.config} >&2 2>{log}
        '''

##=============================================step 8: xenium 的质控报告 ================================================
if config['module']['QCreport']:
    rule xenium_qcReport:
        """
        xenium quality control report
        """
        input:
            config = "config/config.yaml"
            "result/1.Xenium/summary_all.csv"
        output:
            html = f"result/report/{project}_QCReport/Report.html"
        log:
            "logs/QCReport.log"
        shell:
            '''
            python scripts/report.py -i result -c {input.config} >&2 2>{log}
            '''
