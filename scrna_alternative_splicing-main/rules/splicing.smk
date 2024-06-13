if filter_hard:
    rule MARVEL_analysis:
        input:
            MARVEL_data = "result/Splicing/MARVEL.RData"
        output:
            marvelfile = "result/Splicing/{diff_group}/MARVEL.RData",
            DEtable = "result/Splicing/{diff_group}/diff_SJ_exp/splicing_junction_differect_expression.total.csv"
        params:
            treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
            control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
            diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
            outpath = "result/Splicing",
        threads:
            config['MARVEL_analysis']['cpu']
        resources:
            qsub_mem = config['MARVEL_analysis']['mem']
        log:
            "logs/{diff_group}/MARVEL_analysis.log"
        benchmark:
            "benchmarks/{diff_group}/MARVEL_analysis.txt"
        shell:
            """
            {mp}
            {Rscript_splicing}
            Rscript {marvel_analysis} \
                --marvelproject {input.MARVEL_data} \
                --diff_type {params.diff_type} \
                --diff_case \"{params.treatment}\" \
                --diff_control \"{params.control}\" \
                --outpath {params.outpath} >{log} 2>&1
            """
else:
    rule MARVEL_analysis:
        input:
            MARVEL_data = "result/Splicing/MARVEL.RData"
        output:
            gene_expression_dis = "result/Splicing/{diff_group}/diff_SJ_exp/distribution_gene_expression.png",
            SJ_expression_dis = "result/Splicing/{diff_group}/diff_SJ_exp/distribution_SJ_expression.png",
            marvelfile = "result/Splicing/{diff_group}/MARVEL.RData",
            DEtable = "result/Splicing/{diff_group}/diff_SJ_exp/splicing_junction_differect_expression.total.csv"
        params:
            treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
            control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
            diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
            outpath = "result/Splicing",
        threads:
            config['MARVEL_analysis']['cpu']
        resources:
            qsub_mem = config['MARVEL_analysis']['mem']
        log:
            "logs/{diff_group}/MARVEL_analysis.log"
        benchmark:
            "benchmarks/{diff_group}/MARVEL_analysis.txt"
        shell:
            """
            {mp}
            {Rscript_splicing}
            Rscript {marvel_analysis} \
                --marvelproject {input.MARVEL_data} \
                --diff_type {params.diff_type} \
                --diff_case \"{params.treatment}\" \
                --diff_control \"{params.control}\" \
                --outpath {params.outpath} >{log} 2>&1
            """

rule volcano_plot:
    input:
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData"
    output:
        volcano = expand("result/Splicing/{{diff_group}}/volcano/diff_splicing_volcano.{volcan_long_name}.{format}", volcan_long_name = volcan_long_name, format = ["png", "pdf"])
    log:
        "logs/{diff_group}/volcano_plot.log"
    benchmark:
        "benchmarks/{diff_group}/volcano_plot.txt"
    threads:
        config["volcano_plot"]["cpu"]
    resources:
        qsub_mem = config['volcano_plot']['mem']
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing",
    shell:
        """
        {mp}
        {Rscript_splicing}
        Rscript {marvel_plot} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --analysis_type volcano \
            --outpath {params.outpath}  >{log} 2>&1
        """

rule diff_SJ_exp:
    input:
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData"
    output:
        DEtable_sig = expand("result/Splicing/{{diff_group}}/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.csv", diff_long_name = diff_long_name),
        figures = expand("result/Splicing/{{diff_group}}/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.{format}", diff_long_name = diff_long_name, format = ["png", "pdf"]),
    log:
        "logs/{diff_group}/diff_SJ_exp.log"
    benchmark:
        "benchmarks/{diff_group}/diff_SJ_exp.txt"
    threads:
        config["diff_SJ_exp"]["cpu"]
    resources:
        qsub_mem = config['diff_SJ_exp']['mem']
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing"
    shell:
        """
        {mp}
        {Rscript_splicing}
        Rscript {marvel_plot} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --analysis_type diff_SJ_exp \
            --outpath {params.outpath}  >{log} 2>&1
        """

rule diff_SJ_enrich:
    input:
        DEtable_sig = expand("result/Splicing/{{diff_group}}/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.csv", diff_long_name = diff_long_name)
    output:
        genelist_temp = temp(expand("result/Splicing/{{diff_group}}/{diff_long_name}.txt", diff_long_name = diff_long_name)),
        GO = expand("result/Splicing/{{diff_group}}/enrich/1.GO_enrichment/{diff_long_name}/GO.top.Total.png", diff_long_name = diff_long_name),
        KEGG = expand("result/Splicing/{{diff_group}}/enrich/2.KEGG_enrichment/{diff_long_name}/KEGG.top.Total.png", diff_long_name = diff_long_name)
    log:
        "logs/{diff_group}/diff_SJ_enrich.log"
    benchmark:
        "benchmarks/{diff_group}/diff_SJ_enrich.txt"
    params:
        anno_path = f"{genomeDir}/annotation",
        outpath = "result/Splicing/{diff_group}/enrich"
    threads:
        config["diff_SJ_enrich"]["cpu"]
    resources:
        qsub_mem = config["diff_SJ_enrich"]["mem"]
    shell:
        """
        {mp}
        {Rscript_splicing}
        if [[ -d result/Splicing/{wildcards.diff_group}/enrich ]]; then rm -rf result/Splicing/{wildcards.diff_group}/enrich ; fi
        awk '{{print $2}}' {input.DEtable_sig} | sed '1d' | sort | uniq | sed '1igene' > {output.genelist_temp}
        # {enrichwrap} -i {input.DEtable_sig} -g {params.anno_path} -o {params.outpath} --nomap --minListHits 3 >{log} 2>&1
        {enrichwrap} -i {output.genelist_temp} -g {params.anno_path} -o {params.outpath} --nomap --minListHits 3 >{log} 2>&1
        """

rule scatter_plot:
    input:
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData"
    output:
        donefile = (expand("result/Splicing/{{diff_group}}/scatter/{diff_long_name}/scatter.done", diff_long_name = diff_long_name))
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing"
    threads:
        config["scatter_plot"]["cpu"]
    resources: 
        qsub_mem = config["scatter_plot"]["mem"]
    log:
        "logs/{diff_group}/scatter_plot.log"
    benchmark:
        "benchmarks/{diff_group}/scatter_plot.txt"
    shell:
        """
        {mp}
        {Rscript_splicing}
        if [[ -d result/Splicing/{wildcards.diff_group}/scatter ]]; then rm -rf result/Splicing/{wildcards.diff_group}/scatter ; fi
        Rscript {marvel_plot} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --analysis_type scatter \
            --outpath {params.outpath}  >{log} 2>&1
        """


rule tabulate_gene_junction_exp:
    input:
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData",
    output:
        exptable = expand("result/Splicing/{{diff_group}}/tabulate/{diff_long_name}/gene_SJ_expression.csv", diff_long_name = diff_long_name)
    log:
        "logs/{diff_group}/tabulate_gene_junction_exp.log"
    benchmark:
        "benchmarks/{diff_group}/tabulate_gene_junction_exp.txt"
    threads:
        config["tabulate_gene_junction_exp"]["cpu"]
    resources: 
        qsub_mem = config["tabulate_gene_junction_exp"]["mem"]
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing"
    shell:
        """
        {mp}
        {Rscript_splicing}
        if [[ -d result/Splicing/{wildcards.diff_group}/tabulate ]]; then rm -rf result/Splicing/{wildcards.diff_group}/tabulate ; fi
        Rscript {marvel_plot} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --analysis_type tabulate \
            --outpath {params.outpath}  >{log} 2>&1
        """


rule structure_plot:
    input:
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData",
    output:
        transcripts_table = expand("result/Splicing/{{diff_group}}/structure/{diff_long_name}/structure_transtripts.csv",diff_long_name = diff_long_name ),
        junction_table = expand("result/Splicing/{{diff_group}}/structure/{diff_long_name}/structure_junctions.csv",diff_long_name = diff_long_name )
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing"
    threads:
        config["structure_plot"]["cpu"]
    resources:
        qsub_mem = config["structure_plot"]["mem"]
    log:
        "logs/{diff_group}/structure_plot.log"
    benchmark:
        "benchmarks/{diff_group}/structure_plot.txt"
    shell:
        """
        {mp}
        {Rscript_splicing}
        if [[ -d result/Splicing/{wildcards.diff_group}/structure ]]; then rm -rf result/Splicing/{wildcards.diff_group}/structure ; fi
        Rscript {marvel_plot} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --analysis_type structure \
            --outpath {params.outpath}  >{log} 2>&1
        """

rule reintegrateBam:
    input:
        bam = expand("result/Splicing/STARSolo/{sample}/Bam/possorted_genome_bam.bam", sample = samples.index),
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData"
    output:
        bam = expand("result/Splicing/{{diff_group}}/sashimi/{{diff_group}}.{diff}.bam", diff = ['g1', 'g2'])
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing"
    threads:
        config['reintegrateBam']['cpu']
    log:
        "logs/{diff_group}/reintegrateBam.log"
    shell:
        """
        {mp}
        {Rscript_splicing}
        if [[ -d result/Splicing/{wildcards.diff_group}/sashimi ]]; then rm -rf result/Splicing/{wildcards.diff_group}/sashimi ; fi
        Rscript {marvel_subsetBam} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --outpath {params.outpath}  >{log} 2>&1
        """

rule sashimi:
    input:
        marvelfile = "result/Splicing/{diff_group}/MARVEL.RData",
        bam = expand("result/Splicing/{{diff_group}}/sashimi/{{diff_group}}.{diff}.bam", diff = ['g1', 'g2'])
    output:
        misodone = "result/Splicing/{diff_group}/sashimi/sashimi.done"
    params:
        treatment = lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control = lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        diff_type = lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"],
        outpath = "result/Splicing",
        intron_s = config["params"]["splicing"]["plot"]["intron_s"],
        min_count = config["params"]["splicing"]["plot"]["min_count"]
    log:
        "logs/{diff_group}/sashimi.log"
    threads:
        config["sashimi"]["cpu"]
    shell:
        """
        {mp}
        {Rscript_splicing} 
        Rscript {marvel_plot} \
            --marvelproject {input.marvelfile} \
            --diff_type {params.diff_type} \
            --diff_case \"{params.treatment}\" \
            --diff_control \"{params.control}\" \
            --analysis_type sashimi \
            --intron_s {params.intron_s} \
            --min_count {params.min_count} \
            --outpath {params.outpath}  >{log} 2>&1
        """



localrules: report
rule report:
    input:
        marvelfile = expand("result/Splicing/{diff_group}/MARVEL.RData", diff_group=diff_groups.index),
        DEtable = expand("result/Splicing/{diff_group}/diff_SJ_exp/splicing_junction_differect_expression.total.csv", diff_group=diff_groups.index),
        volcano = expand("result/Splicing/{diff_group}/volcano/diff_splicing_volcano.{volcan_long_name}.png", volcan_long_name = volcan_long_name, diff_group=diff_groups.index),
        DEtable_sig = expand("result/Splicing/{diff_group}/diff_SJ_exp/diff_SJ_exp.{diff_long_name}.csv", diff_long_name = diff_long_name, diff_group=diff_groups.index),
        donefile = expand("result/Splicing/{diff_group}/scatter/{diff_long_name}/scatter.done",diff_long_name = diff_long_name, diff_group=diff_groups.index),
        exptable = expand("result/Splicing/{diff_group}/tabulate/{diff_long_name}/gene_SJ_expression.csv", diff_long_name = diff_long_name, diff_group=diff_groups.index),
        transcripts_table = expand("result/Splicing/{diff_group}/structure/{diff_long_name}/structure_transtripts.csv",diff_long_name = diff_long_name, diff_group=diff_groups.index),
        junction_table = expand("result/Splicing/{diff_group}/structure/{diff_long_name}/structure_junctions.csv",diff_long_name = diff_long_name , diff_group=diff_groups.index ),
        misodone = expand("result/Splicing/{diff_group}/sashimi/sashimi.done", diff_group=diff_groups.index)
    output:
        report = f"report/{report_dir}/Report.html"
    params:
        result_path = "result",
        configfile = "config/config.yaml"
    log:
        "logs/report.log"
    shell:
        """
        {mp}
        /public/scRNA_works/works/luyao/test/miniconda3/envs/cloudreport/bin/python {report_script} -i {params.result_path} -c {params.configfile} >{log} 2>&1
        """
