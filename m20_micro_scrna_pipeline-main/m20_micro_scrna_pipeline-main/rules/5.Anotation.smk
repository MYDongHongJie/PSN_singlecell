##==================================================step1: KEGG anno ===================================================
rule kegg_anno:
    '''
    KEGG annotation
    '''
    input:
        cdna=f"{config['database']['ref_genome']}/gene.fa"
    output:
        kegg_blast="result/4.Annotation/KEGG/KEGG.blast.xml",
        kegg_blastm8="result/4.Annotation/KEGG/KEGG.blastm8.xls",
        kegg_best="result/4.Annotation/KEGG/KEGG.blast.best.xls",
        kegg_anno="result/4.Annotation/KEGG/Unigene.KEGG.gene.anno.xls"
    benchmark:
        "benchmarks/anno/kegg.benchmark.txt"
    log:
        "logs/anno/kegg.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        #ref_dir = config['database']['ref_genome'],
        kegg_database=config['database']["kegg_database"],
        qsub_mem=config['kegg_anno']['mem'],
        qsub_p=config['kegg_anno']['cpu']
    shell:
        """(
        diamond blastx \
            -q {input.cdna} \
            -d {params.kegg_database} \
            -p 10 \
            -k 10 \
            --more-sensitive \
            --salltitles \
            --outfmt 5 \
            -e 1e-5 \
            -o {output.kegg_blast}
        perl  scripts/anno/alignxml2tab.pl\
            -nohead \
            -eval 1e-5 \
            -tophit 1 \
            -topmatch 1 \
            -m8 {output.kegg_blastm8}\
             {output.kegg_blast} > {output.kegg_best}
        sed -i "1i #Query_id\\tQuery_length\\tQuery_start\\tQuery_end\\tSubject_id\\tSubject_length\\tSubject_start\\tSubject_end\\tIdentity\\tGap\\tAlign_length\\tScore\\tE_value\\tDescription\" {output.kegg_best}
        python3 scripts/anno/oeanno.py \
            -i {output.kegg_best}\
            -m micro \
            -n Unigene \
            -o result/4.Annotation/KEGG/
        ) >{log} 2>&1"""

##==================================================step2: swissprot anno ==============================================
rule swissprot_anno:
    '''
        SWISSPROT annotation
    '''
    input:
        #kegg_log="logs/anno/kegg.log",
        cdna=f"{config['database']['ref_genome']}/gene.fa"
    output:
        swissprot_blast="result/4.Annotation/SWISSPROT/swissprot.blast.xml",
        swissprot_best="result/4.Annotation/SWISSPROT/swissprot.blast.best.xls",
        swissprot_anno="result/4.Annotation/SWISSPROT/swissprot.blast.anno.xls"
    benchmark:
        "benchmarks/anno/swissprot.benchmark.txt"
    log:
        "logs/anno/swissprot.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        #ref_dir=config['database']['ref_genome'],
        swissprot_database=config['database']["swissprot_database"],
        qsub_mem=config['swissprot_anno']['mem'],
        qsub_p=config['swissprot_anno']['cpu']
    shell:
        """(
        diamond blastx -q {input.cdna}  -d {params.swissprot_database} -p 10  -k 10 --more-sensitive --outfmt 5 -e 1e-5 -o {output.swissprot_blast}
        perl scripts/anno/alignxml2tab.pl -eval 1e-5 -tophit 1 -topmatch 1  {output.swissprot_blast} > {output.swissprot_best}
        perl scripts/anno/blast_tab2anno.pl {output.swissprot_best} {output.swissprot_anno}
        ) >{log} 2>&1"""

##==================================================step3: GO anno =====================================================
rule go_anno:
    '''
    GO annotation
    '''
    input:
        swissprot_log="logs/anno/swissprot.log",
        swissprot_anno="result/4.Annotation/SWISSPROT/swissprot.blast.anno.xls"
    output:
        go_stat="result/4.Annotation/GO/Unigene.GO.classification.stat.pdf",
        go_anno="result/4.Annotation/GO/Unigene.GO.gene.anno.xls"
    benchmark:
        "benchmarks/anno/go.benchmark.txt"
    log:
        "logs/anno/go.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        ref_dir=config['database']['ref_genome'],
        go_database=config['database']["go_database"],
        qsub_mem=config['go_anno']['mem'],
        qsub_p=config['go_anno']['cpu']
    shell:
        """(
        perl scripts/anno/SwissProt2GO.pl -i {input.swissprot_anno} -p Unigene -o result/4.Annotation/GO -g {params.go_database}  
        perl scripts/anno/gene_ontology_graph.pl -i {output.go_anno} -s Unigene -o result/4.Annotation/GO
        ) >{log} 2>&1"""

##==================================================step4: CARD anno =====================================================
rule card_anno:
    '''
    CARD annotation
    '''
    input:
        cdna=f"{config['database']['ref_genome']}/gene.fa"
    output:
        card_anno="result/4.Annotation/CARD/UnigeneCARD.anno.xls"
    benchmark:
        "benchmarks/anno/card.benchmark.txt"
    log:
        "logs/anno/card.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        ref_dir=config['database']['ref_genome'],
        card_database=config['database']["card_database"],
        qsub_mem=config['card_anno']['mem'],
        qsub_p=config['card_anno']['cpu']
    shell:
        """(
        perl scripts/anno/blast_card_local.pl -p {input.cdna}  -s UnigeneCARD -o result/4.Annotation/CARD/ -d {params.card_database} -type nucl -e 1e-10
        perl scripts/anno/statics_ARO_local.pl -i {output.card_anno} -o result/4.Annotation/CARD -r UnigeneCARD
        ) >{log} 2>&1"""

##==================================================step5: CAZy anno ====================================================
rule cazy_anno:
    '''
    CAZy annotation
    '''
    input:
        cdna=f"{config['database']['ref_genome']}/gene.fa"
    output:
        cazy_blast="result/4.Annotation/CAZy/CAZy.blast.xml",
        cazy_best="result/4.Annotation/CAZy/CAZy.blast.best.xls",
        cazt_plot="result/4.Annotation/CAZy/CAZy_class_plot.pdf"
    benchmark:
        "benchmarks/anno/cazy.benchmark.txt"
    log:
        "logs/anno/cazy.log"
    envmodules:
        config['envmodules']['oesinglecell']
    params:
        ref_dir=config['database']['ref_genome'],
        cazy_database=config['database']["cazy_database"],
        cazy_family_info=config['database']["cazy_family_info"],
        cazy_class_info=config['database']["cazy_class_info"],
        qsub_mem=config['cazy_anno']['mem'],
        qsub_p=config['cazy_anno']['cpu']
    shell:
        """(
        diamond blastx -q {input.cdna}  -d {params.cazy_database} -p 10 -k 10 --more-sensitive --outfmt 5 -e 1e-5 -o {output.cazy_blast}
        perl scripts/anno/alignxml2tab.pl -eval 1e-5 -tophit 1 -topmatch 1 {output.cazy_blast} > {output.cazy_best}
        sed -i "s/\#//g\" {output.cazy_best}
        Rscript scripts/anno/CAZy_plot.r -i {output.cazy_best} -f {params.cazy_family_info} -c {params.cazy_class_info} -o result/4.Annotation/CAZy
        ) >{log} 2>&1"""
