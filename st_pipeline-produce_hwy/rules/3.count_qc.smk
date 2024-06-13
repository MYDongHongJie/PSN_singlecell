import os
#=======================================================================================================================
rule  QC_seurat:
    """
    QC_seurat
    """
    input:
         raw_rds = "result/sctools/spatial.rds"
        # mt_gene = os.path.join(config['database']['reference'],"MT_genelist.txt")
    output:
          "result/count_qc/statitics_for_QC.xls",
          "result/count_qc/filtered_seurat.rds",
          "result/count_qc/QC_featureplot_for_nCount_Spatial.pdf",
          "result/count_qc/QC_featureplot_for_nCount_Spatial.png",
          "result/count_qc/QC_featureplot_for_nFeature_Spatial.pdf",
          "result/count_qc/QC_featureplot_for_nFeature_Spatial.png"
    params:
          outdir="result/count_qc",
          method= "nFeature_Spatial,nCount_Spatial,percent.mito" if (config['params']['library_type'] == "fresh" and os.path.isfile(os.path.join(config['database']['reference'],'MT_genelist.gmt'))) or (config['params']['library_type'] == "cytassist" and config['database']['reference'] == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A")  else "nFeature_Spatial,nCount_Spatial",
          low_thresh= "NULL,NULL,NULL" if (config['params']['library_type'] == "fresh" and os.path.isfile(os.path.join(config['database']['reference'],'MT_genelist.gmt'))) or (config['params']['library_type'] == "cytassist" and config['database']['reference'] == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A") else "NULL,NULL",
          high_thresh= "NULL,NULL,NULL" if (config['params']['library_type'] == "fresh" and os.path.isfile(os.path.join(config['database']['reference'],'MT_genelist.gmt'))) or (config['params']['library_type'] == "cytassist" and config['database']['reference'] == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A") else "NULL,NULL",
          vars2regress = "percent.mito" if (config['params']['library_type'] == "fresh" and os.path.isfile(os.path.join(config['database']['reference'],'MT_genelist.gmt'))) or (config['params']['library_type'] == "cytassist" and config['database']['reference'] == "/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A") else "NULL",
          sct_split= "--sct_split  sampleid" if config['params']['sct_split'] and config['report']['Sample_Num']!=1  else " ",
          crop= "TRUE" if config['params']['library_type'] == "cytassist" else "FALSE",
          qsub_mem=config['QC_seurat']['mem'],
          qsub_p=config['QC_seurat']['cpu']
    benchmark:
        "benchmarks/QC_seurat.benchmark.txt"
    log:
        "logs/qc/run.snakemake.output"
    envmodules:
        envmodules["oesinglecell"]
    shell:
        """
        sctool \
              --input {input.raw_rds}\
              --informat rds \
              --output {params.outdir} \
              --assay Spatial \
              --outformat rds \
              --ncores {params.qsub_p} \
              --update FALSE \
              --prefix filtered_seurat\
              qc \
                  --vars2regress  {params.vars2regress} \
                  --filters {params.method} \
                  --lower  {params.low_thresh} \
                  --upper {params.high_thresh} \
                  --normmeth {config[params][normmeth]} \
                   {params.sct_split} \
                  --rmdoublets F \
                  --nvfeatures {config[params][nvfeatures]} \
                  --crop {params.crop} \
                  --pointsize 0  >&2 2>{log}    
        """
