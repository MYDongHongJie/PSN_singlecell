import os
import pandas as pd
localrules: cellranger_aggr_prep
rule cellranger_aggr_prep:
    """
    getting cellranger aggr libraries
    """ 
    input:
        samples_file = config['report_params']['samples_file']
    output:
        #meta = os.path.abspath("config/libraries.csv")
        libraries_file = config['report_params']['libraries_file']
    run:
        samples = pd.read_csv(input.samples_file, dtype=str).set_index("sampleid", drop=False)
        libraries = pd.DataFrame(columns=['library_id', 'molecule_h5'])
        for sample in list(samples.index):
            libraries.loc[sample, 'library_id'] = sample
            libraries.loc[sample, 'molecule_h5'] = f"{sample}/outs/molecule_info.h5"
        #    libraries.loc[sample, 'cloupe_file'] = f"{sample}/outs/cloupe.cloupe"
        #    libraries.loc[sample, 'spatial_folder'] = f"{sample}/outs/spatial"
            libraries.to_csv(output.libraries_file, encoding='utf-8', index=False)

rule cellranger_aggr:
        """
        cellranger_aggr
        """

        input:
            cellranger_outfiles = expand(["result/cellranger/{sample}/outs/web_summary.html",
                                        "result/cellranger/{sample}/outs/metrics_summary.csv",
                                        "result/cellranger/{sample}/outs/cloupe.cloupe",
                                        "result/cellranger/{sample}/outs/raw_feature_bc_matrix.h5",
                                        "result/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5"], sample=samples.index),
            libraries_file = config['report_params']['libraries_file']
        output:
            "result/cellranger/aggr/outs/web_summary.html",
            "result/cellranger/aggr/outs/count/cloupe.cloupe",
            tmpfiles = temp(["result/cellranger/aggr/_cmdline",
                         "result/cellranger/aggr/_filelist",
                         "result/cellranger/aggr/_finalstate",
                         "result/cellranger/aggr/_invocation",
                         "result/cellranger/aggr/_jobmode",
                         "result/cellranger/aggr/_log",
                         "result/cellranger/aggr/_mrosource",
                         "result/cellranger/aggr/_perf",
                         "result/cellranger/aggr/_sitecheck",
                         "result/cellranger/aggr/_tags",
                         "result/cellranger/aggr/_timestamp",
                         "result/cellranger/aggr/_uuid",
                         "result/cellranger/aggr/_vdrkill",
                         "result/cellranger/aggr/_versions",
                         directory("result/cellranger/aggr/SC_RNA_AGGREGATOR_CS"),
                         "result/cellranger/aggr/aggr.mri.tgz"])
        params:
            outdir="result/report/cellranger",
            #normalize=config["params"]["aggr_normalize"]
        benchmark:
            "benchmarks/cellranger_aggr.benchmark.txt"
        resources:
            qsub_mem=64,
            qsub_p=config['cpu']['cellranger_aggr']
        envmodules:
            config['cellranger_params']['envmodules']
        shell:
           """
           cd  result/cellranger/;
           rm -rf aggr;
           cellranger aggr \
               --id=aggr  \
               --csv=../../{input.libraries_file}  --normalize=none  \
           """
