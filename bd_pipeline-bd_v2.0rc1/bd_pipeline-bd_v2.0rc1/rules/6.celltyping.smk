rule celltyping:
    """
    Dimension reduction & clustering  
    """
    input:
        cluster_rds = "result/Count_QC/filtered.h5seurat",
        temp_file="clustering_h5seurat_touched.check" # to make sure running after Clustering
    output:
        expand([
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_statistics.xls",
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_results.xls",
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_simplified_celltyping_results.csv",
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_heatmap.png",
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_heatmap.pdf",
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_plot.png",
            "result/Reference_celltype/{ref_species}ref_{db_name}_main_celltyping_plot.pdf",
            "result/Reference_celltype/{ref_species}ref_{db_name}_top.main_celltyping_plot.png",
            "result/Reference_celltype/{ref_species}ref_{db_name}_top.main_celltyping_plot.pdf"
            ],
            db_name=celltypingdb.split('/')[-1].split('.')[0],
            ref_species=ref_species,
            #level="main"
            #      clustering_method2=config["params"]["clustering_method2"],
            #      clustering_method1=config["params"]["clustering_method1"],
            #      resolution=config["params"]["resolution"]
        ),
        temp_file=temp("celltyping_h5seurat_touched.check")
    params:
        outdir="result/Reference_celltype",
        reduct2=config["report_params"]["reduct2"],
        #  clustering_method1=config["params"]["clustering_method1"],
        pointsize_params = get_pointsize_params,
        resolution=config["report_params"]["resolution"],
        #celltypingdb=config["report_params"]["celltypingdb"],
        #celltypingdb_dir="/data/database/celltype_refdata/logNorm_rds/"+{params.celltypingdb}+".rds",
        #level="main"
    benchmark:
        "benchmarks/celltyping.benchmark.txt"
#    log:  "logs/celltyping/celltyping.log"
    resources:
        qsub_mem=30,
        qsub_p=config['cpu']['celltyping']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript  {script_dir}/sctool \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  -d h5seurat \
  --update T \
  --assay RNA \
  --dataslot counts \
  celltyping \
  -r {celltypingdb} \
  --annolevel main \
  --usecluster F \
  --demethod classic \
  {params.pointsize_params} \
  -n 25 \
  --reduct {params.reduct2} \
  --species {ref_species} 
touch -r benchmarks/qc.benchmark.txt  result/Count_QC/filtered.h5seurat && 
touch celltyping_h5seurat_touched.check && 
touch celltyping_h5seurat_touched.`date +%Y%m%d-%H%M`
        """
