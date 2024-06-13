##==================================================WNN celltyping =======================================================
rule wnn_celltyping:
	'''
	use sctool to wnn celltyping
	'''
	input:
		log_wnnmarker="logs/wnn/wnn_ATAC_vis.log",
		inputfile= "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		expand([
			"result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_heatmap.pdf",
			"result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_plot.pdf",
			"result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_results.xls",
			"result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_simplified_celltyping_results.csv",
			"result/12.WNN_Reference_celltype/{species}ref_{refbuiltin}_top.{annolevel}_celltyping_plot.pdf"
		],
			species=config["report"]["Species"],
			refbuiltin=config["wnn_params"]["refbuiltin"],
			annolevel=config["wnn_params"]["annolevel"])
	benchmark:
		"benchmarks/wnn/wnn_celltyping.benchmark.txt"
	log:
		"logs/wnn/wnn_celltyping.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['wnn_celltyping'],
		qsub_n=1
	params:
		refbuiltin=config["wnn_params"]["refbuiltin"],
		annolevel=config["wnn_params"]["annolevel"],
		demethod=config["wnn_params"]["demethod"],
		RNA_reduct=config["cluster_params"]["RNA_reduct1"],
		ATAC_reduct=config["cluster_params"]["ATAC_reduct1"],
		resolution=config["wnn_params"]["resolution"],
		species=config["report"]["Species"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/12.WNN_Reference_celltype \
              -d h5seurat \
              --assay SCT \
              --dataslot data,counts \
              celltyping    --refbuiltin {params.refbuiltin} \
                            --annolevel {params.annolevel} \
                            --refmarker NULL \
                            --demethod {params.demethod} \
                            --clusterby wnn.res.{params.resolution} \
                            -v 0.8 \
                            -n 25 \
                            --palette blindless \
                            --reduct WNN_{params.RNA_reduct}_{params.ATAC_reduct}_umap \
                            --species {params.species} \
                            --doplot TRUE
		'''
