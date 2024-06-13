##================================================== RNA celltyping=======================================================
rule RNA_celltyping:
	'''
	use sctool to identification of cell types
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		log_marker="logs/marker/ATAC_markervis.log"
	output:
		expand([
			"result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_heatmap.pdf",
			"result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_plot.pdf",
			"result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_celltyping_results.xls",
			"result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_{annolevel}_simplified_celltyping_results.csv",
			"result/5.RNA_Reference_celltype/{species}ref_{refbuiltin}_top.{annolevel}_celltyping_plot.pdf"
		],
			species=config["report"]["Species"],
			refbuiltin=config["celltyping_params"]["refbuiltin"],
			annolevel=config["celltyping_params"]["annolevel"])
	benchmark:
		"benchmarks/celltyping/RNA_celltyping.benchmark.txt"
	log:
		"logs/celltyping/RNA_celltyping.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['RNA_celltyping'],
		qsub_n=1
	params:
		refbuiltin=config["celltyping_params"]["refbuiltin"],
		annolevel=config["celltyping_params"]["annolevel"],
		demethod=config["celltyping_params"]["demethod"],
		reduct=config["cluster_params"]["RNA_reduct2"],
		resolution=config["cluster_params"]["RNA_resolution"],
		species=config["report"]["Species"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/5.RNA_Reference_celltype \
              -d h5seurat \
              --assay SCT \
              --dataslot data,counts \
              celltyping    --refbuiltin {params.refbuiltin} \
                            --annolevel {params.annolevel} \
                            --refmarker NULL \
                            --demethod {params.demethod} \
                            --clusterby SCT.{params.reduct}.res.{params.resolution} \
                            -v 0.8 \
                            -n 25 \
                            --palette blindless \
                            --reduct SCT_{params.reduct} \
                            --species {params.species} \
                            --doplot TRUE
		'''
