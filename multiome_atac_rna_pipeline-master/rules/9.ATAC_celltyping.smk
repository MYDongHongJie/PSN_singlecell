##========================================RNA celltyping mapping to ATAC =======================================================
rule ATAC_celltyping:
	'''
	use sctool to mapping ATAC celltyping
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		log_RNAcelltyping="logs/celltyping/RNA_celltyping.log"
	output:
		"result/8.ATAC_Reference_celltype/visualize_cluster_by_celltype/groupby-celltype_resolution0.4_contrast_plot.png",
		"result/8.ATAC_Reference_celltype/visualize_cluster_by_celltype/groupby-celltype_resolution0.4_contrast_plot.pdf"
	benchmark:
		"benchmarks/celltyping/ATAC_celltyping.benchmark.txt"
	log:
		"logs/celltyping/ATAC_celltyping.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['ATAC_celltyping'],
		qsub_n=1
	params:
		reduct = config["cluster_params"]["ATAC_reduct2"]
	shell:
		'''
		scVis -i {input.inputfile} \
				-f h5seurat \
				--assay ATAC \
				--dataslots data \
				-o result/8.ATAC_Reference_celltype \
				--reduct ATAC_{params.reduct} \
				clusterplot -g celltype \
									-s 0.5
		'''
