##==================================================step1: WNN Dimensionality reduction/ clustring  =======================================================
rule wnn_clustering:
	'''
	use sctool to wnn cluster
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		log_celltyping="logs/celltyping/ATAC_celltyping.log"
	output:
		expand(
			["result/9.WNN_Clustering/wnn_groupby_cluster_resolution{resolution}_plot.pdf",
			 "result/9.WNN_Clustering/wnn_groupby_cluster_resolution{resolution}_plot.png"
			 ],resolution=config["wnn_params"]["resolution"])
	benchmark:
		"benchmarks/wnn/wnn_clustering.benchmark.txt"
	log:
		"logs/wnn/wnn_clustering.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['wnn_bclust'],
		qsub_n=1
	params:
		RNA_reduct=config["cluster_params"]["RNA_reduct1"],
		ATAC_reduct=config["cluster_params"]["ATAC_reduct1"],
		clusteringuse=config["wnn_params"]["clusteringuse"],
		resolution=config["wnn_params"]["resolution"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              --assay SCT \
              --subassay ATAC \
              --dataslot data \
              --update TRUE\
              -o result/9.WNN_Clustering \
              -d h5seurat \
              wnn --reduction SCT_{params.RNA_reduct},ATAC_{params.ATAC_reduct} \
                  --clusteringuse {params.clusteringuse} \
                  --resolution {params.resolution} \
                  --palette blindless
		'''
##==================================================== step 2: WNN cluster visualization ======================================
rule wnn_vis:
	'''
	use sctool to WNN cluster visualization
	'''
	input:
		cluster_log = "logs/wnn/wnn_clustering.log",
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		expand([
			"result/9.WNN_Clustering/visualize_cluster_by_wnn.res.{resolution}/clust_cond_freq_info.xls"
		],
			resolution=config["wnn_params"]["resolution"]
		)
	benchmark:
		"benchmarks/wnn/wnn_vis.benchmark.txt"
	log:
		"logs/wnn/wnn_vis.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['wnn_summarize'],
		qsub_n=1
	params:
		RNA_reduct = config["cluster_params"]["RNA_reduct1"],
		ATAC_reduct = config["cluster_params"]["ATAC_reduct1"],
		resolution=config["wnn_params"]["resolution"],
		groupby=config["wnn_params"]["groupby"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              --assay SCT \
              --subassay ATAC \
              --dataslot data \
              -o result/9.WNN_Clustering \
              summarize --reduct WNN_{params.RNA_reduct}_{params.ATAC_reduct}_umap \
                        --palette blindless \
                        -c wnn.res.{params.resolution} \
                        -b {params.groupby} \
                        --ptsize 0.5 \
                        --dosummary T
		'''
