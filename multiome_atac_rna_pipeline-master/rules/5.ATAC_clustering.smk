##==========================================step 1: ATAC Dimensionality reduction/ clustring ==================================================
rule ATAC_cluster:
	'''
	use sctool to ATAC cluster
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		log_vis="logs/clustering/RNA_vis.log"
	output:
		expand([
			"result/6.ATAC_Clustering/{reduct1}_Dimension_Reduction/{reduct1}_Dimension_Reduction_coordination.csv",
			"result/6.ATAC_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_Dimension_Reduction_coordination.csv",
			"result/6.ATAC_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.png",
			"result/6.ATAC_Clustering/{reduct2}_Dimension_Reduction/{reduct2}_groupby_cluster_resolution{resolution}_plot.pdf"
		],
			reduct1=config["cluster_params"]["ATAC_reduct1"],
			reduct2=config["cluster_params"]["ATAC_reduct2"],
			resolution=config["cluster_params"]["ATAC_resolution"]
		)
	benchmark:
		"benchmarks/clustering/ATAC_cluster.benchmark.txt"
	log:
		"logs/clustering/ATAC_cluster.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['clustering_bclust'],
		qsub_n=1
	params:
		reduct1 = config["cluster_params"]["ATAC_reduct1"],
		reduct2 = config["cluster_params"]["ATAC_reduct2"],
		clusteringuse = config["cluster_params"]["ATAC_clusteringuse"],
		resolution = config["cluster_params"]["ATAC_resolution"],
		palette = config["cluster_params"]["ATAC_palette"],
		component = config["cluster_params"]["ATAC_component"]
	shell:
		'''
		if  [ {params.reduct1} == "harmony" ];then sctool -i {input.inputfile} \
               -f h5seurat \
               --assay ATAC \
               --dataslot data,counts  \
               --update TRUE \
               -o result/6.ATAC_Clustering\
               -d h5seurat \
               bclust --reduct1 lsi \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid batch \
                      --component {params.component} ; fi
		sctool -i {input.inputfile} \
               -f h5seurat \
               --assay ATAC \
               --dataslot data,counts  \
               --update TRUE \
               -o result/6.ATAC_Clustering\
               -d h5seurat \
               bclust --reduct1 {params.reduct1} \
                      --reduct2 {params.reduct2} \
                      --clusteringuse {params.clusteringuse} \
                      --resolution {params.resolution} \
                      --rerun T \
                      --palette {params.palette} \
                      --batchid batch \
                      --component {params.component}
		'''
##==========================================step 2: ATAC cluster visualization ==================================================
rule ATAC_vis:
	'''
	use sctool to ATAC cluster visualization
	'''
	input:
		log_cluster="logs/clustering/ATAC_cluster.log",
		inputfile= "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		expand([
			"result/6.ATAC_Clustering/visualize_cluster_by_ATAC.{reduct2}.res.{resolution}/clust_cond_freq_info.xls"
		],
				reduct2=config["cluster_params"]["ATAC_reduct2"],
				resolution=config["cluster_params"]["ATAC_resolution"]
		)
	benchmark:
		"benchmarks/clustering/ATAC_vis.benchmark.txt"
	log:
		"logs/clustering/ATAC_vis.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['clustering_summarize'],
		qsub_n=1
	params:
		reduct2=config["cluster_params"]["ATAC_reduct2"],
		resolution=config["cluster_params"]["ATAC_resolution"],
		groupby=config["cluster_params"]["ATAC_groupby"]
	shell:
		'''
		sctool -i {input.inputfile} \
               -f h5seurat \
               --assay ATAC \
               --dataslot data,counts  \
               --update  TRUE \
               -o result/6.ATAC_Clustering \
               -d h5seurat \
               summarize --reduct ATAC_{params.reduct2} \
                         --palette blindless \
                         -c ATAC.{params.reduct2}.res.{params.resolution}\
                         -b {params.groupby} \
                         --ptsize 0.5 \
                         --dosummary T
		'''



