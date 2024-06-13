##==================================================step1: find marker =======================================================
rule RNA_findmarker:
	'''
	use sctool to find RNA marker
	'''
	input:
		log_cluster="logs/clustering/ATAC_vis.log",
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		"result/4.RNA_Marker/all_markers_for_each_cluster.xls",
		"result/4.RNA_Marker/top10_markers_for_each_cluster.xls"
	benchmark:
		"benchmarks/marker/RNA_findmarker.benchmark.txt"
	log:
		"logs/marker/RNA_findmarker.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['marker_findallmarkers'],
		qsub_n=1
	params:
		topn_marker=config["marker_params"]["RNA_topn_marker"],
		test_method=config["marker_params"]["RNA_test_method"],
		reduct2=config["cluster_params"]["RNA_reduct2"],
		resolution=config["cluster_params"]["RNA_resolution"]
	shell:
		'''
		sctool -i {input.inputfile} \
               -f h5seurat \
               -o result/4.RNA_Marker \
               --assay SCT \
               --dataslot data,counts \
               -j {resources.qsub_p} \
               findallmarkers --min_pct1 0.5 \
                              --max_pct2 0.5  \
                              --pct_fold 2 \
                              --topn_marker {params.topn_marker} \
                              --avg_log2FC 1 \
                              --pvalue 0.05 \
                              --strict F \
                              --cluster_name  SCT.{params.reduct2}.res.{params.resolution}\
                              --test {params.test_method}
		'''
##==================================================step 2: marker visualization=============================================
rule RNA_markervis:
	'''
	use sctool to RNA marker visualization
	'''
	input:
		log_marker="logs/marker/RNA_findmarker.log",
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		"result/4.RNA_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
		"result/4.RNA_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf"
	benchmark:
		"benchmarks/marker/RNA_markervis.benchmark.txt"
	log:
		"logs/marker/RNA_markervis.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['marker_visualize'],
		qsub_n=1
	params:
		topn=config["marker_params"]["RNA_topn"],
		topby=config["marker_params"]["RNA_topby"],
		vismethod=config["marker_params"]["RNA_vismethod"],
		reduct2=config["cluster_params"]["RNA_reduct2"],
		resolution=config["cluster_params"]["RNA_resolution"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/4.RNA_Marker \
              --assay SCT \
              --dataslot data,scale.data \
              visualize --markers result/4.RNA_Marker/top10_markers_for_each_cluster.xls \
                        --topn {params.topn} \
                        --topby {params.topby} \
                        --vismethod {params.vismethod} \
                        --groupby SCT.{params.reduct2}.res.{params.resolution} \
                        --pointsize 0 \
                        --reduct SCT_{params.reduct2}
		'''
