##==================================================step1:  find marker =======================================================
rule ATAC_findmarker:
	'''
	use sctool to find ATAC marker peak
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		log_cluster="logs/marker/RNA_markervis.log"
	output:
		"result/7.ATAC_Marker/all_markers_for_each_cluster.xls",
		"result/7.ATAC_Marker/top10_markers_for_each_cluster.xls"
	benchmark:
		"benchmarks/marker/ATAC_findmarker.benchmark.txt"
	log:
		"logs/marker/ATAC_findmarker.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['marker_findallmarkers'],
		qsub_n=1
	params:
		topn_marker=config["marker_params"]["ATAC_topn_marker"],
		test_method=config["marker_params"]["ATAC_test_method"],
		reduct2=config["cluster_params"]["ATAC_reduct2"],
		resolution=config["cluster_params"]["ATAC_resolution"]
	shell:
		'''
		sctool -i {input.inputfile} \
               -f h5seurat \
               -o result/7.ATAC_Marker \
               --assay ATAC \
               --dataslot data,counts \
               -j {resources.qsub_p} \
               findallmarkers --min_pct1 0.5 \
                              --max_pct2 0.5  \
                              --pct_fold 2 \
                              --topn_marker {params.topn_marker} \
                              --avg_log2FC 1 \
                              --pvalue 0.05 \
                              --strict F \
                              --cluster_name  ATAC.{params.reduct2}.res.{params.resolution}\
                              --test {params.test_method} 
		'''
##=======================================================step 2: marker visualization ========================================
rule ATAC_markervis:
	'''
	use sctool to ATAC marker visualization
	'''
	input:
		log_marker="logs/marker/ATAC_findmarker.log",
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		"result/7.ATAC_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
		"result/7.ATAC_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf"
	benchmark:
		"benchmarks/marker/ATAC_markervis.benchmark.txt"
	log:
		"logs/marker/ATAC_markervis.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['marker_visualize'],
		qsub_n=1
	params:
		topn=config["marker_params"]["ATAC_topn"],
		topby=config["marker_params"]["ATAC_topby"],
		vismethod=config["marker_params"]["ATAC_vismethod"],
		reduct2=config["cluster_params"]["ATAC_reduct2"],
		resolution=config["cluster_params"]["ATAC_resolution"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/7.ATAC_Marker \
              --assay ATAC \
              --dataslot data,counts \
              visualize --markers result/7.ATAC_Marker/top10_markers_for_each_cluster.xls \
                        --topn {params.topn} \
                        --topby {params.topby} \
                        --vismethod {params.vismethod} \
                        --groupby ATAC.{params.reduct2}.res.{params.resolution} \
                        --pointsize 0 \
                        --reduct SCT_{params.reduct2}
		'''