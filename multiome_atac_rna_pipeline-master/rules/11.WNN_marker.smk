##==================================================step1: WNN find RNA marker=======================================================
rule wnn_RNA_marker:
	'''
	use sctool to find wnn RNA marker
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		wnn_rna_marker="logs/wnn/wnn_vis.log"
	output:
		"result/10.WNN_RNA_Marker/all_markers_for_each_cluster.xls",
		"result/10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls"
	benchmark:
		"benchmarks/wnn/wnn_RNA_marker.benchmark.txt"
	log:
		"logs/wnn/wnn_RNA_marker.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['wnn_findallmarkers'],
		qsub_n=1
	params:
		topn_marker=config["wnn_params"]["wnn_topn_marker"],
		resolution=config["wnn_params"]["resolution"],
		test_method=config["wnn_params"]["wnn_test_method"]
	shell:
		'''
		sctool -i {input.inputfile} \
               -f h5seurat \
               -o result/10.WNN_RNA_Marker \
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
                              --cluster_name  wnn.res.{params.resolution}\
                              --test {params.test_method}
		'''
##================================================ step2: WNN RNA marker visualization ======================================
rule wnn_RNA_vis:
	'''
	use sctool to wnn RNA marker visualization
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		wnn_marker_vis="logs/wnn/wnn_RNA_marker.log"
	output:
		"result/10.WNN_RNA_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
		"result/10.WNN_RNA_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf"
	benchmark:
		"benchmarks/wnn/wnn_RNA_vis.benchmark.txt"
	log:
		"logs/wnn/wnn_RNA_vis.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['wnn_visualize'],
		qsub_n=1
	params:
		topn=config["wnn_params"]["wnn_topn"],
		topby=config["wnn_params"]["wnn_topby"],
		vismethod=config["wnn_params"]["wnn_vismethod"],
		resolution=config["wnn_params"]["resolution"],
		RNA_reduct=config["cluster_params"]["RNA_reduct1"],
		ATAC_reduct=config["cluster_params"]["ATAC_reduct1"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/10.WNN_RNA_Marker \
              --assay SCT \
              --dataslot data,scale.data \
              visualize --markers result/10.WNN_RNA_Marker/top10_markers_for_each_cluster.xls \
                       --topn {params.topn} \
                        --topby {params.topby} \
                        --vismethod {params.vismethod} \
                        --groupby wnn.res.{params.resolution} \
                        --pointsize 0 \
                        --reduct WNN_{params.RNA_reduct}_{params.ATAC_reduct}_umap
		'''
##=============================================== step 3: WNN find ATAC marker ============================================
rule wnn_ATAC_marker:
	'''
	use sctool to find wnn ATAC marker
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		wnn_atac_marker="logs/wnn/wnn_RNA_vis.log"
	output:
		"result/11.WNN_ATAC_Marker/all_markers_for_each_cluster.xls",
		"result/11.WNN_ATAC_Marker/top10_markers_for_each_cluster.xls"
	benchmark:
		"benchmarks/wnn/wnn_ATAC_marker.benchmark.txt"
	log:
		"logs/wnn/wnn_ATAC_marker.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['wnn_findallmarkers'],
		qsub_n=1
	params:
		topn_marker=config["wnn_params"]["wnn_topn_marker"],
		resolution=config["wnn_params"]["resolution"],
		test_method=config["wnn_params"]["wnn_test_method"]
	shell:
		'''
		sctool -i {input.inputfile} \
               -f h5seurat \
               -o result/11.WNN_ATAC_Marker \
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
                              --cluster_name  wnn.res.{params.resolution}\
                              --test {params.test_method}
		'''
##================================================= step4: WNN ATAC marker visualization ====================================
rule wnn_ATAC_vis:
	'''
	use sctool to wnn ATAC marker visualization
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		wnn_marker_vis = "logs/wnn/wnn_ATAC_marker.log"
	output:
		"result/11.WNN_ATAC_Marker/markers_vis4cluster1/marker_gene_featureplot.pdf",
		"result/11.WNN_ATAC_Marker/markers_vis4cluster1/marker_gene_violin_plot.pdf"
	benchmark:
		"benchmarks/wnn/wnn_ATAC_vis.benchmark.txt"
	log:
		"logs/wnn/wnn_ATAC_vis.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem = 50,
		qsub_p = config['cpu']['wnn_visualize'],
		qsub_n = 1
	params:
		topn = config["wnn_params"]["wnn_topn"],
		topby = config["wnn_params"]["wnn_topby"],
		vismethod = config["wnn_params"]["wnn_vismethod"],
		resolution = config["wnn_params"]["resolution"],
		RNA_reduct = config["cluster_params"]["RNA_reduct1"],
		ATAC_reduct = config["cluster_params"]["ATAC_reduct1"]
	shell:
		'''
		sctool -i {input.inputfile} \
			  -f h5seurat \
			  -o result/11.WNN_ATAC_Marker \
			  --assay ATAC \
			  --dataslot data,counts \
			  visualize --markers result/11.WNN_ATAC_Marker/top10_markers_for_each_cluster.xls \
					   --topn {params.topn} \
						--topby {params.topby} \
						--vismethod {params.vismethod} \
						--groupby wnn.res.{params.resolution} \
						--pointsize 0 \
						--reduct WNN_{params.RNA_reduct}_{params.ATAC_reduct}_umap
		'''


