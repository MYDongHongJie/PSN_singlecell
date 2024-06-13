##============================================= ATAC differential expression analysis =======================================================
rule diffexp_ATAC:
	'''
	use sctool to run differential expression analysis with ATAC assay
	'''
	input:
		log_RNAd_diffexp="logs/diffexp/diff_enrich.log",
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		"result/14.Diffexp_ATAC/{diff_group}-all_diffexp_genes.xls",
		"result/14.Diffexp_ATAC/{diff_group}-diff-pval-0.05-FC-1.5.xls"
	params:
		#FC=config["diffexp_params"]["FC"],
		#pvalue=config["diffexp_params"]["pvalue"],
		test=config["diffexp_params"]["test"],
		treatment=lambda w:diff_groups.loc["{}".format(w.diff_group)]["treatment"],
		control=lambda w:diff_groups.loc["{}".format(w.diff_group)]["control"],
		type=lambda w:diff_groups.loc["{}".format(w.diff_group)]["type"]
	benchmark:
		"benchmarks/diffexp/diff_gene_for_{groupname}-pval-{pvalue}-FC-{foldchange}_ATAC.benchmark.txt"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['diffexp'],
		qsub_n=1
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/14.Diffexp_ATAC \
              --assay ATAC \
              -d h5seurat \
              --dataslot data,counts  \
              -j {resources.qsub_p} \
              diffexp --contrast {params.type}:{params.treatment}:{params.control} \
                     --FC 1.5 \
                     --pvalue 0.05  \
                     --test {params.test} \
                     -v {params.type} \
                     -u {params.treatment},{params.control}
		'''
##==================================================step2: plot diff_heatmap ==============================================
rule diffexp_ATAC_heatmap:
	'''
	use sctool to plot diff_heatmap with ATAC assay
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		diff_DEG= "result/14.Diffexp_ATAC/{diff_group}-diff-pval-0.05-FC-1.5.xls"
	output:
		"result/14.Diffexp_ATAC/top25_{diff_group}_heatmap.pdf",
		"result/14.Diffexp_ATAC/top25_{diff_group}_heatmap.png",
		"result/14.Diffexp_ATAC/top25_{diff_group}_genes.xls"
	benchmark:
		"benchmarks/diffexp/diffexp_ATAC_heatmap_{diff_group}.benchmark.txt"
	log:
		"logs/diffexp/diffexp_ATAC_heatmap_{diff_group}.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['diffexp_heatmap'],
		qsub_n=1
	params:
		type = lambda w:diff_groups.loc["{}".format(w.diff_group)]["type"],
		treatment=lambda w:diff_groups.loc["{}".format(w.diff_group)]["treatment"],
		control=lambda w:diff_groups.loc["{}".format(w.diff_group)]["control"]
	shell:
		'''
		sctool -i {input.inputfile} \
              -f h5seurat \
              -o result/14.Diffexp_ATAC \
              --assay ATAC \
              --dataslot data \
              visualize --diffGene {input.diff_DEG} \
                        --topn 25 \
                        --vismethod diff_heatmap \
                        --groupby {params.type} \
                        -v {params.type} \
                        -u {params.treatment},{params.control}
		'''
