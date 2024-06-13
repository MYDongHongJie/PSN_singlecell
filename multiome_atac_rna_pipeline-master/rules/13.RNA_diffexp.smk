##============================================= step1: Differential expression analysis =======================================================
rule diffexp_RNA:
	'''
	use sctool to run differential expression analysis with RNA assay
	'''
	input:
		log_wnn="logs/wnn/wnn_celltyping.log",
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat"
	output:
		"result/13.Diffexp_RNA/{diff_group}-all_diffexp_genes.xls",
		"result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5.xls"
	params:
		#FC=config["diffexp_params"]["FC"],
		#pvalue=config["diffexp_params"]["pvalue"],
		test=config["diffexp_params"]["test"],
		treatment=lambda w:diff_groups.loc["{}".format(w.diff_group)]["treatment"],
		control=lambda w:diff_groups.loc["{}".format(w.diff_group)]["control"],
		type=lambda w:diff_groups.loc["{}".format(w.diff_group)]["type"]
	benchmark:
		"benchmarks/diffexp_RNA_{diff_group}.benchmark"
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
              -o result/13.Diffexp_RNA \
              --assay SCT \
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
rule diffexp_RNA_heatmap:
	'''
	use sctool to plot diff_heatmap with RNA assay
	'''
	input:
		inputfile = "result/2.Count_QC/seurat_ob.h5seurat",
		diff_DEG= "result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5.xls"
	output:
		"result/13.Diffexp_RNA/top25_{diff_group}_heatmap.pdf",
		"result/13.Diffexp_RNA/top25_{diff_group}_heatmap.png",
		"result/13.Diffexp_RNA/top25_{diff_group}_genes.xls"
	benchmark:
		"benchmarks/diffexp_RNA_heatmap_{diff_group}.benchmark.txt"
	log:
		"logs/diffexp/diffexp_RNA_heatmap_{diff_group}.log"
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
              -o result/13.Diffexp_RNA \
              --assay SCT \
              --dataslot data,scale.data \
              visualize --diffGene {input.diff_DEG} \
                        --topn 25 \
                        --vismethod diff_heatmap \
                        --groupby {params.type} \
                        -v {params.type} \
                        -u {params.treatment},{params.control}
		'''
##=================================================step 3: Annotation for diff expression genes ==================================
rule diffgene_annotation:
	'''
	Add annotation for diff expression genes
	'''
	input:
		all_DEG = "result/13.Diffexp_RNA/{diff_group}-all_diffexp_genes.xls",
		diff_DEG = "result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5.xls",
		log_diffexp="logs/diffexp_RNA_heatmap_{diff_group}.log"
	output:
		"result/13.Diffexp_RNA/{diff_group}-all_diffexp_genes_anno.xls",
		"result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5_anno.xls"
	benchmark:
		"benchmarks/diff_gene_annotation_for_{diff_group}.benchmark.txt"
	log:
		"logs/diffexp/diffgene_annotation_{groupname}.log"
	envmodules:
		config['envmodules']['oesinglecell']
	params:
		species=config["report"]["Species"]
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['diffexp_annotation'],
		qsub_n=1
	shell:
		'''
		Rscript scripts/add_annotation.r -i {input.diff_DEG}  -e {params.species}
		Rscript scripts/add_annotation.r -i {input.all_DEG}  -e {params.species}
		'''
##================================================step4: GO and KEGG enrich===================================================
localrules: diff_gene_enrich
rule diff_gene_enrich:
	'''
	Add GO and KEGG enrich
	'''
	input:
		diff_DEG=expand(["result/13.Diffexp_RNA/{diff_group}-diff-pval-0.05-FC-1.5.xls"],diff_group=diff_groups.index),
		log_annotation=expand(["logs/diffgene_annotation_{diff_group}.log"],diff_group=diff_groups.index)
	output:
		"result/15.enrichment/GO_enrichment/enrichment_go.xls",
		"result/15.enrichment/KEGG_enrichment/enrichment_kegg.xls"
	benchmark:
		"benchmarks/diffexp/diff_gene_enrich.benchmark.txt"
	log:
		"logs/diffexp/diff_enrich.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=50,
		qsub_p=config['cpu']['enrich'],
		qsub_n=1
	params:
		reference=config["report"]["Reference"]
	shell:
		'''
		perl  scripts/enrichment/5.2.enrich_go_kegg.pl -infile result/13.Diffexp_RNA/*-vs-*-diff-*.xls \
			-go_bg  {params.reference}/annotation/gene_go.backgroud.xls \
			-category scripts/enrichment/enrich_background/category.xls   \
			-kegg_bg {params.reference}/annotation/gene_kegg.backgroud.xls \
			-outdir result/15.enrichment  \
			-shelldir result/15.enrichment/enrichment_sh \
			-thread {resources.qsub_p} \
			-queue big
		'''

