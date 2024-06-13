##==================================================create seurat object =======================================================
rule create:
	'''
	run sctools to create multimodal seurat file
	'''
	input:
		metadatafile="config/metadata.csv",
		cellranger_out=expand(
			[
				"result/cellranger/{sample}/outs/gex_molecule_info.h5",
				"result/cellranger/{sample}/outs/web_summary.html",
				"result/cellranger/{sample}/outs/summary.csv"
			],sample=samples.index)
	output:
		"result/seurat_ob.h5seurat"
	benchmark:
		"benchmarks/create/create.benchmark.txt"
	log:
		"logs/create/create.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=60,
		qsub_p=config['cpu']['create'],
		qsub_n=1
	params:
		outformat="h5seurat",
		metrics=config['create_params']['metrics'],
		reference=config["report"]["Reference"],
		gset=get_gset_param
	shell:
		'''
		sctool --input  result/cellranger/ \
               --output result/ \
               --outformat {params.outformat} \
               --assay RNA \
               --subassay ATAC \
               --prefix  seurat_ob \
               --ncores {resources.qsub_p} \
               multimodal --source  h5 \
                          --metadata {input.metadatafile} \
                          --metrics {params.metrics} \
                          {params.gset}  \
                          --gcolumn 2 \
                          --feature.meta {params.reference}/genes/genes.gtf.gz \
                          --cell.meta TRUE \
                          --fragment TRUE
		'''



