##==================================================Quality Control =======================================================
rule qc:
	'''
	use sctool to quality control
	'''
	input:
		inputfile="result/seurat_ob.h5seurat"
	output:
		"result/2.Count_QC/seurat_ob.h5seurat",
		"result/2.Count_QC/statitics_before_after_QC.xls"
	benchmark:
		"benchmarks/qc/qc.benchmark.txt"
	log:
		"logs/qc/qc.log"
	envmodules:
		config['envmodules']['oesinglecell']
	resources:
		qsub_mem=100,
		qsub_p=config['cpu']['qc'],
		qsub_n=1
	params:
		vars2regress=config["qc_params"]["vars2regress"],
		filters=get_qc_filter_parameters,
		features2filter=config["qc_params"]["features2filter"],
		mincell4gene=config["qc_params"]["mincell4gene"],
		rmdoublets=config["qc_params"]["rmdoublets"],
		method=config["qc_params"]["method"],
		normmeth=config["qc_params"]["normmeth"]
	shell:
		'''
		sctool --input result/seurat_ob.h5seurat \
              --informat h5seurat \
              --output result/2.Count_QC \
              --outformat h5seurat \
              --assay RNA \
              --subassay ATAC \
              --dataslot counts \
              --prefix  seurat_ob \
              --ncores {resources.qsub_p} \
              --update FALSE \
              qc --vars2regress "{params.vars2regress}" \
                --filters {params.filters[0]} \
                --lower {params.filters[1]} \
                --upper {params.filters[2]} \
                --cut.1 median \
                --cut.2 mad \
                --nfold 2 \
                --features2filter {params.features2filter} \
                --mincell4gene  {params.mincell4gene} \
                --rmdoublets {params.rmdoublets} \
                --method {params.method}\
                --ident NULL \
                --normmeth "{params.normmeth}" \
                --nvfeatures 2000 \
                --pointsize 0.1
		'''
