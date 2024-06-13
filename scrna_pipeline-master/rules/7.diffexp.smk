rule  diffexp:
    """
    diff expressed gene select
    """
    input:
        filtered_rds = "result/Count_QC/filtered.h5seurat"
    output:
        "result/Diffexp/{diff_group}-all_diffexp_genes.xls",
        "result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5.xls"
    params:
        outdir="result/Diffexp/",
        foldchange=1.5,
        pvalue=0.05,
        treatment=lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control=lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        type=lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"]
    benchmark:
        "benchmarks/diffexpr_{diff_group}.benchmark.txt"
#    log:  "logs/diffexp/diffexp_{diff_group}.log"
    resources:
        qsub_mem=50,
        qsub_p=config['cpu']['diffexp']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript {script_dir}/sctool \
    -i {input.filtered_rds} \
    -f h5seurat \
    -o {params.outdir} \
    --assay RNA \
    --dataslot data,counts \
    -j 10 \
    diffexp \
    -c {params.type}:{params.treatment}:{params.control} \
    -k {params.foldchange} \
    -p {params.pvalue} \
    -e presto
        """
# ==================================================================================================
rule diffexp_heatmap:
    """
    Visualizing top DEG in heatmap
    """
    input:
        cluster_rds = "result/Count_QC/filtered.h5seurat",
        #temp_file="clustering_h5seurat_touched.check",
        DEG = "result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5.xls"
    output:
        "result/Diffexp/top20_{diff_group}_genes.xls",
        "result/Diffexp/top20_{diff_group}_heatmap.png",
        "result/Diffexp/top20_{diff_group}_heatmap.pdf"
    params:
        outdir="result/Diffexp/",
        treatment=lambda w: diff_groups.loc["{}".format(w.diff_group)]["treatment"],
        control=lambda w: diff_groups.loc["{}".format(w.diff_group)]["control"],
        type=lambda w: diff_groups.loc["{}".format(w.diff_group)]["type"]
    benchmark:
        "benchmarks/diffexp_heatmap_{diff_group}.benchmark.txt"
    #log:  "logs/diffexp/diffexp_heatmap_{diff_group}.log"
    resources:
        qsub_mem=50,
        qsub_p=config['cpu']['diffexp_heatmap']
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript {script_dir}/scVis \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o {params.outdir} \
  -t 10 \
  --assay RNA \
  --slot data,scale.data \
  --predicate "{params.type} %in% c(\\"{params.treatment}\\",\\"{params.control}\\")" \
  diff_heatmap \
  -d {input.DEG} \
  -n 20 \
  -g {params.type} \
  --group_colors customecol2 \
  --sample_ratio 0.8
        """

# ==================================================================================================
import os
#localrules: diffexp_annotation
rule diffexp_annotation:
    """
    add Annotation for diff expression genes
    """
    input:
        all_DEG="result/Diffexp/{diff_group}-all_diffexp_genes.xls",
        diff_DEG="result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5.xls",
        gene_annotation= os.path.join(reference,"annotation/gene_annotation.xls")
    output:
        "result/Diffexp/{diff_group}-all_diffexp_genes_anno.xls",
        "result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5_anno.xls"
    params:
        outdir="result/Diffexp",
    benchmark:
        "benchmarks/diffexp_annotation_{diff_group}.benchmark.txt"
    resources:
        qsub_mem=10,
        qsub_p=1
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
Rscript  {script_dir}/sctool annotation \
  -g {input.all_DEG} \
  --anno {input.gene_annotation}
Rscript  {script_dir}/sctool annotation \
  -g {input.diff_DEG} \
  --anno {input.gene_annotation}
        """
# =============================================================================================
#localrules: diffexp_enrich
rule diffexp_enrich:
    """
    add Annotation for diff expression genes
    """
    input:
        diff_DEG = expand(["result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5.xls"],diff_group=diff_groups.index)
    output:
        #expand([
            #"result/Diffexp/enrichment/GO_enrichment/{diff_group}/enrichment-go-{diff_group}-Total.xls",
            "result/Diffexp/enrichment/GO_enrichment/enrichment_go.xls",
            #"result/Diffexp/enrichment/KEGG_enrichment/{diff_group}/enrichment-kegg-{diff_group}-Total.xls"
            "result/Diffexp/enrichment/KEGG_enrichment/enrichment_kegg.xls"
        #],diff_group=diff_groups.index)

    params:
        outdir="result/Diffexp/enrichment",
    benchmark:
        "benchmarks/diffexp_enrich.benchmark.txt"
    resources:
        qsub_mem=10,
        qsub_p=1
    log:  "logs/diffexp_enrich.log"
   # envmodules:
   #     config['report_params']['envmodules']
    shell:
        """
rm -r result/Diffexp/enrichment/ ;
mkdir -p result/Diffexp/enrichment/background_files/;
cp {reference}/annotation/gene_kegg.backgroud.xls result/Diffexp/enrichment/background_files/gene_kegg.background.xls;
cp {reference}/annotation/gene_go.backgroud.xls result/Diffexp/enrichment/background_files/gene_go.background.xls;
cp {reference}/annotation/gene_annotation.xls result/Diffexp/enrichment/background_files/gene_annotation.xls;
perl  scripts/enrichment/enrich_go_kegg.pl  \
    -infile result/Diffexp/*-vs-*-diff-pval-0.05-FC-1.5.xls \
    -go_bg {reference}/annotation/gene_go.backgroud.xls  \
    -category scripts/enrichment/category.xls \
    -kegg_bg {reference}/annotation/gene_kegg.backgroud.xls  \
    -outdir  {params.outdir} \
    -shelldir {params.outdir}/enrichment_sh \
    -thread 4  \
    -queue big |& tee {log}
        """
# =============================================================================================
with open("%s/config.yaml" % (reference) ) as f:
    reference_config=yaml.safe_load(f)
stringid = reference_config["stringid"]

#localrules: diffexp_ppi
rule diffexp_ppi:
    """
    add Annotation for diff expression genes
    """
    input:
        diff_DEG="result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5.xls"
    output:
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.new_colors.pdf",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.new_colors.png",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.new_colors.svg",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.pdf",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.png",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.svg",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction.tsv",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.mapping_results.tsv",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.mapping_results-top25.tsv",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.string_protein-protein-interaction-top25.tsv",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.ppi_network.tsv",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.top_25_ppi_network.pdf",
        "result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.top_25_ppi_network.png"
        
    params:
        outdir="result/Diffexp/ppi",
        prefix="{diff_group}-diff-pval-0.05-FC-1.5",
        network_tsv="result/Diffexp/ppi/{diff_group}-diff-pval-0.05-FC-1.5.ppi_network.tsv"
        #stringid=config['params']["stringid"]
    log:   "logs/diffexp_ppi_{diff_group}.log"
    benchmark:
        "benchmarks/diffexp_ppi_{diff_group}.benchmark.txt"
    resources:
        qsub_mem=10,
        qsub_p=1
    envmodules:
        config['report_params']['envmodules']
    shell:
        """
python scripts/ppi/ppi.py \
   -i {input.diff_DEG}  \
   -p {params.prefix} \
   -o {params.outdir} \
   -s {stringid}  |& tee {log}
Rscript  scripts/ppi/ppi_circle.r  \
   -i {params.network_tsv} \
   -d {input.diff_DEG} \
   -o {params.outdir}
        """
#localrules: diffexp_fc12
rule diffexp_fc12:
    """
    check if need to run FC1.2.
    """
    input:
        diff_DEG=expand("result/Diffexp/{diff_group}-diff-pval-0.05-FC-1.5_anno.xls",diff_group=diff_groups.index),
        gene_annotation= os.path.join(reference,"annotation/gene_annotation.xls"),
        cluster_rds = "result/Count_QC/filtered.h5seurat"
    output:
        "logs/diffexp_fc12.check"
    params:
        diff_stat = "result/Diffexp/diffexp_results_stat.xls",
    benchmark:
        "benchmarks/diffexp_fc12.benchmark.txt"
    resources:
        qsub_mem=50,
        qsub_p=10
    run:
        module = config['report_params']['envmodules']
        import pandas as pd
        from numpy import where
        from math import log
        diff_stat = pd.read_table( params.diff_stat , sep="\t" )
        if any(diff_stat["Total_diff(pvalue<0.05&FoldChange>1.5)"] < 20): 
            shell("mkdir -p result/Diffexp/FC1.2/ ")
            #diff_stat.index = [re.sub(r'(^.*\(|\)$)','',x)+"_"+re.sub(r'\(.*\)','',x)+"-vs-"+re.sub(r'\(.*\)','',y) for x in diff_stat.Case for y in diff_stat.Control]
            diff_stat.index = [re.sub(r'(^.*\(|\)$)','',getattr(row,"Case"))+"_"+re.sub(r'\(.*\)','',getattr(row,"Case"))+"-vs-"+re.sub(r'\(.*\)','',getattr(row,"Control")) for row in diff_stat.itertuples()]
            diff_stat = diff_stat.rename(columns={'Total_diff(pvalue<0.05&FoldChange>1.5)':'Total_diff(pvalue<0.05&FoldChange>1.2)'})
            stat_file = "result/Diffexp/FC1.2/diffexp_results_stat.xls"
            for diff_group in diff_groups.index :
                print("\033[32m"+diff_group)
                # params
                sig_file = f"result/Diffexp/FC1.2/{diff_group}-diff-pval-0.05-FC-1.2_anno.xls"
                sig_network_file = f"result/Diffexp/FC1.2/{diff_group}-diff-pval-0.05-FC-1.2.ppi_network.tsv"
                prefix=f"{diff_group}-diff-pval-0.05-FC-1.2"
                # ========== step1. write new file ========== 
                all= pd.read_table(f"result/Diffexp/{diff_group}-all_diffexp_genes_anno.xls")
                sig = all[(all["p-value"] <0.05 )& ((all["FoldChange"]> 1.2  ) | (all["FoldChange"]< 1/1.2 ))]
                sig["Regulation"]= where(sig["FoldChange"] > 1 , "Up","Down")
                sig.to_csv(sig_file,sep='\t',index=False)
                # diffexp_results_stat.xls calculation
                diff_stat["Up_diff"][diff_group] = len(where(sig["Regulation"] == "Up")[0])
                diff_stat["Down_diff"][diff_group] = len(where(sig["Regulation"] == "Down")[0])
                diff_stat["Total_diff(pvalue<0.05&FoldChange>1.2)"][diff_group] = len(where(sig["Regulation"] == "Up")[0])+len(where(sig["Regulation"] == "Down")[0])
                # ========== step2. heatmap ========== 
                case_name = re.sub(r'\(.*\)','',diff_stat.loc[diff_group,"Case"]) 
                ctrl_name = re.sub(r'\(.*\)','',diff_stat.loc[diff_group,"Control"])
                type_name = re.sub(r'(^.*\(|\)$)','',diff_stat.loc[diff_group,"Case"])
                shell("""
module load {module} &&
Rscript {script_dir}/scVis \
  -i {input.cluster_rds} \
  -f h5seurat \
  -o result/Diffexp/FC1.2/ \
  -t 10 \
  --assay RNA \
  --slot data,scale.data \
  --predicate "{type_name} %in% c(\\"{case_name}\\",\\"{ctrl_name}\\")" \
  diff_heatmap \
  -d {sig_file} \
  -n 20 \
  -g {type_name} \
  --group_colors customecol2 \
  --sample_ratio 0.8
  """)
                # ========== step3. ppi ========== 
                shell("module load {module} && python scripts/ppi/ppi.py -i {sig_file} -p {prefix} -o result/Diffexp/FC1.2/ppi -s {stringid} && Rscript  scripts/ppi/ppi_circle.r     -i {sig_network_file}    -d {sig_file}    -o result/Diffexp/FC1.2/ppi") 
            # write diffexp_results_stat.xls
            diff_stat.to_csv(stat_file,sep='\t',index=False)
            # ========== step4. enrich ========== 
            shell("""
module load {module} && perl scripts/enrichment/enrich_go_kegg.pl -infile result/Diffexp/FC1.2/*-vs-*-diff-*.xls \
-go_bg {reference}/annotation/gene_go.backgroud.xls  -kegg_bg {reference}/annotation/gene_kegg.backgroud.xls  \
-category /home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/enrich_background/category.xls  \
-outdir result/Diffexp/FC1.2/enrichment -shelldir result/Diffexp/FC1.2/enrichment/enrichment_sh    \
-thread 2 -queue big """)
            shell ("echo See result in result/Diffexp/FC1.2 > {output}")
        else : shell("echo All diff_group have number of DEGs greater than 20.  > {output}")

