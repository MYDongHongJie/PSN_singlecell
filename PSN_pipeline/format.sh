if [ -d "summary/04_Clustering" ];then
		#echo "1"
    sample=$(sed -n '2p' summary/sample_info.txt |cut -f 1)
    sample_num=$(sed '1d' summary/sample_info.txt |wc -l)
    cp -r /PERSONALBIO/work/singlecell/s00/software/3.StdPipe/Report .
    /PERSONALBIO/work/singlecell/s00/software/miniconda3/envs/stdpipe/bin/perl /PERSONALBIO/work/singlecell/s00/software/script/1.source/10x_webSummary.pl -i CellrangerOut/ -o summary/02_cellranger/
    cp -r  summary/02_cellranger/table*txt Report/table/
    #head -n 100  summary/04_Cluster/cells_GeneCounts.xls |cut -f 1-50 | sed '1s/^/Gene/'> Report/table/table4.txt
    sed '1s/^/Sample/' summary/03_QC_Filter/cells_filter_stat.xls >Report/table/cells_filter_stat.txt
    cp summary/03_QC_Filter/cells_qc_filter.png Report/pictures/seurat/
    cp summary/04_Clustering/1.preprocess/pca.ElbowPlot.png Report/pictures/seurat/pca.ElbowPlot.png

    #cp summary/04_Cluster/cluster_number.png Report/pictures/seurat/
		cp summary/04_Clustering/2.cluster_overview/cluster_summary.xls Report/table/cluster_summary.txt
    cp summary/04_Clustering/4.umap/cluster_umap.png Report/pictures/seurat/umap.png
    cp summary/04_Clustering/4.umap/cluster_umap.group.png Report/pictures/seurat/umap_group.png
		cp summary/04_Clustering/2.cluster_overview/cluster_group_cellcounts.png Report/pictures/seurat/cluster_group_cellcounts.png
		cp summary/04_Clustering/2.cluster_overview/cluster_sample_cellcounts.png Report/pictures/seurat/cluster_sample_cellcounts.png

    #cp summary/04_Cluster/cluster_tsne.label.png Report/pictures/seurat/tsne1.png
    #cp summary/04_Cluster/cluster_tsne_splitgroup.png Report/pictures/seurat/tsne2.png
		echo "marker基因"
		#mkdir Report/pictures/marker
    cp summary/05_Marker/marker_number.png Report/pictures/marker/marker_number.png
    head -n 20 summary/05_Marker/allmarkers.xls | sed '1s/^/Gene/'> Report/table/marker_anno.txt
		cp summary/05_Marker/cluster_top10_markers_heatmap.png Report/pictures/marker/cluster_top10_markers_heatmap.png
		cp summary/05_Marker/cluster_top5_dotplot.png Report/pictures/marker/cluster_top5_dotplot.png
		mkdir -p Report/pictures/marker/Each_celltype_marker
		cp summary/05_Marker/Each_celltype_marker/*/c*_top10_umap.png Report/pictures/marker/Each_celltype_marker
		cp summary/05_Marker/Each_celltype_marker/*/c*_top10_vlnplot.png Report/pictures/marker/Each_celltype_marker
		echo "细胞特征"
		cp summary/04_Clustering/1.preprocess/nCount_RNA_umap.png Report/pictures/seurat/nCount_RNA_umap.png
		cp summary/04_Clustering/1.preprocess/nFeature_RNA_umap.png Report/pictures/seurat/nFeature_RNA_umap.png
		cp summary/04_Clustering/1.preprocess/percent.mt_umap.png Report/pictures/seurat/percent.mt_umap.png
		echo "细胞周期"
		cp summary/04_Clustering/1.preprocess/cell_cycle_per_cluster_umap.png Report/pictures/seurat/cell_cycle_per_cluster_umap.png
		cp summary/04_Clustering/1.preprocess/cell_cycle_per_sample.png Report/pictures/seurat/cell_cycle_per_sample.png
		echo "高变基因"
		cp summary/04_Clustering/1.preprocess/VariableFeatures_distribution.png Report/pictures/seurat/VariableFeatures_distribution.png

		cp summary/06_Enrichment/cluster_0/up/GO_enrichment.xls Report/table/GO_enrichment.txt
		cp summary/06_Enrichment/cluster_0/up/KEGG_enrichment.xls Report/table/KEGG_enrichment.txt
		for i in `ls summary/05_Marker/Each_celltype_marker/` ;do mkdir -p Report/pictures/seurat/enrichment/$i;  cp summary/06_Enrichment/$i/up/*.png Report/pictures/seurat/enrichment/$i ;done

   mkdir log
   path=`pwd`
   report=$(basename "$PWD")
   contrast=$(basename $(dirname "$PWD"))
    mkdir -p ${contrast}/${report}/summary
    mv *.err log
    mv *.out log
    mv cp* log
    mv cr* log
    mv mv* log
    mv ag.slurm log
    mv merge* log
    bash /PERSONALBIO/work/singlecell/s02/software/script/std/rinfo.sh $path/Report/table/project.txt  ${contrast} ${report} $1 $2 $3
		rm $path/Report/table/oldproject.txt
    #echo $2
    cd $path/Report
    #if [ $2 == "1" ];then
       # python /PERSONALBIO/work/singlecell/s00/software/3.StdPipe/md2html.py single_sample.md 
				# sed "s/XXXXXXXX/$3/g" report.html
				# time=`date "+%Y-%m-%d"`
				# sed "s/YYYYYYYY/$time/g" report.html
    #else
        #echo "sample multi"
     python /PERSONALBIO/work/singlecell/s00/software/3.StdPipe/10XRNA/md2html.py report.md  $3
			#fi
    /PERSONALBIO/work/singlecell/s04//Test/donghongjie/Miniconda/envs/scvelo/bin/wkhtmltopdf --print-media-type --enable-javascript --javascript-delay 3000 report.html report.pdf
    cd $path
    if [ ! -s "report.html" ]; then
        echo "文件存在且大小为0KB。"
    else
        rm report.md
    fi
    if [ ! -f 1.raw_data/md5.txt ]; then
        md5sum 1.raw_data/*/* > 1.raw_data/md5.txt
    fi
    ln -s $path/1.raw_data $path/${contrast}/${report}
    ln -s $path/QC $path/${contrast}/${report}/summary
    ln -s $path/summary/0* $path/${contrast}/${report}/summary
    if [ "`ls -A $path/Report/pictures/seurat`" != "" ];then
        ln -s $path/Report $path/${contrast}/${report}/summary
        ln -s $path/summary/07_Rds $path/${contrast}/${report}/summary
    fi

else
		#echo "第一步"
    path=`pwd`
    report=$(basename "$PWD")
    contrast=$(basename $(dirname "$PWD"))
    mkdir -p ${contrast}/${report}/summary
    cd $path
    if [ ! -f 1.raw_data/md5.txt ]; then
        md5sum 1.raw_data/*/* > 1.raw_data/md5.txt
    fi
    for sample in `awk 'NR>1{print $1}' sample_info.txt`;do
        if [ ! -d CellrangerOut/$sample ]; then
            mkdir CellrangerOut
            mv $sample CellrangerOut
        fi

        if [ ! -d summary/02_cellranger/$sample ]; then
            mkdir -p summary/02_cellranger/$sample
        cp -r  \
               CellrangerOut/$sample/outs/cloupe.cloupe  CellrangerOut/$sample/outs/web_summary.html  \
               CellrangerOut/$sample/outs/filtered_feature_bc_matrix  \
               CellrangerOut/$sample/outs/filtered_feature_bc_matrix.h5 summary/02_cellranger/$sample
        fi
    done
    ln -s $path/1.raw_data $path/${contrast}/${report}
    ln -s $path/QC $path/${contrast}/${report}/summary
    ln -s $path/summary/02_cellranger $path/${contrast}/${report}/summary
fi
