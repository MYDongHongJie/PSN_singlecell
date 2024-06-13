#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Cwd;
use FindBin qw($Bin);

########################### software/env/database #####################################
my $env="module load OESingleCell/3.0.d  && module load gcc/6.4.0";
#my $env="export PATH=/data/software/gcc/gcc-v6.4.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
my $qsub_pbs="/data/software/qsub/qsub-sge.pl";
#######################################################################################
my (@diff_infile, $go_bg, $category, $kegg_bg, $anno_kegg, $go_lv2, $kegg_lv2, $html_png, $outdir, $work_sh, $queue, $thread, $help);
GetOptions(
	"infile:s{1,}"  => \@diff_infile,
	"go_bg:s"       => \$go_bg,
	"category:s"    => \$category,
	"kegg_bg:s"     => \$kegg_bg,
	"anno_kegg:s"   => \$anno_kegg,
	"outdir:s"      => \$outdir,
	"shelldir:s"    => \$work_sh,
	"go_lv2:s"      => \$go_lv2,
	"kegg_lv2:s"    => \$kegg_lv2,
	"html_png:s"    => \$html_png,
	"queue:s"       => \$queue,
	"thread:i"      => \$thread,
	"h|help!"       => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: difference gene GO and KEGG enrichment analysis
Options:
	-infile        <file>      The inputfile:eg Quantification/*-vs-*-diff-*.xls      [Required]
	-go_bg         <file>      go backgroud file                                      [Required]
	-category      <file>      The input category.xls                                 [Required]
	-kegg_bg       <file>      kegg backgroud file                                    [Required]
	-anno_kegg     <file>      anno kegg backgroud                                    [Optional]
	-outdir        <dir>       The output directory of result                         [Required]
	-shelldir      <dir>       The output directory of shell scripts.[default: ./]    [Optional]
	-go_lv2        <file>      Unigene.GO.classification.stat.xls                     [Optional]
	-kegg_lv2      <file>      Unigene.KEGG.Classification.xls                        [Optional]
	-html_png      <dir>       html_png directory                                     [Optional]
	                           [default: /public/land/database/kegg/pathway_ko]
	-thread        <num>       max thread. [default: 5]                               [Optional]
	-queue         <str>       queue:all,cu,big. [default: all]                       [Optional]
Example:
	1: perl 5.2.enrich_go_kegg.pl -infile *-vs-*-diff-*.xls -go_bg go.backgroud.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/ -go_lv2 Unigene.GO.classification.stat.xls -kegg_lv2 Unigene.KEGG.Classification.xls
	2: perl 5.2.enrich_go_kegg.pl -infile *-vs-*-diff-*.xls -go_bg go.backgroud.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/
	3: perl 5.2.enrich_go_kegg.pl -infile A-vs-B-diff-pval-0.05-FC-2.gene.xls A-vs-C-diff-pval-0.05-FC-2.gene.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/

USAGE

die $usage if(!@diff_infile || !$outdir || !$go_bg || !$category || !$kegg_bg || $help);
#die $usage if(!@diff_infile || !$outdir || !$go_bg || !$category || !$kegg_bg || !$anno_kegg || $help);
#$html_png ||="/public/land/database/kegg/pathway_ko";
$queue ||= "all";
$thread ||= 5;
$work_sh ||= getcwd;

(-d $outdir) || mkdir $outdir;
(-d $work_sh) || mkdir $work_sh;
(-s $go_bg) || die "Error: don't find go_bg: $go_bg !\n";
(-s $category) || die "Error: don't find category: $category !\n";
(-s $kegg_bg) || die "Error: don't find kegg_bg: $kegg_bg !\n";
# (-s $anno_kegg) || die "Error: don't find anno_kegg: $anno_kegg !\n";
#(-d $html_png) || die "Error: don't find html_png: $html_png !\n";
$outdir=File::Spec->rel2abs($outdir);
$work_sh=File::Spec->rel2abs($work_sh);
$go_bg=File::Spec->rel2abs($go_bg);
$category=File::Spec->rel2abs($category);
$kegg_bg=File::Spec->rel2abs($kegg_bg);
$anno_kegg=File::Spec->rel2abs($anno_kegg);
#$html_png=File::Spec->rel2abs($html_png);
if(defined $go_lv2){
	(-s $go_lv2) || die "Error: don't find go_lv2: $go_lv2 !\n";
	$go_lv2=File::Spec->rel2abs($go_lv2);
}
if(defined $kegg_lv2){
	(-s $kegg_lv2) || die "Error: don't find kegg_lv2: $kegg_lv2 !\n";
	$kegg_lv2=File::Spec->rel2abs($kegg_lv2);
}

(@diff_infile==0) && die "Error: don't find infile !\n";
for(my $i=0;$i<=$#diff_infile;$i++){
	(-s $diff_infile[$i]) || die "Error: don't find file $diff_infile[$i] !\n";
	$diff_infile[$i]=File::Spec->rel2abs($diff_infile[$i]);
}
my $diff_infile=join(" ",@diff_infile);

my @group_name;
open A,">$work_sh/a.enrichment.sh" || die $!;
for my $diff (@diff_infile){
	my $name=basename($diff);
	$name=~s/\.xls$//g;
	$name=(split /-diff-/,$name)[0];
	push(@group_name, $name);
	print A "perl $Bin/diff_enrichment.pl -infile $diff -go_bg $go_bg -category $category -kegg_bg $kegg_bg -outdir $outdir\n";
}
close A;
system("bash  $work_sh/a.enrichment.sh");

open B,">$work_sh/b.stati_enrichment.sh" || die $!;
print B "$env && Rscript $Bin/stati_enrichment.r -j $outdir/GO_enrichment -k $outdir/KEGG_enrichment\n";
close B;
system("bash $work_sh/b.stati_enrichment.sh");

open C,">$work_sh/c.kegg_go_graph.sh" || die $!;
for my $i (@group_name){
	for my $j ("Total", "Up", "Down"){
		print C "$env && Rscript $Bin/top20_KEGG.r -i $outdir/KEGG_enrichment/$i/enrichment-kegg-$i-$j.xls -m $j -o $outdir/KEGG_enrichment/$i/";
		print C " && perl $Bin/diff_kegg_level2.pl -i $outdir/KEGG_enrichment/$i/enrichment-kegg-$i-$j.xls -o $outdir/KEGG_enrichment/$i/";
		print C " ; Rscript $Bin/top10X3_GO.r -i $outdir/GO_enrichment/$i/enrichment-go-$i-$j.xls -m $j -o $outdir/GO_enrichment/$i/";
		# print C " && Rscript $Bin/topGO.r -d $outdir/GO_enrichment/$i/diff-$i-$j.xls -j $go_bg -m $j -o $outdir/GO_enrichment/$i/ && rm $outdir/GO_enrichment/$i/diff-$i-$j.xls";
		print C " && perl $Bin/diff_go_level2.pl -i $outdir/GO_enrichment/$i/enrichment-go-$i-$j.xls -o $outdir/GO_enrichment/$i/\n";
	}
}
close C;
system("bash $work_sh/c.kegg_go_graph.sh");

open D,">$work_sh/d.kegg_go_compare_graph.sh" || die $!;
for my $i (@group_name){
	print D "perl $Bin/compare_kegg_level2.pl -i $outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Up.xls,$outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Down.xls -o $outdir/KEGG_enrichment/$i -m Up,Down -n $i";
	if(defined $kegg_lv2){
		print D " ; perl $Bin/compare_kegg_level2.pl -i $kegg_lv2,$outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Total.xls -o $outdir/KEGG_enrichment/$i -m ALL,DEG -n $i";
	}else{
		print D " ; perl $Bin/format_keggbackgroud.pl -i $kegg_bg -o $outdir/KEGG_enrichment/$i/Unigene.KEGG.Classification.xls";
		print D " && perl $Bin/compare_kegg_level2.pl -i $outdir/KEGG_enrichment/$i/Unigene.KEGG.Classification.xls,$outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Total.xls -o $outdir/KEGG_enrichment/$i -m ALL,DEG -n $i";
		print D " && rm -fr $outdir/KEGG_enrichment/$i/Unigene.KEGG.Classification.xls";
	}

	print D " ; perl $Bin/compare_go_level2.pl -i $outdir/GO_enrichment/$i/GO.level2.stat.Up.xls,$outdir/GO_enrichment/$i/GO.level2.stat.Down.xls -o $outdir/GO_enrichment/$i -m Up,Down -n $i";
	if(defined $go_lv2){
		print D " ; perl $Bin/compare_go_level2.pl -i $go_lv2,$outdir/GO_enrichment/$i/GO.level2.stat.Total.xls -o $outdir/GO_enrichment/$i -m ALL,DEG -n $i";
	}else{
		print D " ; perl $Bin/format_gobackgroud.pl -i $go_bg -s Unigene -c $category -o $outdir/GO_enrichment/$i/";
		print D " && perl $Bin/compare_go_level2.pl -i $outdir/GO_enrichment/$i/Unigene.GO.classification.stat.xls,$outdir/GO_enrichment/$i/GO.level2.stat.Total.xls -o $outdir/GO_enrichment/$i -m ALL,DEG -n $i";
		print D " && rm -fr $outdir/GO_enrichment/$i/Unigene.GO.classification.stat.xls";
	}
	print D "\n";
}
close D;
system("bash  $work_sh/d.kegg_go_compare_graph.sh");

#open E,">$work_sh/e.kegg_map.sh" || die $!;
#print E "[ -d $outdir/KEGG_map ] && rm -rf $outdir/KEGG_map\n";
#print E "mkdir -p $outdir/KEGG_map && cp $anno_kegg $outdir/KEGG_map && cd $outdir/KEGG_map && ln -s $html_png html_png && cp $diff_infile . && sh $Bin/path_mapper_a1.0.sh\n";
#close E;
#system("perl $qsub_pbs --queue $queue --lines 1 --num_proc 4 --interval 50 $work_sh/e.kegg_map.sh");
