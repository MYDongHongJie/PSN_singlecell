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
my $qsub_pbs="/data/software/qsub/qsub-sge.pl";
#######################################################################################
my (@diff_infile, $go_bg, $category, $kegg_bg, $anno_kegg, $html_png, $outdir, $work_sh, $help);
GetOptions(
	"infile:s{1,}"  => \@diff_infile,
	"go_bg:s"       => \$go_bg,
	"category:s"    => \$category,
	"kegg_bg:s"     => \$kegg_bg,
	"anno_kegg:s"   => \$anno_kegg,
	"outdir:s"      => \$outdir,
	"shelldir:s"    => \$work_sh,
	"html_png:s"    => \$html_png,
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
	-anno_kegg     <file>      anno kegg backgroud                                    [Required]
	-outdir        <dir>       The output directory of result                         [Required]
	-shelldir      <dir>       The output directory of shell scripts.[default: ./]    [Optional]
	-go_lv2        <file>      Unigene.GO.classification.stat.xls                     [Optional]
	-kegg_lv2      <file>      Unigene.KEGG.Classification.xls                        [Optional]
	                           [default: /public/land/database/kegg/pathway_ko]
	-thread        <num>       max thread. [default: 5]                               [Optional]
	-queue         <str>       queue:all,cu,big. [default: all]                       [Optional]
Example:
	1: perl 5.2.enrich_go_kegg.pl -infile *-vs-*-diff-*.xls -go_bg go.backgroud.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/ -go_lv2 Unigene.GO.classification.stat.xls -kegg_lv2 Unigene.KEGG.Classification.xls
	2: perl 5.2.enrich_go_kegg.pl -infile *-vs-*-diff-*.xls -go_bg go.backgroud.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/
	3: perl 5.2.enrich_go_kegg.pl -infile A-vs-B-diff-pval-0.05-FC-2.gene.xls A-vs-C-diff-pval-0.05-FC-2.gene.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/

USAGE

die $usage if(!@diff_infile || !$outdir || !$go_bg || !$kegg_bg  || $help);
$work_sh ||= getcwd;
$category ||= "$Bin/category.xls";

(-d $outdir) || mkdir $outdir;
(-d $work_sh) || mkdir $work_sh;
(-s $go_bg) || die "Error: don't find go_bg: $go_bg !\n";
(-s $category) || die "Error: don't find category: $category !\n";
(-s $kegg_bg) || die "Error: don't find kegg_bg: $kegg_bg !\n";
# (-s $anno_kegg) || die "Error: don't find anno_kegg: $anno_kegg !\n";
#(-d $html_png) || die "Error: don't find html_png: $html_png !\n";
$outdir=File::Spec->rel2abs($outdir);
#$work_sh=File::Spec->rel2abs($work_sh);
$go_bg=File::Spec->rel2abs($go_bg);
#$category=File::Spec->rel2abs($category);
$kegg_bg=File::Spec->rel2abs($kegg_bg);
$anno_kegg=File::Spec->rel2abs($anno_kegg);
# $html_png=File::Spec->rel2abs($html_png);
#if(defined $go_lv2){#
#	(-s $go_lv2) || die "Error: don't find go_lv2: $go_lv2 !\n";
#	$go_lv2=File::Spec->rel2abs($go_lv2);
#}
#if(defined $kegg_lv2){
#	(-s $kegg_lv2) || die "Error: don't find kegg_lv2: $kegg_lv2 !\n";
#	$kegg_lv2=File::Spec->rel2abs($kegg_lv2);
#}



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
	#$name=(split /-diff-/,$name)[0];
	$name=(split /-diff-/,$name)[0];
	push(@group_name, $name);
	print A "$env && perl $Bin/diff_enrichment.pl -infile $diff -go_bg $go_bg -category $category -kegg_bg $kegg_bg -outdir $outdir\n";
	#print A "sh $Bin/2.diff_enrichment.sh $diff $go_bg $kegg_bg $outdir\n";
}
close A;
system("sh $work_sh/a.enrichment.sh");

open B,">$work_sh/b.stati_enrichment.sh" || die $!;
print B "$env && Rscript $Bin/stati_enrichment.r -j $outdir/GO_enrichment -k $outdir/KEGG_enrichment\n";
close B;
system("sh $work_sh/b.stati_enrichment.sh");

################GO###############
open C1,">$work_sh/c.go_graph_top.sh" || die $!;
for my $i (@group_name){
	for my $j ("Total", "Up", "Down"){
	print C1 "$env && Rscript $Bin/top10X3_GO_gai.r -i $outdir/GO_enrichment/$i/enrichment-go-$i-$j.xls -m $j -o $outdir/GO_enrichment/$i/\n";
	# print C1 "$env && Rscript $Bin/top10X3_GO.r -i $outdir/GO_enrichment/$i/enrichment-go-$i-$j.xls -m $j -o $outdir/GO_enrichment/$i/\n";	
	}
}
close C1;
system("sh $work_sh/c.go_graph_top.sh");

# open C2,">$work_sh/c.go_graph_level2.sh" || die $!;
# for my $i (@group_name){
# 	for my $j ("Total", "Up", "Down"){
# 	print C2 "perl $Bin/diff_go_level2.pl -i $outdir/GO_enrichment/$i/enrichment-go-$i-$j.xls -o $outdir/GO_enrichment/$i/\n";
# 	}
# }
# close C2;
# system("sh $work_sh/c.go_graph_level2.sh");

# open C3,">$work_sh/c.go_compare_graph.sh" || die $!;
# for my $i (@group_name){
# 	print C3 "sh $Bin/2.compare_go_level2.sh $outdir/GO_enrichment/$i/GO.level2.stat.Up.xls $outdir/1.GO_enrichment/$i/GO.level2.stat.Down.xls $outdir/GO_enrichment/$i/GO.level2.stat.Total.xls $i $outdir/GO_enrichment/$i $go_bg";
# 	print C3 " && rm $outdir/1.GO_enrichment/$i/diff-$i-Up.xls $outdir/GO_enrichment/$i/diff-$i-Down.xls $outdir/GO_enrichment/$i/diff-$i-Total.xls\n";
# }
# close C3;
# system("sh $work_sh/c.go_compare_graph.sh");

open C4,">$work_sh/c.go_chord.sh" || die $!;
for my $diff (@diff_infile){
	my $name=basename($diff);
	$name=~s/\.xls$//g;
	# $name=(split /-diff-/,$name)[0];
	$name=(split /-diff-/,$name)[0];
	for my $j ("Total", "Up", "Down"){
	print C4 "$env && Rscript $Bin/enrich_chord.r -i $outdir/GO_enrichment/$name/GO.top.$j.xls -d $diff -s p-value -o $outdir/GO_enrichment/$name\n";
	}
}
close C4;
system("sh $work_sh/c.go_chord.sh");

open C5,">$work_sh/c.go_circbar.sh" || die $!;
for my $diff (@diff_infile){
	my $name=basename($diff);
	$name=~s/\.xls$//g;
	#$name=(split /-diff-/,$name)[0];
	$name=(split /-diff-/,$name)[0];
	print C5 "$env && Rscript $Bin/enrich_circbarplot.r -i $diff -e $outdir/GO_enrichment/$name/enrichment-go-$name-Total.xls -s p-value -o $outdir/GO_enrichment/$name\n";
}
close C5;
system("sh $work_sh/c.go_circbar.sh");

################kegg###############
open D1,">$work_sh/d.kegg_graph_top.sh" || die $!;
for my $i (@group_name){
	for my $j ("Total", "Up", "Down"){
	print D1 "$env && Rscript $Bin/top20_KEGG_hypelink.r -i $outdir/KEGG_enrichment/$i/enrichment-kegg-$i-$j.xls -m $j -o $outdir/KEGG_enrichment/$i/\n";
	}
}
close D1;
system("sh $work_sh/d.kegg_graph_top.sh");

# open D2,">$work_sh/d.kegg_graph_level2.sh" || die $!;
# for my $i (@group_name){
# 	for my $j ("Total", "Up", "Down"){
# 	print D2 "perl $Bin/diff_kegg_level2.pl -i $outdir/KEGG_enrichment/$i/enrichment-kegg-$i-$j.xls -o $outdir/KEGG_enrichment/$i/\n";
# 	}
# }
# close D2;
# system("sh $work_sh/d.kegg_graph_level2.sh");

# open D3,">$work_sh/d.kegg_compare_graph.sh" || die $!;
# for my $i (@group_name){
# 	print D3 "sh $Bin/1.compare_kegg_level2.sh $outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Up.xls $outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Down.xls $outdir/KEGG_enrichment/$i/diff-KEGG-Classification-Total.xls $i $outdir/KEGG_enrichment/$i $kegg_bg\n";
# }
# close D3;
# system("sh $work_sh/d.kegg_compare_graph.sh");

open D4,">$work_sh/d.kegg_chord.sh" || die $!;
for my $diff (@diff_infile){
	my $name=basename($diff);
	$name=~s/\.xls$//g;
	#$name=(split /-diff-/,$name)[0];
	$name=(split /-diff-/,$name)[0];
	for my $j ("Total", "Up", "Down"){
	print D4 "$env && Rscript $Bin/enrich_chord.r -i $outdir/KEGG_enrichment/$name/KEGG.top.$j.xls -d $diff -o $outdir/KEGG_enrichment/$name\n";
	}
}
close D4;
system("sh $work_sh/d.kegg_chord.sh");

open D5,">$work_sh/d.kegg_circbar.sh" || die $!;
for my $diff (@diff_infile){
	my $name=basename($diff);
	$name=~s/\.xls$//g;
	#$name=(split /-diff-/,$name)[0];
	$name=(split /-diff-/,$name)[0];
	print D5 "$env && Rscript $Bin/enrich_circbarplot.r -i $diff -e $outdir/KEGG_enrichment/$name/enrichment-kegg-$name-Total.xls -o $outdir/KEGG_enrichment/$name\n";
}
close D5;
system("sh $work_sh/d.kegg_circbar.sh");


#open E,">$work_sh/e.kegg_map.sh" || die $!;
#print E "[ -d $outdir/3.KEGG_map ] && rm -rf $outdir/3.KEGG_map\n";
#print E "mkdir -p $outdir/3.KEGG_map && cp $anno_kegg $outdir/3.KEGG_map && cd $outdir/3.KEGG_map && ln -s $html_png html_png && cp $diff_infile . && sh $Bin/path_mapper_a2.1.sh\n";
#close E;
#system("sh $work_sh/e.kegg_map.sh");


