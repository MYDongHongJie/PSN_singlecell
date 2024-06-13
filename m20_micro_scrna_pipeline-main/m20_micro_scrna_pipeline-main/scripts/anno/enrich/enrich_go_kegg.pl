#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Cwd;
use FindBin qw($Bin);

########################### software/env/database #####################################
my $env="module load OESingleCell/2.0.0";
#my $env="export PATH=/data/software/gcc/gcc-v6.4.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
my $qsub_pbs="/data/software/qsub/qsub-sge.pl";
#######################################################################################
my (@diff_infile, $go_bg, $category, $kegg_bg, $outdir, $work_sh, $queue, $thread, $help);
GetOptions(
	"infile:s{1,}"  => \@diff_infile,
	"go_bg:s"       => \$go_bg,
	"category:s"    => \$category,
	"kegg_bg:s"     => \$kegg_bg,
	"outdir:s"      => \$outdir,
	"shelldir:s"    => \$work_sh,
	"queue:s"       => \$queue,
	"thread:i"      => \$thread,
	"h|help!"       => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: difference gene GO and KEGG enrichment analysis
Options:
	-infile        <file>      The inputfile:eg markergene/*.xls      [Required]
	-go_bg         <file>      go backgroud file                                      [Required]
	-category      <file>      The input category.xls                                 [Required]
	-kegg_bg       <file>      kegg backgroud file                                    [Required]
	-outdir        <dir>       The output directory of result                         [Required]
	-shelldir      <dir>       The output directory of shell scripts.[default: ./]    [Optional]
	-thread        <num>       max thread. [default: 5]                               [Optional]
	-queue         <str>       queue:all,cu,big. [default: all]                       [Optional]
Example:
	1: perl 5.2.enrich_go_kegg.pl -infile *.xls -go_bg go.backgroud.xls -category category.xls -kegg_bg kegg.backgroud.xls -anno_kegg anno-kegg.backgroud.xls -outdir enrich/

USAGE

die $usage if(!@diff_infile || !$outdir || !$go_bg || !$category || !$kegg_bg || $help);
$queue ||= "all";
$thread ||= 5;
$work_sh ||= getcwd;

(-d $outdir) || mkdir $outdir;
(-d $work_sh) || mkdir $work_sh;
(-s $go_bg) || die "Error: don't find go_bg: $go_bg !\n";
(-s $category) || die "Error: don't find category: $category !\n";
(-s $kegg_bg) || die "Error: don't find kegg_bg: $kegg_bg !\n";

$outdir=File::Spec->rel2abs($outdir);
$work_sh=File::Spec->rel2abs($work_sh);
$go_bg=File::Spec->rel2abs($go_bg);
$category=File::Spec->rel2abs($category);
$kegg_bg=File::Spec->rel2abs($kegg_bg);

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
	push(@group_name, $name);
	print A "perl $Bin/diff_enrichment.pl -infile $diff -go_bg $go_bg -category $category -kegg_bg $kegg_bg -outdir $outdir\n";
}
close A;
system("perl $qsub_pbs --queue $queue --lines 1  --num_proc ${thread} --interval 50  $work_sh/a.enrichment.sh");

open B,">$work_sh/b.stati_enrichment.sh" || die $!;
print B "$env && Rscript $Bin/stati_enrichment.r -j $outdir/GO_enrichment -k $outdir/KEGG_enrichment\n";
close B;
system("perl $qsub_pbs --queue $queue --lines 1  --num_proc 1 --interval 50  $work_sh/b.stati_enrichment.sh");

open C,">$work_sh/c.kegg_go_graph.sh" || die $!;
for my $i (@group_name){
	for my $j ("Total"){
		print C "$env && Rscript $Bin/top20_KEGG.r -i $outdir/KEGG_enrichment/$i/enrichment-kegg-$i-$j.xls -m $j -o $outdir/KEGG_enrichment/$i/";
		print C " && Rscript $Bin/top10X3_GO.r -i $outdir/GO_enrichment/$i/enrichment-go-$i-$j.xls -m $j -o $outdir/GO_enrichment/$i/\n";
	}
}
close C;
system("perl $qsub_pbs --queue $queue --lines 1  --num_proc ${thread}  --interval 50 $work_sh/c.kegg_go_graph.sh");
