#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;
use FindBin qw($Bin);

########################### software/env/database #####################################
my $env="module load OESingleCell/3.0.d && module load gcc/6.4.0";
#my $env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
#######################################################################################

my ($diff_infile, $go_bg, $category, $kegg_bg, $outdir, $help);
GetOptions(
	"infile:s"     => \$diff_infile,
	"go_bg:s"      => \$go_bg,
	"category:s"   => \$category,
	"kegg_bg:s"    => \$kegg_bg,
	"outdir:s"     => \$outdir,
	"h|help!"      => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: difference gene GO and KEGG enrichment analysis
Options:
	-infile     <file>      The input file of *-vs-*-diff-*.xls            [Required]
	-go_bg      <file>      go backgroud file                              [Required]
	-category   <file>      The input category.xls                         [Required]
	-kegg_bg    <file>      kegg backgroud file                            [Required]
	-outdir     <dir>       The output directory of result                 [Required]
	-h|help                 print help info
Example:
	perl diff_enrichment.pl -infile Group3-vs-Group1-diff-padj-0.05-FC-2.gene.xls -go_bg go.backgroud.xls -category category.xls -kegg_bg kegg.backgroud.xls -outdir enrich

USAGE

die $usage if(!$diff_infile || !$outdir || !$go_bg || !$category || !$kegg_bg || $help);
(-s $diff_infile) || die "Error: don't find infile: $diff_infile !\n";
(-d $outdir) || mkdir $outdir;
(-s $go_bg) || die "Error: don't find go_bg: $go_bg !\n";
(-s $category) || die "Error: don't find category: $category !\n";
(-s $kegg_bg) || die "Error: don't find kegg_bg: $kegg_bg !\n";
$diff_infile=File::Spec->rel2abs($diff_infile);
$outdir=File::Spec->rel2abs($outdir);
$go_bg=File::Spec->rel2abs($go_bg);
$category=File::Spec->rel2abs($category);
$kegg_bg=File::Spec->rel2abs($kegg_bg);

my $name=basename($diff_infile);
$name=~s/\.xls$//g;
$name=(split /-diff-/,$name)[0];
system("mkdir -p $outdir/{GO_enrichment,KEGG_enrichment}/$name/");
#system("awk -F \"\\t\" '{if(NR==1 || \$9==\"Down\"){print \$0}}' $diff_infile > $outdir/$name-Down.xls");
#system("awk -F \"\\t\" '{if(NR==1 || \$9==\"Up\"){print \$0}}' $diff_infile > $outdir/$name-Up.xls");
system("awk -F \"\\t\" 'NR==1 {for (f=1; f <= NF; f+=1) {  if (\$f ~ /up_down/) flag=f } ;print \$0 } ;{if(\$flag ==\"Down\") print \$0} ' $diff_infile > $outdir/$name-Down.xls");
system("awk -F \"\\t\" 'NR==1 {for (f=1; f <= NF; f+=1) {  if (\$f ~ /up_down/) flag=f } ;print \$0 } ;{if(\$flag ==\"Up\") print \$0} ' $diff_infile > $outdir/$name-Up.xls");
system("cp $diff_infile $outdir/$name-Total.xls");
#system("$env && Rscript $Bin/enrichment.r -i $outdir/$name-Down.xls,$outdir/$name-Total.xls,$outdir/$name-Up.xls -j $go_bg -c $category -k $kegg_bg -o $outdir");
system("$env && Rscript $Bin/enrichment_clusterprofiler.r -i $outdir/$name-Down.xls,$outdir/$name-Total.xls,$outdir/$name-Up.xls -j $go_bg -c $category -k $kegg_bg  -p $name-Down,$name-Total,$name-Up  -o $outdir -n Gene -m 2 ");
print "$env && Rscript $Bin/enrichment_clusterprofiler.r -i $outdir/$name-Down.xls,$outdir/$name-Total.xls,$outdir/$name-Up.xls -j $go_bg -c $category -k $kegg_bg  -p $name-Down,$name-Total,$name-Up  -o $outdir -n Gene -m 2 ";
system("mv $outdir/enrichment-go-$name-Up.xls $outdir/GO_enrichment/$name/");
system("mv $outdir/enrichment-go-$name-Total.xls $outdir/GO_enrichment/$name/");
system("mv $outdir/enrichment-go-$name-Down.xls $outdir/GO_enrichment/$name/");
system("mv $outdir/enrichment-kegg-$name-Up.xls $outdir/KEGG_enrichment/$name/");
system("mv $outdir/enrichment-kegg-$name-Total.xls $outdir/KEGG_enrichment/$name/");
system("mv $outdir/enrichment-kegg-$name-Down.xls $outdir/KEGG_enrichment/$name/");
system("mv $outdir/$name-Up.xls $outdir/GO_enrichment/$name/diff-$name-Up.xls");
system("mv $outdir/$name-Total.xls $outdir/GO_enrichment/$name/diff-$name-Total.xls");
system("mv $outdir/$name-Down.xls $outdir/GO_enrichment/$name/diff-$name-Down.xls");

