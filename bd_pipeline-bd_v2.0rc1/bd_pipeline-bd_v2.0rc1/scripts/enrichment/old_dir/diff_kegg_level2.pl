#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);

########################### software/env/database #####################################
my $pathway_type="/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/KEGGpathway_three_levels_v2.xls";
my $env="module load OESingleCell/3.0.d && module load gcc/6.4.0";
#my $env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.2.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
#######################################################################################

my ($enrich_file, $outdir, $help);
GetOptions(
	"i:s"     => \$enrich_file,
	"o:s"     => \$outdir,
	"h|help!" => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: difference genes kegg level2 distribution
Options:
	-i     <infile>    The infile name(enrichment-kegg-*.xls)     [Required]
	-o     <outdir>    The output directory of result             [Required]
	-h|help            print help info
Example:
	perl diff_kegg_level2.pl -i KEGG_enrichment/Group2-vs-Group1/Total/enrichment-kegg-Group2-vs-Group1-Total.xls -o KEGG_enrichment/Group2-vs-Group1/Total/

USAGE

die $usage if(!$enrich_file || !$outdir || $help);
(-s $enrich_file) || die "Error: don't find infile: $enrich_file !\n";
(-d $outdir) || mkdir $outdir;
$enrich_file=File::Spec->rel2abs($enrich_file);
$outdir=File::Spec->rel2abs($outdir); $outdir=~s/\/$//g;
my $groupname=basename($enrich_file);
my $mark=$1 if($groupname=~/-(Total|Down|Up)/);
$groupname=~s/^enrichment-kegg-//;
$groupname=~s/\.(xls|txt)$//;
$groupname=~s/-(Total|Down|Up)/\\\($1\\\)/;

(-s $pathway_type) || die "Error: don't open file: $pathway_type !\n";
open PATH,"<$pathway_type" || die $!;
my %path_class;
while(<PATH>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$path_class{$l[0]}="$l[1]--$l[2]";
}
close PATH;

open EN,"<$enrich_file" || die $!;
my %class_gene;my %totalgene;
while(<EN>){
	chomp;
	next if(/^#/ || /^\s*$/ || /id\tterm/);
	my @t=split /\t/;
	my @gene=split /;/,$t[-1];
	for my $g (@gene){
		$t[0]=~s/[A-Za-z:]+(\d+)/ko$1/;
		$class_gene{$path_class{$t[0]}}{$g}=1;
		$totalgene{$g}=1;
	}
}
close EN;

my $total_genenum=keys %totalgene;
open OUT,">$outdir/diff-KEGG-Classification-$mark.xls" || die $!;
print OUT "Classification_level2\tClassification_level1\tgene_number\tpercentage\tGenes\n";
for my $c (sort keys %class_gene){
	my ($level1, $level2)=split /--/,$c;
	my @gene=sort keys %{$class_gene{$c}};
	my $gene_num=@gene;
	my $per=sprintf("%.2f",$gene_num*100/$total_genenum);
	print OUT "$level2\t$level1\t$gene_num\t$per\t".join("; ",@gene)."\n";
}
close OUT;

system("$env && Rscript $Bin/diff_KEGG_Classification.r -i $outdir/diff-KEGG-Classification-$mark.xls -o $outdir/diff-KEGG-Classification-$mark -m $groupname");
