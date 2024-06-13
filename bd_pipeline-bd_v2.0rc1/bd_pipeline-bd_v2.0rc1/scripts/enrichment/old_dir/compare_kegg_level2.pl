#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin $Script);

#######################################################################################
my $env="module load OESingleCell/2.0.0 && module load gcc/6.4.0";
#my $env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
#######################################################################################

my ($infile, $outdir, $mark, $groupname, $help);
GetOptions(
	"i:s"       =>\$infile,
	"o:s"       =>\$outdir,
	"m:s"       =>\$mark,
	"n:s"       =>\$groupname,
	"h|help!"   =>\$help,
);
my $usage=<< "USAGE";
Program: $0
Description: KEGG level2 distribution compare
Options:
	-i    <infile>    The input file by comma-separated      [Required]
	-o    <outdir>    The output directory of result         [Required]
	-m    <prefix>    The infile prefix by comma-separated   [Required]
	-n    <group>     The group name, eg:Group2-vs-Group1    [Required]
	-h|help           print help info
Example:
	perl compare_kegg_level2.pl -i Up/diff-KEGG_Classification.xls,Down/diff-KEGG_Classification.xls -o outdir/ -m Up,Down -n Group2-vs-Group1

USAGE

die $usage if(!$infile || !$outdir || !$mark || !$groupname || $help);
my @infile=split /,/,$infile;
for(my $i=0; $i<@infile; $i++){
	(-s $infile[$i]) || die "Error: don't find infile: $infile[$i] !\n";
	$infile[$i]=File::Spec->rel2abs($infile[$i]);
}
(-d $outdir) || mkdir $outdir;
$outdir=File::Spec->rel2abs($outdir);
my $mark_str=$mark; $mark_str=~s/,/_vs_/g;

my %hash; my @mark_arr=split /,/,$mark;
for(my $i=0; $i<@infile; $i++){
	open IN,"<$infile[$i]" || die $!;
	while(<IN>){
		chomp;
		next if(/^\s*$/ || /^#/ || /^Classification_level/);
		my @l=split /\t/;
		#$l[1]="$l[1]\t$mark_arr[$i]";
		${$hash{$l[1]}{$l[0]}}[$i]="$l[2]\t$l[3]\t$l[4]";
	}
	close IN;
}

open OUT,">$outdir/$mark_str\.KEGG_Classification.xls" || die $!;
print OUT "Classification_level2\tClassification_level1\tType\tGene_number\tPercentage\tGene\n";
for my $k(reverse(sort keys %hash)){
	for my $j(sort keys %{$hash{$k}}){
		for(my $i=0; $i<@infile; $i++){
			if(defined ${$hash{$k}{$j}}[$i] && ${$hash{$k}{$j}}[$i] ne ""){
				print OUT "$j\t$k\t$mark_arr[$i]\t${$hash{$k}{$j}}[$i]\n";
			}
		}
	}
}
close OUT;

system("$env && Rscript $Bin/diff-KEGG_Classification_all_DEG_up_down.r -i $outdir/$mark_str\.KEGG_Classification.xls -n $groupname -o $outdir/$mark_str\.KEGG_Classification");
