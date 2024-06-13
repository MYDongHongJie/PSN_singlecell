#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin $Script);

#######################################################################################
my $env="module load OESingleCell/3.0.d  && module load gcc/6.4.0";
#my $env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
my $level2="$Bin/endNode3.csv";
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
Description: GO level2 distribution compare
Options:
	-i    <infile>    The input file by comma-separated      [Required]
	-o    <outdir>    The output directory of result         [Required]
	-m    <prefix>    The infile prefix by comma-separated   [Required]
	-n    <group>     The group name, eg:Group2-vs-Group1    [Required]
	-h|help           print help info
Example:
	perl compare_go_level2.pl -i Up/GO.level2.stat.xls,Down/GO.level2.stat.xls -o outdir/ -m Up,Down -n Group2-vs-Group1

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

my %level2_64;
(-s $level2) || die "Error: don't open file: $level2 !\n";
open LE,"<$level2" || die $!;
while(<LE>){
	chomp;
	my @l=split /,/;
	$level2_64{"$l[2]\t$l[1]"}=1;
}
close LE;
my @level2_64=sort keys %level2_64;
unshift(@level2_64,"#Total_gene\t");
unshift(@level2_64,"#GO_classify1\tGO_classify2");

my %hash;my @Total_gene; my %gene;
for(my $i=0; $i<@infile; $i++){
	open IN,"<$infile[$i]" || die $!;
	while(<IN>){
		chomp;
		next if(/^\s*$/);
		my @l=split /\t/,$_,4;
		if(/^#GO_classify1/){
			my $key="$l[0]\t$l[1]";
			${$hash{$key}}[$i]=$l[2];
			${$gene{$key}}[$i]="$l[3]($l[2])";
		}elsif(/^#Total_gene/){
			$Total_gene[$i]=$l[2];
			my $key="$l[0]\t$l[1]";
			${$hash{$key}}[$i]=$l[2];
			${$gene{$key}}[$i]="";
		}else{
			my $key="$l[0]\t$l[1]";
			${$hash{$key}}[$i]=$l[2]*100/$Total_gene[$i];
			${$gene{$key}}[$i]=$l[3];
		}
	}
	close IN;
}

open TMP,">$outdir/$mark_str\.GO.level2.stat.xls_tmp" || die $!;
open OUT,">$outdir/$mark_str\.GO.level2.stat.xls" || die $!;
#for my $k(sort keys %hash){
for my $k(@level2_64){
	print TMP "$k";
	print OUT "$k";
	for(my $i=0; $i<@infile; $i++){
		if(defined ${$hash{$k}}[$i] && ${$hash{$k}}[$i] ne ""){
			print TMP "\t${$hash{$k}}[$i]";
			print OUT "\t${$hash{$k}}[$i]";
		}else{
			print TMP "\t0";
			print OUT "\t0";
		}
	}
	for(my $i=0; $i<@infile; $i++){
		if(defined ${$hash{$k}}[$i] && ${$hash{$k}}[$i] ne ""){
			print OUT "\t${$gene{$k}}[$i]";
		}else{
			print OUT "\t--";
		}
	}
	print TMP "\n";
	print OUT "\n";
}
close TMP;
close OUT;

system("$env && Rscript $Bin/GO.level2_plot_compare.r -i $outdir/$mark_str\.GO.level2.stat.xls_tmp -o $outdir -m $mark -n $groupname -p $mark_str\.GO.level2.stat && rm $outdir/$mark_str\.GO.level2.stat.xls_tmp");
