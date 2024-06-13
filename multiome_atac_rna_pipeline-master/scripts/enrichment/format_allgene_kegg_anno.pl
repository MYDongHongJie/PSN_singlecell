#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use FindBin qw($Bin);

my ($infile, $outfile, $help);
GetOptions(
	"i:s"     => \$infile,
	"o:s"     => \$outfile,
	"h|help!" => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: difference genes kegg level2 distribution
Options:
	-i     <infile>    The infile name(prefix.KEGG.pathway.xls)   [Required]
	-o     <outfile>   The output file                            [Required]
	-h|help            print help info
Example:
	perl diff_kegg_level2.pl -i Unigene.KEGG.pathway.xls -o Unigene.KEGG_Classification.xls

USAGE

die $usage if(!$infile || !$outfile || $help);
(-s $infile) || die "Error: don't find infile: $infile !\n";
$infile=File::Spec->rel2abs($infile);
$outfile=File::Spec->rel2abs($outfile); 

open EN,"<$infile" || die $!;
my %class_gene;my %totalgene;
while(<EN>){
	chomp;
	next if(/^#/ || /^\s*$/ || /id\tterm/);
	my @t=split /\t/;
	my @gene=split /;/,$t[4];
	for my $g (@gene){
		$class_gene{$t[0]}{$g}=1;
		$totalgene{$g}=1;
	}
}
close EN;

my $total_genenum=keys %totalgene;
open OUT,">$outfile" || die $!;
print OUT "Classification_level2\tClassification_level1\tgene_number\tpercentage\tGenes\n";
for my $c (sort keys %class_gene){
	my ($level1, $level2)=split /--/,$c;
	my @gene=sort keys %{$class_gene{$c}};
	my $gene_num=@gene;
	my $per=sprintf("%.2f",$gene_num*100/$total_genenum);
	print OUT "$level2\t$level1\t$gene_num\t$per\t".join("; ",@gene)."\n";
}
close OUT;
