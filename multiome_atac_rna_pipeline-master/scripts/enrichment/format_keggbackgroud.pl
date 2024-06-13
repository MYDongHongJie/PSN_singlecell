#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use FindBin qw($Bin);

############################################################################
my $pathwaydb="/home/luyao/10X_scRNAseq_v3/enrichment_of_different_expressed_gene/KEGGpathway_three_levels_v2.xls";
############################################################################

my ($infile, $outfile, $help);
GetOptions(
	"i:s"     => \$infile,
	"o:s"     => \$outfile,
	"h|help!" => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: total genes kegg level2 distribution
Options:
	-i     <infile>    The infile name(kegg.backgroud.xls)        [Required]
	-o     <outfile>   The output file                            [Required]
	-h|help            print help info
Example:
	perl format_keggbackgroud.pl -i kegg.backgroud.xls -o Unigene.KEGG.Classification.xls

USAGE

die $usage if(!$infile || !$outfile || $help);
(-s $infile) || die "Error: don't find infile: $infile !\n";
$infile=File::Spec->rel2abs($infile);
$outfile=File::Spec->rel2abs($outfile); 

my %pathway_class;
open PY,"<$pathwaydb" || die $!;
while(<PY>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$pathway_class{$l[0]}="$l[1]--$l[2]";
}
close PY;

open EN,"<$infile" || die $!;
my %class_gene;my %totalgene;
while(<EN>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	my @pathway=split /,/,$l[1];
	$totalgene{$l[0]}=1;
	for my $id (@pathway){
		$id=~s/[A-Za-z:]+(\d+)/ko$1/;
		$class_gene{$pathway_class{$id}}{$l[0]}=1;
	}
}
close EN;
undef %pathway_class;

my $total_genenum=keys %totalgene;
undef %totalgene;
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
