#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin);

my ($infile, $odir, $prefix, $Uniprot2GO, $help);
GetOptions(
	"i:s"      => \$infile,
	"o:s"      => \$odir,
	"p:s"      => \$prefix,
	"g:s"      => \$Uniprot2GO,
	"help|h!"  => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: The conversion of SwissProt result to GO
Options:
	-i     <file>      The infile of prefix.Swissprot.blast.best.anno.xls
	-o     <odir>      The output directory of result
	-p     <str>       The prefix of outfile
	-g     <file>      The infile of Uniprot2GO,
	                   [default: $Bin/GO/Uniprot2GO.txt]
	-h --help          print help info
Example:
	perl SwissProt2GO.pl -i Unigene.Swissprot.blast.best.anno.xls -o GO/ -p Unigene

USAGE

die $usage if(!$infile || !$odir || !$prefix || $help);
$Uniprot2GO ||= "$Bin/Uniprot2GO.txt";

(-s $infile) || die "Error: don't open infile: $infile !\n";
(-d $odir) || die "Error: don't find odir: $odir !\n";
(-s $Uniprot2GO) || die "Error: don't open Uniprot2GO: $Uniprot2GO !\n";
$infile=File::Spec->rel2abs($infile);
$odir=File::Spec->rel2abs($odir);
$Uniprot2GO=File::Spec->rel2abs($Uniprot2GO);

open IN,"<$infile" || die $!;
my %swissprot_gene; my @gene;
while(<IN>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$gene[@gene]=$l[0];
	$l[1]=~s/(\.\w+)$//;
	${$swissprot_gene{$l[1]}}[@{$swissprot_gene{$l[1]}}]=$l[0];
}
close IN;

my %class=("C"=>"Cellular Component", "F"=>"Molecular Function", "P"=>"Biological Process");
open GO,"<$Uniprot2GO" || die $!;
my %gene_go;my %class_gene;
while(<GO>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/,$_,6;
	next if($l[2] eq "");
	my @go_id=split /,/,$l[2];
	my @category=split /\|/,$l[3];
	my @go_term=split /\|/,$l[4];
	if(exists $swissprot_gene{$l[0]}){
		for my $g(@{$swissprot_gene{$l[0]}}){
			for(my $i=0; $i<=$#go_id; $i++){
				$gene_go{$g}{$go_id[$i]}="$class{$category[$i]}: $go_term[$i] ($go_id[$i]);";
				$class_gene{$class{$category[$i]}}{$g}=1;
			}
		}
		delete $swissprot_gene{$l[0]};
	}
	undef @go_id; undef @category; undef @go_term;
}
close GO;
undef %swissprot_gene;
undef %class;

open GENE2GO,">$odir/$prefix.GO.gene.anno.xls" || die $!;
print GENE2GO "#Gene\tNum\tGO_Anno\n";
my %go_gene;
for my $g(@gene){
	if(exists $gene_go{$g}){
		my $anno_num=scalar(keys %{$gene_go{$g}});
		my $anno;
		for my $k(sort keys %{$gene_go{$g}}){
			$anno.="$gene_go{$g}{$k}|";
			$go_gene{$gene_go{$g}{$k}}{$g}++;
		}
		$anno=~s/\|$//g;
		print GENE2GO "$g\t$anno_num\t$anno\n";
	}
}
close GENE2GO;
undef %gene_go;

open GO2GENE,">$odir/$prefix.GO.class.stat.xls" || die $!;
print GO2GENE "#GO_Function\tgene_number\tgene_ID\n";
my %num;
for my $i (sort keys %go_gene){
	if($i=~/(Biological Process|Cellular Component|Molecular Function)/){
		my $classname=$1;
		$num{$classname}++;
		if($num{$classname}==1){
			my $total_num=scalar(keys %{$class_gene{$classname}});
			print GO2GENE "$classname\t$total_num\n";
		}
	}
	my @anno_gene=sort keys %{$go_gene{$i}};
	my $gene_num=@anno_gene;
	my $gene=join(";",@anno_gene);
	print GO2GENE "$i\t$gene_num\t$gene\n";
}
close GO2GENE;
