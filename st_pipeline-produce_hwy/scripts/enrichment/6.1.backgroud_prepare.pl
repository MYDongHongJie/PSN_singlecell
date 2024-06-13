#!/usr/bin/env perl
# by The Coder, 20180110
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Cwd;

my ($gofile, $outdir, $help, $keggfile);
GetOptions(
	"g:s"       =>\$gofile,
	"k:s"       =>\$keggfile,
	"o:s"       =>\$outdir,
	"h|help!"   =>\$help,
);

my $usage=<< "USAGE";
Program: $0
Description: generate go.backgroud.xls,category.xls and kegg.backgroud.xls,anno-kegg.backgroud.xls
Options:
	-g    <infile>    The input file(prefix.GO.gene.anno.xls)   [Optional, -g -k at least one]
	-k    <infile>    The input file(prefix.KEGG.gene.anno.xls) [Optional, -g -k at least one]
	-o    <outdir>    The output directory of result            [Optional, default is ./]
	-h|help           print help info
Example:
	perl backgroud_prepare.pl -g unigene.GO.gene.anno.xls -k unigene.KEGG.gene.anno.xls -o ./

USAGE

die $usage if($help);
if(!$gofile && !$keggfile){
	print "Warn: -g -k at least one !\n";
	die $usage;
}
$outdir ||= getcwd;
if(defined $gofile){
	(-s $gofile) || die "Error: don't find infile: $gofile !\n";
	$gofile=File::Spec->rel2abs($gofile);
}
if(defined $keggfile){
	(-s $keggfile) || die "Error: don't find infile: $keggfile !\n";
	$keggfile=File::Spec->rel2abs($keggfile);
}
(-d $outdir) || mkdir $outdir;
$outdir=File::Spec->rel2abs($outdir);

if(defined $gofile){
	open INGO,"<$gofile" || die $!;
	open OUTGO,">$outdir/go.backgroud.xls" || die $!;
	my %category;
	while(<INGO>){
		chomp;
		next if(/^#/ || /^\s*$/);
		my @l=split /\t/;
		my @go=split /;\|/,$l[2];
		my $go_id_str; my $go_def_str;
		for my $i (@go){
			$i=~/^(\w+ \w+): (.+) \((GO:\d+)\)/;
			my $class=$1; my $go_def=$2; my $go_id=$3;
			$go_id_str.="$go_id,";
			$go_def_str.="$go_def|";
			$class=~s/ /_/g;
			$category{$go_id}=lc($class);
		}
		$go_id_str=~s/,$//g;
		$go_def_str=~s/\|$//g;
		print OUTGO "$l[0]\t$go_id_str\t$go_def_str\n";
	}
	close INGO;
	close OUTGO;
	open CATE,">$outdir/category.xls" || die $!;
	for my $i (sort keys %category){
		print CATE "$i\t$category{$i}\n";
	}
	close CATE;
}

if(defined $keggfile){
	open INKG,"<$keggfile" || die $!;
	open OUTKG,">$outdir/kegg.backgroud.xls" || die $!;
	open OUTAKG,">$outdir/anno-kegg.backgroud.xls" || die $!;
	while(<INKG>){
		chomp;
		next if(/^#/ || /^\s*$/);
		my @l=split /\t/;
		next if($l[7] eq "--" || $l[7] eq "");
		$l[7]=~s/\|/,/g;
		print OUTKG "$l[0]\t$l[7]\t$l[8]\n";
		print OUTAKG "$l[0]\t$l[4]\t$l[7]\n";
	}
	close INKG;
	close OUTKG;
	close OUTAKG;
}
