#!/usr/bin/perl
use warnings;
use warnings;
use File::Spec;

(@ARGV==2) || die "Discription: Conversion *blast.best.xls into *blast.best.anno.xls\nUsage: perl $0 [infile] [outfile]\n\n";
my ($infile,$outfile)=@ARGV;

(-s $infile) || die "Error: don\'t open file:$infile line $.!\n";
$infile=File::Spec->rel2abs($infile);
$outfile=File::Spec->rel2abs($outfile);

open IN,"<$infile" || die $!;
open OUT,">$outfile" || die $!;
print OUT "Cluster_id\tDatabase_ID\tE_value\tIdentity\tScore\tAnnotation\n";
while(<IN>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @anno = split /\t/, $_;
	print OUT "$anno[0]\t$anno[4]\t$anno[12]\t$anno[8]\t$anno[11]\t$anno[13]\n";
}
close IN;
close OUT;

