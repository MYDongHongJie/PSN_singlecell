#!/usr/bin/perl
#xiaoyue.wang\@oebiotech.com    2017.09.27
my $version="v1.0.0";
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Spec;

my %opts;
GetOptions(\%opts,"aa=s","eval=s","pre=s","odir=s","help!");

my $usage =<<"USAGE";
	Program: $0 $version
	Discription: PHI database annotation;
	Usage: perl $0 options
		-aa    <file>      input file (protein fasta), required
		-eval  <str>       blast E-value, default:1e-10
		-pre   <str>       output prefix, required
		-odir  <dir>       output file dir, default: ./
		-help              print this help
	Example: perl $0 -aa sample.faa -pre sample -odir PHI/

USAGE

die $usage if(!$opts{aa} || !$opts{pre} || $opts{help});
$opts{eval} ||= "1e-10";
$opts{odir} ||= getcwd;
$opts{aa}=File::Spec->rel2abs($opts{aa});
$opts{odir}=File::Spec->rel2abs($opts{odir});
mkdir $opts{odir} unless( -d $opts{odir});
$opts{odir}=~s/\/$//g;

######################################################################################
my $diamond="/home/sunkun/Miniconda3/bin/diamond";
my $PHI_DB="/home/sunkun/database/scMicro_refdata/PHI";
my $alignxml2tab="/home/sunkun/data/alignxml2tab.pl";
######################################################################################
(-s $opts{aa}) || die "Error: don\'t open file:$opts{aa}!\n";
system("$diamond blastp -q $opts{aa} -d $PHI_DB/phi_accessions.dmnd -p 5 --more-sensitive --max-target-seqs 10 --outfmt 5 -e $opts{eval} -o $opts{odir}/$opts{pre}.blastp.xml >/dev/null 2>&1");

(-s "$opts{odir}/$opts{pre}.blastp.xml") || die "Error: blastp failure, blast result file($opts{odir}/$opts{pre}.PHI.blastp.xml) don\'t find!\n";
system("perl $alignxml2tab -tophit 1 -topmatch 1 -eval $opts{eval} $opts{odir}/$opts{pre}.blastp.xml > $opts{odir}/$opts{pre}.blastp.best.xls");

open XLS,"<$PHI_DB/phi_accessions_abridgment.txt" || die "Error: don\'t open file:$PHI_DB/phi_accessions_abridgment.txt!\n";
my %disc;
while(<XLS>){
	chomp;
	my @a=split /\t/;
	my $str=$a[0]."\t".join("\t",@a[2..$#a]);
	$a[1]=~s/^\s//g;
	$a[1]=~s/\s$//g;
	$disc{$a[1]}.=$str."\n";
}
close XLS;

open BL,"<$opts{odir}/$opts{pre}.blastp.best.xls" || die "Error: don\'t open file:$opts{odir}/$opts{pre}.PHI.blastp.best.xls, may be alignxml2tab.pl run err!\n";
open OUT,">$opts{odir}/$opts{pre}.anno.xls" || die $!;
print OUT "Cluster_id\tIdentity\tScore\tE_value\tProtein_ID\t$disc{Protein_ID}";
while(<BL>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$l[0]=(split /\s/,$l[0])[0];
	if(exists $disc{$l[4]}){
		my $anno=$disc{$l[4]};
		$anno=~s/\n$//g;
		$anno=~s/\n/\n\t\t\t\t\t/g;
		print OUT "$l[0]\t$l[8]\t$l[11]\t$l[12]\t$l[4]\t$anno\n";
	}else{
		#print OUT "$l[0]\t$l[8]\t$l[11]\t$l[12]\t$l[4]\tnot find anno\n";
	}
}
close BL;
close OUT;

