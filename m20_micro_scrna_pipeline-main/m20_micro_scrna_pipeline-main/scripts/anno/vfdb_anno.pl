#!/usr/bin/perl
#xiaoyue.wang\@oebiotech.com     2017.10.09
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
	Discription: VFDB database annotation;
	Usage: perl $0 options
		-aa    <file>      input file (protein fasta), required
		-eval  <str>       blast E-value, default:1e-10
		-pre   <str>       output prefix, required
		-odir  <dir>       output file dir, default: ./
		-help              print this help
	Example: perl $0 -aa sample.faa -pre sample -odir VFDB/

USAGE

die $usage if(!$opts{aa} || !$opts{pre} || $opts{help});
$opts{eval} ||= "1e-10";
$opts{odir} ||= getcwd;
$opts{aa}=File::Spec->rel2abs($opts{aa});
$opts{odir}=File::Spec->rel2abs($opts{odir});
mkdir $opts{odir} unless( -d $opts{odir});
$opts{odir}=~s/\/$//g;

######################################################################################
my $diamond="/data/software/diamond/diamond-v0.9.7/diamond";
my $VFDB_DB="/data/database/micro/VFDB";
my $alignxml2tab="/home/zjming/metagenome/Anno/alignxml2tab.pl";
######################################################################################
(-s $opts{aa}) || die "Error: don\'t open file:$opts{aa}!\n";
(-s "$VFDB_DB/VFDB_setB_pro.dmnd") || die "Error: Please create a diamond index(VFDB_setB_pro.dmnd)!\n";
system("$diamond blastp -q $opts{aa} -d $VFDB_DB/VFDB_setB_pro.dmnd -p 5 --more-sensitive --max-target-seqs 10 --outfmt 5 -e $opts{eval} -o $opts{odir}/$opts{pre}.blastp.xml >/dev/null 2>&1");

(-s "$opts{odir}/$opts{pre}.blastp.xml") || die "Error: blastp failure, blast result file($opts{odir}/$opts{pre}.VFDB.blastp.xml) don\'t find!\n";
system("perl $alignxml2tab -tophit 1 -topmatch 1 -eval $opts{eval} $opts{odir}/$opts{pre}.blastp.xml > $opts{odir}/$opts{pre}.blastp.best.xls");

open XLS,"<$VFDB_DB/pro_geneList.xls" || die "Error: don\'t open $VFDB_DB/pro_geneList.xls!\n";
my %anno;
while(<XLS>){
	chomp;
	my @l=split /\t/;
	$anno{$l[1]}=$l[0]."\t".join("\t",@l[2..$#l]);
}
close XLS;

open BL,"<$opts{odir}/$opts{pre}.blastp.best.xls" || die "Error: don\'t open $opts{odir}/$opts{pre}.VFDB.blastp.best.xls!\n";
open OUT,">$opts{odir}/$opts{pre}.anno.xls" || die $!;
print OUT "Cluster_id\tIdentity\tScore\tE_value\tVFDB_type\tVFDB_Gene_ID\tGene_ID\tVF_type_ID\tVF_Description\tVF_Name\tVFDB_VF_ID\tVF_species\n";
while(<BL>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$l[0]=(split /\s/,$l[0])[0];
	print OUT "$l[0]\t$l[8]\t$l[11]\t$l[12]\t$anno{$l[4]}\n";
}
close BL;
close OUT;
system("cp $VFDB_DB/VFs.xls $opts{odir}/VFs_reference.xls");

#_end_
