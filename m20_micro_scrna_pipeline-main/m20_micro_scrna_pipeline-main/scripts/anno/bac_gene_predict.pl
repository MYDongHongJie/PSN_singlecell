#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Spec;
use Cwd;

my $prodigal="/home/zjming/softwares/prodigal_v2.6.3/prodigal";
my $qsub="/data/software/qsub/qsub-sge.pl";
my $get_Gene_Information="/home/zhouxingya/Denovo_Pacbio/bac_pipeline/get_Gene_Information.pl";

my %opts;
GetOptions(\%opts, "fa=s", "od=s", "p=s", "h!");

$opts{od} ||= getcwd;

my $usage =<< "USAGE";
Program: $0
Descriptions: This script is gene predict workflow(software:Prodigal).
Contact: xiaoyue.wang\@oebiotech.com
Usage : perl $0 [options]
	-fa     <file>       input genome fasta file
	-od     <dir>        output directory(default:./)
	-p      <prefix>     output file prefix
	-h                   print help document
Example: perl $0 -fa query.fasta -od ./ -p AO-01
USAGE

die $usage if($opts{h} || ! $opts{fa} || ! $opts{p});

my $outdir=File::Spec->rel2abs($opts{od});
my $genome=File::Spec->rel2abs($opts{fa});
$outdir=~s/\/$//g;
(-d $outdir) || mkdir $outdir;

(-s $genome) || die "don\'t open $opts{fa} file!";
open SH,">$outdir/$opts{p}\.gene_predict.sh" || die $!;
print SH "$prodigal -p single -i $genome -d $outdir/$opts{p}\.temp.ffn -a $outdir/$opts{p}\.temp.faa -f gff -o $outdir/$opts{p}\.temp.gff && rm -rf $outdir/$opts{p}\.temp.faa $outdir/$opts{p}\.temp.ffn\n";
print SH "perl $Bin/format_GenePredict_result.pl -fa $genome -gff $outdir/$opts{p}\.temp.gff -od $outdir -p $opts{p} && rm -rf $outdir/$opts{p}\.temp.gff\n";
print SH "perl $get_Gene_Information $genome $outdir/$opts{p}.ffn $opts{p} $outdir/$opts{p}.gene.summary.xls\n";
close SH;
system("sh $outdir/$opts{p}\.gene_predict.sh");
