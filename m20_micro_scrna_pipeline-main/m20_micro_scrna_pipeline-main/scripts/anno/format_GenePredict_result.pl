#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Cwd;
use File::Spec;

my $transeq="/home/sunkun/Miniconda3/bin/transeq";

my %opts;
GetOptions(\%opts, "fa=s", "gff=s", "od=s", "p=s", "h!");

$opts{od} ||= getcwd;

my $usage =<< "USAGE";
Program: $0
Descriptions: format results from the Prodigal to predict
Contact: xiaoyue.wang\@oebiotech.com
Usage:perl $0 [options]
	-fa   <file>    input genome fasta file
	-gff  <file>    input gene gff file
	-od   <dir>     output directory(dafault:./)
	-p    <prefix>  output file prefix(produce: prefix.ffn prefix.faa prefix.gff)
	-h              print help document
USAGE

die $usage if($opts{h} || ! $opts{fa} || ! $opts{gff} || ! $opts{p});

my $fa=File::Spec->rel2abs($opts{fa});
my $gff=File::Spec->rel2abs($opts{gff});
my $outdir=File::Spec->rel2abs($opts{od});

(-d $outdir) || mkdir $outdir;

open FA,"<$fa" || die "don\'t open $opts{fa} file !\n";
$/=">";my %genome_len;my %genome_seq;my $total_len=0;
while(<FA>){
	chomp;
	next if(/^\s*$/);
	my ($head, $seq)=split /\n/,$_,2;
	$seq=~s/\s//g;
	my $id =(split /\s/,$head)[0];
	$genome_seq{$id}=$seq;
	$total_len+=length($seq);
	$genome_len{$id}=$total_len-length($seq);
}
close FA;
$/="\n";	

open GFF,"<$gff" || die "don\'t open $opts{gff} file!\n";
open OG,">$outdir/$opts{p}.gene.gff" || die $!;
open ONT,">$outdir/$opts{p}.ffn" || die $!;
my $num=1;
while(<GFF>){
	chomp;
	next if(/^#/);
	my @l=split /\t/;
	my $ctg_num=$1 if($l[8]=~/ID=\d+_(\d+);/i);
	my $start_code=$1 if($l[8]=~/(start_type=[A-Z]+);/i);
	my $id=&addprefix($num);
	print OG (join("\t",@l[0..7]))."\tID=$id;$start_code;\n";
	my $gene_seq;
	if($l[6]=~/-/){
		$gene_seq=reverse(substr($genome_seq{$l[0]},$l[3]-1,$l[4]-$l[3]+1));
		$gene_seq=~tr/ATCGatcg/TAGCtagc/;
		print ONT ">$id ".($l[4]+$genome_len{$l[0]})." ".($l[3]+$genome_len{$l[0]})." $l[0]\_$ctg_num $l[4] $l[3]\n";
	}else{
		$gene_seq=substr($genome_seq{$l[0]},$l[3]-1,$l[4]-$l[3]+1);
		print ONT ">$id ".($l[3]+$genome_len{$l[0]})." ".($l[4]+$genome_len{$l[0]})." $l[0]\_$ctg_num $l[3] $l[4]\n";
	}
	print ONT "$gene_seq\n";
	$num++;
}
close GFF;
close OG;
close ONT;

(-s "$outdir/$opts{p}.ffn") || die "Error: don't open file $outdir/$opts{p}.ffn!\n";
system("$transeq -sequence $outdir/$opts{p}.ffn -outseq $outdir/$opts{p}.faa -table 11");
system("sed -i 's/_1 / /g' $outdir/$opts{p}.faa");
system("sed -i 's/\*\$//g' $outdir/$opts{p}.faa");

sub addprefix(){
	my $number=shift;
	my $index="0"x(4-length($number));
	my $id="gene".$index.$number;
	return $id;
}

