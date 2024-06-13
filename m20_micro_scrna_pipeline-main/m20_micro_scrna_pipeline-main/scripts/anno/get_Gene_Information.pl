#!/usr/bin/perl
use strict;
use warnings;
#The author: xiaoyue.wang\@majorbio.com

(@ARGV==4) || die "Usage: perl $0 [fasta] [ffn] [sample name] [outfile]\nDiscription: Statistical gene predict information\n";

my($fna,$ffn,$sample_name,$outfile)=@ARGV;

my $contig_num=0;my $leng_string=0;my ($num_g,$num_c,$num_a,$num_t,$num_n)=(0)x5;
(-s $fna) || die "$fna is not exists!\n";
(-s $ffn) || die "$ffn is not exists!\n";
open(FNA,"<$fna") or die $!;
while(<FNA>){
	chomp $_;
	if($_ =~ /^\>(\S+)/){
		$contig_num++;
	}else{
		$leng_string += length($_);
		$num_g += tr/gG/gG/;
		$num_c += tr/Cc/Cc/;
		$num_a += tr/Aa/Aa/;
		$num_t += tr/Tt/Tt/;
		$num_n += tr/Nn/Nn/;
	}
}
close(FNA);
my $total_len = $num_a+$num_t+$num_g+$num_c+$num_n;
if($leng_string != $total_len){
       die "The sequence file contains orther chars expect AGCTN\n";
}

my %orf;my $name;my $gene_num=0;my($gene_a,$gene_t,$gene_g,$gene_c,$gene_n)=(0)x5;
open(FFN,"<$ffn") or die $!;
while(<FFN>){
	chomp;
	if(/^>/){
		$_ =~ /^>(\S+)/;
		$name = $1;
		$orf{$name} = "";
		$gene_num++;
	}else{
		$_ =~ s/[^a-zA-z]//g;
		if($name ne ""){
			$orf{$name} .= $_;
		}
		$gene_a += tr/Aa/Aa/;
		$gene_t += tr/Tt/Tt/;
		$gene_g += tr/gG/gG/;
		$gene_c += tr/Cc/Cc/;
		$gene_n += tr/Nn/Nn/;
	}
}
my $gene_total = $gene_g +$gene_c+$gene_a+$gene_t+$gene_n;
my $gene_percent = sprintf("%.2f",($gene_total/$total_len)*100);
my $gene_average = sprintf("%.2f",$gene_total/$gene_num);
my $gc_gene = ($gene_g+$gene_c)/($gene_a+$gene_t+$gene_g+$gene_c)*100;
my $gc_gene_percent= sprintf("%.2f",$gc_gene);

#(-s $outfile) || system("echo -e \"Sample\\tGenome size(bp)\\tGene_num(#)\\ttotal_len(bp)\\taverage_len(bp)\\tGC% (gene region)\\tGene_len/Genome(%)\" > $outfile");
open OUT,">$outfile" || die $!;
print OUT "Sample\tGenome_size(bp)\tGene_num(#)\tTotal_len(bp)\tAverage_len(bp)\tGC%(gene region)\tGene_len/Genome(%)\n";
print OUT "$sample_name\t$total_len\t$gene_num\t$gene_total\t$gene_average\t$gc_gene_percent\t$gene_percent\n";
close OUT;
