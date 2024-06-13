#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use File::Spec;
use Cwd;

my $rfam_family="/home/sunkun/database/Rfam/family.txt";
my($outdir,$fasta,$prefix);

$outdir||=getcwd;
$genomesize||=0;
GetOptions(
    "fa:s"=>\$fasta,       
    "o:s"=>\$outdir,
    "p:s"=>\$prefix,
           );
sub usage{
    print qq{
This script is ncRNA workflow.
usage:
perl $0 -fa contig.fa -o outdir -g 4.5
options:
-fa         input contig fasta sequence
-o          output directory(default:$outdir)
-p          output file prefix
};
    exit;
}
if(!$fasta || !$prefix)
{
    &usage();
}

(-s $fasta) || die "Error: don't open $fasta !\n";
$fasta=File::Spec->rel2abs($fasta);
$outdir=File::Spec->rel2abs($outdir);
(-d $outdir) || mkdir $outdir;

open GM,"<$fasta" || die $!;
my %seq;my $seq_name;
while(<GM>){
	chomp;
	next if(/^\s*$/);
	if(/^>/){
		$seq_name=$1 if(/^>(\S+)/);
		$seq{$seq_name}="";
	}else{
		s/\s+//g;
		$seq{$seq_name}.=uc($_);
	}
}
close GM;

open(FA,"$rfam_family");
my %hash;
while(<FA>)
{
    chomp;
    my @array=split(/\t/,$_);
    if($array[8] !~/rRNA/ && $array[18]!~/tRNA/)
    {
         $array[18]=~s/; /:/g;
         $hash{$array[0]}=$array[18];
    }
}
system("less $outdir/$prefix\.Rfam_blast.out|awk '{if(\$0!~\"#\") print \$4,\"Rfam\",\$2,\$10,\$11,\$17,\$12,\".\",\$3}'>$outdir/$prefix\.rfam.out");

open(RF,"$outdir/$prefix\.rfam.out") || die $!;
open(RF2,">$outdir/$prefix\.sRNA.gff") || die $!;
open SA,">$outdir/$prefix\.sRNA.fna" || die $!;
my $mark="0001";
while(<RF>)
{
    chomp;
    if($_!~ /tRNA/ && $_!~ /rRNA/)
    {
        my @array=split;
	my ($start,$end)= ($array[3]<$array[4])?($array[3],$array[4]):($array[4],$array[3]);
	$array[3]=$start;
	$array[4]=$end;
        if(exists($hash{$array[$#array]})){
           $hash{$array[$#array]}=~s/\s+//;
           print RF2 "$array[0]\t$array[1]\tsRNA\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\tID=sRNA_$mark;Name=$array[2];class=$hash{$array[$#array]}\n";
	   my $str;
	   if($array[6]=~/\+/){
		$str=substr($seq{$array[0]},$array[3]-1,$array[4]-$array[3]+1);
	   }elsif($array[6]=~/-/){
		my $str_tmp=substr($seq{$array[0]},$array[3]-1,$array[4]-$array[3]+1);
		$str=reverse($str_tmp);
		$str=~tr/ATCG/TAGC/;
	   }
 	   print SA ">sRNA_$mark Name=$array[2];class=$hash{$array[$#array]}\n$str\n";
           $mark++;
        }
    }
}
close SA;
close RF;
close RF2;

system("cat $outdir/rRNA/$prefix\.rRNA.gff $outdir/$prefix\.sRNA.gff $outdir/tRNA/$prefix\.tRNA.gff >$outdir/$prefix\.ncRNA.gff");

