#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Cwd;

my $pilercr="/home/sunkun/Miniconda3/bin/pilercr";
my $crt="java -cp /home/sunkun/softwore/CRT1.2-CLI.jar crt";

my %opts;
GetOptions(\%opts, "fa=s", "od=s", "p=s", "h!");

my $usage =<< "USAGE";
Program: $0S
Descriptions: This script is CRISPR predict workflow.
Contact: xiaoyue.wang\@oebiotech.com
Usage : perl $0 [options]
        -fa     <file>       input genome fasta file [requred]
        -p      <prefix>     output file prefix [requred]
        -od     <dir>        output directory [default:./]
        -h                   print help document
Example: perl $0 -fa query.fasta -od ./ -p AO-01

USAGE

die $usage if($opts{h} || ! $opts{fa} || ! $opts{p});
$opts{od} ||= getcwd;
(-s $opts{fa}) || die "Error: don\'t open $opts{fa} file!";
my $outdir=File::Spec->rel2abs($opts{od});
my $genome=File::Spec->rel2abs($opts{fa});
$outdir=~s/\/$//g;(-d $outdir) || mkdir $outdir;

#####################pilercr predict
system("$pilercr -in $genome -out $outdir/$opts{p}\.pilercr.txt -noinfo");
open AA,"<$outdir/$opts{p}\.pilercr.txt" || die "Error: don't open file $outdir/$opts{p}\.pilercr.txt!\n";
my $open="close";my %pilercr;my $crisprid;my $equal=0;my $copy_num=0;my $chr;my $start;my $end;my $pilercr_num=0;
while(<AA>){
	chomp;
	if(/^DETAIL REPORT/){
		$open="open";
		next;
	}
	$open="close" if(/^SUMMARY BY SIMILARITY/);
	next if($open eq "close");
	next if(/^\s*$/);
	if(/^Array\s(\d+)/){
		$crisprid=$1;
		$pilercr_num++;
		next;
	}
	if(/^>(\S+)/){
		$chr=$1;
		next;
	}
	if(/^====/){
		$equal++;
		next;
	}
	s/^\s+//g;
	next if(/^Pos  Repeat/);
	my @l=split /\s+/,$_;
	if($equal==2){
		die "Error: Please check Copies(#) of Array $crisprid!\n" if($copy_num != $l[0]);
		$pilercr{$crisprid}{total}="crispr$crisprid\t$chr\t$l[0]\t$l[1]\t$l[2]\t-\t$start-$end\t$l[3]\t-\n";
		$equal=0;$copy_num=0;$chr=undef;$start=undef;$end=undef;$crisprid=undef;
	}else{
		my $str;
		$copy_num++;
		if(@l==6){
			my $spacer_len=length($l[-1]);
			$end=$l[0]+$l[1]+$spacer_len-1;
			$str="\t\t\t$l[1]\t$spacer_len\t$l[2]\t$l[0]\t$l[4]\t$l[5]\n";
		}else{
			$start=$l[0] if($copy_num==1);
			$str="\t\t\t$l[1]\t$l[3]\t$l[2]\t$l[0]\t$l[5]\t$l[6]\n";
		}
		$pilercr{$crisprid}{detail}.=$str;
	}
}
close AA;

open OUTA,">$outdir/$opts{p}\.pilercr.xls" || die $!;
print OUTA "CRISPR ID\tChr\tCopies(#)\tRepeat length(bp)\tSpacer length(bp)\t%id\tPosition\tRepeat\tSpacer\n";
for my $k(sort {$a<=>$b} keys %pilercr){
	print OUTA "$pilercr{$k}{total}";
	print OUTA "$pilercr{$k}{detail}";
}
close OUTA;

#####################crt predict
open FA,"<$genome" || die $!;
$/=">";my $len_sum=0;my %length;
open OUTFA,">$outdir/$opts{p}\.genome.tmp.fasta" || die $!;
print OUTFA ">$opts{p}\n";
while(<FA>){
	chomp;
	next if(/^\s*$/);
	my @b=split /\n/;
	my $id=(split /\s/,$b[0])[0];
	my $seq=join("",@b[1..$#b]);
	print OUTFA "$seq";
	my $len=length($seq);
	my $s=$len_sum+1;
	$len_sum+=$len;
	$length{$id}="$s-$len_sum";
}
close FA;
print OUTFA "\n";
close OUTFA;

system("$crt $genome $outdir/$opts{p}\.crt.txt");
system("$crt $outdir/$opts{p}\.genome.tmp.fasta $outdir/$opts{p}\.crt.txt");
$/="\n\n\n";
open BB,"<$outdir/$opts{p}\.crt.txt" || die "Error: don't open file $outdir/$opts{p}\.crt.txt!\n";
my %crt;my $crt_num=0;
while(<BB>){
	chomp;
	next if($_!~/CRISPR\s*(\d+)\s*Range:/);
	my($crt_id,$crt_start,$crt_end,$jian,$chr);
	my @line=split /\n/,$_;
	for my $l(@line){
		next if($l=~/^\s*$/ || $l=~/^POSITION/ || $l=~/^---/);
		if($l=~/^CRISPR\s*(\d+)\s*Range:\s*(\d+)\s*-\s*(\d+)/){
			$crt_id=$1;$crt_start=$2;$crt_end=$3;
			for my $i(sort keys %length){
				my ($tmp_s,$tmp_e)=split /-/,$length{$i};
				if($crt_start>=$tmp_s && $crt_end<=$tmp_e){
					$jian=$tmp_s-1;
					$chr=$i;
				}
			}
			last if(!defined $jian || !defined $chr);
			$crt_start=$crt_start-$jian;
			$crt_end=$crt_end-$jian;
			$crt_num++;
		}elsif($l=~/^Repeats:/){
			$l=~s/ //g;
			$l=~s/Repeats://g;
			$l=~s/AverageLength://g;
			$l=~s/\t\t/\t/g;
			$crt{$crt_id}{total}="crispr$crt_id\t$chr\t$l\t$crt_start-$crt_end\t-\t-\n";
		}else{
			my @tmp=split /\t+/,$l;
			$tmp[0]=$tmp[0]-$jian;
			if(@tmp==2){
				my $repeat_len=length($tmp[1]);
				$crt{$crt_id}{detail}.="\t\t\t$repeat_len\t0\t$tmp[0]\t$tmp[1]\t-\n";
			}else{
				$tmp[-1]=~s/[\[\] ]//g;
				my ($repeat_len,$spacer_len)=split /,/,$tmp[-1];
				$crt{$crt_id}{detail}.="\t\t\t$repeat_len\t$spacer_len\t$tmp[0]\t$tmp[1]\t$tmp[2]\n";
			}
		}
	}
}
close BB;
$/="\n";

open OUTB,">$outdir/$opts{p}\.crt.xls" || die $!;
print OUTB "CRISPR ID\tChr\tCopies(#)\tRepeat length(bp)\tSpacer length(bp)\tPosition\tRepeat\tSpacer\n";
for my $k(sort {$a<=>$b} keys %crt){
	print OUTB "$crt{$k}{total}";
	print OUTB "$crt{$k}{detail}";
}
close OUTB;

#total tab
(-s "$outdir/CRISPR.summary.xls") || system("echo -e \"Sample\\tpilercr(#)\\tCRT(#)\" > $outdir/CRISPR.summary.xls");
open TA,">>$outdir/CRISPR.summary.xls" || die $!;
print TA "$opts{p}\t$pilercr_num\t$crt_num\n";
close TA;

#per tab
open PER,">$outdir/$opts{p}.CRISPR.summary.xls" || die $!;
print PER "Sample\tPILER-CR(#)\tCRT(#)\n";
print PER "$opts{p}\t$pilercr_num\t$crt_num\n";
close PER;

system("rm -rf $outdir/$opts{p}\.pilercr.txt $outdir/$opts{p}\.crt.txt $outdir/$opts{p}\.genome.tmp.fasta");
#_end_
