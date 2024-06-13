#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Spec;
use File::Basename;
my $repeatmasker="/home/sunkun/Miniconda3/bin/RepeatMasker";

my ($genome,$odir,$prefix);
GetOptions(
    "i:s"=>\$genome,
    "odir:s"=>\$odir,
    "pre:s"=>\$prefix,
);

my $usage=<<"USAGE";
Program: $0
Discription: predict repeat sequence by RepeatMasker software
Usage:
	-i      <file>        genome fasta [required]
	-pre    <str>         the prefix of outfile [required]
	-odir   <dir>         output dir [default: ./ ]

USAGE

die $usage if(!$genome || !$prefix);
$genome=File::Spec->rel2abs($genome);
$odir ||= getcwd;
$odir=File::Spec->rel2abs($odir);$odir=~s/\/$//g;
(-d $odir) || mkdir $odir;

my $dir_out="$odir/$prefix\_repeat_tmp";
( -e "$dir_out" ) || system("mkdir -p $dir_out");
my $gff="$odir/$prefix.repeat.gff";
my $file_name=basename($genome);
my $masker_out="$dir_out/$file_name\.out";
my $tbl="$dir_out/$file_name\.tbl";
my $xls="$odir/$prefix.repeat.xls";

&Repeat_masker($genome,$dir_out);
&Repeat_prepare($masker_out,$gff);
&tbl2xls($tbl,$xls);
system("cp $dir_out/$file_name\.masked $odir/$prefix\.masked.fna");
system("rm -rf $dir_out/* && rm -rf $dir_out/");

##################################################
#-------------------+
#    sub program    +
#-------------------+
sub Repeat_masker{
	my $genome=shift;
	my $dir=shift;
	my $species="AAA";
	my $pwd=getcwd;
	chdir("$dir");
	system("$repeatmasker $genome -dir $dir");
	chdir("$pwd");
}

sub Repeat_prepare{
	my $masker=shift;  ##$rmout
	my $gff=shift;  ##$gff="$dir/Bisulfite_Gff/Repeat.gff";
	open I,"<$masker" || die "Error: don\'t open file $masker!\n";
	open OUT,">$gff" || die $!;
	$_=<I>;
	$_=<I>;
	$_=<I>; #ǰ��������ȥ
	my $N = 1;
	while (<I>) {
		chomp;
		my @a = split;
		my $Chr = $a[4];
		my $Class = $a[10];
		$Class=~s/\?//g; #���ַ�����Ϣ������?������������滻��
		my $Repeat = $a[9];
		my @U = split/\//,$Class;
		my $Family = $U[0];
		my $Start = $a[5];
		my $End = $a[6];
		my $Strand = $a[8] ;
		if($a[8] eq "C"){
			$Strand = "-";
		};
		print OUT "$Chr\tRepeatMasker-v4.0.7\trepeat\t$Start\t$End\t.\t$Strand\t.\tID=Repeat$N;repeat=$Repeat;class=$Class\n";
		$N++;
		# body...
	};
	close I;
	close OUT;
};

sub tbl2xls{
	my $tbl_file=shift;
	my $xls_file=shift;
	open TBL,"<$tbl_file" || die "Error: don\'t open file $tbl_file!\n";
	open XLS,">$xls_file" || die $!;
	my $len;my $interspersed=0;
	print XLS "##sample name: $prefix\n";
	while(<TBL>){
		chomp;
		if(/^sequences:\s*(\d+)/){
			print XLS "##sequences num: $1\n";
		}elsif(/^total length:\s*(\d+)\s*bp/){
			$len=$1;
			print XLS "##total length: $len bp\n";
		}elsif(/^GC level:\s*([\d\.]+)\s*\%/){
			print XLS "##GC%: $1 %\n";
		}elsif(/^bases masked:\s*(\d+)\s*bp/){
			my $masked=$1;
			my $count=sprintf("%.2f",$masked*100/$len);
			print XLS "##bases masked: $masked bp ($count %)\n\n#class\tfamily\tnumber(#)\tlength(bp)\tpercentage(%)\n------\t------\t---------\t----------\t-------------\n";
		}elsif(/^[=-]+/ || /^\s+(number of|elements)/ || /^\s*$/ || /^file name:/){
			next;
		}elsif(/^\* most repeats/){
			last;
		}else{
			my $line=$_;
			$line=~s/ \s+/\t/g;
			$line=~s/ bp//g;$line=~s/ \%//g;
			my @t=split /\t/,$line;
			if(/^\S/){
				if(/^Total interspersed repeats:/){
					print XLS "$t[0]\t\t$interspersed\t$t[1]\t$t[2]\n------\t------\t---------\t----------\t-------------\n";
				}else{
					$interspersed=$interspersed+$t[1];
					print XLS "$t[0]\t\t".join("\t",@t[1..$#t])."\n";
				}
			}else{
				print XLS "$line\n";
			}
		}
	}
}

