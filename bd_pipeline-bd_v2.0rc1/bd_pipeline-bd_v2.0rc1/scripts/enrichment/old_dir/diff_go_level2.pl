#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#######################################################################################
my $env="module load OESingleCell/3.0.d && module load gcc/6.4.0";
#my $env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
#$env.=" && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
#######################################################################################

my ($infile, $outdir, $help);
GetOptions(
	"i:s"       =>\$infile,
	"o:s"       =>\$outdir,
	"h|help!"   =>\$help,
);

my $usage=<< "USAGE";
Program: $0
Description: difference genes GO level2 distribution
Options:
	-i    <infile>    The input file(enrichment-go-*.xls)       [Required]
	-o    <outdir>    The output directory of result            [Required]
	-h|help           print help info
Example:
	perl diff_go_level2.pl -i enrichment-go-Group2-vs-Group1-Total.xls -o ./

USAGE

die $usage if(!$infile || !$outdir || $help);
(-s $infile) || die "Error: don't find infile: $infile !\n";
(-d $outdir) || mkdir $outdir;
$infile=File::Spec->rel2abs($infile);
$outdir=File::Spec->rel2abs($outdir);
my $prefix;
if($infile=~/enrichment-go\S+(Total|Down|Up)\.\S+/i){
	$prefix=$1;
};
my $groupname=basename($infile);
my $mark=$1 if($groupname=~/-(Total|Down|Up)/);
$groupname=~s/^enrichment-go-//;
$groupname=~s/\.(xls|txt)$//;
$groupname=~s/-(Total|Down|Up)/\\\($1\\\)/;

my %handle;
my $BP; open $BP,">$outdir/BP.go2gene.$mark.txt" || die $!; $handle{"biological process"}=$BP; print $BP "go_id\n";
my $MF; open $MF,">$outdir/MF.go2gene.$mark.txt" || die $!; $handle{"molecular function"}=$MF; print $MF "go_id\n";
my $CC; open $CC,">$outdir/CC.go2gene.$mark.txt" || die $!; $handle{"cellular component"}=$CC; print $CC "go_id\n";

open IN,"<$infile" || die $!;
my %gene2go;my %class;my %go_class;
while(<IN>){
	chomp;
	next if(/^#/ || /^\s*$/ || /id\tterm/);
	my @l=split /\t/;
	my @gene=split /;/,$l[-1];
	for my $i (@gene){
		push(@{$gene2go{$i}},$l[0]);
	}
	$l[2]=~s/_/ /g;
	$class{$l[0]}=$l[2];
	my $out =$handle{$l[2]};
	print $out "$l[0]\n";
}
close IN;

my $total_gene=keys %gene2go;

system("$env && Rscript $Bin/GO.level2.r --bpGOterm $outdir/BP.go2gene.$mark.txt --mfGOterm $outdir/MF.go2gene.$mark.txt --ccGOterm $outdir/CC.go2gene.$mark.txt --level2 $Bin/endNode3.csv");

my %level2_def;
open EN,"<$Bin/endNode3.csv" || die $!;
while(<EN>){
	chomp;
	my @t=split /,/;
	$level2_def{$t[0]}=$t[1];
}
close EN;

my @file=("$outdir/BP.go2gene.$mark.ancestor.xls", "$outdir/MF.go2gene.$mark.ancestor.xls", "$outdir/CC.go2gene.$mark.ancestor.xls");
my %go_anno;
for my $f (@file){
	open TX,"<$f" || die $!;<TX>;
	while(<TX>){
		chomp;
		my @l=split /\t/,$_,2;
		if($l[1] eq ""){
			$l[1]=$l[0] if(exists $level2_def{$l[0]});
		}
		next if($l[0] eq "NA" || $l[1] eq "");
		my @ancestor=split /,/,$l[1];
		for my $ance (@ancestor){
			$go_anno{$l[0]}{$ance}="$l[0]\t$ance\t$level2_def{$ance}\t$class{$l[0]}";
		}
	}
	close TX;
}
undef %level2_def;

open OUT,">$outdir/gene2go.$mark.ancestor.xls" || die $!;
print OUT "#Gene\tGO_id\tGO_classify2\tGO_classify2_definition\tGO_classify1\n";
for my $i (sort keys %gene2go){
	for my $j (@{$gene2go{$i}}){
		for my $z (sort keys %{$go_anno{$j}}){
			print OUT "$i\t$go_anno{$j}{$z}\n";
		}
	}
}
close OUT;
undef %gene2go;
undef %go_anno;

open RE,"<$outdir/gene2go.$mark.ancestor.xls" || die $!;
my %level2_gene;my %level2_des;
while(<RE>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$level2_gene{$l[2]}{$l[0]}=1;
	$level2_des{$l[2]}="$l[4]\t$l[3]";
}
close RE;

open STAT,">$outdir/GO.level2.stat.$mark.xls" || die $!;
open TMP,">$outdir/GO.level2.stat.$mark.xls_tmp" || die $!;
print STAT "#GO_classify1\tGO_classify2\tdiff-$prefix\tGene\n#Total_gene\t\t$total_gene\t\n";
print TMP "#GO_classify1\tGO_classify2\tdiff-$prefix\n#Total_gene\t\t$total_gene\n";
for my $i(sort{$level2_des{$a} cmp $level2_des{$b}} keys %level2_des){
	my @lv2gene=sort keys %{$level2_gene{$i}};
	my $lv2_num=@lv2gene;
	print STAT "$level2_des{$i}\t$lv2_num\t".join("; ",@lv2gene)."\n";
	print TMP "$level2_des{$i}\t$lv2_num\n";
}
close STAT;
close TMP;
system("$env && Rscript $Bin/GOClassificationMap.r --infile $outdir/GO.level2.stat.$mark.xls_tmp --outpath $outdir --fileName GO.level2.stat.$mark --mark $groupname && rm $outdir/GO.level2.stat.$mark.xls_tmp");
system("rm -fr $outdir/BP.go2gene.$mark.txt $outdir/MF.go2gene.$mark.txt $outdir/CC.go2gene.$mark.txt $outdir/BP.go2gene.$mark.ancestor.xls $outdir/MF.go2gene.$mark.ancestor.xls $outdir/CC.go2gene.$mark.ancestor.xls $outdir/gene2go.$mark.ancestor.xls");
