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

my ($infile, $category, $outdir, $help, $sample);
GetOptions(
	"i:s"       =>\$infile,
	"c:s"       =>\$category,
	"s:s"       =>\$sample,
	"o:s"       =>\$outdir,
	"h|help!"   =>\$help,
);

my $usage=<< "USAGE";
Program: $0
Description: draw GO level2 distribution picture
Options:
	-i    <infile>    The input file(go.backgroud.xls)          [Required]
	-c    <infile>    The input file(category.xls)              [Required]
	-s    <prefix>    The sample name                           [Required]
	-o    <outdir>    The output directory of result            [Required]
	-h|help           print help info
Example:
	perl format_gobackgroud.pl -i go.backgroud.xls -c category.xls -s sample -o out/

USAGE

die $usage if(!$infile || !$category || !$sample || !$outdir || $help);
(-s $infile) || die "Error: don't find infile: $infile !\n";
(-s $category) || die "Error: don't find infile: $category !\n";
(-d $outdir) || mkdir $outdir;
$infile=File::Spec->rel2abs($infile);
$category=File::Spec->rel2abs($category);
$outdir=File::Spec->rel2abs($outdir);

my %handle;
my $BP; open $BP,">$outdir/$sample.BP.go2gene.txt" || die $!; $handle{"biological process"}=$BP; print $BP "go_id\n";
my $MF; open $MF,">$outdir/$sample.MF.go2gene.txt" || die $!; $handle{"molecular function"}=$MF; print $MF "go_id\n";
my $CC; open $CC,">$outdir/$sample.CC.go2gene.txt" || die $!; $handle{"cellular component"}=$CC; print $CC "go_id\n";

open CY,"<$category" || die $!;
my %class;
while(<CY>){
	chomp;
	my @l=split /\t/;
	$l[1]=~s/_/ /g;
	#$l[1]=~s/\b\w/\u$&/g;
	$class{$l[0]}=lc($l[1]);
}
close CY;

open IN,"<$infile" || die $!;
my %gene2go;my %total_gene;my %go_class;
while(<IN>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	my @go=split /,/,$l[1];
	for my $go_id (@go){
		push(@{$gene2go{$l[0]}},$go_id);
		$go_class{$go_id}=1;
	}
	$total_gene{$l[0]}=1;
}
close IN;
my $total_gene=keys %total_gene; undef %total_gene;

for my $i (sort keys %go_class){
	my $out =$handle{$class{$i}};
	print $out "$i\n";
}
undef %go_class;

system("$env && Rscript $Bin/GO.level2.r --bpGOterm $outdir/$sample.BP.go2gene.txt --mfGOterm $outdir/$sample.MF.go2gene.txt --ccGOterm $outdir/$sample.CC.go2gene.txt --level2 $Bin/endNode3.csv");

my %level2_def;
open EN,"<$Bin/endNode3.csv" || die $!;
while(<EN>){
	chomp;
	my @t=split /,/;
	$level2_def{$t[0]}=$t[1];
}
close EN;

my @file=("$outdir/$sample.BP.go2gene.ancestor.xls", "$outdir/$sample.MF.go2gene.ancestor.xls", "$outdir/$sample.CC.go2gene.ancestor.xls");
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

open OUT,">$outdir/$sample.gene2go.ancestor.xls" || die $!;
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

open RE,"<$outdir/$sample.gene2go.ancestor.xls" || die $!;
my %level2_gene;my %level2_des;
while(<RE>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	$level2_gene{$l[2]}{$l[0]}=1;
	$level2_des{$l[2]}="$l[4]\t$l[3]";
}
close RE;

open STAT,">$outdir/$sample.GO.classification.stat.xls" || die $!;
print STAT "#GO_classify1\tGO_classify2\t$sample\tGene\n#Total_gene\t\t$total_gene\t\n";
for my $i(sort{$level2_des{$a} cmp $level2_des{$b}} keys %level2_des){
	my @lv2gene=sort keys %{$level2_gene{$i}};
	my $lv2_num=@lv2gene;
	print STAT "$level2_des{$i}\t$lv2_num\t".join("; ",@lv2gene)."\n";
}
close STAT;
system("rm -fr $outdir/$sample.BP.go2gene.txt $outdir/$sample.MF.go2gene.txt $outdir/$sample.CC.go2gene.txt $outdir/$sample.BP.go2gene.ancestor.xls $outdir/$sample.MF.go2gene.ancestor.xls $outdir/$sample.CC.go2gene.ancestor.xls $outdir/$sample.gene2go.ancestor.xls");
