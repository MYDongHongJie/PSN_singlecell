#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Data::Dumper;
#use autodie;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.1";
#######################################################################################
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($ftab,$foutdir,$fkey,$kobas,$plant);
GetOptions(
				"help|h|?" =>\&USAGE,
				"tab:s"=>\$ftab,
				"od:s"=>\$foutdir,
				"key:s"=>\$fkey,
				"plant!"=>\$plant,
				"kos:s"=>\$kobas,
				) or &USAGE;
&USAGE unless ($ftab and $foutdir and $fkey and $kobas);
$ftab = Cwd::realpath($ftab);
$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);

########database file############
my $class7="/public/land/database/kegg/KEGGpathway_three_levels_v2.xls";
my $Rscript="/data/software/Anaconda3/envs/mro-v3.4.3/bin/Rscript";
my $anno="/public/land/database/kegg/KO2definition_nonredundancy.txt";
#######################given a seven class number
open MAP,"$class7" || die $!;
my %path_class;
while(<MAP>){
	chomp;
	my @array=split(/\t/,$_);
	next if($array[1] eq "Human Diseases");
	if((defined $plant) && ($array[1] eq "Organismal Systems")){
		next;
	}else{
		$path_class{$array[0]}="$array[1]--$array[2]";
	}
}
close MAP;

my %annotation;
open I,"$anno"  || die $!;
while(<I>){
	chomp;
	next if(/^\s*$/ || /^#/);
	my ($id,$tmp,$des)=split /\t/;
	$annotation{$id}=$des;
}
close I;

my %path_title;my %genes_ko;
my %result_path;my %ko_genename;my %query_path;
open I,"$kobas" || die $!;
local $/="\/\/\/\/";<I>;
while(<I>){
	chomp;
	next if (/^\s*$/);
	my $gName;my $path_num=0;my $path_filter=0;
	for my $line (split /\n/,$_){
		next if $line=~/^\s*$/;
		my @tmp = split /\t/, $line;
		$tmp[0]=~s/\s*//g;
		if($tmp[0] eq "Query:"){
			$gName=$tmp[1];
		}elsif($tmp[0] eq "KO:"){
			last if(! exists $annotation{$tmp[1]});
			$genes_ko{$gName}=$tmp[1];
			$ko_genename{$gName}{$tmp[1]}=$tmp[2];
		}else{
			$path_num++;
			if(exists $path_class{$tmp[3]}){
				$path_title{$tmp[3]}=$tmp[1];
				$result_path{$tmp[3]}{$gName}=$genes_ko{$gName};
				${$query_path{$gName}}[0].="$tmp[3]|";
				${$query_path{$gName}}[1].="$tmp[1]|";
			}else{
				$path_filter++;
			}
		}
	}
	if($path_num==$path_filter && $path_num!=0){
		delete $genes_ko{$gName};
		delete $ko_genename{$gName};
	}
}
close I;
$/="\n";

open TAB,"$ftab" or die "cannot open $ftab \n";
open OUTKO,">$foutdir/$fkey.gene.anno.xls" or die "cannot open $foutdir/$fkey.gene.anno.xls ,$!\n";
print OUTKO "#Gene_id\te_value\tidentity\tDatabase_GeneID\tKO\tgene_name\tDescription\tPathway\tPathway_definition\n";
while(<TAB>){
	chomp;
	next if (/^#/ || /^\s*$/);
	my @temp = split /\t/, $_;
	if(defined $genes_ko{$temp[0]}){
		if(defined ${$query_path{$temp[0]}}[0]){
			${$query_path{$temp[0]}}[0]=~s/\|$//g;
			${$query_path{$temp[0]}}[1]=~s/\|$//g;
		}else{
			${$query_path{$temp[0]}}[0]="--";
			${$query_path{$temp[0]}}[1]="--";
		}
		#$annotation{$temp[4]}="--" if(! exists $annotation{$temp[4]});
		#print OUTKO "$temp[0]\t$temp[-2]\t$temp[8]\t$temp[4]\t$annotation{$temp[4]}\t$genes_ko{$temp[0]}\t$ko_genename{$temp[0]}{$genes_ko{$temp[0]}}\t${$query_path{$temp[0]}}[0]\t${$query_path{$temp[0]}}[1]\n";
		print OUTKO "$temp[0]\t$temp[-2]\t$temp[8]\t$temp[4]\t$genes_ko{$temp[0]}\t$ko_genename{$temp[0]}{$genes_ko{$temp[0]}}\t$annotation{$genes_ko{$temp[0]}}\t${$query_path{$temp[0]}}[0]\t${$query_path{$temp[0]}}[1]\n";
	}
}
close TAB;
close OUTKO;

my %total_gene;
my %class_num;
open OUTPATH,">$foutdir/$fkey.pathway.xls" or die "cannot open $foutdir/$fkey.pathway.xls ,$!\n";
print OUTPATH "#class\tPathway_definition\tPathway\tGene_number\tGene_id\tKOs\n";
foreach my $key (sort keys %result_path){
	print OUTPATH "$path_class{$key}\t$path_title{$key}\t$key\t";
	my ($num, $realgenes,$kos);
	$num = 0;
	foreach my $ge (sort keys %{$result_path{$key}}){
		$num++;
		$realgenes .= "$ge;";
		$kos .= "$result_path{$key}{$ge};";
		$class_num{$path_class{$key}}{$ge}=1;
		$total_gene{$ge}=1;
	}
	$kos=~s/;$//g; $realgenes=~s/;$//g;
	print OUTPATH "$num\t$realgenes\t$kos\n";
	#$class_num{$path_class{$key}}+=$num;
}
close OUTPATH;
my $total_gene=keys %total_gene; undef %total_gene;


open STATS,">$foutdir/$fkey.classification.xls";
#print STATS "Group\tOrder\tPathway\tGene_number\n";
print STATS "Classification_level2\tClassification_level1\tgene_number\tpercentage\tGenes\n";
foreach my $key1 (sort keys %class_num){
	my($clas1,$clas2)=split /--/,$key1;
	my @gnum=sort keys %{$class_num{$key1}};
	my $gnums=@gnum;
	my $per=sprintf("%.2f",$gnums*100/$total_gene);
	#print STATS "$clas1\t$num\t$clas2\t$gnums\n";
	print STATS "$clas2\t$clas1\t$gnums\t$per\t".join("; ",@gnum)."\n";
}
close STATS;

open KOBARPLOT,">$foutdir/$fkey.KObarplot.R";
print KOBARPLOT "#!$Rscript
library\(ggplot2\)
library\(reshape2\)
dat <- read.table\(\"$foutdir/$fkey.classification.xls\",header=T,sep=\"\\t\",check.names=F\)
dat\$Classification_level2=factor(dat\$Classification_level2, levels=dat\$Classification_level2)
p <- ggplot\(dat\, aes(x=Classification_level2,y=gene_number,fill=Classification_level1)) + ylab\(\"Counts of Genes\"\) +xlab\(\"\"\) + labs(title=\"KEGG Classification\"\)
p <- p + geom_bar(stat=\"identity\",position=\"dodge\",width=0.7) +geom_text\(aes\(label=gene_number\),vjust=0.4,hjust=-0.2,size=2.5\) +coord_flip\(\) 
p <- p + theme\(panel.border=element_rect\(fill=NA,colour=\"black\"\)\)
p <- p + theme\(panel.background = element_rect\(fill=\"transparent\",colour=NA\),panel.grid.minor = element_blank\(\),panel.grid.major = element_blank\(\),plot.background = element_rect\(fill=\"transparent\",colour=NA\) \)
p <- p + theme(legend.title=element_blank())+ theme(legend.key.height=unit(0.5,\"cm\"),legend.text=element_text(size=8))
p <- p + theme(axis.text.y=element_text(size=8,color=\"black\")) + theme(axis.text.x=element_text(hjust=1, size=8,color=\"black\"))
ggsave\(\"$foutdir/$fkey.classification.pdf\",width=12,height=6,plot=p\)
ggsave\(\"$foutdir/$fkey.classification.png\",type=\"cairo-png\",width=12,height=6,plot=p\)\n";
close KOBARPLOT;
system("export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH && $Rscript $foutdir/$fkey.KObarplot.R");
system("rm -rf $foutdir/$fkey.KObarplot.R");

#######################################################################################
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub USAGE {
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version
	Description:	Extract KEGG Pathway and KOs file From Blast tab file;
	Usage:
		-tab         The blast_tab Out of seq with KEGG           must be given
		-od          Output dir                                   must be given
		-key         Prefix of OUT files (eg. key.path key.ko)    must be given
		-kos         The  out of KOBAS                            must be given
		-plant       If species is plant, please add this parameter
		-h           Help document
USAGE
	print $usage;
	exit;
}
