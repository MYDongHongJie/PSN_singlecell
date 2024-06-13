#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Spec;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.1.0";
my $R="/data/software/Anaconda3/envs/mro-v3.4.3/bin/R";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$limit_max);

GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"m:s"=>\$limit_max,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
$limit_max||=10;
(-s $fIn) || die "Error: don\'t open file:$fIn!\n";
$fIn=File::Spec->rel2abs($fIn);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
&nr_pie_stat($fIn,$fOut);
my $pre=$fOut;
$pre=~s/\.[a-zA-Z]+$//;
open R,">$pre\.r" || die $!;
print R "library(ggplot2)
library(reshape2)
par(mar=c(2,2,2,2))
df<-read.table(\"$fOut\",header=T,sep=\"\\t\",check.names=F)
name = as.vector(df\$Species_Name)
pct=round(df\$Gene_Number/sum(df\$Gene_Number) * 100, 2)
myLabel = paste(name, \"(\", pct, \"%)\", sep = \"\")
p=ggplot(df, aes(x = \"\", y = Gene_Number, fill = Species_Name)) +
geom_bar(stat = \"identity\", width = 1) +
coord_polar(theta = \"y\") +
theme_bw() +
theme(axis.ticks = element_blank()) +
theme(axis.text.x = element_blank()) +
labs(x = \"\", y = \"\", title = \"Top $limit_max species distribution\") +
theme(plot.title = element_text(hjust = 0.5)) +
theme(legend.title = element_blank(), legend.position = \"right\") + 
scale_fill_discrete(breaks = df\$Species_Name, labels = myLabel) +
theme(panel.grid=element_blank()) + 
theme(panel.border=element_blank())
ggsave(\"$pre.pdf\",width=8,height=7,plot=p)
ggsave(\"$pre.png\",type=\"cairo-png\",width=8,height=7,plot=p)\n";

system("export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH && $R --restore --no-save < $pre\.r && rm -rf $pre\.r");

#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

################################################################################################################
sub nr_pie_stat{#&nr_pie_stat($in,$out)
	my ($in,$out)=@_;
	my %H; 
	open (IN,$in) or die $!;
	while(<IN>){
		chomp;
		next if(/^#/ || /^\s*$/);
		my @A=split/\t/,$_;;
		my $name=$A[-1];
		$name=~s/\s*$//;
		next unless($name=~/\[([^\]]+)\]$/);
		my $full_name=$1;
		$H{$full_name}++;
	}
	close (IN) ;
	
	my @name=sort {$H{$b} <=> $H{$a}} keys %H;
	my %count;my @name_sort;
	for(my $i=0;$i<@name;$i++){
		if($i<$limit_max){
			$count{$name[$i]}=$H{$name[$i]};
			push(@name_sort,$name[$i]);
		}else{
			$count{"other species"}+=$H{$name[$i]};
		}
	}
	push(@name_sort,"other species");

	open (OUT,">$out") or die $!;
	print OUT "Species_Name\tGene_Number\n";
	for my $k(@name_sort){
		print OUT "$k\t$count{$k}\n";
	}
	close OUT;
}
################################################################################################################
sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:   $0
     Version:	$version
 Description:	This program is used to statistics of the number of top ten species
       Usage:
		-i <infile>	NR/*.NR.blast.best.anno.xls
		-o <outfile>	NR/*.NR.species.topnum.xls
		-m <int>	maximum species number to show [10]
		-h		help

USAGE
	print $usage;
	exit;
}
