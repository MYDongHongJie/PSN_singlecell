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
my ($ftab,$foutdir,$fkey,$kobas,$mode);
GetOptions(
				"help|h|?" =>\&USAGE,
				"tab:s"=>\$ftab,
				"od:s"=>\$foutdir,
				"key:s"=>\$fkey,
				) or &USAGE;
&USAGE unless ($ftab and $foutdir and $fkey );
$ftab = Cwd::realpath($ftab);
$foutdir = Cwd::realpath($foutdir);
mkdir $foutdir if (!-d $foutdir);

#######################given a seven class number
open KOBARPLOT,">$foutdir/$fkey.KObarplot.R";
print KOBARPLOT "#!Rscript
library\(ggplot2\)
library\(reshape2\)
dat <- read.table\(\"$foutdir/$fkey.KEGG.classification.xls\",header=T,sep=\"\\t\",check.names=F\)
dat\$Classification_level2=factor(dat\$Classification_level2, levels=dat\$Classification_level2)
p <- ggplot\(dat\, aes(x=Classification_level2,y=gene_number,fill=Classification_level1)) + ylab\(\"Counts of Genes\"\) +xlab\(\"\"\) + labs(title=\"KEGG Classification\"\)
p <- p + geom_bar(stat=\"identity\",position=\"dodge\",width=0.7) +geom_text\(aes\(label=gene_number\),vjust=0.4,hjust=-0.2,size=2.5\) +coord_flip\(\) 
p <- p + theme\(panel.border=element_rect\(fill=NA,colour=\"black\"\)\)
p <- p + theme\(panel.background = element_rect\(fill=\"transparent\",colour=NA\),panel.grid.minor = element_blank\(\),panel.grid.major = element_blank\(\),plot.background = element_rect\(fill=\"transparent\",colour=NA\) \)
p <- p + theme(legend.title=element_blank())+ theme(legend.key.height=unit(0.5,\"cm\"),legend.text=element_text(size=8))
p <- p + theme(axis.text.y=element_text(size=8,color=\"black\")) + theme(axis.text.x=element_text(hjust=1, size=8,color=\"black\"))
ggsave\(\"$foutdir/$fkey.classification.pdf\",width=12,height=6,plot=p,bg=\"white\"\)
ggsave\(\"$foutdir/$fkey.classification.png\",type=\"cairo-png\",width=12,height=6,plot=p,bg=\"white\"\)\n";
close KOBARPLOT;
system("Rscript $foutdir/$fkey.KObarplot.R");
#system("rm -rf $foutdir/$fkey.KObarplot.R");

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
		-h           Help document
USAGE
	print $usage;
	exit;
}
